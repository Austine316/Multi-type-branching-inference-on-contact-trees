# Ensure to run this script covid_mcmc.R inside run_covid_mcmc.R 
# Do not RUN covid_mcmc.R as Standalone !!!

# =============================================================================
# covid_mcmc.R  —  Bayesian MCMC for R0 and k from Karnataka COVID-19 data
#                  Latent birth times + NegBin prior on k
# =============================================================================
# Sourced by run_covid_mcmc.R.
# Expects in global environment:
#   K_MAX, K_INIT, STEP_K,
#   NBINOM_MU, NBINOM_PHI,          (NegBin prior parameters for k)
#   POBS_FIX, LAM_FIX, INIT_R0,
#   N_ITER, N_BURN, THIN,
#   INIT_RW, TARGET_ACC, ADAPT_INT,
#   ODE_DT, SEED,
#   covid_edges (data frame),
#   neg_log_lik, precompute_ED_for_R0, ll_from_ED.
#
# THREE-STEP MCMC
# ---------------
# Step A: Gaussian random walk on log(R0)
# Step B: Uniform proposal for one latent time (birth time OR internal
#         branching time), chosen uniformly at random across all latent times
# Step C: Random walk on k with NegBin prior correction
#
# LATENT TIMES — TWO KINDS
# ------------------------
# Kind 1 — birth times (tau_born of each individual):
#   We observe detection time t_conf(i) but NOT when individual i was
#   infected. The birth time is latent:
#     tau_born(i) in [t_conf(parent(i)), t_conf(i)]   non-root individuals
#     tau_born(i) in [0,                 t_conf(i)]   root individuals
#   Updating tau_born(i) changes:
#     tau_a of individual i's first segment (seg_idx = 0)
#     tau_b of the parent's last segment (which ends at the child's birth)
#     All delta values on those two rows
#
# Kind 2 — internal branching times within a multi-transmission individual:
#   Individual with j_obs = m has m internal branching events at unknown times.
#   tau_branch(i,j) in [tau_born(i), t_conf(i)]
#   Updating one such time changes the two adjacent rows of individual i.
#
# PRIOR ON k: NEGATIVE BINOMIAL
# ------------------------------
# k ~ NegBin(mu = NBINOM_MU, size = NBINOM_PHI)
# This is overdispersed relative to Poisson, appropriate for heterogeneous
# contact networks and superspreading events in COVID-19.
# The NegBin prior does NOT cancel in the Hastings ratio even with a
# symmetric proposal — the prior ratio log[pi(k*)] - log[pi(k_cur)] must
# be included in the acceptance probability.
#
# ODE CACHING
# -----------
# The ODE depends on (R0, k) via beta = R0 * lam / k.
# It is solved once per proposal in Steps A and C, then cached.
# Step B reuses the current cache: no ODE solve needed.
# =============================================================================

suppressPackageStartupMessages(library(deSolve))

cat("=== COVID-19 Bayesian MCMC: R0, k, and latent times ===\n")
cat(sprintf("k prior: NegBin(mu=%.1f, phi=%.2f)  |  k in {1,...,%d}  |  k init = %d\n",
            NBINOM_MU, NBINOM_PHI, K_MAX, K_INIT))
cat(sprintf("lam = %.3f  |  pobs = %.3f\n", LAM_FIX, POBS_FIX))
cat(sprintf("Init R0 = %.3f  |  %d iter, %d burn, thin = %d\n\n",
            INIT_R0, N_ITER, N_BURN, THIN))

# =============================================================================
# 1. Priors
# =============================================================================
log_prior_R0 <- function(R0) {
  if (!is.finite(R0) || R0 <= 0) return(-Inf)
  dnorm(log(R0), mean = log(3), sd = 1, log = TRUE) - log(R0)
}

# Negative binomial prior on k: k ~ NegBin(mu=NBINOM_MU, size=NBINOM_PHI)
# Support: {1, 2, ..., K_MAX}  (truncated at K_MAX; k=0 excluded)
# The truncation normalisation constant is computed once at startup.
# log_prior_k <- function(k) {
#   if (k < 1L || k > K_MAX) return(-Inf)
#   dnbinom(k, size = NBINOM_PHI, mu = NBINOM_MU, log = TRUE) -
#     LOG_NBINOM_NORM
# }

log_prior_k <- function(k) {
  if (k < 1L || k > K_MAX) return(-Inf)
  -log(K_MAX)   # constant; cancels in MH ratio
}

# Precompute log normalisation constant for the truncated NegBin
nbinom_probs <- dnbinom(1L:K_MAX, size = NBINOM_PHI, mu = NBINOM_MU)
LOG_NBINOM_NORM <- log(sum(nbinom_probs))
rm(nbinom_probs)

# =============================================================================
# 2. Build per-replicate data structures (birth times + internal times)
# =============================================================================
raw <- covid_edges

build_rep_data_covid <- function(rid) {
  re   <- raw[raw$rep_id == rid, , drop = FALSE]
  tips <- re[!is.na(re$j_obs), ]
  if (nrow(tips) < 1L) return(NULL)
  
  # ── Kind 2: internal branching times within multi-transmission individuals ──
  # These are non-last segments of any individual with j_obs >= 1
  int_rows   <- which(re$seg_idx < (re$n_segs - 1L))
  c_rows_int <- integer(length(int_rows))
  p_rows_int <- integer(length(int_rows))
  lo_int     <- numeric(length(int_rows))
  hi_int     <- numeric(length(int_rows))
  valid_int  <- logical(length(int_rows))
  
  for (ii in seq_along(int_rows)) {
    pr   <- int_rows[ii]
    xid  <- re$ind_id[pr]
    nseg <- re$seg_idx[pr] + 1L
    cr   <- which(re$ind_id == xid & re$seg_idx == nseg)
    if (length(cr) != 1L) next
    lo <- re$tau_a[pr]; hi <- re$tau_b[cr]
    if (!is.finite(lo) || !is.finite(hi) || hi - lo < 1e-6) next
    p_rows_int[ii] <- pr; c_rows_int[ii] <- cr
    lo_int[ii]     <- lo; hi_int[ii]     <- hi
    valid_int[ii]  <- TRUE
  }
  keep_int <- which(valid_int)
  
  # ── Kind 1: birth times (tau_born) of each individual ────────────────────
  # For each individual, tau_born is latent in [lo_born, hi_born].
  # Updating tau_born changes:
  #   tau_a of this individual's first segment (first_row)
  #   tau_b of the parent's last segment (parent_last_row) — if it exists
  # lo_born = t_conf(parent) for non-root, or 0 for root
  # hi_born = t_conf(self) - small buffer so delta > 0
  all_ids  <- unique(re$ind_id)
  born_first_row <- integer(0L)   # first segment row of this individual
  born_par_row   <- integer(0L)   # parent's last segment row (NA if root)
  born_lo        <- numeric(0L)
  born_hi        <- numeric(0L)
  
  for (xid in all_ids) {
    rows_x    <- which(re$ind_id == xid)
    first_row <- rows_x[which.min(re$seg_idx[rows_x])]
    tau_born_cur <- re$tau_a[first_row]
    t_det_x      <- max(re$tau_b[rows_x], na.rm = TRUE)
    
    # Find parent's last segment row
    # Parent's last segment ends at the birth of this individual
    # We identify it as the row whose tau_b is closest to tau_born_cur
    # and whose ind_id != xid
    tol <- 1e-4
    par_candidates <- which(abs(re$tau_b - tau_born_cur) < tol &
                              re$ind_id != xid)
    par_last_row <- if (length(par_candidates) >= 1L) par_candidates[1L] else NA_integer_
    
    # Feasibility bounds
    if (is.na(par_last_row)) {
      # Root individual: lo = 0 (or a small positive number)
      lo_b <- 0.0
    } else {
      # lo = tau_a of parent's last segment (parent must exist before child born)
      lo_b <- re$tau_a[par_last_row]
    }
    hi_b <- t_det_x - 1e-4   # born strictly before own detection
    
    if (!is.finite(lo_b) || !is.finite(hi_b) || hi_b - lo_b < 1e-4) next
    
    born_first_row <- c(born_first_row, first_row)
    born_par_row   <- c(born_par_row,   if (is.na(par_last_row)) NA_integer_ else par_last_row)
    born_lo        <- c(born_lo,        lo_b)
    born_hi        <- c(born_hi,        hi_b)
  }
  
  list(
    rep_id       = rid,
    edges        = re,
    n_tips       = nrow(tips),
    maxTime      = max(re$tau_b, na.rm = TRUE) + 0.5,
    # Kind 2: internal branching times
    tau_int_init = if (length(keep_int) > 0L) re$tau_b[int_rows[keep_int]] else numeric(0L),
    c_row_int    = c_rows_int[keep_int],
    p_row_int    = p_rows_int[keep_int],
    lo_int       = lo_int[keep_int],
    hi_int       = hi_int[keep_int],
    n_int        = length(keep_int),
    lat_rows_int = p_rows_int[keep_int],
    # Kind 1: birth times
    tau_born_init = re$tau_a[born_first_row],
    first_row     = born_first_row,
    par_last_row  = born_par_row,
    lo_born       = born_lo,
    hi_born       = born_hi,
    n_born        = length(born_first_row)
  )
}

cat("Building replicate data structures...\n")
all_rids <- sort(unique(raw$rep_id))
rep_data <- Filter(Negate(is.null), lapply(all_rids, build_rep_data_covid))
N_rep    <- length(rep_data)

n_with_int  <- sum(sapply(rep_data, function(pr) pr$n_int  > 0L))
n_with_born <- sum(sapply(rep_data, function(pr) pr$n_born > 0L))
total_int   <- sum(sapply(rep_data, `[[`, "n_int"))
total_born  <- sum(sapply(rep_data, `[[`, "n_born"))

cat(sprintf("  %d replicates total\n", N_rep))
cat(sprintf("  %d with latent internal branching times (%d times total)\n",
            n_with_int, total_int))
cat(sprintf("  %d with latent birth times (%d times total)\n",
            n_with_born, total_born))
cat(sprintf("  mean n_tips = %.1f\n\n",
            mean(sapply(rep_data, `[[`, "n_tips"))))

if (N_rep == 0L) stop("No valid replicates. Check covid_data_prep.R output.")

GLOBAL_MAXTIME <- max(sapply(rep_data, `[[`, "maxTime"))

# =============================================================================
# 3. Edge mutation helpers
# =============================================================================

# Update one internal branching time within a multi-transmission individual
mutate_int <- function(edges, c_row, p_row, tau_new) {
  edges$tau_b[p_row] <- tau_new
  edges$delta[p_row] <- tau_new - edges$tau_a[p_row]
  edges$tau_a[c_row] <- tau_new
  edges$delta[c_row] <- edges$tau_b[c_row] - tau_new
  edges
}

# Update the birth time (tau_born) of one individual
# Changes: tau_a of first_row, tau_b of par_last_row (if not NA)
mutate_born <- function(edges, first_row, par_last_row, tau_new) {
  old_born           <- edges$tau_a[first_row]
  edges$tau_a[first_row] <- tau_new
  edges$delta[first_row] <- edges$tau_b[first_row] - tau_new
  if (!is.na(par_last_row)) {
    edges$tau_b[par_last_row] <- tau_new
    edges$delta[par_last_row] <- tau_new - edges$tau_a[par_last_row]
  }
  edges
}

# Build edges from current latent state for one replicate
build_edges_covid <- function(pr, tau_int, tau_born) {
  e <- pr$edges
  for (i in seq_along(tau_int))
    e <- mutate_int(e, pr$c_row_int[i], pr$p_row_int[i], tau_int[i])
  for (i in seq_along(tau_born))
    e <- mutate_born(e, pr$first_row[i], pr$par_last_row[i], tau_born[i])
  e
}

# Log-likelihood for one replicate given cached ED
ll_rep_cached <- function(pr, tau_int, tau_born, ED_cache) {
  edges    <- build_edges_covid(pr, tau_int, tau_born)
  lat_rows <- if (length(pr$lat_rows_int) > 0L) pr$lat_rows_int else NULL
  ll_from_ED(edges, ED_cache, latent_rows = lat_rows)
}

init_latent_random <- function(pr) {
  tau_int  <- if (pr$n_int > 0L)
    mapply(function(lo, hi) runif(1L, lo, hi), pr$lo_int,  pr$hi_int)
  else numeric(0L)
  tau_born <- if (pr$n_born > 0L)
    mapply(function(lo, hi) runif(1L, lo, hi), pr$lo_born, pr$hi_born)
  else numeric(0L)
  list(tau_int = tau_int, tau_born = tau_born)
}

# =============================================================================
# 4. Initialise
# =============================================================================
set.seed(SEED)

k_cur    <- K_INIT
R0_cur   <- INIT_R0
tau_int_cur  <- lapply(rep_data, function(pr) pr$tau_int_init)
tau_born_cur <- lapply(rep_data, function(pr) pr$tau_born_init)

cat("Precomputing initial ED cache...\n")
ED_cur <- precompute_ED_for_R0(R0_cur, k = k_cur, lam = LAM_FIX,
                               pobs = POBS_FIX,
                               maxTime = GLOBAL_MAXTIME, dt = ODE_DT)
if (is.null(ED_cur)) stop("Failed to compute initial ED. Check parameters.")

ll_cache <- sapply(seq_len(N_rep), function(i)
  ll_rep_cached(rep_data[[i]], tau_int_cur[[i]], tau_born_cur[[i]], ED_cur))
ll_sum <- sum(ll_cache)
lp_cur <- ll_sum + log_prior_R0(R0_cur) + log_prior_k(k_cur)

cat(sprintf("Initial log-posterior: %.2f (finite: %s)\n",
            lp_cur, is.finite(lp_cur)))

if (!is.finite(lp_cur)) {
  cat("Trying random tau initialisation...\n")
  for (.a in seq_len(300L)) {
    init_list    <- lapply(rep_data, init_latent_random)
    tau_int_cur  <- lapply(init_list, `[[`, "tau_int")
    tau_born_cur <- lapply(init_list, `[[`, "tau_born")
    ll_cache <- sapply(seq_len(N_rep), function(i)
      ll_rep_cached(rep_data[[i]], tau_int_cur[[i]], tau_born_cur[[i]], ED_cur))
    ll_sum <- sum(ll_cache)
    lp_cur <- ll_sum + log_prior_R0(R0_cur) + log_prior_k(k_cur)
    if (is.finite(lp_cur)) { cat(sprintf("  Finite lp at attempt %d\n", .a)); break }
  }
}
if (!is.finite(lp_cur))
  stop("Cannot initialise MCMC. Check parameters and edge table.")

cat(sprintf("Starting MCMC: %d iter, %d burn, thin = %d\n\n",
            N_ITER, N_BURN, THIN))

# =============================================================================
# 5. MCMC loop
# =============================================================================
sigma_rw <- INIT_RW
n_post   <- (N_ITER - N_BURN) %/% THIN
R0_chain <- numeric(n_post)
k_chain  <- integer(n_post)
post_idx <- 0L

acc_R0 <- 0L; prop_R0 <- 0L
acc_k  <- 0L; prop_k  <- 0L
acc_bt <- 0L; prop_bt <- 0L
t0 <- proc.time()

# Pre-identify replicates with each type of latent time
has_int  <- which(sapply(rep_data, function(pr) pr$n_int  > 0L))
has_born <- which(sapply(rep_data, function(pr) pr$n_born > 0L))

for (iter in seq_len(N_ITER)) {
  
  # ── Step A: propose new R0 ─────────────────────────────────────────────────
  log_R0_p <- log(R0_cur) + rnorm(1L, 0, sigma_rw)
  R0_p     <- exp(log_R0_p)
  prop_R0  <- prop_R0 + 1L
  
  if (is.finite(R0_p) && R0_p > 0.1 && R0_p < 20.0) {
    ED_p <- precompute_ED_for_R0(R0_p, k = k_cur, lam = LAM_FIX,
                                 pobs = POBS_FIX,
                                 maxTime = GLOBAL_MAXTIME, dt = ODE_DT)
    if (!is.null(ED_p)) {
      ll_p_vec <- sapply(seq_len(N_rep), function(i)
        ll_rep_cached(rep_data[[i]], tau_int_cur[[i]], tau_born_cur[[i]], ED_p))
      lp_p <- sum(ll_p_vec) + log_prior_R0(R0_p) + log_prior_k(k_cur)
      
      if (is.finite(lp_p) && log(runif(1L)) < lp_p - lp_cur) {
        R0_cur   <- R0_p;  ED_cur <- ED_p
        ll_cache <- ll_p_vec;  ll_sum <- sum(ll_cache)
        lp_cur   <- lp_p
        if (iter <= N_BURN) acc_R0 <- acc_R0 + 1L
      }
    }
  }
  
  # ── Step B: propose one latent time (birth time or internal time) ──────────
  # Choose uniformly over all latent times across both kinds and all replicates
  # Total latent times = total_born + total_int
  # We draw a type (born vs int) proportional to their counts, then a replicate
  any_latent <- (total_born + total_int) > 0L
  if (any_latent) {
    prop_bt <- prop_bt + 1L
    
    # Pick type proportional to count of each kind
    if (total_born == 0L) {
      pick_born <- FALSE
    } else if (total_int == 0L) {
      pick_born <- TRUE
    } else {
      pick_born <- (runif(1L) < total_born / (total_born + total_int))
    }
    
    if (pick_born && length(has_born) > 0L) {
      # Propose a new birth time for one individual in one replicate
      r      <- has_born[sample.int(length(has_born), 1L)]
      pr_r   <- rep_data[[r]]
      i_born <- sample.int(pr_r$n_born, 1L)
      tau_new <- runif(1L, pr_r$lo_born[i_born], pr_r$hi_born[i_born])
      
      tau_born_prop      <- tau_born_cur
      tau_born_prop[[r]][i_born] <- tau_new
      ll_r_new <- ll_rep_cached(pr_r, tau_int_cur[[r]], tau_born_prop[[r]], ED_cur)
      lp_prop  <- lp_cur - ll_cache[r] + ll_r_new
      
      if (is.finite(lp_prop) && log(runif(1L)) < lp_prop - lp_cur) {
        tau_born_cur[[r]] <- tau_born_prop[[r]]
        ll_cache[r]       <- ll_r_new
        ll_sum            <- sum(ll_cache)
        lp_cur            <- lp_prop
        if (iter <= N_BURN) acc_bt <- acc_bt + 1L
      }
      
    } else if (!pick_born && length(has_int) > 0L) {
      # Propose a new internal branching time
      r      <- has_int[sample.int(length(has_int), 1L)]
      pr_r   <- rep_data[[r]]
      i_node <- sample.int(pr_r$n_int, 1L)
      tau_new <- runif(1L, pr_r$lo_int[i_node], pr_r$hi_int[i_node])
      
      tau_int_prop             <- tau_int_cur
      tau_int_prop[[r]][i_node] <- tau_new
      ll_r_new <- ll_rep_cached(pr_r, tau_int_prop[[r]], tau_born_cur[[r]], ED_cur)
      lp_prop  <- lp_cur - ll_cache[r] + ll_r_new
      
      if (is.finite(lp_prop) && log(runif(1L)) < lp_prop - lp_cur) {
        tau_int_cur[[r]] <- tau_int_prop[[r]]
        ll_cache[r]      <- ll_r_new
        ll_sum           <- sum(ll_cache)
        lp_cur           <- lp_prop
        if (iter <= N_BURN) acc_bt <- acc_bt + 1L
      }
    }
  }
  
  
  # ── Step C: propose new k (random walk + NegBin prior correction) ──────────
  # k* = k_cur + delta, delta ~ Uniform{-STEP_K,...,STEP_K} \ {0}
  # Valid: k* in {1,...,K_MAX}.  Rows with j_obs > k* are penalised by
  # log(TINY) in log_lik_rep — no hard rejection based on data values.
  # Hastings ratio includes NegBin prior ratio since prior != proposal.
  prop_k <- prop_k + 1L
  deltas  <- setdiff(-STEP_K:STEP_K, 0L)
  k_star  <- k_cur + sample(deltas, 1L)
  
  if (k_star >= 1L && k_star <= K_MAX) {
    n_valid_cur  <- sum((k_cur  + deltas) >= 1L & (k_cur  + deltas) <= K_MAX)
    n_valid_star <- sum((k_star + deltas) >= 1L & (k_star + deltas) <= K_MAX)
    
    ED_k <- precompute_ED_for_R0(R0_cur, k = k_star, lam = LAM_FIX,
                                 pobs = POBS_FIX,
                                 maxTime = GLOBAL_MAXTIME, dt = ODE_DT)
    if (!is.null(ED_k) && n_valid_star > 0L) {
      ll_k_vec <- sapply(seq_len(N_rep), function(i)
        ll_rep_cached(rep_data[[i]], tau_int_cur[[i]], tau_born_cur[[i]], ED_k))
      ll_k_sum <- sum(ll_k_vec)
      
      # MH acceptance: likelihood ratio + Hastings correction + NegBin prior ratio
      log_alpha <- (ll_k_sum - ll_sum) +
        (log(n_valid_cur) - log(n_valid_star)) +
        (log_prior_k(k_star) - log_prior_k(k_cur))
      
      if (is.finite(log_alpha) && log(runif(1L)) < log_alpha) {
        k_cur    <- k_star;  ED_cur <- ED_k
        ll_cache <- ll_k_vec;  ll_sum <- ll_k_sum
        lp_cur   <- ll_sum + log_prior_R0(R0_cur) + log_prior_k(k_cur)
        if (iter <= N_BURN) acc_k <- acc_k + 1L
      }
    }
  }
  
  # ── Adaptive tuning of R0 step (burn-in only) ──────────────────────────────
  if (iter <= N_BURN && iter %% ADAPT_INT == 0L) {
    rate     <- acc_R0 / max(prop_R0, 1L)
    sigma_rw <- max(0.01, min(sigma_rw * exp(rate - TARGET_ACC), 1.0))
  }
  
  # ── Store ──────────────────────────────────────────────────────────────────
  if (iter > N_BURN && (iter - N_BURN) %% THIN == 0L) {
    post_idx           <- post_idx + 1L
    R0_chain[post_idx] <- R0_cur
    k_chain[post_idx]  <- k_cur
  }
  
  # ── Progress ───────────────────────────────────────────────────────────────
  if (iter %% max(1L, N_ITER %/% 10L) == 0L) {
    el <- (proc.time() - t0)[3L]
    cat(sprintf(
      "  iter %5d | R0 = %.3f  k = %2d | acc_R0 = %2.0f%%  acc_k = %2.0f%%  acc_bt = %2.0f%%  sigma_rw = %.3f | %.1fs\n",
      iter, R0_cur, k_cur,
      100 * acc_R0 / max(prop_R0, 1L),
      100 * acc_k  / max(prop_k,  1L),
      100 * acc_bt / max(prop_bt, 1L),
      sigma_rw, el))
  }
}

R0_chain <- R0_chain[seq_len(post_idx)]
k_chain  <- k_chain[seq_len(post_idx)]

# =============================================================================
# 6. Convergence diagnostics
# =============================================================================
hpdi <- function(x, p = 0.95) {
  x <- sort(x); n <- length(x)
  g <- max(1L, floor(p * n))
  w <- x[(g + 1L):n] - x[1L:(n - g)]
  lo <- which.min(w)
  c(lo = x[lo], hi = x[lo + g])
}

n1 <- max(1L, floor(0.10 * length(R0_chain)))
n2 <- max(1L, floor(0.50 * length(R0_chain)))
s1 <- R0_chain[1:n1]
s2 <- R0_chain[(length(R0_chain) - n2 + 1L):length(R0_chain)]
geweke_z <- (mean(s1) - mean(s2)) / sqrt(var(s1)/n1 + var(s2)/n2)

batch_sz <- max(10L, length(R0_chain) %/% 20L)
n_bat    <- length(R0_chain) %/% batch_sz
batches  <- sapply(1:n_bat, function(b)
  mean(R0_chain[((b-1L)*batch_sz + 1L):(b*batch_sz)]))
ess_R0   <- length(R0_chain) * var(R0_chain) / (batch_sz * var(batches))

ci_R0  <- hpdi(R0_chain)
k_tab  <- table(k_chain)
k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
k_top3 <- sort(k_tab, decreasing = TRUE)[1:min(3L, length(k_tab))]

# =============================================================================
# 7. Results summary
# =============================================================================
cat("\n")
cat(strrep("=", 65), "\n")
cat("POSTERIOR RESULTS — Karnataka COVID-19\n")
cat(strrep("=", 65), "\n")
cat(sprintf("R0 posterior mean:    %.3f\n",  mean(R0_chain)))
cat(sprintf("R0 posterior median:  %.3f\n",  median(R0_chain)))
cat(sprintf("R0 posterior SD:      %.3f\n",  sd(R0_chain)))
cat(sprintf("R0 95%% HPDI:         [%.3f,  %.3f]  (width = %.3f)\n",
            ci_R0["lo"], ci_R0["hi"], ci_R0["hi"] - ci_R0["lo"]))
cat(sprintf("Geweke z (R0):        %.3f  (|z| < 2 => converged)\n", geweke_z))
cat(sprintf("ESS (R0):             %.0f / %d  (%.0f%%)\n",
            ess_R0, length(R0_chain), 100 * ess_R0 / length(R0_chain)))
cat(strrep("-", 65), "\n")
cat(sprintf("k posterior mean:     %.2f\n",  mean(k_chain)))
cat(sprintf("k posterior SD:       %.2f\n",  sd(k_chain)))
cat(sprintf("k posterior mode:     %d\n",    k_mode))
cat(sprintf("k top-3 values:       %s\n",
            paste(sprintf("k=%s (%.1f%%)", names(k_top3),
                          100 * as.numeric(k_top3) / length(k_chain)), collapse = ", ")))
cat(strrep("-", 65), "\n")
cat(sprintf("pobs (fixed):         %.3f  (%d / %d directly detected)\n",
            POBS_FIX, 538L, 956L))
cat(sprintf("gamma (fixed):        1.0 model unit = %.2f days\n", 5.65))
cat(sprintf("k prior:              NegBin(mu=%.1f, phi=%.2f) on {1,...,%d}\n",
            NBINOM_MU, NBINOM_PHI, K_MAX))
cat(sprintf("Latent times:         %d birth times + %d internal times\n",
            total_born, total_int))
cat(sprintf("Implied beta:         %.4f  (= R0_mean * lam / k_mode)\n",
            mean(R0_chain) * LAM_FIX / k_mode))
cat(strrep("=", 65), "\n\n")

# =============================================================================
# 8. Save
# =============================================================================
saveRDS(list(
  R0_chain  = R0_chain,   k_chain    = k_chain,
  ci_R0     = ci_R0,      geweke_z   = geweke_z,  ess_R0 = ess_R0,
  k_mode    = k_mode,
  N_ITER    = N_ITER,     N_BURN     = N_BURN,    THIN   = THIN,
  K_MAX     = K_MAX,      K_INIT     = K_INIT,    STEP_K = STEP_K,
  NBINOM_MU = NBINOM_MU,  NBINOM_PHI = NBINOM_PHI,
  LAM_FIX   = LAM_FIX,   POBS_FIX   = POBS_FIX,
  INIT_R0   = INIT_R0,   ODE_DT     = ODE_DT,
  N_rep     = N_rep,      total_born = total_born, total_int = total_int,
  MEAN_SI   = 5.65
), "covid_mcmc_results.rds")
cat("Saved: covid_mcmc_results.rds\n")

write.csv(data.frame(iter = seq_along(R0_chain), R0 = R0_chain, k = k_chain),
          "covid_mcmc_chain.csv", row.names = FALSE)
cat("Saved: covid_mcmc_chain.csv\n\n")

# =============================================================================
# To reload without re-running
# =============================================================================

cat("  res <- readRDS('covid_mcmc_results.rds')\n")
cat("  R0_chain <- res$R0_chain;  k_chain <- res$k_chain\n\n")

res <- readRDS('covid_mcmc_results.rds')
R0_chain <- res$R0_chain;  k_chain <- res$k_chain
# =============================================================================
# 9. Diagnostic plots
# =============================================================================
C_BLACK <- "black"; C_DARK <- "gray20"; C_MID <- "gray50"

par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)

# ============= PLOT 1 Trace: R0 =================================
plot(R0_chain, type = "l", col = C_DARK, lwd = 0.5,
     xlab = "Post-burn-in sample", ylab = expression(R[0]), main = "")
abline(h = mean(R0_chain), col = C_BLACK, lty = 2, lwd = 2)
legend("topleft", bty = "n", cex = 1.0,
       legend = c(expression(R[0]~chain), 
                  sprintf("Post. mean")),
       lty = c(1, 2), col = c(C_DARK, C_BLACK), lwd = c(0.5, 2))

# ============= PLOT 2 Trace: muk =================================
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
plot(k_chain, type = "l", col = C_DARK, lwd = 0.5,
     xlab = "Post-burn-in sample", ylab = expression(mu[K]), main = "")
abline(h = mean(k_chain), col = C_BLACK, lty = 2, lwd = 2)
# abline(h = 11.8, col = C_BLACK, lty = 3, lwd = 2)
# abline(h = k_mode, col = C_BLACK, lty = 3, lwd = 2)

legend("top", bty = "n", cex = 0.8,
       legend = c(expression(mu[K]~chain), 
                  sprintf("Post. mean = %.1f", mean(k_chain))),
       lty = c(1, 2), col = c(C_DARK, C_BLACK), lwd = c(0.5, 2))

# ============= PLOT 3 Posterior Ro =================================
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
hist(R0_chain, xlab = expression(R[0]), 
     breaks = 40, 
     col = "gray70",
     probability = TRUE, main = " ", ylim = c(0, 0.5))
# Posterior mean
abline(v = mean(R0_chain), col = C_BLACK, lty = 2, lwd = 2)

# LogNormal(log 11, 0.5) prior for R0
x_grid <- seq(0.01, max(R0_chain) * 1.1, length.out = 500)
prior_dens <- dlnorm(x_grid, meanlog = log(3), sdlog = 1)
lines(x_grid, prior_dens, col = C_BLACK, lwd = 2, lty = 1)

legend("topright", bty = "n", cex = 1.0,
       legend = c("Post. mean", "Prior: LogNormal(log 3, 1)"),
       lty = c(2, 1), 
       col = c(C_BLACK, C_BLACK), 
       lwd = c(2, 2))

# ============= PLOT 4. Posterior vs prior: k =================================
k_prob <- as.numeric(k_tab) / sum(k_tab)
k_vals <- as.integer(names(k_tab))
k_prior_vals <- rep(1, K_MAX)
k_prior_vals <- k_prior_vals / sum(k_prior_vals)

ci_R0  <- hpdi(R0_chain)
k_tab  <- table(k_chain)
k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
k_top3 <- sort(k_tab, decreasing = TRUE)[1:min(3L, length(k_tab))]

par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
plot(k_vals, k_prob, type = "b", lwd = 2, col = C_DARK,
     xlab = expression(mu[K]), ylab = "Probability", main = "",
     ylim = c(0, max(k_prob) * 1.2)
)

x_grid <- 1L:K_MAX
prior_dens <- dlnorm(x_grid, meanlog = log(11), sdlog = 0.5)
lines(x_grid, prior_dens, type = "l", lwd = 1, col = C_DARK, lty = 1)
abline(v = k_mode, col = C_BLACK, lty = 3, lwd = 2)
abline(v = mean(k_chain), col = C_BLACK, lty = 2, lwd = 2)

legend("topright", bty = "n", cex = 1.0,
       legend = c("Posterior", "Post. mean", "Post. mode"),
       col = c(C_DARK, C_BLACK, C_BLACK), 
       lwd = c(2, 2, 2), 
       lty = c(1, 2, 3),
       pch = c(1, NA, NA)
)

# ============= PLOT 5. Auto correlation Ro =================================
ac_R0 <- acf(R0_chain, lag.max = 40, plot = FALSE)

# Compute ESS (if you haven't already; requires coda package)
if (!exists("ess_R0")) {
  library(coda)
  ess_R0 <- effectiveSize(R0_chain)
}

par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
plot(ac_R0$lag, ac_R0$acf, type = "b", lwd = 2.0, col = C_DARK,
     xlab = "Lag", ylab = "ACF", main = "")

abline(h =  1.96 / sqrt(length(R0_chain)), col = C_MID, lty = 2)
abline(h = -1.96 / sqrt(length(R0_chain)), col = C_MID, lty = 2)
abline(h = 0, col = C_BLACK)
 
legend("topright", 
       legend = bquote(ACF: ~ R[0] ~~ (ESS == .(round(ess_R0, 0)))),
       bty = "n", cex = 1.0)

# ============= PLOT 6. Joint posterior: R0 vs k ===============================
# Scatter plot with legend
idx_sub <- sample(seq_along(R0_chain), min(600L, length(R0_chain)))
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
plot(k_chain[idx_sub], R0_chain[idx_sub],
     pch = "*", col = C_DARK,
     xlab = expression(mu[K]), ylab = expression(R[0]),
     main = "")
abline(h = mean(R0_chain), col = C_BLACK, lty = 2, lwd = 2)
abline(v = mean(k_chain),          col = C_BLACK,   lty = 3, lwd = 2)
abline(v = k_mode,           col = C_BLACK,   lty = 4, lwd = 2)
legend("topright", bty = "n", cex = 1.0,
       legend = c("Posterior samples", 
                  expression(Mean ~ R[0]), 
                  expression(Mean ~ mu[K]),
                  expression(Mode ~ mu[K])
                  ),
       pch = c('*', NA, NA, NA), 
       lty = c(NA, 2, 3, 4), 
       col = c(C_DARK, C_BLACK, C_BLACK, C_BLACK),
       lwd = c(NA, 2, 2, 2))

par(mfrow = c(1, 1))
cat("Plots drawn — save via RStudio Export\n")
cat("=== DONE ===\n")
