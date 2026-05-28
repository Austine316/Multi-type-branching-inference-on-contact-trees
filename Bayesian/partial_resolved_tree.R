# =============================================================================
# partial_resolved_tree.R  —  Joint MCMC over (R0, k, branching times)
#                         with state marginalisation
# =============================================================================
#
# WHAT THIS IMPLEMENTS
# ─────────────────────────────────────────────────────────────────────────────
# Treats internal branching times {tau_m}, their states {s_m}, AND the contact
# degree k as jointly LATENT.  The MCMC has three alternating steps:
#
#   Step A: Gaussian random walk on log R0  [Metropolis-Hastings]
#   Step B: Uniform proposal for one latent branching time  [MH, symmetric]
#   Step C: Discrete MH proposal for k  [MH, independence sampler]
#
# States are analytically marginalised over s in {0,...,k} at each latent
# branching node (equation 5 of the likelihood section), so the MCMC samples
# only continuous and integer-valued unknowns.
#
# PRIORS
# ─────────────────────────────────────────────────────────────────────────────
# R0   ~ LogNormal(log(5), 1)      weakly informative; median 5, 95% in [0.7, 36]
# k    ~ DiscreteUniform{1,...,K_MAX}   flat over plausible contact degrees
# tau_m ~ Uniform[lo_m, u_m]           maximally non-informative within feasible set
#
# WHY ESTIMATE k?
# ─────────────────────────────────────────────────────────────────────────────
# In the full-tree MLE, k is selected by AIC, which fixes it at a point
# estimate and ignores uncertainty.  In the Bayesian framework, k is a model
# parameter whose uncertainty should propagate to the posterior of R0.
# Moreover, the coupling k*beta ~ const means that misspecifying k introduces
# a compensating bias in beta but not in R0 (since R0 = k*beta/gamma).
# Estimating k jointly therefore provides:
#   (a) a posterior distribution over contact degrees (useful when k is unknown)
#   (b) exact propagation of k uncertainty into the R0 posterior
#   (c) a check: if the data are informative about k, the posterior will
#       concentrate near the true value; otherwise it will remain diffuse
#
# ODE CACHING
# ─────────────────────────────────────────────────────────────────────────────
# The backward Kolmogorov ODE depends on (R0, k) via beta = R0*gamma/k.
# It is solved ONCE per proposal in Steps A and C and cached as an ED_cache
# object.  Step B reuses the current cache without any new ODE solve.
# This limits the ODE solve cost to at most 2 per iteration.
#
# OUTPUT
# ─────────────────────────────────────────────────────────────────────────────
# Saves per p_res level:
#   partial_tree_chains_p{p_res}.rds  — list(R0_chain, k_chain, p_res)
#   partial_tree_summary.csv          — posterior summaries across p_res levels
# =============================================================================

suppressPackageStartupMessages(library(deSolve))
source("parameter_estimate_bayesian_unknown_internals_unknown_k.R")

cat('Current working directory', getwd())
setwd(getwd())

# =============================================================================
# 0. Settings
# =============================================================================
LAM_FIX    <- 1.0
POBS_FIX   <- 0.5
TRUE_R0    <- 6.0
TRUE_K     <- 4L        # true value for simulation (not used in inference)
K_INIT     <- 4L        # starting value for k chain
K_MAX      <- 12L       # upper bound for discrete uniform prior on k

if (!exists("N_REPS_USE"))  N_REPS_USE  <- 30L
if (!exists("N_ITER"))      N_ITER      <- 2000L
if (!exists("N_BURN"))      N_BURN      <- 400L
if (!exists("THIN"))        THIN        <- 2L
if (!exists("INIT_RW"))     INIT_RW     <- 0.12
if (!exists("TARGET_ACC"))  TARGET_ACC  <- 0.23
if (!exists("ADAPT_INT"))   ADAPT_INT   <- 50L
if (!exists("ODE_DT"))      ODE_DT      <- 0.02
if (!exists("SEED"))        SEED        <- 42L
if (!exists("SAVE_CHAINS")) SAVE_CHAINS <- TRUE   # save .rds per p_res level

if (!exists("P_RES_VALS"))
  P_RES_VALS <- c(0.0, 0.25, 0.5, 0.75, 1.0)

FULL_TREE_MLE <- 6.031
FULL_TREE_CI  <- c(5.731, 6.331)

C_BLACK <- "black"; C_DARK <- "gray20"; C_MID <- "gray50"; C_FILL <- "gray88"

cat("=== Joint MCMC: R0, k, branching times (with state marginalisation) ===\n")
cat(sprintf("Fixed: lam=%.1f, pobs=%.1f | True R0=%.1f, true k=%d\n",
            LAM_FIX, POBS_FIX, TRUE_R0, TRUE_K))
cat(sprintf("k prior: DiscreteUniform{1,...,%d}  |  k init = %d\n",
            K_MAX, K_INIT))
cat(sprintf("p_res grid: %s\n",
            paste(sprintf("%.2f", P_RES_VALS), collapse = ", ")))
cat(sprintf("MCMC: %d iter (%d burn, thin=%d) | %d reps per level\n\n",
            N_ITER, N_BURN, THIN, N_REPS_USE))

# =============================================================================
# 1. Priors
# =============================================================================
# log prior for R0: LogNormal(log(5), 1)
log_prior_R0 <- function(R0) {
  if (R0 <= 0) return(-Inf)
  dnorm(log(R0), mean = log(5), sd = 1, log = TRUE) - log(R0)  # log-normal
}

# log prior for k: Discrete Uniform{1,...,K_MAX}
log_prior_k <- function(k) {
  if (k < 1L || k > K_MAX) return(-Inf)
  -log(K_MAX)   # log(1/K_MAX)
}

# =============================================================================
# 2. Load data
# =============================================================================
cat("Loading edge table...\n")
raw_full <- read.csv("full_tree_edges.csv", stringsAsFactors = FALSE)
raw_full <- raw_full[is.finite(raw_full$delta) & raw_full$delta >= 0 &
                       !is.nan(raw_full$s) & !is.nan(raw_full$c), ]
raw_full$j_obs[is.nan(raw_full$j_obs)] <- NA_integer_
raw <- raw_full[raw_full$n_tips > 0, ]

use_reps <- sort(unique(raw$rep_id))[
  seq_len(min(N_REPS_USE, length(unique(raw$rep_id))))]
cat(sprintf("  Using %d replicates\n\n", length(use_reps)))

# =============================================================================
# 3. Per-replicate data structures
# =============================================================================
build_rep_data <- function(rid) {
  re      <- raw[raw$rep_id == rid, , drop = FALSE]
  re_full <- raw_full[raw_full$rep_id == rid, , drop = FALSE]
  tips    <- re[re$edge_type == "tip" & !is.na(re$j_obs), ]
  if (nrow(tips) < 2L) return(NULL)
  
  int_rows <- which(re$seg_idx == 0L & re$edge_type == "internal")
  if (length(int_rows) == 0L) return(NULL)
  
  child_rows  <- integer(length(int_rows))
  parent_rows <- integer(length(int_rows))
  lo_bounds   <- numeric(length(int_rows))
  hi_bounds   <- numeric(length(int_rows))
  valid       <- logical(length(int_rows))
  
  for (ii in seq_along(int_rows)) {
    cr       <- int_rows[ii]
    tau_born <- re$tau_a[cr]
    xid      <- re$ind_id[cr]
    
    pr_full <- which(abs(re_full$tau_b - tau_born) < 1e-7 &
                       re_full$ind_id != xid)
    if (length(pr_full) == 0L) next
    
    pf <- re_full[pr_full[1L], ]
    pr <- which(abs(re$tau_b  - pf$tau_b)  < 1e-9 &
                  re$ind_id  == pf$ind_id  &
                  re$seg_idx == pf$seg_idx)
    if (length(pr) != 1L) next
    
    lo <- re$tau_a[pr]; hi <- re$tau_b[cr]
    if (hi - lo < 1e-8) next
    
    child_rows[ii]  <- cr; parent_rows[ii] <- pr
    lo_bounds[ii]   <- lo; hi_bounds[ii]   <- hi
    valid[ii]       <- TRUE
  }
  
  keep <- which(valid)
  if (length(keep) == 0L) return(NULL)
  
  list(
    rep_id     = rid,
    edges      = re,
    n_tips     = nrow(tips),
    maxTime    = max(re$tau_b, na.rm = TRUE) + 0.5,
    tau_true   = re$tau_a[int_rows[keep]],
    child_row  = child_rows[keep],
    parent_row = parent_rows[keep],
    lo         = lo_bounds[keep],
    hi         = hi_bounds[keep],
    n_int      = length(keep)
  )
}

cat("Building replicate data structures...\n")
rep_data <- Filter(Negate(is.null), lapply(use_reps, build_rep_data))
N_rep    <- length(rep_data)
cat(sprintf("  %d valid replicates | mean n_tips=%.1f | mean n_int=%.1f\n\n",
            N_rep,
            mean(sapply(rep_data, `[[`, "n_tips")),
            mean(sapply(rep_data, `[[`, "n_int"))))

# Global maxTime
GLOBAL_MAXTIME <- max(sapply(rep_data, `[[`, "maxTime"))

# =============================================================================
# 4. MCMC helper functions
# =============================================================================
mutate_edge <- function(edges, child_row, parent_row, tau_new) {
  edges$tau_b[parent_row] <- tau_new
  edges$delta[parent_row] <- tau_new - edges$tau_a[parent_row]
  edges$tau_a[child_row]  <- tau_new
  edges$delta[child_row]  <- edges$tau_b[child_row] - tau_new
  edges
}

build_edges <- function(pr, tau_latent, latent_idx) {
  e <- pr$edges
  for (i in seq_along(latent_idx)) {
    m <- latent_idx[i]
    e <- mutate_edge(e, pr$child_row[m], pr$parent_row[m], tau_latent[i])
  }
  e
}

# Log-likelihood for one replicate using cached ED (no ODE solve)
ll_rep_cached <- function(pr, tau_latent, latent_idx, ED_cache) {
  edges    <- build_edges(pr, tau_latent, latent_idx)
  lat_rows <- if (length(latent_idx) > 0L) pr$parent_row[latent_idx] else NULL
  ll_from_ED(edges, ED_cache, latent_rows = lat_rows)
}

# Log-likelihood summed over all replicates
lp_all_cached <- function(tau_list, latent_idx_list, ED_cache) {
  total <- 0.0
  for (i in seq_len(N_rep)) {
    ll <- ll_rep_cached(rep_data[[i]], tau_list[[i]],
                        latent_idx_list[[i]], ED_cache)
    total <- total + if (is.finite(ll)) ll else log(1e-300)
  }
  total
}

# Random tau initialisation
init_latent <- function(pr, latent_idx) {
  if (length(latent_idx) == 0L) return(numeric(0L))
  tau <- numeric(length(latent_idx))
  for (i in seq_along(latent_idx))
    tau[i] <- runif(1L, pr$lo[latent_idx[i]], pr$hi[latent_idx[i]])
  tau
}

# =============================================================================
# 5. Assign known / latent nodes for a given p_res
# =============================================================================
make_partial <- function(rd, p_res) {
  n_known    <- round(p_res * rd$n_int)
  set.seed(SEED + rd$rep_id)
  known_idx  <- if (n_known == 0L) integer(0L) else
    sort(sample.int(rd$n_int, n_known))
  latent_idx <- setdiff(seq_len(rd$n_int), known_idx)
  list(
    rep_id     = rd$rep_id,
    edges_base = rd$edges,   # base edge table (true tau values)
    n_tips     = rd$n_tips,
    maxTime    = rd$maxTime,
    tau_true   = rd$tau_true,
    child_row  = rd$child_row,
    parent_row = rd$parent_row,
    lo         = rd$lo,
    hi         = rd$hi,
    n_int      = rd$n_int,
    known_idx  = known_idx,
    latent_idx = latent_idx
  )
}

# =============================================================================
# 6. Main MCMC: joint over (log R0, k, {tau_latent})
# =============================================================================
run_mcmc <- function(pr_list, p_res) {
  
  NR          <- length(pr_list)
  any_latent  <- any(sapply(pr_list, function(pr) length(pr$latent_idx) > 0L))
  n_latent    <- sum(sapply(pr_list, function(pr) length(pr$latent_idx)))
  latent_idxs <- lapply(pr_list, `[[`, "latent_idx")
  
  # ── Initialise ──────────────────────────────────────────────────────────────
  tau_cur  <- lapply(pr_list, function(pr) pr$tau_true[pr$latent_idx])
  k_cur    <- K_INIT
  R0_cur   <- TRUE_R0     # start at true value for numerical stability
  
  ED_cur <- precompute_ED_for_R0(R0_cur, k = k_cur, lam = LAM_FIX,
                                 pobs = POBS_FIX,
                                 maxTime = GLOBAL_MAXTIME, dt = ODE_DT)
  if (is.null(ED_cur)) {
    cat(sprintf("  [p=%.2f] Cannot compute initial ED — skipping\n", p_res))
    return(NULL)
  }
  
  ll_cache <- sapply(seq_len(NR), function(i)
    ll_rep_cached(pr_list[[i]], tau_cur[[i]], latent_idxs[[i]], ED_cur))
  ll_sum   <- sum(ll_cache)
  lp_cur   <- ll_sum + log_prior_R0(R0_cur) + log_prior_k(k_cur)
  
  # Fallback: random tau if initial point is infeasible
  if (!is.finite(lp_cur)) {
    for (.a in seq_len(100L)) {
      tau_cur  <- lapply(pr_list, function(pr) init_latent(pr, pr$latent_idx))
      ll_cache <- sapply(seq_len(NR), function(i)
        ll_rep_cached(pr_list[[i]], tau_cur[[i]], latent_idxs[[i]], ED_cur))
      ll_sum <- sum(ll_cache)
      lp_cur <- ll_sum + log_prior_R0(R0_cur) + log_prior_k(k_cur)
      if (is.finite(lp_cur)) break
    }
  }
  if (!is.finite(lp_cur)) {
    cat(sprintf("  [p=%.2f] Cannot initialise — skipping\n", p_res))
    return(NULL)
  }
  
  cat(sprintf("  [p=%.2f] init lp=%.1f | k=%d | %d latent nodes\n",
              p_res, lp_cur, k_cur, n_latent))
  
  sigma_rw <- INIT_RW
  n_post   <- (N_ITER - N_BURN) %/% THIN
  R0_chain <- numeric(n_post)
  k_chain  <- integer(n_post)
  post_idx <- 0L
  
  acc_R0 <- 0L; prop_R0 <- 0L
  acc_k  <- 0L; prop_k  <- 0L
  acc_bt <- 0L; prop_bt <- 0L
  t0 <- proc.time()
  
  for (iter in seq_len(N_ITER)) {
    
    # ── Step A: propose new R0 (Gaussian RW on log scale) ───────────────────
    log_R0_p <- log(R0_cur) + rnorm(1L, 0, sigma_rw)
    R0_p     <- exp(log_R0_p)
    prop_R0  <- prop_R0 + 1L
    
    if (is.finite(R0_p) && R0_p > 0.2 && R0_p < 30.0) {
      ED_p <- precompute_ED_for_R0(R0_p, k = k_cur, lam = LAM_FIX,
                                   pobs = POBS_FIX,
                                   maxTime = GLOBAL_MAXTIME, dt = ODE_DT)
      if (!is.null(ED_p)) {
        ll_p_vec <- sapply(seq_len(NR), function(i)
          ll_rep_cached(pr_list[[i]], tau_cur[[i]], latent_idxs[[i]], ED_p))
        lp_p <- sum(ll_p_vec) + log_prior_R0(R0_p) + log_prior_k(k_cur)
        
        if (is.finite(lp_p) && log(runif(1L)) < lp_p - lp_cur) {
          R0_cur   <- R0_p; ED_cur <- ED_p
          ll_cache <- ll_p_vec; ll_sum <- sum(ll_cache)
          lp_cur   <- lp_p
          if (iter <= N_BURN) acc_R0 <- acc_R0 + 1L
        }
      }
    }
    
    # ── Step B: propose new tau for one latent node ──────────────────────────
    # Uses the current ED_cur — no ODE solve needed
    if (any_latent) {
      prop_bt  <- prop_bt + 1L
      eligible <- which(sapply(pr_list, function(pr) length(pr$latent_idx) > 0L))
      r        <- eligible[sample.int(length(eligible), 1L)]
      pr_r     <- pr_list[[r]]
      lat_r    <- pr_r$latent_idx
      i_node   <- sample.int(length(lat_r), 1L)
      m        <- lat_r[i_node]
      tau_new  <- runif(1L, pr_r$lo[m], pr_r$hi[m])
      
      tau_prop      <- tau_cur
      tau_prop[[r]][i_node] <- tau_new
      
      ll_r_new <- ll_rep_cached(pr_r, tau_prop[[r]], lat_r, ED_cur)
      lp_prop  <- lp_cur - ll_cache[r] + ll_r_new
      
      if (is.finite(lp_prop) && log(runif(1L)) < lp_prop - lp_cur) {
        tau_cur[[r]] <- tau_prop[[r]]
        ll_cache[r]  <- ll_r_new
        ll_sum       <- sum(ll_cache)
        lp_cur       <- lp_prop
        if (iter <= N_BURN) acc_bt <- acc_bt + 1L
      }
    }
    
    # ── Step C: propose new k (independence sampler from prior) ──────────────
    # Draw k* ~ DiscreteUniform{1,...,K_MAX}
    # Since proposal = prior, the Hastings ratio is:
    #   [L(D|k*,R0,tau) * p(k*) / q(k*|k)] / [L(D|k, R0,tau) * p(k) / q(k|k*)]
    # = L(D|k*,...) / L(D|k,...) * [p(k*)/p(k)] * [q(k|k*)/q(k*|k)]
    # With independence sampler q(k*|k) = p(k*), so q ratio = p(k)/p(k*) 
    # Simplifies to: L(D|k*,...) / L(D|k,...)  [prior cancels exactly]
    prop_k <- prop_k + 1L
    k_star <- sample.int(K_MAX, 1L)   # k* ~ Uniform{1,...,K_MAX}
    
    if (k_star != k_cur) {
      ED_k <- precompute_ED_for_R0(R0_cur, k = k_star, lam = LAM_FIX,
                                   pobs = POBS_FIX,
                                   maxTime = GLOBAL_MAXTIME, dt = ODE_DT)
      if (!is.null(ED_k)) {
        # Need to check: if k_star > max(j_obs), the likelihood may be undefined
        # because j_obs values can't exceed k in the model.
        # Check: max j_obs across all replicates
        max_jobs <- max(sapply(pr_list, function(pr) {
          jo <- pr$edges_base$j_obs
          if (all(is.na(jo))) 0L else max(jo, na.rm = TRUE)
        }))
        
        if (k_star >= max_jobs) {   # k must be >= max(j_obs)
          ll_k_vec <- sapply(seq_len(NR), function(i)
            ll_rep_cached(pr_list[[i]], tau_cur[[i]], latent_idxs[[i]], ED_k))
          ll_k_sum <- sum(ll_k_vec)
          # Independence sampler: acceptance = L(k*) / L(k)  (prior cancels)
          log_alpha <- ll_k_sum - ll_sum
          
          if (is.finite(log_alpha) && log(runif(1L)) < log_alpha) {
            k_cur    <- k_star; ED_cur <- ED_k
            ll_cache <- ll_k_vec; ll_sum <- ll_k_sum
            lp_cur   <- ll_sum + log_prior_R0(R0_cur) + log_prior_k(k_cur)
            if (iter <= N_BURN) acc_k <- acc_k + 1L
          }
        }
        # If k_star < max_jobs: reject (constraint violation)
      }
    }
    
    # ── Adaptive tuning of R0 step size ─────────────────────────────────────
    if (iter <= N_BURN && iter %% ADAPT_INT == 0L) {
      rate     <- acc_R0 / max(prop_R0, 1L)
      sigma_rw <- max(0.015, min(sigma_rw * exp(rate - TARGET_ACC), 0.8))
    }
    
    # ── Store post-burn-in samples ───────────────────────────────────────────
    if (iter > N_BURN && (iter - N_BURN) %% THIN == 0L) {
      post_idx         <- post_idx + 1L
      R0_chain[post_idx] <- R0_cur
      k_chain[post_idx]  <- k_cur
    }
    
    # ── Progress ─────────────────────────────────────────────────────────────
    if (iter %% max(1L, N_ITER %/% 5L) == 0L) {
      el <- (proc.time() - t0)[3L]
      cat(sprintf(
        "  [p=%.2f] iter%5d | R0=%.3f k=%2d | acc_R0=%2.0f%% acc_k=%2.0f%% acc_bt=%2.0f%% | %.0fs\n",
        p_res, iter, R0_cur, k_cur,
        100 * acc_R0 / max(prop_R0, 1L),
        100 * acc_k  / max(prop_k,  1L),
        100 * acc_bt / max(prop_bt, 1L), el))
    }
  }
  
  list(
    R0_chain = R0_chain[seq_len(post_idx)],
    k_chain  = k_chain[seq_len(post_idx)],
    p_res    = p_res,
    n_post   = post_idx
  )
}

# =============================================================================
# 7. Sweep over p_res values
# =============================================================================
hpdi <- function(x, p = 0.95) {
  x <- sort(x); n <- length(x); g <- max(1L, floor(p * n))
  w <- x[(g + 1L):n] - x[1L:(n - g)]; lo <- which.min(w)
  c(lo = x[lo], hi = x[lo + g])
}

results <- vector("list", length(P_RES_VALS))

for (pi_idx in seq_along(P_RES_VALS)) {
  p_res   <- P_RES_VALS[pi_idx]
  pr_list <- lapply(rep_data, make_partial, p_res = p_res)
  
  n_int_avg   <- mean(sapply(rep_data, `[[`, "n_int"))
  n_known_avg <- mean(sapply(pr_list, function(pr) length(pr$known_idx)))
  
  cat(sprintf("\n=== p_res = %.2f  |  %.0f / %.0f internal nodes fixed ===\n",
              p_res, n_known_avg, n_int_avg))
  
  set.seed(SEED + pi_idx)
  mcmc_out <- run_mcmc(pr_list, p_res)
  
  if (is.null(mcmc_out) || length(mcmc_out$R0_chain) < 10L) {
    results[[pi_idx]] <- list(
      p_res = p_res, R0_mean = NA_real_, R0_sd = NA_real_,
      R0_hpdi_lo = NA_real_, R0_hpdi_hi = NA_real_, R0_bias = NA_real_,
      k_mean = NA_real_, k_sd = NA_real_, k_mode = NA_integer_,
      R0_chain = NULL, k_chain = NULL)
    next
  }
  
  R0_ch <- mcmc_out$R0_chain
  k_ch  <- mcmc_out$k_chain
  ci_R0 <- hpdi(R0_ch)
  k_tab <- table(k_ch)
  k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
  
  res <- list(
    p_res      = p_res,
    R0_mean    = mean(R0_ch),    R0_sd      = sd(R0_ch),
    R0_hpdi_lo = ci_R0["lo"],   R0_hpdi_hi = ci_R0["hi"],
    R0_bias    = mean(R0_ch) - TRUE_R0,
    k_mean     = mean(k_ch),    k_sd       = sd(k_ch),
    k_mode     = k_mode,
    R0_chain   = R0_ch,         k_chain    = k_ch
  )
  results[[pi_idx]] <- res
  
  cat(sprintf("  R0: mean=%.3f  SD=%.4f  HPDI=[%.3f, %.3f]  bias=%+.3f\n",
              res$R0_mean, res$R0_sd, res$R0_hpdi_lo, res$R0_hpdi_hi, res$R0_bias))
  cat(sprintf("   k: mean=%.2f  SD=%.2f  mode=%d  (true k=%d)\n",
              res$k_mean, res$k_sd, res$k_mode, TRUE_K))
  
  # ── Save individual chain to .rds ─────────────────────────────────────────
  if (SAVE_CHAINS) {
    fname <- sprintf("partial_chains_p%.2f.rds", p_res)
    saveRDS(list(R0_chain = R0_ch, k_chain = k_ch,
                 p_res = p_res, N_ITER = N_ITER, N_BURN = N_BURN,
                 THIN = THIN, K_MAX = K_MAX,
                 LAM_FIX = LAM_FIX, POBS_FIX = POBS_FIX,
                 TRUE_R0 = TRUE_R0, TRUE_K = TRUE_K),
            fname)
    cat(sprintf("  Saved: %s\n", fname))
  }
}

# =============================================================================
# 8. Summary table + CSV
# =============================================================================
cat("\n"); cat(strrep("=", 72), "\n")
cat(sprintf("%-8s  %-7s %-7s %-18s %-7s | %-6s %-5s %-5s\n",
            "p_res", "R0_mean", "R0_SD", "95% HPDI R0", "R0_bias",
            "k_mean", "k_SD", "k_mode"))
cat(strrep("-", 72), "\n")
for (res in results) {
  if (is.na(res$R0_mean)) {
    cat(sprintf("p=%.2f  FAILED\n", res$p_res))
  } else {
    cat(sprintf("p=%.2f  %-7.3f %-7.4f [%.3f, %.3f]  %+.3f | %-6.2f %-5.2f %d\n",
                res$p_res, res$R0_mean, res$R0_sd,
                res$R0_hpdi_lo, res$R0_hpdi_hi, res$R0_bias,
                res$k_mean, res$k_sd, res$k_mode))
  }
}
cat(strrep("-", 72), "\n")
cat(sprintf("%-8s  %-7.3f %-7s [%.3f, %.3f]  %+.3f | (MLE reference, k=4)\n",
            "Full-MLE", FULL_TREE_MLE, "--",
            FULL_TREE_CI[1], FULL_TREE_CI[2], FULL_TREE_MLE - TRUE_R0))
cat(strrep("=", 72), "\n")
cat(sprintf("True R0=%.1f, True k=%d | lam=%.1f, pobs=%.1f\n\n",
            TRUE_R0, TRUE_K, LAM_FIX, POBS_FIX))

# Save summary CSV
smry <- do.call(rbind, lapply(results, function(r) {
  data.frame(p_res=r$p_res, R0_mean=r$R0_mean, R0_sd=r$R0_sd,
             R0_hpdi_lo=r$R0_hpdi_lo, R0_hpdi_hi=r$R0_hpdi_hi,
             R0_bias=r$R0_bias, R0_hpdi_width=r$R0_hpdi_hi - r$R0_hpdi_lo,
             k_mean=r$k_mean, k_sd=r$k_sd, k_mode=r$k_mode)
}))
write.csv(smry, "partial_tree_summary.csv", row.names = FALSE)
cat("Saved: partial_tree_summary.csv\n")

# Save all results as RDS for later replotting
saveRDS(list(results=results, P_RES_VALS=P_RES_VALS,
             TRUE_R0=TRUE_R0, TRUE_K=TRUE_K, K_MAX=K_MAX,
             FULL_TREE_MLE=FULL_TREE_MLE, FULL_TREE_CI=FULL_TREE_CI,
             N_ITER=N_ITER, N_BURN=N_BURN, THIN=THIN,
             LAM_FIX=LAM_FIX, POBS_FIX=POBS_FIX),
        "partial_tree_all_results.rds")
cat("Saved: partial_tree_all_results.rds\n\n")

# =============================================================================
# 9. Plots
# =============================================================================
# Helper to reload results later:
saved <- readRDS("partial_tree_all_results.rds")
results <- saved$results; P_RES_VALS <- saved$P_RES_VALS
smry <- read.csv('partial_tree_summary.csv', header = TRUE)

pv  <- smry$p_res
mn  <- smry$R0_mean
lo  <- smry$R0_hpdi_lo
hi  <- smry$R0_hpdi_hi
bv  <- smry$R0_bias
wv  <- smry$R0_hpdi_width
km  <- smry$k_mean
ks  <- smry$k_sd

# ── Panel 1: R0 posterior mean + HPDI vs p_res ───────────────────────────────
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
ylim1 <- range(c(lo, hi, TRUE_R0, FULL_TREE_MLE), na.rm = TRUE) + c(-0.3, 0.5)

plot(pv, mn, type = "n", 
     xlab = "Resolution fraction (p_res)",
     ylab = expression(hat(R)[0]),
     main = "",
     ylim = ylim1, xaxt = "n")
axis(1, at = P_RES_VALS, labels = sprintf("%.2f", P_RES_VALS))

ok <- !is.na(lo)
polygon(c(pv[ok], rev(pv[ok])), c(lo[ok], rev(hi[ok])),
        col = 'gray87', border = NA)

points(pv, mn, type = "b", pch = 19, cex = 1.5, lwd = 2, col = C_DARK)
lines(pv[ok], lo[ok], lty = 2, col = C_MID, lwd = 2)
lines(pv[ok], hi[ok], lty = 2, col = C_MID, lwd = 2)
abline(h = TRUE_R0,       col = C_BLACK, lwd = 2,   lty = 2)
abline(h = FULL_TREE_MLE, col = C_DARK,  lwd = 2.0, lty = 3)

legend("topright", bty = "n", cex = 1.0,
       legend = c(expression(True~R[0]),
                  sprintf("Full tree MLE"),
                  expression(Posterior~mean~R[0]),
                  "95% HPDI"),
       col = c(C_BLACK, C_DARK, C_DARK, C_MID),
       lty = c(2,3,1,2), lwd = c(2,2,2,2), pch = c(NA,NA,19,NA))

# ── Panel 2: k posterior mean ± 1 SD vs p_res ────────────────────────────────
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
ylim_k <- c(0, K_MAX + 1)

plot(pv, km, type = "n",
     xlab = "Resolution fraction (p_res)",
     ylab = "k (contact degree)",
     main = "",
     ylim = ylim_k, xaxt = "n")
axis(1, at = P_RES_VALS, labels = sprintf("%.2f", P_RES_VALS))

ok_k <- !is.na(km)
polygon(c(pv[ok_k], rev(pv[ok_k])),
        c(km[ok_k] + ks[ok_k], rev(km[ok_k] - ks[ok_k])),
        col = 'gray87', border = NA)

points(pv, km, type = "b", pch = 19, cex = 1.5, lwd = 2, col = C_DARK)
lines(pv[ok_k], km[ok_k] + ks[ok_k], lty = 2, col = C_MID, lwd = 2)
lines(pv[ok_k], km[ok_k] - ks[ok_k], lty = 2, col = C_MID, lwd = 2)
abline(h = TRUE_K,  col = C_BLACK, lwd = 2,   lty = 2)
abline(h = K_MAX/2, col = C_MID,   lwd = 2,   lty = 3)

legend("topright", bty = "n", cex = 1.0,
       legend = c(sprintf("True k"),
                  "Post. mean k", "Mean ± 1 SD"),
       col = c(C_BLACK, C_DARK, C_MID),
       lty = c(2, 1, 2), lwd = c(2, 2, 2), pch = c(NA, 19, NA))

# ── Panel 3: k posterior histograms at selected p_res ────────────────────────
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
sel_idx <- c(1L, ceiling(length(P_RES_VALS)/2L), length(P_RES_VALS))
cols3   <- c(C_MID, C_BLACK, C_BLACK)
ltys3   <- c(1, 2, 3)
plot(NULL, xlim = c(1, K_MAX), ylim = c(0, 1),
     xlab = "k", ylab = "Posterior probability",
     main = "")

for (ii in seq_along(sel_idx)) {
  ri <- results[[sel_idx[ii]]]
  if (is.null(ri$k_chain)) next
  k_tab <- table(factor(ri$k_chain, levels = 1:K_MAX))
  k_vals <- as.integer(names(k_tab))
  k_prob <- as.numeric(k_tab) / sum(k_tab)
  k_mode <- as.integer(names(k_tab)[which.max(k_tab)])
  lines(1:K_MAX, k_prob, col = cols3[ii], lwd = 2, 
        lty = ltys3[ii], type = "b")
  points(k_mode, k_prob[k_vals == k_mode], pch = 16, col = C_BLACK, cex = 1.5)
}
abline(v = TRUE_K, col = C_MID, lwd = 2, lty = 2)
legend("topright", bty = "n", cex = 1.0,
       legend = c(sprintf("p_res=%.2f", P_RES_VALS[sel_idx]),
                  "True k",   # fixed
                  "k mode"),
       col = c(cols3, C_MID, C_BLACK),   # now length 5
       lty = c(ltys3, 2, NA),            # line type for first 4, none for Mode
       lwd = c(2, 2, 2, 2, NA),          # line width for first 4
       pch = c(NA, NA, NA, NA, 16))      # point symbol only for Mode

# ── Panel 4: Joint R0 posterior densities at selected p_res ──────────────────
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
dens_list <- lapply(sel_idx, function(si) {
  ri <- results[[si]]
  if (is.null(ri$R0_chain)) return(NULL)
  density(ri$R0_chain, adjust = 1.2)
})
valid_d <- !sapply(dens_list, is.null)
if (any(valid_d)) {
  xlim_d  <- range(sapply(dens_list[valid_d], function(d) range(d$x)))
  ylim_d  <- c(0, max(sapply(dens_list[valid_d], function(d) max(d$y))) * 1.2)
  plot(NULL, xlim = xlim_d, ylim = ylim_d,
       xlab = expression(R[0]), ylab = "Density", main = "")
  for (ii in seq_along(sel_idx)) {
    if (!is.null(dens_list[[ii]]))
      lines(dens_list[[ii]], col = cols3[ii], lwd = 2, lty = ltys3[ii])
  }
  abline(v = TRUE_R0, col = C_MID, lwd = 2.0, lty = 2)
  legend("topright", bty = "n", cex = 1.0,
         legend = c(sprintf("p_res=%.2f", P_RES_VALS[sel_idx]),
                    bquote(True~R[0])),
         col = c(cols3, C_MID), lty = c(ltys3, 2), lwd = c(2,2,2,2.0))
}

par(mfrow = c(1, 1))
cat("Plots drawn — save via RStudio Export\n")
cat("=== DONE ===\n")
cat("\nTo reload results later without re-running:\n")
cat("  saved   <- readRDS('partial_tree_all_results.rds')\n")
cat("  results <- saved$results\n")
cat("  # individual chains:\n")
cat("  ch0 <- readRDS('partial_chains_p0.00.rds')\n")
cat("  plot(density(ch0$R0_chain))\n")
