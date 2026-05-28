# =============================================================================
# parameter_estimation_mle.R
# =============================================================================
#
# ESTIMAND AND PARAMETRISATION
# ─────────────────────────────────────────────────────────────────────────────
# This model has two natural reproduction numbers:
#
#   R0   = k * beta / lam
#          Expected infections from an individual starting in state i = 0
#          (no contacts infected yet).  THIS IS THE ESTIMAND.
#          Inversion: beta = R0 * lam / k
#
#   Rbar = (beta / lam) * S_pi(k),   S_pi(k) = sum_i pi_i * (k - i)
#          Equilibrium-weighted average reproduction number.
#          For our parameter range S_pi(k) = k - 1 exactly, so
#          Rbar = R0 * (k-1) / k  < R0  always.
# 
#
# The code estimates R0 and derives beta = R0 * lam / k.
# Rbar can be reported as a secondary quantity but must not replace R0.
#
#
# FIXED PARAMETERS:
#   lam  = mu + sigma = 1.0
#   pobs = sigma / lam = 0.5
#   Free: R0 (continuous), k (integer grid)
#
# TRUE VALUES: k=4, beta=1.5, mu=0.5, sigma=0.5  =>  R0 = 6
# =============================================================================
#
# LIKELIHOOD FORMULA — ONE ROW AT A TIME
# =============================================================================
# is_last = (seg_idx == n_segs - 1)
# is_obs  = !is.na(j_obs)
#
# (A) is_last & is_obs     ->  D^{s}_{j_obs}(delta)
#       Tip density: probability density of being sampled in state j_obs
#       after an elapsed time delta from the previous event.
#
# (B) is_last & !is_obs    ->  E_{s}(delta)
#       Extinction probability: probability the entire unobserved sub-clade
#       goes undetected over the remaining time delta.
#
# (C) !is_last             ->  exp( -(lam + (k-s)*beta) * delta )
#       Survival over [tau_a, tau_b]: no removal and no new infection
#       in the inter-branch interval of length delta.
#       NOTE: D^{s}_{s+1}(delta) must NOT be used here — that would
#       double-count the newborn lineage which has its own rows.
#
# (D) !is.na(branch_rate)  ->  (k - s) * beta   [recomputed, not from table]
#       Instantaneous infection rate at a branching event.
#       branch_rate column is used only as a MASK (it stores the true
#       simulation beta, not the current proposal).
#
# (E) Root conditioning:
#       log(pi_{i0}) - log(1 - E_{i0}(Delta))
#       Conditions on the tagged root having at least one observed tip.
#       pi_{i0} = equilibrium frequency of starting state i0.
#       Delta = length of the full observation window.
# =============================================================================

suppressPackageStartupMessages(library(deSolve))

LAM  <- 1.0   # fixed: mu + sigma
POBS <- 0.5   # fixed: sigma / (mu + sigma)


# =============================================================================
# 1.  ODE: E_i(t) and D^i_{j_obs}(t)
# =============================================================================

solve_ED <- function(j_obs, beta, mu, sigma, k, maxTime, dt = 0.005) {
  stopifnot(j_obs >= 0L, j_obs <= k)
  times <- seq(0, maxTime, by = dt)
  n     <- k + 1L

  y0                 <- c(rep(1.0, n), rep(0.0, n))
  y0[n + j_obs + 1L] <- sigma          # D^{j_obs}_{j_obs}(0) = sigma

  rhs <- function(t, y, parms) {
    E  <- y[seq_len(n)]
    D  <- y[(n + 1L):(2L * n)]
    dE <- numeric(n)
    dD <- numeric(n)

    dE[n] <- mu - (mu + sigma) * E[n]
    dD[n] <- -(mu + sigma) * D[n]

    for (i in 0L:(k - 1L)) {
      r          <- (k - i) * beta
      dE[i + 1L] <- mu - (mu + sigma + r) * E[i + 1L] +
        r * E[1L] * E[i + 2L]
      dD[i + 1L] <- -(mu + sigma + r) * D[i + 1L] +
        r * E[1L]     * D[i + 2L] +
        r * E[i + 2L] * D[1L]
    }
    list(c(dE, dD))
  }

  sol <- ode(y0, times, rhs, parms = NULL,
             method = "lsoda", rtol = 1e-9, atol = 1e-11)
  list(
    E_mat = sol[, 2L:(n + 1L),          drop = FALSE],
    D_mat = sol[, (n + 2L):(2L*n + 1L), drop = FALSE],
    times = times
  )
}


# =============================================================================
# 2.  Precompute interpolators for all j_obs = 0..k
# =============================================================================

precompute_ED <- function(beta, mu, sigma, k, maxTime, dt = 0.005) {
  lapply(0L:k, function(jj) {
    sol <- solve_ED(jj, beta, mu, sigma, k, maxTime, dt)
    list(
      approxE = lapply(0L:k, function(ii)
        approxfun(sol$times, sol$E_mat[, ii + 1L], rule = 2L)),
      approxD = lapply(0L:k, function(ii)
        approxfun(sol$times, sol$D_mat[, ii + 1L], rule = 2L))
    )
  })
}


# =============================================================================
# 3.  Equilibrium frequencies pi_i
# =============================================================================

compute_pi <- function(beta, mu, sigma, k) {
  gam     <- mu + sigma
  char_eq <- function(r) {
    lhs <- r + gam + k * beta
    rhs <- k * beta; pv <- 1.0
    for (j in seq_len(k - 1L)) {
      pv  <- pv * (k - j + 1L) * beta / (gam + (k - j) * beta + r)
      rhs <- rhs + (k - j) * beta * pv
    }
    lhs - rhs
  }
  r <- tryCatch(uniroot(char_eq, c(1e-8, 1e4), tol = 1e-12)$root,
                error = function(e) NA_real_)
  if (is.na(r)) return(rep(1.0/(k + 1L), k + 1L))
  pi    <- numeric(k + 1L); pi[1L] <- 1.0
  for (i in seq_len(k))
    pi[i + 1L] <- pi[i] * (k - i + 1L) * beta /
    (gam + (k - i) * beta + r)
  pi / sum(pi)
}


# =============================================================================
# 4.  Log-likelihood for one replicate
#
#  The five product terms of the full likelihood are computed here.
#
#  latent_rows : integer vector of row indices (within rep_edges) whose
#                branching-event state is UNKNOWN.  For these rows the
#                fixed-state terms (C) and (D) are replaced by the
#                pi_s-weighted marginalised contribution
#                  log( sum_{s=0}^{k} pi_s*(k-s)*beta*
#                                     exp(-(lam+(k-s)*beta)*delta) )
#                NULL (default) = all states known; standard full-tree.
# =============================================================================

log_lik_rep <- function(rep_edges, ED, pi_vec, beta, lam, k,
                        latent_rows = NULL) {

  TINY    <- 1e-300
  log_lik <- 0.0

  # ── (E) Root conditioning ─────────────────────────────────────────────────
  i0     <- as.integer(rep_edges$root_state[1L])
  if (i0 < 0L || i0 > k) return(log(TINY))

  burn_t     <- min(rep_edges$tau_a)
  Delta      <- max(rep_edges$tau_b) - burn_t
  E_i0_Delta <- ED[[1L]]$approxE[[i0 + 1L]](Delta)
  denom      <- 1.0 - E_i0_Delta
  if (!is.finite(denom) || denom <= 0.0) return(log(TINY))
  log_lik <- log_lik + log(pi_vec[i0 + 1L]) - log(denom)

  # ── Pre-classify all rows ──────────────────────────────────────────────────
  is_last    <- rep_edges$seg_idx == rep_edges$n_segs - 1L
  is_obs     <- !is.na(rep_edges$j_obs)
  has_branch <- !is.na(rep_edges$branch_rate)

  s_vec     <- as.integer(rep_edges$s)
  c_vec     <- as.integer(rep_edges$c)
  delta_vec <- rep_edges$delta

  # Mark latent rows (parent rows of branching nodes with unknown state)
  is_latent <- logical(nrow(rep_edges))
  if (!is.null(latent_rows) && length(latent_rows) > 0L)
    is_latent[latent_rows] <- TRUE

  # ── (A) Observed tip edges: D^{s}_{j_obs}(delta) ─────────────────────────
  idx_A <- which(is_last & is_obs)
  for (ii in idx_A) {
    s     <- s_vec[ii]; cc <- c_vec[ii]; delta <- delta_vec[ii]
    if (s < 0L || s > k || cc < 0L || cc > k || delta < 0.0) {
      log_lik <- log_lik + log(TINY); next
    }
    v       <- ED[[cc + 1L]]$approxD[[s + 1L]](delta)
    log_lik <- log_lik + log(max(v, TINY))
  }

  # ── (B) Unobserved last edges: E_{s}(delta) ──────────────────────────────
  idx_B <- which(is_last & !is_obs)
  for (ii in idx_B) {
    s     <- s_vec[ii]; delta <- delta_vec[ii]
    if (s < 0L || s > k || delta < 0.0) {
      log_lik <- log_lik + log(TINY); next
    }
    v       <- ED[[1L]]$approxE[[s + 1L]](delta)
    log_lik <- log_lik + log(max(v, TINY))
  }

  # ── (C) Non-last edges: survival (known-state rows only) ─────────────────
  idx_C_known <- which(!is_last & !is_latent)
  if (length(idx_C_known) > 0L) {
    rates   <- lam + (k - s_vec[idx_C_known]) * beta
    log_lik <- log_lik + sum(-rates * delta_vec[idx_C_known])
  }

  # ── (D) Branching factors (known-state rows only) ─────────────────────────
  idx_D_known <- which(has_branch & !is_latent)
  if (length(idx_D_known) > 0L) {
    br <- (k - s_vec[idx_D_known]) * beta
    ok <- br > 0.0 & is.finite(br)
    if (any(!ok)) log_lik <- log_lik + sum(!ok) * log(TINY)
    if (any(ok))  log_lik <- log_lik + sum(log(br[ok]))
  }

  # ── (C+D) Latent branching nodes: marginalise over s in {0,...,k} ─────────
  # For each latent parent row replace the fixed-state C+D contribution with:
  #   log( sum_{s=0}^{k} pi_s * (k-s)*beta * exp(-(lam+(k-s)*beta)*delta) )
  idx_D_latent <- which(has_branch & is_latent)
  if (length(idx_D_latent) > 0L) {
    for (ii in idx_D_latent) {
      delta <- delta_vec[ii]
      if (!is.finite(delta) || delta < 0) next
      marg <- sum(sapply(0L:k, function(s) {
        rate <- (k - s) * beta
        if (rate <= 0) return(0.0)
        pi_vec[s + 1L] * rate * exp(-(lam + rate) * delta)
      }))
      log_lik <- log_lik + log(max(marg, TINY))
    }
  }

  # Latent survival-only rows (non-last & latent but NOT a branching event)
  idx_C_latent_only <- which(!is_last & is_latent & !has_branch)
  if (length(idx_C_latent_only) > 0L) {
    for (ii in idx_C_latent_only) {
      delta <- delta_vec[ii]
      if (!is.finite(delta) || delta < 0) next
      marg <- sum(sapply(0L:k, function(s)
        pi_vec[s + 1L] * exp(-(lam + (k - s) * beta) * delta)))
      log_lik <- log_lik + log(max(marg, TINY))
    }
  }

  log_lik
}


# =============================================================================
# 5.  Negative log-likelihood over all replicates
#
#  latent_nodes : NULL (default) -> all states known; standard full-tree.
#                 named list     -> each element named by rep_id (character),
#                                   containing integer vector of parent-row
#                                   indices with unknown branching state.
#
#  Example for a partially resolved contact-tracing tree:
#    neg_log_lik(R0, dat, k=4, lam=1, pobs=0.5, maxTime=10,
#                latent_nodes = list("rep1" = c(12L, 47L), "rep2" = c(5L)))
# =============================================================================

neg_log_lik <- function(R0, dat, k, lam = LAM, pobs = POBS,
                        maxTime, dt = 0.005, verbose = FALSE,
                        latent_nodes = NULL) {
  if (!is.finite(R0) || R0 <= 0.0 || R0 > 40.0) return(1e15)

  beta  <- R0 * lam / k
  sigma <- pobs * lam
  mu    <- (1.0 - pobs) * lam

  pi_vec <- tryCatch(compute_pi(beta, mu, sigma, k),
                     error = function(e) rep(1.0/(k+1L), k+1L))
  ED     <- tryCatch(precompute_ED(beta, mu, sigma, k, maxTime, dt),
                     error = function(e) NULL)
  if (is.null(ED)) return(1e15)

  total <- 0.0
  for (rid in unique(dat$rep_id)) {
    re <- dat[dat$rep_id == rid, , drop = FALSE]

    # Retrieve latent row indices for this replicate (NULL if all states known)
    lat_rows <- if (!is.null(latent_nodes)) latent_nodes[[as.character(rid)]] else NULL

    ll <- tryCatch(
      log_lik_rep(re, ED, pi_vec, beta, lam, k, latent_rows = lat_rows),
      error = function(e) log(1e-300)
    )
    total <- total + if (is.finite(ll)) ll else log(1e-300)
  }

  nll <- -total
  if (verbose)
    cat(sprintf("  R0=%.4f  k=%d  beta=%.4f  NLL=%.4f\n", R0, k, beta, nll))
  nll
}



# =============================================================================
# 5b. Precomputed-ED variants for high-performance MCMC
#
#  When the MCMC proposes a new R0, the ODE system (E and D functions)
#  must be re-solved because beta = R0*lam/k changes.  In a naive loop,
#  this happens once per replicate per iteration.  With 538 replicates and
#  k=30, that is 538 ODE solves per iteration — prohibitively slow.
#
#  The solution: solve the ODE ONCE per R0 proposal, cache the result as
#  an (ED, pi_vec, beta) triple, and evaluate log_lik_rep for all replicates
#  using the cached triple.  This reduces ODE calls from N_rep to 1 per
#  iteration, yielding an ~N_rep-fold speedup.
#
#  Two new functions:
#    precompute_ED_for_R0(R0, k, lam, pobs, maxTime, dt)
#      -> returns list(ED, pi_vec, beta) or NULL on failure
#
#    ll_from_ED(re, ED_cache, latent_rows)
#      -> evaluates log_lik_rep for one replicate using a cached ED_cache
# =============================================================================

precompute_ED_for_R0 <- function(R0, k, lam = LAM, pobs = POBS,
                                  maxTime, dt = 0.005) {
  # Returns NULL on any failure; caller should check.
  if (!is.finite(R0) || R0 <= 0.0 || R0 > 40.0) return(NULL)
  beta  <- R0 * lam / k
  sigma <- pobs * lam
  mu    <- (1.0 - pobs) * lam
  pi_vec <- tryCatch(compute_pi(beta, mu, sigma, k),
                     error = function(e) rep(1.0 / (k + 1L), k + 1L))
  ED     <- tryCatch(precompute_ED(beta, mu, sigma, k, maxTime, dt),
                     error = function(e) NULL)
  if (is.null(ED)) return(NULL)
  list(ED = ED, pi_vec = pi_vec, beta = beta, lam = lam, k = k)
}


ll_from_ED <- function(re, ED_cache, latent_rows = NULL) {
  # Evaluate log_lik_rep for replicate re using a pre-cached ED object.
  # ED_cache: output of precompute_ED_for_R0()
  # latent_rows: integer vector of parent-row indices with unknown state (or NULL)
  tryCatch(
    log_lik_rep(re, ED_cache$ED, ED_cache$pi_vec,
                ED_cache$beta, ED_cache$lam, ED_cache$k,
                latent_rows = latent_rows),
    error = function(e) log(1e-300)
  )
}


# =============================================================================
# 6.  Estimation: grid over k, 1-D optimise over R0
# =============================================================================

estimate_full_tree <- function(dat, k_grid = 3L:6L,
                               lam = LAM, pobs = POBS,
                               maxTime, dt = 0.005,
                               R0_lo = 0.5, R0_hi = 30.0,
                               verbose = FALSE) {
  all_res <- list()

  for (kk in k_grid) {
    cat(sprintf("\n=== k = %d ===\n", kk))

    res <- tryCatch(
      optimize(neg_log_lik, interval = c(R0_lo, R0_hi),
               dat = dat, k = kk, lam = lam, pobs = pobs,
               maxTime = maxTime, dt = dt, verbose = verbose,
               tol = 1e-4),
      error = function(e) list(minimum = NA_real_, objective = Inf)
    )

    R0_hat <- res$minimum
    nll    <- res$objective

    cat(sprintf("  R0_MLE = %.4f  NLL = %.4f\n", R0_hat, nll))
    if (!is.finite(nll)) next

    beta_h  <- R0_hat * lam / kk
    sigma_h <- pobs * lam
    mu_h    <- (1.0 - pobs) * lam
    nr      <- length(unique(dat$rep_id))
    np      <- 2L

    all_res[[as.character(kk)]] <- list(
      k     = kk,
      R0    = R0_hat,
      beta  = beta_h,
      mu    = mu_h,
      sigma = sigma_h,
      lam   = lam,
      pobs  = pobs,
      nll   = nll,
      aic   = 2L * np + 2.0 * nll,
      bic   = log(nr) * np + 2.0 * nll
    )
  }

  if (length(all_res) == 0L) stop("All k values failed.")

  tab <- do.call(rbind, lapply(all_res, as.data.frame))
  rownames(tab) <- NULL
  tab <- tab[order(tab$nll), ]

  cat("\n\n===== RESULTS TABLE (sorted by NLL) =====\n")
  print(round(tab, 4L))

  w <- tab[1L, ]
  cat("\n===== BEST MODEL =====\n")
  cat(sprintf("  k     = %d\n",               as.integer(w$k)))
  cat(sprintf("  R0    = %.4f  [true 6.0]\n", w$R0))
  cat(sprintf("  beta  = %.4f  [true 1.5]  (= R0*lam/k)\n", w$beta))
  cat(sprintf("  mu    = %.4f  [true 0.5]\n", w$mu))
  cat(sprintf("  sigma = %.4f  [true 0.5]\n", w$sigma))
  cat(sprintf("  lam   = %.4f  [fixed]\n",    w$lam))
  cat(sprintf("  pobs  = %.4f  [fixed]\n",    w$pobs))
  cat(sprintf("  NLL   = %.4f\n",             w$nll))

  invisible(list(table = tab, best_k = as.integer(w$k),
                 best_R0 = w$R0, results = all_res))
}


# =============================================================================
# 7.  Profile likelihood plot
# =============================================================================

plot_profile <- function(R0_hat, dat, k, lo, hi, n = 50L,
                         lam = LAM, pobs = POBS, maxTime, dt = 0.005) {
  grid  <- seq(lo, hi, length.out = n)
  nll_v <- sapply(grid, function(r0)
    neg_log_lik(r0, dat, k, lam, pobs, maxTime, dt))

  ci_h <- min(nll_v) + qchisq(0.95, df = 1L) / 2.0
  plot(grid, nll_v, type = "l", lwd = 2L,
       ylim = c(min(nll_v), min(nll_v) + qchisq(0.99, df = 1L)),
       xlab = "R0", ylab = "NLL",
       main = sprintf("Profile: R0  (k=%d, lam=%.1f, pobs=%.1f fixed)",
                      k, lam, pobs))
  abline(v = R0_hat, col = "blue",   lty = 2L, lwd = 2L)
  abline(v = 6.0,    col = "red",    lty = 3L, lwd = 2L)
  abline(h = ci_h,   col = "gray60", lty = 3L)
  legend("topright",
         legend = c(sprintf("MLE = %.3f", R0_hat), "True = 6.000", "95% CI"),
         col = c("blue","red","gray60"), lty = c(2L,3L,3L),
         lwd = c(2L,2L,1L), bty = "n", cex = 0.85)
  invisible(data.frame(R0 = grid, nll = nll_v))
}


# =============================================================================
# 8.  Entry point
# =============================================================================

if (sys.nframe() == 0L) {

  f <- "full_tree_edges.csv"
  if (!file.exists(f)) stop("Run full_tree_sim.py first to produce ", f)

  raw <- read.csv(f, stringsAsFactors = FALSE)
  dat <- raw[is.finite(raw$delta) & raw$delta >= 0 &
               !is.nan(raw$s) & !is.nan(raw$c), ]
  dat$j_obs[is.nan(dat$j_obs)] <- NA_integer_
  rownames(dat) <- NULL

  cat(sprintf("Loaded %d rows from %d replicates.\n",
              nrow(dat), length(unique(dat$rep_id))))
  maxT <- max(dat$delta[is.finite(dat$delta)], na.rm = TRUE) + 0.5

  res <- estimate_full_tree(dat     = dat,
                            k_grid  = 1L:12L,
                            lam     = LAM,
                            pobs    = POBS,
                            maxTime = maxT,
                            dt      = 0.002)

  bk     <- res$best_k
  R0_hat <- res$best_R0
  plot_profile(R0_hat, dat, bk,
               lo = max(0.5, R0_hat - 5), hi = R0_hat + 5,
               maxTime = maxT)
}
