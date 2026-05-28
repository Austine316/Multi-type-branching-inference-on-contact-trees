# =============================================================================
# parameter_estimation.R
# =============================================================================
#
# ESTIMAND AND PARAMETRISATION
# ─────────────────────────────────────────────────────────────────────────────
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


cat('Current working directory', getwd())
setwd(getwd())

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
#  The five product terms of the full likelihood are computed here:
#
#  log L = (E) root conditioning
#        + sum over OBSERVED last segments    (A): log D^s_{j_obs}(delta)
#        + sum over UNOBSERVED last segments  (B): log E_s(delta)
#        + sum over NON-LAST segments         (C): -(lam+(k-s)*beta)*delta
#        + sum over BRANCHING EVENTS          (D): log((k-s)*beta)
# =============================================================================

log_lik_rep <- function(rep_edges, ED, pi_vec, beta, lam, k) {
  
  TINY    <- 1e-300
  log_lik <- 0.0
  
  # ── (E) Root conditioning: log pi_{i0} - log(1 - E_{i0}(Delta)) ──────────
  # Conditions on observing at least one tip from the root.
  # pi_{i0}      : equilibrium prob of root starting in state i0
  # E_{i0}(Delta): prob the entire clade goes undetected over window Delta
  # Denominator 1 - E_{i0}(Delta) = prob at least one tip is observed.
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
  has_branch <- !is.na(rep_edges$branch_rate)   # MASK only — not for the rate value
  
  s_vec     <- as.integer(rep_edges$s)
  c_vec     <- as.integer(rep_edges$c)           # = j_obs for tip rows
  delta_vec <- rep_edges$delta
  
  # ── (A) Observed tip edges: D^{s}_{j_obs}(delta) ─────────────────────────
  # Last segment of an OBSERVED individual.
  # D^s_{j}(delta) = density of being sampled in state j after time delta,
  # given the individual started segment in state s.
  # j_obs is stored in column c (end-state) for tip rows.
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
  # Last segment of an UNOBSERVED individual.
  # E_s(delta) = probability the entire sub-clade below this individual
  # produces no observed tips over the remaining time delta.
  idx_B <- which(is_last & !is_obs)
  for (ii in idx_B) {
    s     <- s_vec[ii]; delta <- delta_vec[ii]
    if (s < 0L || s > k || delta < 0.0) {
      log_lik <- log_lik + log(TINY); next
    }
    v       <- ED[[1L]]$approxE[[s + 1L]](delta)
    log_lik <- log_lik + log(max(v, TINY))
  }
  
  # ── (C) Non-last (inter-branch) edges: exp(-(lam + (k-s)*beta) * delta) ──
  # Survival probability over [tau_a, tau_b]: individual in state s survives
  # without removal (rate mu+sigma = lam) and without infecting (rate (k-s)*beta)
  # for the entire inter-branch duration delta.
  # Total hazard = lam + (k-s)*beta.
  idx_C <- which(!is_last)
  if (length(idx_C) > 0L) {
    rates   <- lam + (k - s_vec[idx_C]) * beta
    log_lik <- log_lik + sum(-rates * delta_vec[idx_C])
  }
  
  # ── (D) Branching factors: (k-s)*beta ─────────────────────────────────────
  # Instantaneous rate of infecting the next contact when in state s.
  # One factor per observed branching event.
  # IMPORTANT: beta is recomputed from the current R0 proposal — the
  # branch_rate column stores the simulation's true beta and is used
  # only as a mask to identify which rows correspond to branching events.
  idx_D <- which(has_branch)
  if (length(idx_D) > 0L) {
    br <- (k - s_vec[idx_D]) * beta
    ok <- br > 0.0 & is.finite(br)
    if (any(!ok)) log_lik <- log_lik + sum(!ok) * log(TINY)
    if (any(ok))  log_lik <- log_lik + sum(log(br[ok]))
  }
  
  log_lik
}


# =============================================================================
# 5.  Negative log-likelihood over all replicates
#
#  ESTIMAND: R0 = k * beta / lam  (reproduction number at state i=0)
#  BETA FORMULA: beta = R0 * lam / k
#
# =============================================================================

neg_log_lik <- function(R0, dat, k, lam = LAM, pobs = POBS,
                        maxTime, dt = 0.005, verbose = FALSE) {
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
    ll <- tryCatch(
      log_lik_rep(re, ED, pi_vec, beta, lam, k),
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
    
    # ── Derived quantities from R0_hat ───────
    beta_h  <- R0_hat * lam / kk   # FIXED
    sigma_h <- pobs * lam
    mu_h    <- (1.0 - pobs) * lam
    nr      <- length(unique(dat$rep_id))
    np      <- 2L   # R0 and k counted as two parameters
    
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
  
  raw <- read.csv(f, stringsAsFactors = FALSE)
  dat <- raw[is.finite(raw$delta) & raw$delta >= 0 &
               !is.nan(raw$s) & !is.nan(raw$c), ]
  dat$j_obs[is.nan(dat$j_obs)] <- NA_integer_
  rownames(dat) <- NULL
  
  cat(sprintf("Loaded %d rows from %d replicates.\n",
              nrow(dat), length(unique(dat$rep_id))))
  
  n_total <- nrow(dat)
  n_tips <- sum(!is.na(dat$j_obs))
  
  cat(sprintf("Total edge rows:  %d\n", n_total))
  cat(sprintf("Observed tips:    %d\n", n_tips))
  cat(sprintf("Empirical p_obs:  %d / %d = %.4f\n",
              n_tips, n_total, n_tips / n_total))
  
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

