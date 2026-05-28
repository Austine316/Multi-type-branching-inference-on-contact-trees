# =============================================================================
# paramater_estimate_mle_sensitivity_pobs.R
# =============================================================================
# 
# This script runs a sensitivity analysis: it estimates R0 for a range of
# fixed pobs values (0, 0.1, 0.2, ..., 1.0) using the same simulated dataset
# (true pobs = 0.5, true R0 = 6, true k = 4).  The results are plotted as
# R0_hat vs pobs.
# 
# =============================================================================

cat('Current working directory', getwd())
setwd(getwd())

suppressPackageStartupMessages(library(deSolve))

LAM  <- 1.0   # fixed: mu + sigma

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
# =============================================================================

log_lik_rep <- function(rep_edges, ED, pi_vec, beta, lam, k) {
  
  TINY    <- 1e-300
  log_lik <- 0.0
  
  # Root conditioning
  i0     <- as.integer(rep_edges$root_state[1L])
  if (i0 < 0L || i0 > k) return(log(TINY))
  
  burn_t     <- min(rep_edges$tau_a)
  Delta      <- max(rep_edges$tau_b) - burn_t
  E_i0_Delta <- ED[[1L]]$approxE[[i0 + 1L]](Delta)
  denom      <- 1.0 - E_i0_Delta
  if (!is.finite(denom) || denom <= 0.0) return(log(TINY))
  
  log_lik <- log_lik + log(pi_vec[i0 + 1L]) - log(denom)
  
  # Pre‑classify rows
  is_last    <- rep_edges$seg_idx == rep_edges$n_segs - 1L
  is_obs     <- !is.na(rep_edges$j_obs)
  has_branch <- !is.na(rep_edges$branch_rate)   # MASK only
  
  s_vec     <- as.integer(rep_edges$s)
  c_vec     <- as.integer(rep_edges$c)           # = j_obs for tip rows
  delta_vec <- rep_edges$delta
  
  # (A) Observed tip edges: D^{s}_{j_obs}(delta)
  idx_A <- which(is_last & is_obs)
  for (ii in idx_A) {
    s     <- s_vec[ii]; cc <- c_vec[ii]; delta <- delta_vec[ii]
    if (s < 0L || s > k || cc < 0L || cc > k || delta < 0.0) {
      log_lik <- log_lik + log(TINY); next
    }
    v       <- ED[[cc + 1L]]$approxD[[s + 1L]](delta)
    log_lik <- log_lik + log(max(v, TINY))
  }
  
  # (B) Unobserved last edges: E_{s}(delta)
  idx_B <- which(is_last & !is_obs)
  for (ii in idx_B) {
    s     <- s_vec[ii]; delta <- delta_vec[ii]
    if (s < 0L || s > k || delta < 0.0) {
      log_lik <- log_lik + log(TINY); next
    }
    v       <- ED[[1L]]$approxE[[s + 1L]](delta)
    log_lik <- log_lik + log(max(v, TINY))
  }
  
  # (C) Non-last edges: survival probability
  idx_C <- which(!is_last)
  if (length(idx_C) > 0L) {
    rates   <- lam + (k - s_vec[idx_C]) * beta
    log_lik <- log_lik + sum(-rates * delta_vec[idx_C])
  }
  
  # (D) Branching factors: (k-s)*beta
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

estimate_full_tree <- function(dat, k_grid = 4L:5L,
                               lam = LAM, pobs = POBS,
                               maxTime, dt = 0.005,
                               R0_lo = 0.5, R0_hi = 30.0,
                               verbose = FALSE) {
  all_res <- list()
  
  for (kk in k_grid) {
    if (verbose) cat(sprintf("\n=== k = %d ===\n", kk))
    
    res <- tryCatch(
      optimize(neg_log_lik, interval = c(R0_lo, R0_hi),
               dat = dat, k = kk, lam = lam, pobs = pobs,
               maxTime = maxTime, dt = dt, verbose = verbose,
               tol = 1e-4),
      error = function(e) list(minimum = NA_real_, objective = Inf)
    )
    
    R0_hat <- res$minimum
    nll    <- res$objective
    
    if (verbose) cat(sprintf("  R0_MLE = %.4f  NLL = %.4f\n", R0_hat, nll))
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
  
  if (length(all_res) == 0L) return(NULL)
  
  tab <- do.call(rbind, lapply(all_res, as.data.frame))
  rownames(tab) <- NULL
  tab <- tab[order(tab$nll), ]
  
  if (verbose) {
    cat("\n\n===== RESULTS TABLE (sorted by NLL) =====\n")
    print(round(tab, 4L))
    w <- tab[1L, ]
    cat("\n===== BEST MODEL =====\n")
    cat(sprintf("  k     = %d\n",               as.integer(w$k)))
    cat(sprintf("  R0    = %.4f\n", w$R0))
    cat(sprintf("  beta  = %.4f  (= R0*lam/k)\n", w$beta))
    cat(sprintf("  mu    = %.4f\n", w$mu))
    cat(sprintf("  sigma = %.4f\n", w$sigma))
    cat(sprintf("  lam   = %.4f  [fixed]\n",    w$lam))
    cat(sprintf("  pobs  = %.4f  [fixed]\n",    w$pobs))
    cat(sprintf("  NLL   = %.4f\n",             w$nll))
  }
  
  list(table = tab, best_k = as.integer(tab[1L, "k"]),
       best_R0 = tab[1L, "R0"], results = all_res)
}


# =============================================================================
# 7.  Profile likelihood plot (kept for reference, not used in this analysis)
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
# 8.  Entry point: Sensitivity analysis over pobs
# =============================================================================

if (sys.nframe() == 0L) {
  
  # --- Load data ------------------------------------------------------------
  f <- "full_tree_edges.csv"
  if (!file.exists(f)) stop("Run full_tree_sim.py first to produce ", f)
  
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
  
  # --- pobs grid ------------------------------------------------------------
  pobs_seq <- seq(0.1, 0.9, by = 0.1)
  results <- data.frame(pobs = pobs_seq, R0_hat = NA, k_hat = NA)
  
  cat("\n=== Sensitivity analysis: fix pobs and estimate R0 ===\n")
  
  for (i in seq_along(pobs_seq)) {
    p <- pobs_seq[i]
    cat(sprintf("\n--- pobs = %.2f ---\n", p))
    
    # Run estimation with fixed pobs = p
    est <- tryCatch(
      estimate_full_tree(dat, k_grid = 4L:5L, lam = LAM, pobs = p,
                         maxTime = maxT, dt = 0.002, verbose = FALSE),
      error = function(e) NULL
    )
    
    if (!is.null(est) && !is.na(est$best_R0) && is.finite(est$best_R0)) {
      results[i, "R0_hat"] <- est$best_R0
      results[i, "k_hat"]  <- est$best_k
      cat(sprintf("  Best: k = %d, R0 = %.4f\n", est$best_k, est$best_R0))
    } else {
      cat("  Estimation failed for this pobs.\n")
    }
  }
  
  # --- Plot results ---------------------------------------------------------
  # Remove rows where estimation failed
  results_ok <- results[!is.na(results$R0_hat), ]
  
  if (nrow(results_ok) > 0) {
    # Print summary
    cat("\n\n===== Sensitivity analysis summary =====\n")
    print(results_ok)

    # Save results for later use
    write.csv(results_ok, "ro_vs_pobs_sensitivity_results.csv", row.names = FALSE)
    cat("\nResults saved to 'sensitivity_results.csv'\n")
  } else {
    cat("\nNo successful estimations.\n")
  }
} else {
  # If script is sourced (not run as main), check if saved results exist
  if (file.exists("sensitivity_results.csv")) {
    results_ok <- read.csv("sensitivity_results.csv")
    cat("Loaded saved results from 'sensitivity_results.csv'\n")
    cat(sprintf("Loaded %d sensitivity analysis points.\n", nrow(results_ok)))
  } else {
    stop("No saved results found. Run the full analysis first.\n")
  }
}


# Set up plotting area
par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)

# Plot R0_hat vs pobs
plot(results_ok$pobs, results_ok$R0_hat,
     type = "b", pch = 17, col = "black", cex=1.5, lwd = 2, lty = 1,
     xlab = expression(p[obs]), 
     ylab = expression(hat(R)[0]),
     main = "",
     xlim = c(0, 1), ylim = c(0, max(results_ok$R0_hat, 6) + 1))

# Add true R0 line
abline(h = 6, col = "gray50", lty = 2, lwd = 2)

# Optionally add k_hat as text
text(results_ok$pobs, results_ok$R0_hat + 0.3,
     labels = paste("k=", results_ok$k_hat),
     cex = 0.8, col = "black")

legend("topright", cex=1.0,
       legend = c("Estimated Ro", "True Ro = 6", "k (best)"),
       col = c("black", "gray50", "black"),
       lty = c(1, 2, NA), pch = c(17, NA, NA),
       lwd = c(2, 2, NA), bty = "n")
