# =============================================================================
# mle_results_summary.R
# =============================================================================
# Produces all outputs for the full-tree MLE section:
#   1. MLE across k: estimates, CIs, AIC weights, CSV, LaTeX table
#   2. R0 vs k plot with 95% CIs
#   3. Equilibrium frequencies: theory vs simulation

# Usage: source("mle_results_summary.R")
# =============================================================================

cat('Current working directory', getwd())
setwd(getwd())

suppressPackageStartupMessages(library(deSolve))
if (!exists("neg_log_lik")) source("parameter_estimate_mle.R")

# =============================================================================
# 0.  True simulation parameters
# =============================================================================
TRUE_K    <- 4L
TRUE_BETA <- 1.5
TRUE_LAM  <- 1.0
TRUE_POBS <- 0.5
TRUE_MU   <- (1 - TRUE_POBS) * TRUE_LAM
TRUE_SIG  <- TRUE_POBS  * TRUE_LAM
TRUE_R0   <- TRUE_K * TRUE_BETA / TRUE_LAM

LAM_FIX  <- TRUE_LAM
POBS_FIX <- TRUE_POBS
K_GRID   <- 1L:12L

# Greyscale palette
C_BLACK <- "black"
C_DARK  <- "gray20"
C_MID   <- "gray50"
C_LIGHT <- "gray75"
C_FILL  <- "gray88"

cat(sprintf("True parameters: k=%d, beta=%.2f, gamma=%.1f, pobs=%.1f => R0=%.1f\n\n",
            TRUE_K, TRUE_BETA, TRUE_LAM, TRUE_POBS, TRUE_R0))

# =============================================================================
# 1.  Load data
# =============================================================================
cat("Loading edge table...\n")
raw <- read.csv("full_tree_edges.csv", stringsAsFactors = FALSE)
dat <- raw[is.finite(raw$delta) & raw$delta >= 0 &
             !is.nan(raw$s) & !is.nan(raw$c), ]
dat$j_obs[is.nan(dat$j_obs)] <- NA_integer_
rownames(dat) <- NULL

n_total <- nrow(dat)
n_tips <- sum(!is.na(dat$j_obs))

cat(sprintf("Total edge rows:  %d\n", n_total))
cat(sprintf("Observed tips:    %d\n", n_tips))
cat(sprintf("Empirical p_obs:  %d / %d = %.4f\n",
            n_tips, n_total, n_tips / n_total))
# =============================================================================
# 1.  MLE for each k
# =============================================================================
cat("Running MLE across k =", paste(K_GRID, collapse = ","), "...\n")
results <- list()
N_rep = dim(dat)[1]
for (kk in K_GRID) {
  res_k <- tryCatch(
    optimize(neg_log_lik, interval = c(0.5, 30),
             dat = dat, k = kk, lam = LAM_FIX, pobs = POBS_FIX,
             maxTime = maxT, dt = 0.002, tol = 1e-4),
    error = function(e) list(minimum = NA_real_, objective = Inf)
  )
  R0h <- res_k$minimum;  nll <- res_k$objective
  if (!is.finite(nll)) next
  beta_h <- R0h * LAM_FIX / kk
  np     <- 2L
  results[[as.character(kk)]] <- list(
    k = kk, R0 = R0h, beta = beta_h, kbeta = kk * beta_h,
    mu = (1 - POBS_FIX)*LAM_FIX, sigma = POBS_FIX*LAM_FIX,
    lam = LAM_FIX, pobs = POBS_FIX, nll = nll,
    aic = 2*np + 2*nll, bic = log(N_rep)*np + 2*nll
  )
  cat(sprintf("  k=%2d:  R0=%.4f  beta=%.4f  k*beta=%.4f  NLL=%.2f\n",
              kk, R0h, beta_h, kk*beta_h, nll))
}

tab  <- do.call(rbind, lapply(results, as.data.frame))
rownames(tab) <- NULL
tab  <- tab[order(tab$nll), ]
best <- tab[1L, ]
cat(sprintf("\nBest model: k=%d  R0=%.4f  beta=%.4f\n\n",
            as.integer(best$k), best$R0, best$beta))


# =============================================================================
# 2.  Profile-likelihood 95% CIs
# =============================================================================
cat("Computing 95% profile-likelihood CIs...\n")
tab$ci_lo <- NA_real_;  tab$ci_hi <- NA_real_

for (i in seq_len(nrow(tab))) {
  kk      <- as.integer(tab$k[i])
  r0h     <- tab$R0[i];  nll_min <- tab$nll[i]
  ci_cut  <- nll_min + qchisq(0.95, df = 1L) / 2.0
  grid    <- seq(max(0.5, r0h - 3), r0h + 3, by = 0.05)
  nlls    <- sapply(grid, function(r0)
    neg_log_lik(r0, dat, k = kk, lam = LAM_FIX, pobs = POBS_FIX,
                maxTime = maxT, dt = 0.005))
  ok <- nlls <= ci_cut
  tab$ci_lo[i] <- if (any(ok)) round(min(grid[ok]), 3) else NA_real_
  tab$ci_hi[i] <- if (any(ok)) round(max(grid[ok]), 3) else NA_real_
  cat(sprintf("  k=%2d:  CI = [%.3f, %.3f]\n", kk, tab$ci_lo[i], tab$ci_hi[i]))
}


# =============================================================================
# 3.  AIC weights and model-averaged R0
# =============================================================================
tab_valid       <- tab[tab$k >= 1L & tab$k <= 12L, ]
aic_min         <- min(tab_valid$aic)
aic_w           <- exp(-0.5*(tab_valid$aic - aic_min)); aic_w <- aic_w/sum(aic_w)
R0_avg          <- sum(aic_w * tab_valid$R0)
tab_valid$dAIC  <- round(tab_valid$aic - aic_min, 2)
tab_valid$aic_w <- round(aic_w, 5)
tab_valid$bias  <- round(tab_valid$R0 - TRUE_R0, 4)
tab_valid$beta_ci_lo <- round(tab_valid$ci_lo / tab_valid$k, 3)
tab_valid$beta_ci_hi <- round(tab_valid$ci_hi / tab_valid$k, 3)

best_k  <- as.integer(best$k)
best_ci <- tab_valid[tab_valid$k == best_k,
                     c("ci_lo","ci_hi","beta_ci_lo","beta_ci_hi")]

cat(sprintf("\nAIC-weighted R0 average = %.4f  (true = %.1f, bias = %+.4f)\n\n",
            R0_avg, TRUE_R0, R0_avg - TRUE_R0))
cat("=== BEST MODEL SUMMARY ===\n")
cat(sprintf("  k               = %d\n",   best_k))
cat(sprintf("  R0 MLE          = %.4f  [true %.1f]\n", best$R0, TRUE_R0))
cat(sprintf("  95%% CI (R0)    = [%.3f, %.3f]\n", best_ci$ci_lo,      best_ci$ci_hi))
cat(sprintf("  beta MLE        = %.4f  [true %.1f]\n", best$beta, TRUE_BETA))
cat(sprintf("  95%% CI (beta)  = [%.3f, %.3f]\n", best_ci$beta_ci_lo, best_ci$beta_ci_hi))
cat(sprintf("  True R0 in CI?  = %s\n\n",
            if (!is.na(best_ci$ci_lo) &&
                TRUE_R0 >= best_ci$ci_lo && TRUE_R0 <= best_ci$ci_hi) "YES" else "NO"))

write.csv(tab_valid, "mle_results_table.csv", row.names = FALSE)
cat("Saved mle_results_table.csv\n")
print(round(tab_valid[, c("k","R0","beta","kbeta","nll","aic","dAIC",
                          "ci_lo","ci_hi","bias","aic_w")], 4))


# ==================== LOAD RESULTS IF ALREADY EXIST ==========================
tab_valid <- read.csv("mle_results_table.csv", header = TRUE)
# =============================================================================

# =============================================================================
# 4.  Plot 1: R0 vs k with 95% CIs
# =============================================================================
tab_p <- tab_valid[order(tab_valid$k), ]
ok_ci <- !is.na(tab_p$ci_lo)
ylim  <- c(min(tab_p$ci_lo[ok_ci], na.rm = TRUE) - 0.4,
           max(tab_p$ci_hi[ok_ci], na.rm = TRUE) + 0.6)

par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
plot(tab_p$k, tab_p$R0,
     type = "b", pch = 19, cex = 1.2, col = C_DARK, lwd = 2.0,
     xlab = "k",
     ylab = expression(hat(R)[0]),
     ylim = ylim, xaxt = "n")
axis(1, at = tab_p$k)

arrows(tab_p$k[ok_ci], tab_p$ci_lo[ok_ci],
       tab_p$k[ok_ci], tab_p$ci_hi[ok_ci],
       length = 0.06, angle = 90, code = 3, col = C_DARK, lwd = 1.2)

abline(h = TRUE_R0, col = C_BLACK, lwd = 2.0, lty = 2)
abline(v = TRUE_K,  col = C_MID,   lwd = 2.0, lty = 3)

best_row <- tab_p[tab_p$k == best_k, ]
points(best_row$k, best_row$R0, pch = 8, cex = 1.8, col = C_BLACK, lwd = 2)

legend("bottomright", bty = "n", cex = 0.80,
       legend = c(expression(hat(R)[0] * " with 95% CI"),
                  sprintf("True Ro = %.1f", TRUE_R0),
                  sprintf("True k = %d", TRUE_K),
                  sprintf("Best k = %d (AIC)", best_k)),
       col = c(C_DARK, C_BLACK, C_MID, C_BLACK),
       lty = c(1, 2, 3, NA), lwd = c(2.0, 2.0, 2.0, NA),
       pch = c(19, NA, NA, 8))


# =============================================================================
# 5.  Plot 2: Equilibrium frequencies — theory vs simulation
# =============================================================================
states    <- 0L:TRUE_K
MLE_BETA  <- best$R0 / best_k
pi_theory <- compute_pi(MLE_BETA, TRUE_MU, TRUE_SIG, best_k)

root_per_rep <- tapply(dat$root_state, dat$rep_id, function(x) x[1])
emp_counts   <- table(factor(root_per_rep, levels = states))
emp_pi       <- as.numeric(emp_counts) / sum(emp_counts)

chi2 <- sum((as.numeric(emp_counts) - N_rep*pi_theory)^2 / (N_rep*pi_theory))
pval <- pchisq(chi2, df = TRUE_K, lower.tail = FALSE)
cat(sprintf("\nEquilibrium pi goodness-of-fit: X2 = %.3f, df = %d, p = %.4f\n",
            chi2, TRUE_K, pval))
cat("Theory:   ", round(pi_theory, 4), "\n")
cat("Empirical:", round(emp_pi, 4), "\n")

par(mar = c(5, 4.5, 3.5, 2), cex.main = 0.95, cex.lab = 1.0, cex.axis = 1.0)
bp <- barplot(
  rbind(pi_theory, emp_pi),
  beside    = TRUE,
  col       = c(C_FILL, C_MID),
  border    = c(C_DARK, C_BLACK),
  names.arg = paste0("i = ", states),
  xlab      = "i",
  ylab      = expression(pi[i]),
  ylim      = c(0, max(pi_theory, emp_pi) * 1.35),
  cex.names = 0.90, las = 1
)

text(bp[1,], pi_theory + max(pi_theory)*0.04,
     labels = round(pi_theory, 3), col = C_DARK,  cex = 0.75)
text(bp[2,], emp_pi    + max(pi_theory)*0.04,
     labels = round(emp_pi,    3), col = C_BLACK, cex = 0.75)

legend("topright", bty = "n", cex = 0.85,
       legend = c("Theory", "Empirical"),
       fill   = c(C_FILL, C_MID),
       border = c(C_DARK, C_BLACK),
       )
       
mtext(sprintf("Chi-squared goodness-of-fit:  X2(%d) = %.2f,  p = %.3f",
              TRUE_K, chi2, pval),
      side = 1, line = 3.8, cex = 0.78, col = C_MID)


# =============================================================================
# 6.  LaTeX table
# =============================================================================
cat("\n\n% ============================================================\n")
cat("% LaTeX table\n")
cat("% ============================================================\n\n")

writeLines("\\begin{table}[ht]")
writeLines("\\centering")
writeLines(paste0(
  "\\caption{Maximum-likelihood estimates of $R_0$, $\\hat{\\beta}$,",
  " and $k\\hat{\\beta}$ across degrees $k=1,\\ldots,12$ from $N=",
  N_rep, "$ simulated trees.",
  " Fixed: $\\gamma=", LAM_FIX, "$, $p_{\\mathrm{obs}}=", POBS_FIX, "$.",
  " True values: $k=", TRUE_K, "$, $R_0=", TRUE_R0, "$, $\\beta=", TRUE_BETA, "$.",
  " Profile-likelihood 95\\% confidence intervals shown throughout.",
  " $\\Delta$AIC relative to the best model ($k=", TRUE_K, "$, starred).}"
))
writeLines("\\label{tab:mle_full}")
writeLines("\\begin{tabular}{crrrrrrr}")
writeLines("\\toprule")
writeLines(paste(
  "$k$ & $\\hat{R}_0$ & 95\\% CI & $\\hat{\\beta}$ & 95\\% CI",
  "& $k\\hat{\\beta}$ & $\\Delta$AIC & AIC wt \\\\"
))
writeLines("\\midrule")

for (i in seq_len(nrow(tab_valid))) {
  r    <- tab_valid[i, ]
  star <- if (r$k == best_k) "$^*$" else ""
  ci_r <- if (!is.na(r$ci_lo))
    paste0("[", format(round(r$ci_lo,3), nsmall=3), ",\\,",
           format(round(r$ci_hi,3), nsmall=3), "]") else "---"
  ci_b <- if (!is.na(r$beta_ci_lo))
    paste0("[", format(round(r$beta_ci_lo,3), nsmall=3), ",\\,",
           format(round(r$beta_ci_hi,3), nsmall=3), "]") else "---"
  writeLines(paste0(
    r$k, star, " & ", round(r$R0,3), " & ", ci_r,
    " & ", round(r$beta,4), " & ", ci_b,
    " & ", round(r$kbeta,3),
    " & ", round(r$dAIC,1),
    " & ", round(r$aic_w,4), " \\\\"
  ))
}
writeLines("\\midrule")
writeLines(paste0(
  "\\multicolumn{8}{l}{AIC-weighted model-average: $\\hat{R}_0^{\\mathrm{avg}}=",
  round(R0_avg,3), "$;\\enspace bias $=", sprintf("%+.3f", R0_avg-TRUE_R0),
  "$;\\enspace $p_{\\mathrm{obs}}=", POBS_FIX, "$ fixed.}\\\\"))
writeLines("\\bottomrule")
writeLines("\\end{tabular}")
writeLines("\\end{table}")

cat("\n\n=== DONE ===\n")
cat("mle_results_table.csv saved.\n")
cat("Three plots drawn — save each via RStudio Export.\n")
cat("LaTeX table printed above.\n")