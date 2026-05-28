# =============================================================================
# run_covid_mcmc.R  —  Entry point for COVID-19 joint R0 and k estimation
# =============================================================================
# USAGE:  source("run_covid_mcmc.R")
#
# Required files in working directory:
#   traced.csv, contacts.csv        (raw Karnataka COVID-19 data)
#   parameter_estimation_mle.R      (likelihood with latent_nodes + ED cache)
#   covid_data_prep.R               (builds covid_edges.csv)
#   covid_mcmc.R                    (three-step MCMC sampler)
# =============================================================================

# ── FAST TEST (uncomment to use) ─────────────────────────────────────────────
# N_ITER <- 300L;  N_BURN <- 60L;  THIN <- 2L;  ODE_DT <- 0.20

# ── PRODUCTION ────────────────────────────────────────────────────────────────
N_ITER <- 2000L
N_BURN <- 400L
THIN   <- 2L
ODE_DT <- 0.02

# ── Fixed epidemiological parameters ─────────────────────────────────────────
LAM_FIX    <- 1.0     # gamma = 1 (time unit = 5.65 days mean serial interval)
POBS_FIX   <- 0.7454  # empirical: 2401 / 3221 (March-May 2020, parentid=0 / total)
INIT_RW    <- 0.15
TARGET_ACC <- 0.23
ADAPT_INT  <- 50L
SEED       <- 1234L

# ── k prior: Negative Binomial ───────────────────────────────────────────────
# k ~ NegBin(mu = NBINOM_MU, size = NBINOM_PHI), truncated to {1,...,K_MAX}
# NBINOM_MU:  prior mean contact degree
# NBINOM_PHI: dispersion (smaller = more overdispersed; 0.5 reflects
#             superspreading heterogeneity typical of COVID-19)
# K_MAX:      upper truncation point; set to cover plausible range
# K_INIT:     starting value for k chain
# STEP_K:     random walk step size for k proposals
NBINOM_MU  <- 11.8    # prior mean
NBINOM_PHI <- 0.29    # dispersion (phi -> 0: very overdispersed; phi -> Inf: Poisson)
K_MAX      <- 60L     # upper bound
K_INIT     <- 11L     # starting value
STEP_K     <- 2L      # step size for random walk on k

# ── Initial R0 ────────────────────────────────────────────────────────────────
INIT_R0 <- 1.843   # updated below by quick MLE

# =============================================================================
# Step 1: Build edge table
# =============================================================================
if (!file.exists("covid_edges.csv")) {
  cat("Building edge table from raw data...\n")
  source("covid_data_prep.R")
} else {
  cat("Using existing covid_edges.csv\n\n")
}

# =============================================================================
# Step 2: Load likelihood functions
# =============================================================================
source("parameter_estimation_mle.R")

# =============================================================================
# Step 3: Load and validate edge table
# =============================================================================
covid_edges <- read.csv("covid_edges.csv", stringsAsFactors = FALSE)
covid_edges$j_obs[is.nan(covid_edges$j_obs)] <- NA_integer_
covid_edges <- covid_edges[
  is.finite(covid_edges$delta) & covid_edges$delta >= 0 &
    !is.nan(covid_edges$s)       & !is.nan(covid_edges$c), ]

cat(sprintf("Edge table: %d rows, %d clusters\n", nrow(covid_edges),
            length(unique(covid_edges$rep_id))))
tips_rows     <- covid_edges[!is.na(covid_edges$j_obs), ]
max_j_obs_inf <- as.integer(max(tips_rows$j_obs, na.rm = TRUE))
cat(sprintf("Max j_obs in data: %d  (informational; k < this value is penalised by likelihood)\n\n",
            max_j_obs_inf))

# =============================================================================
# Step 4: Quick MLE at K_INIT for R0 starting value
# =============================================================================
cat("=== Quick MLE for R0 starting value ===\n")
MAXTIME_MLE <- max(covid_edges$tau_b, na.rm = TRUE) + 0.5
res_mle <- tryCatch(
  optimize(neg_log_lik, interval = c(0.3, 8.0),
           dat = covid_edges, k = K_INIT, lam = LAM_FIX, pobs = POBS_FIX,
           maxTime = MAXTIME_MLE, dt = 0.20, tol = 1e-3),
  error = function(e) list(minimum = INIT_R0, objective = Inf)
)
if (is.finite(res_mle$objective)) {
  INIT_R0 <- res_mle$minimum
  cat(sprintf("MLE R0 at k=%d: %.3f  (NLL = %.2f)\n\n", K_INIT, INIT_R0, res_mle$objective))
} else {
  cat(sprintf("MLE failed — using default INIT_R0 = %.3f\n\n", INIT_R0))
}

# =============================================================================
# Step 5: Run MCMC
# =============================================================================
source("covid_mcmc.R")

