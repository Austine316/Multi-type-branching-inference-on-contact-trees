# =============================================================================
# run_partial_resolved_tree  —  Joint MCMC over (R0, k, branching times)
#                   with state marginalisation over s in {0,...,k}
# =============================================================================
# FAST TEST (~5 min):  uncomment block 1
# PRODUCTION (~40 min): block 2 active
# =============================================================================

cat('Current working directory', getwd())
setwd(getwd())

## FAST TEST
# N_REPS_USE <- 10L; N_ITER <- 300L; N_BURN <- 60L; THIN <- 2L
# P_RES_VALS <- c(0.0, 0.5, 1.0); ODE_DT <- 0.05

## PRODUCTION
N_REPS_USE <- 30L
N_ITER     <- 2000L
N_BURN     <- 400L
THIN       <- 2L
P_RES_VALS <- c(0.0, 0.25, 0.5, 0.75, 1.0)
ODE_DT     <- 0.02
SEED       <- 234L
SAVE_CHAINS <- TRUE     # saves partial_chains_p{p_res}.rds per level

source("partial_resolved_tree.R")

# =============================================================================
# To reload and replot without re-running:
# =============================================================================
# saved   <- readRDS("partial_tree_all_results.rds")
# results <- saved$results
# P_RES_VALS <- saved$P_RES_VALS
#
# # Load individual chain for p_res=0:
# ch <- readRDS("partial_chains_p0.00.rds")
# plot(density(ch$R0_chain), main="R0 posterior (p_res=0)")
# plot(table(ch$k_chain)/length(ch$k_chain), main="k posterior (p_res=0)")
