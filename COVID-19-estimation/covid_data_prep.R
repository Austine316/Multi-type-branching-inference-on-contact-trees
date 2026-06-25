# =============================================================================
# covid_data_prep.R  —  Build edge table from Karnataka COVID-19 data
#                       (ann2.csv + contacts.csv)
# =============================================================================
#
# DATA SOURCES
# ------------
# ann2.csv     : 71,068 individuals. Columns used:
#                  ID, parentid, children_pri, dt_conf, cluster_id, cat_4
#                Restricted to March-May 2020 (early Karnataka outbreak).
# contacts.csv : 5,105 directed links (from -> to) = who infected whom.
#                These are exactly the parentid > 0 individuals.
#
# MODEL MAPPING
# -------------
# State i = number of secondary transmissions already made by individual i.
#
# OBSERVED per individual:
#   t_conf = confirmation date in model time (detection = sigma event)
#   j_obs  = children_pri = contacts already infected BY detection time
#            = individual's STATE at the moment of being detected
#
# LATENT per individual:
#   tau_born         = infection time (UNKNOWN; only t_conf is observed)
#   tau_1,...,tau_j  = times of the j_obs within-individual transmissions
#
# TIPS (enter the sigma likelihood):
#   parentid == 0 cases: detected by surveillance (sigma process).
#   j_obs = children_pri is their observed state.
#   Their tip density term is D^{j_obs}_0(delta) where delta = t_conf - tau_born.
#
# NON-TIPS (parentid > 0, Local Traced):
#   Detected via contact tracing, NOT the sigma process.
#   They contribute an E (extinction) term to the likelihood.
#   Their t_conf gives a FEASIBILITY UPPER BOUND for their parent's
#   transmission time to them (= their own tau_born).
#
# TIME IS LATENT
# --------------
# We observe t_conf for everyone; we do NOT observe tau_born.
# tau_b of the last segment of each individual = t_conf (observed, fixed).
# tau_a of the first segment                  = tau_born (LATENT, proposed by MCMC).
#
# MCMC proposes tau_born within:
#   For index cases (parentid=0):
#     lo = max(0, t_conf - MAX_INFECTIOUS)
#     hi = t_conf - MIN_INCUB
#   For secondary cases (parentid>0):
#     lo = max(0, t_conf(parent) - MAX_SI)   [parent born before transmitting]
#     hi = t_conf(self) - MIN_INCUB           [self born before own detection]
#
# Within-individual branching times tau_1,...,tau_j are also latent and
# proposed by the MCMC (Step B) within (tau_born, t_conf).
#
# TIME SCALE
# ----------
# Mean serial interval = 5.88 days (from 54_serial_interval_data.xlsx).
# Fix gamma = 1 and scale all calendar times by 1/MEAN_SI.
# 1 model unit = 5.88 days.  MIN_INCUB = 1 day = 0.17 model units.
# =============================================================================

MEAN_SI        <- 5.88
SCALE          <- 1.0 / MEAN_SI
REF_DATE       <- as.Date("2020-03-09")
MIN_INCUB      <- 1.0  * SCALE   # ~1 day minimum incubation (model units)
MAX_SI         <- 18.0 * SCALE   # 18 days max serial interval (model units)
MAX_INFECTIOUS <- 21.0 * SCALE   # generous maximum infectious period
MIN_DELTA      <- 0.001          # minimum edge length (model units)

cat("=== COVID-19 Data Preparation ===\n")
cat(sprintf("Source: ann2.csv + contacts.csv\n"))
cat(sprintf("Time scale: 1 model unit = %.2f days (mean serial interval)\n", MEAN_SI))
cat(sprintf("Feasibility bounds: MIN_INCUB=%.3f  MAX_SI=%.3f  MAX_INFECT=%.3f\n\n",
            MIN_INCUB, MAX_SI, MAX_INFECTIOUS))

# =============================================================================
# 1. Load
# =============================================================================
cat("Loading data...\n")
ann <- read.csv("covid-19-data/ann2.csv",      stringsAsFactors = FALSE, fileEncoding = "UTF-8")
mean(ann$cluster_size)
table(ann$cluster_size)
table(ann$children_pri)

hist(ann$cluster_size, breaks = 40)
hist(ann$children_pri, breaks = 40)

contacts <- read.csv("covid-19-data/contacts.csv",  stringsAsFactors = FALSE)

# =============================================================================
# 2. Restrict to March-May 2020
# =============================================================================
ann_early <- ann[grepl("Mar|Apr|May", ann$dt_conf, ignore.case = TRUE), ]

# Parse confirmation dates
parse_dtconf <- function(x) {
  d <- suppressWarnings(as.Date(paste0(x, "-2020"), format = "%d-%b-%Y"))
  bad <- is.na(d)
  if (any(bad))
    d[bad] <- suppressWarnings(as.Date(x[bad], format = "%d.%m.%Y"))
  d
}
ann_early$t_conf_date <- parse_dtconf(ann_early$dt_conf)
ann_early$t_conf      <- as.numeric(ann_early$t_conf_date - REF_DATE) * SCALE
ann_early             <- ann_early[!is.na(ann_early$t_conf), ]

POBS_EMP <- 0.7454

cat(sprintf("March-May 2020: %d individuals\n", n_total))
cat(sprintf("  Empirical p_obs = %.4f\n\n", POBS_EMP))

# Build parent lookup: ID -> t_conf
parent_tconf <- setNames(ann_early$t_conf, ann_early$ID)
parent_jpri  <- setNames(ann_early$children_pri, ann_early$ID)

# =============================================================================
# 3. Build edge table
# =============================================================================
# Every cluster with at least one index case becomes a replicate.
# Single-individual clusters (one index case, no traced) still contribute
# a single tip likelihood term.

usable_clusters <- sort(unique(ann_early$cluster_id[ann_early$parentid == 0L]))
cat(sprintf("Building edges for %d clusters...\n", length(usable_clusters)))

build_cluster_edges <- function(cid) {
  cl <- ann_early[ann_early$cluster_id == cid, ]
  if (!any(cl$parentid == 0L)) return(NULL)
  cl <- cl[order(cl$parentid, cl$t_conf), ]
  
  row_list <- list()
  
  for (i in seq_len(nrow(cl))) {
    ind    <- cl[i, ]
    id     <- ind$ID
    par    <- as.integer(ind$parentid)
    t_det  <- ind$t_conf
    j_gen  <- max(0L, as.integer(ind$children_pri))
    is_tip <- (par == 0L)
    
    # ── Feasibility bounds for tau_born (LATENT infection time) ────────────
    if (par == 0L) {
      # Index case: no parent in data
      lo_born <- max(0.0, t_det - MAX_INFECTIOUS)
      hi_born <- t_det - MIN_INCUB
    } else {
      # Secondary case: parent's detection time upper-bounds transmission
      t_par <- parent_tconf[as.character(par)]
      if (is.na(t_par)) t_par <- max(0.0, t_det - MAX_SI)
      lo_born <- max(0.0, t_par - MAX_SI)   # generous lower bound
      hi_born <- t_det - MIN_INCUB           # must be infected before own detection
    }
    if (!is.finite(lo_born) || !is.finite(hi_born) || hi_born <= lo_born)
      hi_born <- lo_born + MIN_DELTA
    
    # Initial tau_born: one model unit before detection (rough guess for MCMC start)
    tau_born_init <- max(lo_born, min(hi_born - MIN_DELTA, t_det - 1.0))
    
    # ── Within-individual branching times (also LATENT) ────────────────────
    # Initial values: evenly spaced in (tau_born_init, t_det)
    if (j_gen > 0L) {
      bt <- seq(tau_born_init, t_det, length.out = j_gen + 2L)[2L:(j_gen + 1L)]
      bt <- pmin(pmax(bt, tau_born_init + MIN_DELTA), t_det - MIN_DELTA)
    } else {
      bt <- numeric(0L)
    }
    
    n_segs    <- j_gen + 1L
    j_obs_val <- if (is_tip) as.integer(j_gen) else NA_integer_
    
    # ── One row per segment ────────────────────────────────────────────────
    for (seg in 0L:j_gen) {
      tau_a_s <- if (seg == 0L)   tau_born_init else bt[seg]
      tau_b_s <- if (seg < j_gen) bt[seg + 1L]  else t_det
      
      if (!is.finite(tau_a_s) || !is.finite(tau_b_s)) next
      if (tau_b_s <= tau_a_s) tau_b_s <- tau_a_s + MIN_DELTA
      
      is_last    <- (seg == j_gen)
      s_val      <- as.integer(seg)
      c_val      <- if (is_last) {
        if (is_tip) as.integer(j_gen) else s_val
      } else as.integer(s_val + 1L)
      has_branch <- if (!is_last) 1.0 else NA_real_
      
      row_list[[length(row_list) + 1L]] <- data.frame(
        rep_id      = as.integer(cid),
        ind_id      = as.integer(id),
        par_id      = as.integer(par),
        is_tip      = as.integer(is_tip),
        s           = s_val,
        c           = c_val,
        j_obs       = j_obs_val,
        tau_a       = tau_a_s,
        tau_b       = tau_b_s,
        delta       = tau_b_s - tau_a_s,
        seg_idx     = as.integer(seg),
        n_segs      = as.integer(n_segs),
        branch_rate = has_branch,
        n_tips      = as.integer(sum(cl$parentid == 0L)),
        root_state  = 0L,
        # Feasibility bounds for tau_born — used by MCMC for birth time proposals
        lo_born     = lo_born,
        hi_born     = hi_born,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(row_list) == 0L) return(NULL)
  do.call(rbind, row_list)
}

edge_list   <- lapply(usable_clusters, function(cid)
  tryCatch(build_cluster_edges(cid), error = function(e) NULL))
edge_list   <- Filter(Negate(is.null), edge_list)
covid_edges <- do.call(rbind, edge_list)
rownames(covid_edges) <- NULL

# Remove degenerate rows
covid_edges <- covid_edges[
  is.finite(covid_edges$delta) & covid_edges$delta > 0 &
    !is.nan(covid_edges$s) & !is.nan(covid_edges$c), ]

# =============================================================================
# 4. Summary
# =============================================================================
tips_rows <- covid_edges[!is.na(covid_edges$j_obs) &
                           covid_edges$seg_idx == covid_edges$n_segs - 1L, ]
n_clusters   <- length(unique(covid_edges$rep_id))
n_with_trace <- sum(tapply(covid_edges$par_id, covid_edges$rep_id,
                           function(x) any(x > 0L)))

cat(sprintf("\nEdge table complete:\n"))
cat(sprintf("  Rows:          %d\n", nrow(covid_edges)))
cat(sprintf("  Clusters:      %d\n", n_clusters))
cat(sprintf("  Tips (j_obs):  %d\n", nrow(tips_rows)))
cat(sprintf("  With traced:   %d clusters have >= 1 contact-traced member\n",
            n_with_trace))
cat(sprintf("  Mean n_tips per cluster: %.2f\n",
            mean(tapply(covid_edges$is_tip, covid_edges$rep_id, sum))))

cat("\nj_obs distribution:\n")
print(table(tips_rows$j_obs))

write.csv(covid_edges, "covid_edges.csv", row.names = FALSE)
cat(sprintf("\nSaved: covid_edges.csv\n"))

# =============================================================================
# 5. Recommended MCMC settings
# =============================================================================
cat("\n=== RECOMMENDED SETTINGS FOR run_covid_mcmc.R ===\n")
cat(sprintf("LAM_FIX   <- 1.0       # gamma = 1 (model time; 1 unit = %.2f days)\n", MEAN_SI))
cat(sprintf("POBS_FIX  <- %.4f   # empirical: %d / %d\n", POBS_EMP, n_index, n_total))
cat(sprintf("K_MAX     <- 60L       # upper bound for k prior\n"))
cat(sprintf("K_INIT    <- 4L        # starting value for k chain\n"))
cat(sprintf("STEP_K    <- 2L        # random walk step for k\n"))
cat(sprintf("GLOBAL_MAXTIME to set from: max(covid_edges$tau_b) + 1 = %.2f\n",
            max(covid_edges$tau_b, na.rm = TRUE) + 1.0))

