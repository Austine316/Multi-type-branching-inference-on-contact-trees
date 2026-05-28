library(deSolve)
library(splines)

cat('Current working directory', getwd())
setwd(getwd())

##########################################
## READ DATA
##########################################
readTable <- function(i_foc) {
  filename <- "Simulated data/phylo-epi-sim-data-fixed.csv"
  dat <- read.table(file = filename, sep = ',', header = TRUE)
  dat <- dat[dat$i_foc == i_foc, ]
  dat$t_obs <- dat$t_obs - dat$t_foc
  N <- nrow(dat)
  list(dat = dat, N = N, i_foc = i_foc)
}

##########################################
## PARAMATERS
##########################################

# Global Param
k       <- 4
maxTime <- 3.0

# SIM Param
hh      <- 0.05
nn      <- maxTime / hh
ttStep  <- (0:nn) * hh         # length nn+1
tt      <- (0:(nn - 1)) * hh   # length nn

# ODE Param
beta    <- 1.5
mu      <- 0.5
sigma   <- 0.5
dt      <- 0.01

##########################################
## SIMULATION SOLVER SECTION
##########################################

# Bin-by-time helper
solveSim <- function(df, N, obsType) {
  frac_D <- frac_E <- rep(NA_real_, nn)
  for (ii in 1:nn) {
    in_bin <- (df$t_obs >= ttStep[ii]) & (df$t_obs < ttStep[ii + 1])
    res_D  <- sum(in_bin & (df$j_obs == obsType))
    res_E  <- sum(df$t_obs >= ttStep[ii])
    frac_D[ii] <- res_D / (N * hh)  # density per unit time
    frac_E[ii] <- res_E / N         # survival fraction
  }
  list(frac_D = frac_D, frac_E = frac_E)
}

# Containers
simDRes <- vector("list", length = k + 1)
simERes <- vector("list", length = k + 1)
names(simDRes) <- paste0("j_", 0:k)
names(simERes) <- paste0("j_", 0:k)

# Main loop: for each observed state j, collect columns i_0..i_k + time
for (obs in 0:k) {
  D_states <- vector("list", length = k + 1)
  E_states <- vector("list", length = k + 1)
  for (foc in 0:k) {
    si  <- readTable(foc)
    out <- solveSim(si$dat, si$N, obs)
    D_states[[foc + 1]] <- out$frac_D
    E_states[[foc + 1]] <- out$frac_E
  }
  names(D_states) <- paste0("i_", 0:k)
  names(E_states) <- paste0("i_", 0:k)
  
  df_D <- data.frame(D_states, time = tt, check.names = FALSE)
  df_E <- data.frame(E_states, time = tt, check.names = FALSE)
  
  simDRes[[obs + 1]] <- df_D
  simERes[[obs + 1]] <- df_E
}

head(simDRes[["j_0"]])
head(simERes[["j_0"]])

# save(simDRes, simERes, file = "simulation_results.RData")
Dij_sim <- simDRes[["j_0"]]
plot(Dij_sim$time, Dij_sim$i_0, xlab = 't', 
     ylim = c(0, max(Dij_sim$i_0)*1.15),
     ylab = expression('D'[0]^0*(t))) # TEST Plot D_0^0 SIM

#############################################################################################
## THEORY (ODE) SOLVER SECTION
#############################################################################################

compute_ED_matrices <- function(obs, beta, k, mu, sigma, maxTime, dt) {
  # Define the ODE system for E and D
  ode_system <- function(t, y, parms) {
    E <- y[1:(k + 1)]  # First (k+1) entries correspond to E
    D <- y[(k + 2):(2 * (k + 1))]  # Next (k+1) entries correspond to D
    dE_dt <- rep(0, k + 1)
    dD_dt <- rep(0, k + 1)
    
    # For maximum state
    dE_dt[k + 1] <- mu - (mu + sigma) * E[k + 1]
    dD_dt[k + 1] <- - (mu + sigma) * D[k + 1]
    
    # For all other states
    for (jj in 1:k) {
      dE_dt[jj] <- mu - (mu + sigma + (k + 1 - jj) * beta) * E[jj] +
        (k + 1 - jj) * beta * E[jj + 1] * E[1]
      
      dD_dt[jj] <- - (mu + sigma + (k + 1 - jj) * beta) * D[jj] +
        (k + 1 - jj) * beta * D[jj + 1] * E[1] +
        (k + 1 - jj) * beta * E[jj + 1] * D[1]
    }
    
    return(list(c(dE_dt, dD_dt)))
  }
  
  # Initial conditions for E and D matrices
  initial_conditions <- c(rep(1, k + 1), rep(0, k + 1))  # E values first, then D values
  initial_conditions[k + 1 + obs + 1] <- sigma  # Set sigma for the observed state
  
  # Time steps
  times <- seq(0, maxTime, by = dt)
  
  # Solve the ODE system using ode from deSolve package
  ode_result <- ode(y = initial_conditions, times = times, func = ode_system, parms = NULL)
  
  # Extract E and D matrices from the ODE result
  E_matrix <- ode_result[, 2:(k + 2)]  # First (k+1) columns after time column for E
  D_matrix <- ode_result[, (k + 3):(2 * (k + 1) + 1)]  # Next (k+1) columns for D
  
  # Return the computed matrices as a list
  return(list(E_matrix = E_matrix, D_matrix = D_matrix, time = times))
}

# Initialize list to store D matrices for each observed individual state
theoryDRes <- vector("list", length = k + 1)  # List to store D results for each obs
theoryERes <- vector("list", length = k + 1)  # List to store E results for each obs

names(theoryDRes) <- paste0("j_", 0:k)
names(theoryERes) <- paste0("j_", 0:k)

# Loop through observed individual states and store the corresponding D matrix for each
for (obs in 0:k) {
  cat("Computing D_matrix for observed state", obs, "\n")
  
  # Call the external function to compute the E and D matrices using ODE solver
  matrices <- compute_ED_matrices(obs, beta, k, mu, sigma, maxTime, dt)
  
  # Store each calculated D_matrix in the D_all list
  D_matrix <- matrices$D_matrix
  E_matrix <- matrices$E_matrix
  
  theoryDRes[[obs + 1]] <- as.data.frame(D_matrix)
  colnames(theoryDRes[[obs + 1]]) <- paste0("i_", 0:k)
  theoryDRes[[obs + 1]]$time <- matrices$time
  
  theoryERes[[obs + 1]] <- as.data.frame(E_matrix)
  colnames(theoryERes[[obs + 1]]) <- paste0("i_", 0:k)
  theoryERes[[obs + 1]]$time <- matrices$time
  
  
  cat("Finished computing D_matrix for observed state", obs, "\n")
}

head(theoryDRes[["j_0"]])
head(theoryERes[["j_0"]])

# D_0_0 <- theoryDRes[["j_0"]]['i_2']
# D_0_1 <- theoryDRes[["j_1"]]['i_2']
# D_0_2 <- theoryDRes[["j_2"]]['i_2']
# D_0_3 <- theoryDRes[["j_3"]]['i_2']
# D_0_4 <- theoryDRes[["j_4"]]['i_2']
# 
# ss <- (D_0_0 + D_0_1 + D_0_2 + D_0_3 + D_0_4) * dt
# (sum(D_0_4 * dt))/sum(ss)


lines(theoryDRes[["j_0"]]$time, theoryDRes[["j_0"]]$'i_0', lwd=2) # FIT TEST WITH D_0^0 ODE THEORY


# ---- Minimal plotting D_j^i(t) and E_i(t) from results ----

Dij <- function(i_foc, j_obs) {
  j_name <- paste0("j_", j_obs)
  i_name <- paste0("i_", i_foc)
  
  sim    <- simDRes[[j_name]]
  theory <- theoryDRes[[j_name]]
  if (is.null(sim) || is.null(theory)) stop("Missing results for ", j_name)
  if (!(i_name %in% names(sim)) || !(i_name %in% names(theory)))
    stop("Column ", i_name, " not found for ", j_name)
  
  # pull numeric vectors
  t_sim   <- sim$time
  y_sim   <- sim[[i_name]]
  t_theo  <- theory$time
  y_theo  <- theory[[i_name]]
  
  maxY <- max(c(y_sim, y_theo), na.rm = TRUE) * 1.15
  if (!is.finite(maxY)) maxY <- 1
  
  ylab_expr <- substitute(D[j]^i * "(" * t * ")", list(j = j_obs, i = i_foc))
  # ylab_expr <- substitute(D[j]^{"(" * i * "," * k *")"} * "(" * t * ")", list(j = j_obs, i = i_foc, k=k))
  par(mar=c(4,5,4,1)+.1)
  
  # define the box first
  # plot(c(-100, -100),
  #      xlim = c(0, maxTime), ylim = c(0, maxY),
  #      xlab = "t", ylab = ylab_expr)
  
  # overlay theory and sim
  lines(t_theo, y_theo, lwd = 2, col = "gray10")
  points(t_sim,  y_sim,  pch = 1, col = "gray10", lwd=2)
  
  yLeg = maxY / 1.1
  if(maxY == 0){
    yLeg = 0.6
  }
  
  # legend (with the same color/linetype)
  legend(2, yLeg,
         legend = c("ODE", "Simulation"),
         lty    = c(1, NA),
         pch    = c(NA, 1),
         lwd    = c(2, 2),
         col    = c("gray10", "gray10"),
         cex=1.0)
  
  invisible(list(t_sim=t_sim, y_sim=y_sim, t_theo=t_theo, y_theo=y_theo))
}


Ei <- function(i_foc) {
  j_name <- paste0("j_", 0) # Choose any j state since E_i is independent of j
  i_name <- paste0("i_", i_foc)
  
  sim    <- simERes[[j_name]]
  theory <- theoryERes[[j_name]]
  if (!(i_name %in% names(sim)) || !(i_name %in% names(theory)))
    stop("Column ", i_name, " not found")
  
  # pull numeric vectors
  t_sim   <- sim$time
  y_sim   <- sim[[i_name]]
  t_theo  <- theory$time
  y_theo  <- theory[[i_name]]
  
  maxY <- max(c(y_sim, y_theo), na.rm = TRUE) * 1.15
  if (!is.finite(maxY)) maxY <- 1
  
  ylab_expr <- substitute(E[i] * "(" * t * ")", list(i = i_foc))
  par(mar=c(4,5,4,1)+.1)
  
  # # define the box first
  # plot(c(-100, -100),
  #      xlim = c(0, maxTime), ylim = c(0, maxY),
  #      xlab = "t", ylab = ylab_expr)
  
  # overlay theory and sim
  lines(t_theo, y_theo, lwd = 2, col = "gray10")
  points(t_sim,  y_sim,  pch = 1, col = "gray10", lwd=2)
  
  yLeg = maxY / 1.1
  if(maxY == 0){
    yLeg = 0.6
  }
  
  # legend (with the same color/linetype)
  legend(2, yLeg,
         legend = c("ODE", "Simulation"),
         lty    = c(1, NA),
         pch    = c(NA, 1),
         lwd    = c(2, 2),
         col    = c("gray10", "gray10"),
         cex=1.0)
  
  invisible(list(t_sim=t_sim, y_sim=y_sim, t_theo=t_theo, y_theo=y_theo))
}

i_foc = 4
j_obs = 0
Ei(i_foc = i_foc) # E_{i}(t)
Dij(i_foc = i_foc, j_obs = j_obs) # D_{i}^{j}(t)


# ----------------- PLOT ALL E_i in one plot --------------- #
plot_all_Ei <- function(i_vals = 0:4, j_val = 0, k_val = k, maxTime = NULL) {
  
  j_name <- paste0("j_", j_val)
  sim    <- simERes[[j_name]]
  theory <- theoryERes[[j_name]]
  
  if (is.null(sim) || is.null(theory))
    stop("No data for ", j_name)
  
  # Determine limits
  y_max <- -Inf; y_min <- Inf; x_max <- -Inf
  for (i in i_vals) {
    i_name <- paste0("i_", i)
    if (!(i_name %in% names(sim)) || !(i_name %in% names(theory))) next
    y_sim   <- sim[[i_name]]
    y_theo  <- theory[[i_name]]
    y_max   <- max(y_max, y_sim, y_theo, na.rm = TRUE)
    y_min   <- min(y_min, y_sim, y_theo, na.rm = TRUE)
    x_max   <- max(x_max, sim$time, theory$time, na.rm = TRUE)
  }
  
  if (!is.finite(y_max)) y_max <- 1
  if (y_min > 0) y_min <- 0
  if (is.null(maxTime)) maxTime <- x_max
  
  ylab_expr <- expression(E["(" * i * "," * k * ")"](t))
  
  # par(mar = c(4,5,4,1)+0.1)
  plot(NA, NA,
       xlim = c(0, maxTime), ylim = c(y_min, y_max * 1.15),
       xlab = "t", ylab = ylab_expr,  cex.lab = 1.2,
       main = "")
  
  # Point shapes and line types
  pch_vec <- c(1, 2, 3, 4, 5)
  lty_vec <- rep(1, length(i_vals))
  col_all <- "gray12"
  
  legend_labels <- c()
  legend_pch    <- c()
  legend_lty    <- c()
  legend_lwd    <- c()
  
  for (idx in seq_along(i_vals)) {
    i <- i_vals[idx]
    i_name <- paste0("i_", i)
    if (!(i_name %in% names(sim)) || !(i_name %in% names(theory))) {
      warning("Column ", i_name, " not found – skipping i = ", i)
      next
    }
    
    lines(theory$time, theory[[i_name]],
          lty = lty_vec[idx], lwd = 2, col = col_all)
    points(sim$time, sim[[i_name]],
           pch = pch_vec[idx], col = col_all, cex = 0.9, lwd=1.0)
    
    leg_expr <- substitute(E["(" *i * "," * k_val * ")"](t), list(i = i, k_val = k_val))
    legend_labels <- c(legend_labels, leg_expr)
    legend_pch    <- c(legend_pch, pch_vec[idx])
    legend_lty    <- c(legend_lty, lty_vec[idx])
    legend_lwd    <- c(legend_lwd, 2)
  }
  
  maxY <- max(c(y_sim, y_theo), na.rm = TRUE) * 1.15
  if (!is.finite(maxY)) maxY <- 1
  
  yLeg = maxY / 1.1
  if(maxY == 0){
    yLeg = 0.6
  }
  
  legend(2, yLeg,
         legend = legend_labels,
         lty    = legend_lty,
         pch    = legend_pch,
         lwd    = legend_lwd,
         col    = col_all,
         cex    = 1.0,
         # bg     = "white"
         )
  
  invisible(list(sim = sim, theory = theory))
}

plot_all_Ei()


# -------------------- PLOT ALL D's in one plot --------------------- #
plotDijRow <- function(j_obs, k = 4, maxTime = NULL, showLegend = TRUE) {
  
  j_name <- paste0("j_", j_obs)
  sim    <- simDRes[[j_name]]
  theory <- theoryDRes[[j_name]]
  if (is.null(sim) || is.null(theory)) stop("Missing results for ", j_name)
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  
  par(mfrow    = c(1, k + 1),
      oma      = c(0, 0, 0, 0),
      mgp      = c(2.6, 0.7, 0),
      cex.lab  = 1.6,
      cex.axis = 1.3,
      cex.main = 1.5)
  
  for (i_foc in 0:k) {
    
    # widen left margin ONLY on the first panel so the y-label fits
    if (i_foc == 0) {
      par(mar = c(4.0, 5.1, 2.0, 0.0) + 0.1)
    } else {
      par(mar = c(4.0, 3.0, 2.0, 0.0) + 0.1)
    }
    
    i_name <- paste0("i_", i_foc)
    if (!(i_name %in% names(sim)) || !(i_name %in% names(theory))) {
      plot.new(); next
    }
    
    t_sim  <- sim$time;    y_sim  <- sim[[i_name]]
    t_theo <- theory$time; y_theo <- theory[[i_name]]
    
    maxY <- max(c(y_sim, y_theo), na.rm = TRUE) * 1.15
    if (!is.finite(maxY) || maxY <= 0) maxY <- 1
    mT <- if (is.null(maxTime)) max(c(t_sim, t_theo), na.rm = TRUE) else maxTime
    
    # panel identifier — was the title, now goes into the legend
    panel_label <- substitute(D[j]^i * "(" * t * ")",
                              list(j = j_obs, i = i_foc))
    
    ylab_expr <- if (i_foc == 0)
      expression(D[j]^i * "(" * t * ")")
    else ""
    
    plot(NA, NA,
         xlim = c(0, mT), ylim = c(0, maxY),
         xlab = "t", ylab = ylab_expr)   # no main = ... anymore
    
    lines(t_theo, y_theo, lwd = 2, col = "gray10")
    points(t_sim, y_sim,  pch = 1, col = "gray10", lwd = 2)
    
    yLeg <- if (maxY == 0) 0.6 else maxY / 1.1
    
    if (showLegend) {
      legend(1.5, yLeg,
             title  = as.expression(panel_label),
             legend = c("ODE", "Simulation"),
             lty    = c(1, NA),
             pch    = c(NA, 1),
             lwd    = c(2, 2),
             col    = c("gray10", "gray10"),
             cex    = 1.1,
             bty    = "n")
    } else {
      legend(1.5, yLeg,
             legend = as.expression(panel_label), 
             bty = "n", cex = 1.1)
    }
  }
  invisible(NULL)
}

j_obs <- 4 # Adjust j_obs, that is j_obs <- 0, 1, ..k
plotDijRow(j_obs = j_obs, k = k)
