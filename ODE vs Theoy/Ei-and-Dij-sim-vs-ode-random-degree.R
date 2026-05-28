library(deSolve)
library(splines)

cat('Current working directory', getwd())
setwd(getwd())

# ----------------------------
# PARAMETERS (decouple k_foc and lambda!)
# ----------------------------
# Global Param

focal_k  <- 4        # <-- Chosen latent/focal degree for transitions
lambda   <- 4        # <-- Poisson mean for newborn mixing; can be != focal_k
Kmax     <- 25       # truncate mixing weights (ensure tail mass negligible)
maxTime  <- 3.0

# ODE Param
beta     <- 1.5
mu       <- 0.5
sigma    <- 0.5
dt       <- 0.01

# SIM Param
hh       <- 0.05
nn       <- maxTime / hh
ttStep   <- (0:nn) * hh
tt       <- (0:(nn - 1)) * hh

# ----------------------------
# READ SIM DATA (Poisson file, condition on focal_k)
# ----------------------------
readTable <- function(i_foc, k_foc) {
  dat <- read.table("Simulated data/phylo-epi-sim-data-Poisson.csv", sep = ",", header = TRUE)
  dat <- dat[dat$i_foc == i_foc & dat$k_foc == k_foc, ]
  dat$t_obs <- dat$t_obs - dat$t_foc
  N <- nrow(dat)
  list(dat = dat, N = N, i_foc = i_foc, k_foc = k_foc)
}

# ----------------------------
# SIM FRACTIONS
# ----------------------------
solveSim <- function(df, N, obsType) {
  if (N <= 0) {
    return(list(frac_D = rep(0, nn), frac_E = rep(0, nn)))
  }
  frac_D <- frac_E <- numeric(nn)
  for (ii in 1:nn) {
    in_bin <- (df$t_obs >= ttStep[ii]) & (df$t_obs < ttStep[ii + 1])
    res_D  <- sum(in_bin & (df$j_obs == obsType))
    res_E  <- sum(df$t_obs >= ttStep[ii])
    frac_D[ii] <- res_D / (N * hh)
    frac_E[ii] <- res_E / N
  }
  list(frac_D = frac_D, frac_E = frac_E)
}

# ----------------------------
# BUILD simDRes / simERes exactly like fixed case,
# but with i = 0..focal_k
# ----------------------------

# number of i-states for the chosen focal degree, 
# j-states can be more than focal_k as newborn degree is indepedent of parent.
# For convenience we set it to focal_k

k <- focal_k

simDRes <- vector("list", length = k + 1)
simERes <- vector("list", length = k + 1)
names(simDRes) <- names(simERes) <- paste0("j_", 0:k)

for (obs in 0:k) {
  D_states <- vector("list", length = k + 1)
  E_states <- vector("list", length = k + 1)
  for (foc in 0:k) {
    si  <- readTable(foc, focal_k)
    out <- solveSim(si$dat, si$N, obs)
    D_states[[foc + 1]] <- out$frac_D
    E_states[[foc + 1]] <- out$frac_E
  }
  names(D_states) <- names(E_states) <- paste0("i_", 0:k)
  simDRes[[obs + 1]] <- data.frame(D_states, time = tt, check.names = FALSE)
  simERes[[obs + 1]] <- data.frame(E_states, time = tt, check.names = FALSE)
}

head(simDRes[["j_0"]])
head(simERes[["j_0"]])

# save(simDRes, simERes, file = "simulation_results.RData")
Dij_sim <- simDRes[["j_1"]]
plot(Dij_sim$time, Dij_sim$i_0, ylim = c(0, max(Dij_sim$i_0)*1.15),  
     xlab = 't', ylab = expression('D'[1]^0*(t))) # Plot D_{1}^{0} SIM condition of focal k

# ----------------------------
# COUPLED RANDOM-DEGREE ODE (mix newborns via w, transitions at each k)
# ----------------------------
compute_ED_coupled_one_obs <- function(obs, beta, mu, sigma, w, maxTime, dt){
  if (is.null(names(w))) names(w) <- as.character(seq(0, length(w)-1))
  Kmax <- length(w) - 1L
  
  # layout: for each k (0..Kmax): E_{0..k}, then D_{0..k}
  block_offsets <- cumsum(c(0, sapply(0:Kmax, function(k) 2*(k+1))))[-(Kmax+2)]
  idx_E <- function(k,i) block_offsets[k+1] + i + 1
  idx_D <- function(k,i) block_offsets[k+1] + (k+1) + i + 1
  total_len <- tail(block_offsets,1) + 2*(Kmax+1)
  
  ode_system <- function(t, y, parms){
    dy <- numeric(total_len)
    # newborn mixtures
    Ehat0  <- sum(w * sapply(0:Kmax, function(kk) y[idx_E(kk,0)]))
    Dhat0j <- sum(w * sapply(0:Kmax, function(kk) y[idx_D(kk,0)]))
    for (kk in 0:Kmax){
      for (ii in 0:kk){
        avail <- kk - ii
        loss  <- mu + sigma + avail*beta
        Ei_k  <- y[idx_E(kk,ii)]
        Di_k  <- y[idx_D(kk,ii)]
        dEi_k <- mu - loss*Ei_k
        dDi_k <- -loss*Di_k
        if (ii < kk){
          Ei1_k <- y[idx_E(kk,ii+1)]
          Di1_k <- y[idx_D(kk,ii+1)]
          dEi_k <- dEi_k + avail*beta*Ehat0*Ei1_k
          dDi_k <- dDi_k + avail*beta*(Ehat0*Di1_k + Ei1_k*Dhat0j)
        }
        dy[idx_E(kk,ii)] <- dEi_k
        dy[idx_D(kk,ii)] <- dDi_k
      }
    }
    list(dy)
  }
  
  # ICs: E_{i,k}(0)=1 ; D_{i,k}(0)=sigma * 1{i=obs}
  y0 <- numeric(total_len)
  for (kk in 0:Kmax) for (ii in 0:kk){
    y0[idx_E(kk,ii)] <- 1
    y0[idx_D(kk,ii)] <- if (ii == obs) sigma else 0
  }
  
  times <- seq(0, maxTime, by = dt)
  sol   <- ode(y=y0, times=times, func=ode_system, parms=NULL, method="ode45")
  
  # unpack into per-k data.frames with columns i_0..i_k and time
  E_by_k <- vector("list", Kmax+1)
  D_by_k <- vector("list", Kmax+1)
  for (kk in 0:Kmax){
    Ek <- Dk <- matrix(NA_real_, nrow=nrow(sol), ncol=kk+1)
    for (ii in 0:kk){
      Ek[, ii+1] <- sol[, idx_E(kk,ii)+1]  # +1 for time column
      Dk[, ii+1] <- sol[, idx_D(kk,ii)+1]
    }
    colnames(Ek) <- colnames(Dk) <- paste0("i_", 0:kk)
    E_by_k[[kk+1]] <- data.frame(Ek, time=times, check.names = FALSE)
    D_by_k[[kk+1]] <- data.frame(Dk, time=times, check.names = FALSE)
  }
  list(time=times, E_by_k=E_by_k, D_by_k=D_by_k)
}

# newborn weights w_k from Poisson(lambda), truncated at Kmax
w <- dpois(0:Kmax, lambda); w <- w / sum(w)

# ----------------------------
# THEORY CONTAINERS (per j, extract k=focal_k)
# ----------------------------
theoryDRes <- vector("list", length = k + 1)
theoryERes <- vector("list", length = k + 1)
names(theoryDRes) <- names(theoryERes) <- paste0("j_", 0:k)

for (j_obs in 0:k) { # Grab observed state upto focal_k (again we only set j_obs <= focal_k for plotting and code convenience)
  res_cpl <- compute_ED_coupled_one_obs(j_obs, beta, mu, sigma, w, maxTime, dt)
  # extract only the focal degree block (k = focal_k)
  D_lat <- res_cpl$D_by_k[[focal_k + 1]]
  E_lat <- res_cpl$E_by_k[[focal_k + 1]]
  # sure columns include exactly i_0..i_k (k=focal_k) + time
  theoryDRes[[j_obs + 1]] <- D_lat
  theoryERes[[j_obs + 1]] <- E_lat
}

Dij_ode <- theoryDRes[["j_1"]]
lines(Dij_ode$time, Dij_ode$i_0) # Plot D_1^0 ODE

# ----------------------------
# PLOTTER
# ----------------------------
Dij_random <- function(i_foc, j_obs) {
  j_name <- paste0("j_", j_obs)
  i_name <- paste0("i_", i_foc)
  sim    <- simDRes[[j_name]]
  theory <- theoryDRes[[j_name]]
  if (is.null(sim) || is.null(theory)) stop("Missing results for ", j_name)
  if (!(i_name %in% names(sim)) || !(i_name %in% names(theory)))
    stop("Column ", i_name, " not found for ", j_name)
  t_sim  <- sim$time;   y_sim  <- sim[[i_name]]
  t_theo <- theory$time; y_theo <- theory[[i_name]]
  maxY <- max(c(y_sim, y_theo), na.rm=TRUE) * 1.1; if (!is.finite(maxY)) maxY <- 1
  plot(c(-100,-100), xlim=c(0,maxTime), ylim=c(0,maxY),
       xlab="t", ylab=bquote(D[.(j_obs)]^{.(i_foc)}(t)))
  lines(t_theo, y_theo, lwd=2, col="gray54")
  points(t_sim,  y_sim,  pch=19, col="gray54")
  yLeg = maxY / 1.1
  if(maxY == 0){
    yLeg = 0.6
  }
  legend(2, yLeg, legend=c("ODE","Simulation"),
         lty=c(1,NA), pch=c(NA,19), lwd=c(2,NA),
         col=c("gray54","gray54"), cex=0.8)
}

Ei_random <- function(i_foc) {
  j_name <- paste0("j_", 0) # Use any j_obs as E_i is independent of j
  i_name <- paste0("i_", i_foc)
  sim    <- simERes[[j_name]]
  theory <- theoryERes[[j_name]]
  if (is.null(sim) || is.null(theory)) stop("Missing results for ", j_name)
  if (!(i_name %in% names(sim)) || !(i_name %in% names(theory)))
    stop("Column ", i_name, " not found for ", j_name)
  t_sim  <- sim$time;   y_sim  <- sim[[i_name]]
  t_theo <- theory$time; y_theo <- theory[[i_name]]
  maxY <- max(c(y_sim, y_theo), na.rm=TRUE) * 1.1; if (!is.finite(maxY)) maxY <- 1
  plot(c(-100,-100), xlim=c(0,maxTime), ylim=c(0,maxY),
       xlab="t", ylab=bquote(D[.(j_obs)]^{.(i_foc)}(t)))
  lines(t_theo, y_theo, lwd=2, col="gray54")
  points(t_sim,  y_sim,  pch=19, col="gray54")
  yLeg = maxY / 1.1
  if(maxY == 0){
    yLeg = 0.6
  }
  legend(2, yLeg, legend=c("ODE","Simulation"),
         lty=c(1,NA), pch=c(NA,19), lwd=c(2,NA),
         col=c("gray54","gray54"), cex=0.8)
}

# ----------------------------
# CALL Plotter
# ----------------------------
Dij_random(2, 0) # D_{i}^{j}(t)
Dij_random(0, 2) # D_{i}^{j}(t)


Ei_random(1)
Ei_random(2)
Ei_random(3)

# ------------ PLOT ALL D's in one plot -------------- #

plotDijRow_random <- function(j_obs, k = 4, maxTime = NULL, showLegend = TRUE) {
  
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
    
    # widen left margin only on the first panel so y-label fits
    if (i_foc == 0) {
      par(mar = c(4.0, 5.2, 2.0, 0.0) + 0.1)
    } else {
      par(mar = c(4.0, 3.0, 2.0, 0.0) + 0.1)
    }
    
    i_name <- paste0("i_", i_foc)
    if (!(i_name %in% names(sim)) || !(i_name %in% names(theory))) {
      plot.new(); next
    }
    
    t_sim  <- sim$time;    y_sim  <- sim[[i_name]]
    t_theo <- theory$time; y_theo <- theory[[i_name]]
    
    maxY <- max(c(y_sim, y_theo), na.rm = TRUE) * 1.1
    if (!is.finite(maxY) || maxY <= 0) maxY <- 1
    mT <- if (is.null(maxTime)) max(c(t_sim, t_theo), na.rm = TRUE) else maxTime
    
    panel_label <- bquote(D[.(j_obs)]^.(i_foc) * (t))
    
    ylab_expr <- if (i_foc == 0)
      expression(D[j]^i * "(" * t * ")")
    else ""
    
    plot(c(-100, -100),
         xlim = c(0, mT), ylim = c(0, maxY),
         xlab = "t", ylab = ylab_expr)
    
    lines(t_theo, y_theo, lwd = 2, col = "gray10")
    points(t_sim, y_sim,  pch = 1, col = "gray10", lwd=2)
    
    yLeg <- if (maxY == 0) 0.6 else maxY / 1.1
    
    if (showLegend) {
      legend(1.5, yLeg,
             title  = as.expression(panel_label),
             legend = c("ODE", "Simulation"),
             lty    = c(1, NA),
             pch    = c(NA, 1),
             lwd    = c(2, 2),
             col    = c("gray10", "gray10"),
             cex    = 1.0,
             bty    = "n")
    } else {
      legend(1.5, yLeg,
             legend = as.expression(panel_label),
             bty = "n", cex = 1.2)
    }
  }
  invisible(NULL)
}


j_obs = 4 j_obs <- 4 # Adjust j_obs, that is j_obs <- 0, 1, ..focal_k
plotDijRow_random(j_obs = j_obs, k = k)
