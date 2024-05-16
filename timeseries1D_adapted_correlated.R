library("Langevin")
# Load library used for plotting
library("plotrix")

# Set the number of bins to estimate drift and diffusion on 
bins <- 100
# Set the time lags for computing the conditional moments
steps <- c(1:3)
# Set the sampling frequency
sf <- 10

##########################################
# One dimensional example
##########################################

# copy paste
library("Langevin")
set.seed(4711)
Ux <- timeseries1D(N = 400*10, startpoint = 0, d10 = 0, d11 = 1, d12 = 0, d13 = -1, d22 = 0, d21 = 0, d20 = 1, sf = 10, dt = 1/sf)

# bins <- 100
# steps <- c(1:3)
# ests <- Langevin1D(x, bins, steps)

timeseries1D_adapted <- function(N, startpoint, d13, d12, d11, d10, d22, d21, d20, sf, dt) {
  # initialize time series with NA's
  ts <- rep(NA, N)
  ts[1] <- startpoint
  
  # calculate integration time step and related values
  stime <- 1 / sf
  if (stime < dt || dt == 0) {
    dt <- stime
  }
  
  # ratio between sampling time and integration time step
  m <- ceiling(stime / dt)
  dt <- stime / m
  
  # integrate
  x <- ts[1]
  
  # initialize previous gamma terms
  gamma_prev1 <- 0
  gamma_prev2 <- 0
  
  for (i in 1:N) {
    # integrate m steps and save only mth step
    for (j in 1:m) {
      # get a single gaussian random number (original)
      # gamma <- rnorm(1, 0, sqrt(2))
      
      # correlate noise (added)
      # "whereas white noise can be produced by randomly choosing each sample independently, Brown noise can be produced by adding a random offset to each sample to obtain the next one" 
      # corr_gamma <- phi * prevNoise + ...??
      # add phi argument to function 
      # change gamma below to corr_gamma? 
      
      # update the prev_noise
      # prev_noise <- corr_gamma
      
      # correlate gamma (added)
      gamma <- 0.2 * gamma_prev1 + 0.1 * gamma_prev2 + rnorm(1, 0, sqrt(2))
      
      # update previous gamma terms
      gamma_prev2 <- gamma_prev1
      gamma_prev1 <- gamma
      
      # iterate with integration time step dt
      x <- x + (d13 * x^3 + d12 * x^2 + d11 * x + d10) * dt + sqrt((d22 * x^2 + d21 * x + d20) * dt) * gamma
    }
    # save every mth step
    ts[i] <- x
  }
  
  # designate class and tsp attributes for timeseries (ts) object
  class(ts) <- "ts"
  tsp(ts) <- c(1, 1 + (N - 1) / sf, sf)
  
  return(ts)
}

set.seed(4711)
Ux_adapted <- timeseries1D_adapted(N = 400*10, startpoint = 0, d10 = 0, d11 = 1, d12 = 0, d13 = -1, d22 = 0, d21 = 0, d20 = 1, sf = 10, dt = 1/sf)
# Plot time series and probability density function (PDF)
op <- par(no.readonly = TRUE)
par(mar = c(4.2, 5.4, 0.5, 0.5))
layout(matrix(c(1, 1, 2), 1, 3))
plot((1:length(Ux_adapted)/sf), Ux_adapted, xlab = "t [s]", ylab = "x [a.u.]", t = 'l')
# plot(density(Ux), xlab = "x [a.u.]", ylab = "density", t = 'l',main = "")
par(op)

round(Ux, 2) == round(Ux_adapted, 2)

plot(Ux, type = "l", col = "blue", lty = 1, xlab = "Time", ylab = "x", main = "Time Series Comparison")
lines(Ux_adapted, col = "red", lty = 2)

# markov tests
## note, it is not violated 
spgs::markov.test(Ux, type = "lb.test")$p.values
spgs::markov.test(Ux_adapted, type = "lb.test")$p.values
spgs::markov.test(Ux, type = "rank.test")$p.values
spgs::markov.test(Ux_adapted, type = "rank.test")$p.values

# test est_D_Carp
example$D
D_test <- list(d10 = 0, d11 = 1, d12 = 0, d13 = -1, d22 = 0, d21 = 0, d20 = 1)
stabs_test = get_stability(D_test)

# Apply Exit Time Analysis (Carpenter's (2022) approach adapted )
est_Carp = est_D_Carp(
  Ux_adapted,
  sf = 10, 
  D_test,
  stabs = stabs_test,
  # Tstep = 1:length(as.vector(Ux)),
  Tstep = 1:length(as.vector(Ux_adapted)) / 10,
  ntau = c(3),
  bins = 100,
  bw_sd = 0.3,
  interpol_steps = 50
)
