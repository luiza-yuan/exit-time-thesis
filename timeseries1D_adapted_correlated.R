# Load necessary libraries and helper functions
library(bvpSolve)
library(cubature)
library(stats)
library(Langevin)
library(boot)
library(tseries)
library(spgs)
library(dplyr)
library(ggplot2)
library(parallel) 
library(doParallel) 
library(foreach)
library(cowplot)
library(ggnewscale)
library(latex2exp)
library(plotrix)
source("~/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/May/helper_functions_Rinn_May.R")

##########################################
# One dimensional example
##########################################
timeseries1D_adapted_correlated <- function(N, startpoint = 0, d13 = 0, d12 = 0, d11 = -1, d10 = 0, d22 = 0, d21 = 0, d20 = 1, sf = 1000, dt = 0) {
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
      
      # print(c("previous gamma is", gamma_prev1))
      # print(c("previous 2 gamma is", gamma_prev2))
      # print(c("gamma is", gamma))
      
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

### Demo
nr_steps_bif = 3
strength_D2 = 0.2
N_high_res = c(12000)
sf_high_res = c(1000)
N_burned = c(10)
interval = c(250, 125, 100)[2]
timeseries_length = c(60000)
interpol_steps = 50 # c(50, 100, 500),
bins <- 10 # c(10, 50)
ntau = c(2) # c(3, 5, 10),
add_to_x = 0.5 # 0.5
bw_sd = 0.3
noise_iter = c(1:10)

set.seed(4711)
Ux_original <- timeseries1D(N = N_high_res*sf_high_res, d10 = 0, d11 = 3, d12 = 0.25, d13 = -3, d22 = strength_D2, d21 = 0, d20 = strength_D2, sf = sf_high_res)
Ux_original_correlated <- timeseries1D_adapted(N = N_high_res*sf_high_res, d10 = 0, d11 = 3, d12 = 0.25, d13 = -3, d22 = strength_D2, d21 = 0, d20 = strength_D2, sf = sf_high_res)

Ux_original_burned <- window(Ux_original, start = c(N_burned + 1, 1), end = c(12000, 1000))
Ux_original_burned_correlated <- window(Ux_original_correlated, start = c(N_burned + 1, 1), end = c(12000, 1000))

Ux <- downsample_Langevin(Ux = Ux_original_burned, interval = interval, timeseries_length = timeseries_length)
Ux_correlated <- downsample_Langevin(Ux = Ux_original_burned_correlated, interval = interval, timeseries_length = timeseries_length)

### Plot time series and probability density function (PDF)
# timeseries with uncorrelated noise (Markovian)
op <- par(no.readonly = TRUE)
par(mar = c(4.2, 5.4, 0.5, 0.5))
layout(matrix(c(1, 1, 2), 1, 3))
plot(1:length(Ux), Ux, xlab = "t [s]", ylab = "x [a.u.]", t = 'l')
plot(density(Ux), xlab = "x [a.u.]", ylab = "density", t = 'l',main = "")
par(op)

# timeseries with correlated noise (non-Markovian)
op <- par(no.readonly = TRUE)
par(mar = c(4.2, 5.4, 0.5, 0.5))
layout(matrix(c(1, 1, 2), 1, 3))
plot(1:length(Ux_correlated), Ux_correlated, xlab = "t [s]", ylab = "x [a.u.]", t = 'l')
plot(density(Ux_correlated), xlab = "x [a.u.]", ylab = "density", t = 'l',main = "")
par(op)

# round(Ux, 2) == round(Ux_correlated, 2)
# plot(Ux[1:7000], type = "l", col = "blue", lty = 1, xlab = "Time", ylab = "x", main = "Time Series Comparison")
# lines(Ux_correlated[1:7000], col = "red", lty = 2)

# markov tests
## note, it is not violated 
# spgs::markov.test(Ux, type = "lb.test")$p.values
# spgs::markov.test(Ux_correlated, type = "lb.test")$p.values
# spgs::markov.test(Ux, type = "rank.test")$p.values
# spgs::markov.test(Ux_correlated, type = "rank.test")$p.values

### test est_D_Carp
D <- list(d10 = 0, d11 = 3, d12 = 0.25, d13 = -3, d22 = strength_D2, d21 = 0, d20 = strength_D2)
stabs = get_stability(D)

# Get exit times
est_Carp = est_D_Carp(
  Ux,
  sf = frequency(Ux), 
  D,
  stabs = stabs,
  ntau = ntau,
  bins = bins,
  bw_sd = bw_sd,
  interpol_steps = interpol_steps
)

est_Carp = est_D_Carp(
  Ux_correlated,
  sf = frequency(Ux_correlated), 
  D,
  stabs = stabs,
  ntau = ntau,
  bins = bins,
  bw_sd = bw_sd,
  interpol_steps = interpol_steps
)
