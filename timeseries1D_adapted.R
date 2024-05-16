library("Langevin")
# Load library used for plotting
library("plotrix")

# Original timeseries1D function from Langevin package
library("Langevin")
N = 500
sf <- 10
set.seed(4711)
Ux <- timeseries1D(N = N, startpoint = 0, d11 = 1, d13 = -1, d22 = 0, d20 = 1, sf = sf, dt = 1/sf)

# Timeseries function translated from C++ sourcecode to R code
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
  
  # define initial noise term (added)
  # prev_noise <- rnorm(1, 0, sqrt(2))
  
  for (i in 1:N) {
    # integrate m steps and save only mth step
    for (j in 1:m) {
      # get a single gaussian random number
      gamma <- rnorm(1, 0, sqrt(2))
      
      # correlate noise (added)
      # "whereas white noise can be produced by randomly choosing each sample independently, Brown noise can be produced by adding a random offset to each sample to obtain the next one" 
      # corr_gamma <- phi * prevNoise + ...??
      # add phi argument to function 
      # change gamma below to corr_gamma? 
      
      # update the prev_noise
      # prev_noise <- corr_gamma
      
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
Ux_adapted <- timeseries1D_adapted(N = N, startpoint = 0, d10 = 0, d11 = 1, d12 = 0, d13 = -1, d22 = 0, d21 = 0, d20 = 1, sf = sf, dt = 1/sf)

# Check if the adapted timeseries1D function produces the same timeseries as the original timeseries1D function
round(Ux, 2) == round(Ux_adapted, 2)

plot(Ux, type = "l", col = "blue", lty = 1, xlab = "Time", ylab = "x", main = "Time Series Comparison")
lines(Ux_adapted, col = "red", lty = 2)


