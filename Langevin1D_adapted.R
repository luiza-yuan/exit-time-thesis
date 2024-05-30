### Langevin1D_adapted function ###
# This is documentation for the function Langevin1D_adapted. It was was adapted from the C++ sourcecode for the Langevin1D function from Rinn's Langevin Package.  
# There are two main modifications in the Langevin1D_adapted function: (1) Rinn's Langevin::Langevin1D function takes 'steps' as an argument, which takes a vector giving the τ steps to calculate the conditional moments (i.e. D1, D2, D4 coefficients). In other words, 'steps' is a vector of integers representing different lag sizes to be used when calculating the differences (increments) in the timeseries data. In the Langevin1D_adapted function, this is modified to the 'ntau' (τ) argument. Instead of specifying a vector, the user specifies the maximum lag value to calculate the conditional moments (i.e. D1, D2, D4 coefficients). This is similar to the 'ntau' argument in Carpenter's DDbins function (also used to estimate D1 and D2 coefficients) where the user specifies the maximum lag value. (2) the bw_sd argument is added to the Langevin1D_adapted function. The adapted function outputs both the raw estimates (identical to estimates using Rinn's Langevin::Langevin1D function) as well as smoothed estimates. The bandwidth is applied when obtaining the smoothed estimates for D1 and D2

filepath_base <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(filepath_base)

# Load necessary packages
library("Langevin")
library("plotrix")
source(file.path(filepath_base, "helper_functions_Rinn.R"))

## Demo of Langevin1D function from Rinn's Langevin Package
# Set necessary parameters 
N = 1000
sf = 10
noise_iter = 17
bins = 50
steps = 1:3

# Generate timeseries to estimate D1 and D2 coefficients on
Ux = do.call(generate_Langevin, utils::modifyList(D, list(
  N = N,
  sf = sf,
  noise_iter = noise_iter
)))

ts <- timeseries1D(N = N, d13 = -1, d11 = 1, d10 = 0.3, sf = sf)

# Estimate D1 and D2 coefficients using Langevin::Langevin1D function
est_Rinn <- Langevin::Langevin1D(data = Ux, bins = bins, steps = steps, bin_min = 100)
plot(est_Rinn$mean_bin, est_Rinn$D1)

## Langevin1D function adapted from C++ source code (can be found on https://github.com/cran/Langevin/blob/master/src/Langevin1D.cpp)
Langevin1D_adapted <- function(Ux, 
                               bins, 
                               ntau, 
                               sf, 
                               bw_sd,
                               bin_min = 50
) {
  # if (!exists("omp_get_num_procs")) {
  #   stop("OpenMP is not available.")
  # }
  # 
  # if (!is.numeric(reqThreads)) {
  #   stop("reqThreads should be numeric")
  # }
  
  # Create vector of different lags to be used when calculating the differences (increments) in the data
  steps <- 1:ntau 
  
  # Specify bandwidth used for obtaining smoothed D1 and D2 estimates
  bw <- bw_sd * sd(Ux)
  
  # haveCores <- as.integer(Sys.getenv("OMP_NUM_THREADS", 1))
  # if (reqThreads <= 0 || reqThreads > haveCores) {
  #   reqThreads <- haveCores
  # }
  
  U <- seq(min(Ux, na.rm = TRUE), max(Ux, na.rm = TRUE), length.out = bins + 1)
  nsteps <- length(steps)
  
  M1 <- matrix(data = NA, nrow = bins, ncol = nsteps)
  eM1 <- matrix(data = NA, nrow = bins, ncol = nsteps)
  M2 <- matrix(data = NA, nrow = bins, ncol = nsteps)
  eM2 <- matrix(data = NA, nrow = bins, ncol = nsteps)
  M4 <- matrix(data = NA, nrow = bins, ncol = nsteps)
  D1 <- numeric(bins)
  eD1 <- numeric(bins)
  D2 <- numeric(bins)
  eD2 <- numeric(bins)
  D4 <- numeric(bins)
  dens <- numeric(bins)
  mean_bin <- numeric(bins)
  
  for (i in 1:bins) {
    sum_m1 <- numeric(nsteps)
    sum_m2 <- numeric(nsteps)
    sum_m4 <- numeric(nsteps)
    len_step <- numeric(nsteps)
    len_bin <- 0
    
    # For each data point and each lag size, calculate the increment as the difference between the data point at the current position and the data point at the position offset by the lag size (for all lag sizes up to the maximum timelag specified by ntau)
    
    #Iterate over each data point (n)
    for (n in 1:(length(Ux) - max(steps))) {
      if (Ux[n] >= U[i] && Ux[n] < U[i + 1] && is.finite(Ux[n])) {
        # Iterate over each time lag (s)
        for (s in 1:nsteps) {
          # For each time lag, the increment (inc) is calculated as the difference between the data point at position n + steps[s] and the data point at position n. Then, the increments are used to accumulate the first moment (sum_m1), second moment (sum_m2), and fourth moment (sum_m4).
          if (is.finite(Ux[n + steps[s]])) {
            inc <- Ux[n + steps[s]] - Ux[n]
            sum_m1[s] <- sum_m1[s] + inc
            sum_m2[s] <- sum_m2[s] + inc^2
            sum_m4[s] <- sum_m4[s] + inc^4
            len_step[s] <- len_step[s] + 1
          }
        }
        mean_bin[i] <- mean_bin[i] + Ux[n]
        len_bin <- len_bin + 1
      }
    }
    
    mean_bin[i] <- mean_bin[i] / len_bin
    dens[i] <- max(len_step)
    
    if (len_bin >= bin_min) {
      # Average each conditional moment by the number of valid increments to obtain the average moments
      M1[i, ] <- sum_m1 / len_step
      M2[i, ] <- sum_m2 / len_step
      M4[i, ] <- sum_m4 / len_step
      
      # Calculate errors for the first and second conditional moments
      eM1[i, ] <- sqrt((M2[i, ] - M1[i, ]^2) / len_step)
      eM2[i, ] <- sqrt((M4[i, ] - M2[i, ]^2) / len_step)
      
      # Use linear regression to estimate drift coefficient (first conditional moment)
      coef <- lm(M1[i, ] ~ steps, weights = 1 / eM1[i, ])$coefficients
      D1[i] <- sf * coef[2]
      
      # Use linear regression to estimate diffusion coefficient (second conditional moment)
      y <- M2[i, ] - (coef[2] * steps)^2
      coef <- lm(y ~ steps, weights = 1 / eM2[i, ])$coefficients
      D2[i] <- sf * coef[2] / 2
      
      # Use linear regression to estimate fourth conditional moment
      coef <- lm(M4[i, ] ~ steps)$coefficients
      D4[i] <- sf * coef[2] / 24
      
      eD1[i] <- sqrt((2 * sf * D2[i] - D1[i]^2) / dens[i])
      eD2[i] <- sqrt((2 * sf * D4[i] - D2[i]^2) / dens[i])
    }
  }
  
  # Added: Fill values that couldn't be estimated with NA's instead of 0's
  replace_near_zero_with_na <- function(vec, tol) {
    vec[abs(vec) < tol] <- NA
    return(vec)
  }
  
  tolerance <- 1e-10 
  
  D1 <- replace_near_zero_with_na(D1, tolerance)
  D2 <- replace_near_zero_with_na(D2, tolerance)
  D4 <- replace_near_zero_with_na(D4, tolerance)
  eD1 <- replace_near_zero_with_na(eD1, tolerance)
  eD2 <- replace_near_zero_with_na(eD2, tolerance)
  
  # Added: Smooth D1 and D2 estimates (not in original Langevin::Langevin1D)
  # Remove NA's (where D1 and D2 couldn't be estimated) before smoothing
  D1_idx = which(!is.na(D1))  
  D2_idx = which(!is.na(D2))
  idx = sort(
    dplyr::intersect(D1_idx,D2_idx)
  )
  D1_no_na = D1[idx]
  D2_no_na = D2[idx]
  
  # X01 = #Xvar lagged by tau???
  D1s = ksmooth(x= mean_bin,y= D1,kernel= 'normal',bandwidth= bw,
                x.points=mean_bin)
  D2s = ksmooth(x= mean_bin ,y= D2,kernel= 'normal',bandwidth=bw,
                x.points=mean_bin)
  
  ret <- list(
    D1 = D1, # raw estimates (same output as from Langevin::Langevin1D function)
    D1s = D1s, # smoothed estimates 
    eD1 = eD1,
    D2 = D2, # raw estimates (same output as from Langevin::Langevin1D function)
    D2s = D2s, # smoothed estimates 
    eD2 = eD2,
    D4 = D4,
    mean_bin = mean_bin, 
    density = dens,
    M1 = M1,
    eM1 = eM1,
    M2 = M2,
    eM2 = eM2,
    M4 = M4,
    U = U
  )
  
  # class(ret) <- "Langevin"
  return(ret)
}

## Demo of Langevin1D_adapted using the same parameter values as the demo of Langevin::Langevin1D above

# Set necessary parameters (same as used in demo of Langevin::Langevin1D)
N = 1000
sf = 10
noise_iter = 17
bins = 50

# Add necessary parameters for Langvin1D_adapted
ntau = 3 # note: the steps argument in Langevin::Langevin1D is 1:ntau in Langevin1D_adapted function
bw_sd = 0.3 

est_adapted <- Langevin1D_adapted(Ux = Ux, bins = bins, ntau = ntau, sf = sf, bw_sd = bw_sd, bin_min = 100)

## Plot D1 and D2 estimates from Langevin::Langevin1D function and Langevin1D_adapted
# The points should exactly match, showing that the two functions produce exactly the same estimates!
plot(est_Rinn$mean_bin, est_Rinn$D1)
points(est_adapted$mean_bin, est_adapted$D1, col = "red", pch = 3)