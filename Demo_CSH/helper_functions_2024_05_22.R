# The following helper functions are developed based on the papers:
# (1) Arani, B. M., Carpenter, S. R., Lahti, L., Van Nes, E. H., & Scheffer, M. (2021). Exit time as a measure of ecological resilience. Science, 372(6547), eaay4895.
# (2) Carpenter, S. R., Arani, B. M., Van Nes, E. H., Scheffer, M., & Pace, M. L. (2022). Resilience of phytoplankton dynamics to trophic cascades and nutrient enrichment. Limnology and Oceanography, 67, S258-S265.
# (3) Rinn, P., Lind, P. G., Wächter, M., & Peinke, J. (2016). The Langevin approach: An R package for modeling Markov processes. arXiv preprint arXiv:1603.02036.

# Generate timeseries using Langevin
generate_Langevin <- function(N, # a scalar denoting the length of the time-series to generate
                              sf, # a scalar denoting the sampling frequency
                              noise_iter, # set seed
                              d13, # scalars denoting the coefficients for the drift polynomial (d13 is the coefficient for x^3, d12 is the coefficient for x^2, etc.)
                              d12,
                              d11,
                              d10,
                              d22, # scalars denoting the coefficients for the diffusion polynomial (d22 is the coefficient for x^2, d21 is the coefficient for x, etc.)
                              d21,
                              d20) {
  # Generate the time series
  set.seed(noise_iter)
  Ux <- Langevin::timeseries1D(
    N = N * sf,
    d13 = d13,
    d12 = d12,
    d11 = d11,
    d10 = d10,
    d22 = d22,
    d21 = d21,
    d20 = d20,
    sf = sf,
    dt = 1 / sf # maximal time step of integration (default dt=0 yields dt=1/sf.)
  )
  
  # Cut the timeseries if it went to infinity
  if (any(is.na(Ux) | is.infinite(Ux))) {
    message("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
    return()
  } else {
    # Return timeseries object
    return(Ux)
  }
}


# Function for downsampling
downsample_Langevin <- function(Ux, # Timeseries (high resolution)
                                n_points # Length of desired final timeseries (downsampled)
) {
  # Calculate the interval for downsampling: This represents the number of data points in the original time series that will be represented by a single data point in the downsampled time series.
  interval <- length(Ux) / (n_points) 
  
  # Generate a sequence of indices for downsampling: This creates a vector of evenly spaced indices that will be used to extract data points from the original time series.
  selected_indices <- seq(1, length(Ux), by = interval) 
  
  # Extract the downsampled time series using the selected indices
  downsampled_Ux <- Ux[selected_indices] 
  
  # Calculate the new "sampling frequency" after downsampling 
  new_frequency <- frequency(Ux) / interval 
  
  # Create a new time-series object with the downsampled data and the new sampling frequency
  downsampled_Ux <- ts(downsampled_Ux, frequency = new_frequency) 
  
  # Return the downsampled time series
  return(downsampled_Ux) 
}


# Polynomial (third-order)
poly3 <-
  function(x, a0, a1, a2, a3) {
    a0 + a1 * x + a2 * x ^ 2 + a3 * x ^ 3
  }


# Derive theoretical drift, diffusion, and potential analytically
get_theoretical_D <- function(D, min_x = -5, max_x = 5, interpol_steps) {
  with(D, {
    # Calculate coefficients for the potential function from the drift coefficients
    # Drift = d10 + d11 * x + d12 * x^2 + d13 * x^3 = negative derivative of potential
    # Potential = constant + a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4
    # Coefficients are derived as follows:
    # d13 = a4 * 4
    # d12 = a3 * 3
    # d11 = a2 * 2
    # d10 = a1 * 1
    # The constant term falls away during differentiation
    
    # Calculate potential coefficients
    a4 = d13 / 4
    a3 = d12 / 3
    a2 = d11 / 2
    a1 = d10
    
    # Generate a sequence of x values from min_x to max_x
    x = seq(min_x, max_x, length.out = interpol_steps)
    
    # Create a data frame with theoretical drift, diffusion, and potential values
    theoretical_df = data.frame(
      x = x,
      # Calculate drift (derivative of potential function)
      drift = d10 + d11 * x + d12 * x^2 + d13 * x^3,
      # Calculate diffusion 
      diffusion = d20 + d21 * x + d22 * x^2,
      # Calculate potential 
      potential = -(a1 * x + a2 * x^2 + a3 * x^3 + a4 * x^4)
    ) %>%
      # Add effective potential field (EPF)
      dplyr::mutate(
        EPF = c(NA, get_effective_potential(
          list(x = x, y = drift), 
          list(x = x, y = diffusion), 
          interpol_steps = interpol_steps, 
          xvec = x
        )$EPF)
      ) %>%
      # Reshape the data frame from wide to long format
      tidyr::gather(variable, value, -x)
    
    # Return the data frame with theoretical values
    return(theoretical_df)
  })
}


get_D <- function(nr_steps_bif, # Number of bifurcation steps 
                  scenario, # Bifurcation scenario 
                  type_D2, # Type of D2 (diffusion) parameter 
                  strength_D2) { # Strength of D2 (diffusion) parameter
  
  # Define d13, d12, d11, d10 based on the bifurcation scenario
  if (scenario == "2fps-balanced-deepening") {
    # For "2fps-balanced-deepening" scenario:
    # d13 decreases linearly from -3 to -1
    d13 = rev(seq(-3,-1, length.out = nr_steps_bif)) 
    # d12 is constant at 0
    d12 = 0 
    # d11 decreases linearly from 3 to 1
    d11 = rev(seq(3, 1, length.out = nr_steps_bif)) 
    # d10 is constant at 0
    d10 = 0 
  } else if (scenario == "right-fp-gains-dominance") {
    # For "right-fp-gains-dominance" scenario:
    # d13 is constant at -3
    d13 = -3 
    # d12 increases linearly from 0.05 to 0.25
    d12 = seq(0.05, 0.25, length.out = nr_steps_bif)  
    # d11 is constant at 3
    d11 = 3 
    # d10 is constant at 0
    d10 = 0 
  } 
  
  # Define d22, d21, d20 based on the type of D2 (diffusion) parameter
  if (type_D2 == "constant") {
    # For "constant" type_D2:
    # d22 is 0
    d22 = 0 
    # d21 is 0
    d21 = 0 
    # d20 is 1
    d20 = 1 
  } else if (type_D2 == "quadratic") {
    # For "quadratic" type_D2:
    # d22 is 1
    d22 = 1 
    # d21 is 0
    d21 = 0 
    # d20 is 1
    d20 = 1 
  }
  
  # Create a data frame with all the d parameters
  Ds = data.frame(
    d13 = d13,
    d12 = d12,
    d11 = d11,
    d10 = d10,
    # Multiply d22, d21, d20 by diffusion strength (noise intensity)
    d22 = d22 * strength_D2,
    d21 = d21 * strength_D2,
    d20 = d20 * strength_D2
  )  %>% 
    # Transpose the data frame
    purrr::transpose() %>% 
    # Keep only unique rows
    unique() %>%
    # Repeat the sequence if there's only one row
    rep(ifelse(length(.) == 1, nr_steps_bif, 1))
  
  # Return the data frame
  return(Ds)
}


# Calculate the stability of the theoretical fixed points for the specified drift function
get_stability <- function(D) {
  with(D, {
    # Find fixed points
    # Calculate the roots of the polynomial defined by d10, d11, d12, and d13
    # These roots represent the fixed points of the system
    fps = sort(unique(round(Re(polyroot(
      c(d10, d11, d12, d13)
    )), 10))) 
    
    # Calculate the derivative of the polynomial
    # This derivative will be used to determine the stability of the fixed points
    f_prime = D(expression(d10 + d11 * x + d12 * x ^ 2 + d13 * x ^ 3), "x") 
    
    # Evaluate the derivative at each fixed point
    # The sign of the derivative at a fixed point determines its stability
    stab = plyr::laply(fps, function(x) {
      eval(f_prime)
    }) 
    
    # Recode stability values to "stable" and "unstable"
    # Convert the stability values from logical (TRUE/FALSE) to character ("stable"/"unstable")
    stab_fps = dplyr::recode(as.numeric(stab < 0), "0" = "unstable", "1" = "stable")
    
    # Print the theoretical equilibria and their stability
    print('Theoretical Equilibria', quote = F)
    print(fps, quote = F)
    print(stab_fps, quote = F)
    
    # Return the fixed points and their stability as a list
    return(list(fps = fps, stab_fps = stab_fps))
  })
}


# Estimate drift and diffusion coefficients as well as resilience metrics
# (adapted partly following approach used by Arani et al. (2021) and Carpenter et al. (2022))
est_D_Carp <- function(Ux,
                       sf,
                       D,
                       stabs,
                       ntau,
                       bins,
                       bw_sd,
                       interpol_steps, 
                       add_to_x = .5) {
  
  # Estimate drift and diffusion
  DD = apply_Langevin1D_adapted(
    Ux = Ux,
    bins = bins,
    ntau = ntau,
    sf = sf, 
    bin_min = 50,
    bw_sd = bw_sd
  )
  
  # Create an interpolation vector specifying points on the x-axis. Exclude the outer edges of the data and make sure the interpolation vector falls within the range of the drift and diffusion function by a margin
  # nearest = .1
  xvec = seq(min(DD$D1s$x, na.rm = TRUE) - add_to_x, 
             max(DD$D1s$x, na.rm = TRUE) + add_to_x,
             length.out = interpol_steps)
  
  # Get potential and effective potential 
  Pot = get_potential(
    D1s = DD$D1s,
    D2s = DD$D2s,
    xeq = DD$xeq,
    interpol_steps = interpol_steps,
    xvec = xvec
  )
  
  # Get theoretical drift and diffusion functions
  Theoretical_df = get_theoretical_D(D, min_x = min(xvec), max_x = max(xvec), interpol_steps = interpol_steps) # KCE
  
  # Only if there are two stable equilibria on the outside and on unstable equilibrium on the inside, compute exit time
  # if (identical(stabs$stab_fps, c("stable", "unstable", "stable"))) {
  # Get exit time using theoretical drift and diffusion
  TheoreticalExitTime = get_exit_time(
    Theoretical_df %>% dplyr::filter(variable == "drift") %>% dplyr::rename(y = value) %>% dplyr::arrange(x) %>% dplyr::select(-variable),
    Theoretical_df %>% dplyr::filter(variable == "diffusion")  %>% dplyr::rename(y = value)  %>% dplyr::arrange(x) %>% dplyr::select(-variable),
    xeq = stabs$fps, # theoretical fixed points
    xvec = xvec,
    interpol_steps = interpol_steps
  ) 
  
  # Only if there are two stable equilibria on the outside and on unstable equilibrium on the inside, compute exit time
  if (identical(Pot$stab_xeq, c("stable", "unstable", "stable"))) {
    
    # Exit time calculated with estimated drift and diffusion
    EstimatedExitTime = get_exit_time(
      D1s = DD$D1s, # smoothed D1 estimates
      D2s = DD$D2s, # smoothed D2 estimates
      xeq = DD$xeq, # fixed points estimated from D1
      xvec = xvec, # interpolation vector specifying points on the x-axis
      interpol_steps = interpol_steps
    )
    
    # Create data frame for fixed points (both theoretical and estimated)
    fp_df = rbind(data.frame(xintercept = stabs$fps,
                             color = stabs$stab_fps, source = "Theoretical"),
                  data.frame(xintercept = DD$xeq,
                             color = Pot$stab_xeq, source = "Estimated")) %>%
      mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
    
  } else if(identical(Pot$stab_pot_xeq, c("stable", "unstable", "stable"))) {
    
    # Exit time calculated with estimated drift and diffusion
    EstimatedExitTime = get_exit_time(
      D1s = DD$D1s, # smoothed D1 estimates
      D2s = DD$D2s, # smoothed D2 estimates
      xeq = Pot$pot_xeq, # fixed points estimated from potential
      xvec = xvec, # interpolation vector specifying points on the x-axis
      interpol_steps = interpol_steps
    )
    
    # Create data frame for fixed points (both theoretical and estimated)
    fp_df = rbind(data.frame(xintercept = stabs$fps,
                             color = stabs$stab_fps, source = "Theoretical"),
                  data.frame(xintercept = Pot$pot_xeq,
                             color = Pot$stab_pot_xeq, source = "Estimated")) %>%
      mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
  } else {
    EstimatedExitTime = list(ET_df = data.frame())
    
    # Create data frame for fixed points (both theoretical and estimated)
    fp_df = rbind(data.frame(xintercept = stabs$fps,
                             color = stabs$stab_fps, source = "Theoretical"),
                  data.frame(xintercept = DD$xeq,
                             color = Pot$stab_xeq, source = "Estimated")) %>%
      mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
  }
  
  # Format results (estimated D1, D2, potential, and effective potential)
  
  # Dataframe with estimated drift and diffusion coefficients (smoothed)
  est_df_Rinn = data.frame(
    x = DD$bin.mid,
    D1 = DD$D1s$y,
    D2 = DD$D2s$y
  ) %>% tidyr::gather(variable, value,-x) %>%
    dplyr::mutate(err = NA)
  
  # Dataframe with potential, and effective potential
  Pot_df = data.frame(
    x = Pot$xvec,
    negPF = Pot$negPF,
    EPF = c(NA, Pot$EPF),
    EPF_err = c(NA, Pot$EPF_err),
    drift = Pot$drift,
    diffusion = Pot$diff
  ) %>% dplyr::rename(potential = negPF) %>% tidyr::gather(variable, value,-c(x, EPF_err)) %>%
    dplyr::rename(err = EPF_err) %>% dplyr::mutate(err = ifelse(variable == "EPF", NA, err))
  
  # Format results of theoretical and estimated exit times and weights from stationary probability distribution
  compl_df = rbind(Theoretical_df %>% mutate(source = "Theoretical", err = NA),
                   rbind(
                     # est_df_Carp,
                     Pot_df) %>% mutate(source = "Estimated"),
                   TheoreticalExitTime$ET_df %>% mutate(source = "Theoretical"),
                   EstimatedExitTime$ET_df %>% mutate(source = "Estimated")
  ) %>% mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
  
  # Get resilience metrics from potential function (basin width and potential depth)
  alt_metrics_potential = get_resilience_potential(fp_df = fp_df, compl_df = compl_df)
  
  # Get resilience metrics from probability distribution of system states (variance around dominant mode)
  alt_metrics_prob_dist = get_resilience_prob_dist(Ux = Ux, bw_sd = bw_sd)
  
  return(
    list(
      DD = DD, # D1 and D2 estimates (raw and smoothed)
      est_df_Rinn = est_df_Rinn, # D1 and D2 estimates (smoothed)
      compl_df = compl_df, # Complete dataframe of theoretical and estimated exit times with weights
      fp_df = fp_df, # Dataframe with theoretical and estimated fixed points
      xeq = DD$xeq, # Estimated fixed points from D1 
      alt_metrics_potential = alt_metrics_potential, # Resilience metrics from potential function
      alt_metrics_prob_dist = alt_metrics_prob_dist, # Resilience metrics from probability distribution
      stab_xeq = Pot$stab_xeq, # Estimated fixed points from potential function
      Dratio = Pot$Dratio,
      D1D2 = Pot$D1D2,
      D1D2adj = Pot$D1D2adj,
      theo_meanETl = TheoreticalExitTime$meanETl, # Theoretical mean exit time (left basin)
      theo_meanETr = TheoreticalExitTime$meanETr, # Theoretical mean exit time (right basin)
      theo_atolL = TheoreticalExitTime$atolL, 
      theo_atolR = TheoreticalExitTime$atolR,
      theo_diagn_ETL = TheoreticalExitTime$diagn_ETL, # Diagnostics
      theo_diagn_ETR = TheoreticalExitTime$diagn_ETR,
      meanETl = EstimatedExitTime$meanETl, # Estimated mean exit time (left basin)
      meanETr = EstimatedExitTime$meanETr, # Estimated mean exit time (right basin)
      atolL = EstimatedExitTime$atolL,
      atolR = EstimatedExitTime$atolR,
      diagn_ETL = EstimatedExitTime$diagn_ETL,
      diagn_ETR = EstimatedExitTime$diagn_ETR
    )
  )
}


# Langevin1D function adapted from C++ source code 
# (can be found on https://github.com/cran/Langevin/blob/master/src/Langevin1D.cpp)
Langevin1D_adapted <- function(Ux, # Time series data
                               bins, # Number of bins to divide the data into
                               ntau, # Maximum time lag to consider
                               sf, # Sampling frequency of the time series
                               bw_sd, # Bandwidth for kernel smoothing
                               bin_min = 50 # Minimum number of data points required in a bin for estimation
) {
  
  # Create vector of different lags to be used when calculating the differences (increments) in the data
  steps <- 1:ntau 
  
  # Specify bandwidth used for obtaining smoothed D1 and D2 estimates
  bw <- bw_sd * sd(Ux) 
  
  # Define bin edges
  U <- seq(min(Ux, na.rm = TRUE), max(Ux, na.rm = TRUE), length.out = bins + 1) 
  nsteps <- length(steps)
  
  # Initialize matrices to store results
  M1 <- matrix(data = NA, nrow = bins, ncol = nsteps) # First moment
  eM1 <- matrix(data = NA, nrow = bins, ncol = nsteps) # Error of the first moment
  M2 <- matrix(data = NA, nrow = bins, ncol = nsteps) # Second moment
  eM2 <- matrix(data = NA, nrow = bins, ncol = nsteps) # Error of the second moment
  M4 <- matrix(data = NA, nrow = bins, ncol = nsteps) # Fourth moment
  D1 <- numeric(bins) # Drift coefficient
  eD1 <- numeric(bins) # Error of the drift coefficient
  D2 <- numeric(bins) # Diffusion coefficient
  eD2 <- numeric(bins) # Error of the diffusion coefficient
  D4 <- numeric(bins) # Fourth moment coefficient
  dens <- numeric(bins) # Density of data points in each bin
  mean_bin <- numeric(bins) # Mean of data points in each bin
  
  # Loop through each bin
  for (i in 1:bins) {
    sum_m1 <- numeric(nsteps) # Initialize sum of first moment for current bin
    sum_m2 <- numeric(nsteps) # Initialize sum of second moment for current bin
    sum_m4 <- numeric(nsteps) # Initialize sum of fourth moment for current bin
    len_step <- numeric(nsteps) # Initialize count of valid increments for each lag in current bin
    len_bin <- 0 # Initialize count of data points in current bin
    
    # For each data point and each lag size, calculate the increment as the difference between the data point at the current position and the data point at the position offset by the lag size (for all lag sizes up to the maximum timelag specified by ntau)
    
    #Iterate over each data point (n)
    for (n in 1:(length(Ux) - max(steps))) { 
      # Check if the current data point falls within the current bin
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
    
    mean_bin[i] <- mean_bin[i] / len_bin # Calculate mean of data points in current bin
    dens[i] <- max(len_step) # Set density to the maximum number of valid increments across all lags in the current bin
    
    # Proceed if the number of data points in the bin is above the minimum threshold
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
  
  D1s = ksmooth(x= mean_bin,y= D1,kernel= 'normal',bandwidth= bw,
                x.points=mean_bin) # Smooth D1 estimates using kernel smoothing
  D2s = ksmooth(x= mean_bin ,y= D2,kernel= 'normal',bandwidth=bw,
                x.points=mean_bin) # Smooth D2 estimates using kernel smoothing
  
  ret <- list(
    D1 = D1, # Raw estimates (same output as from Langevin::Langevin1D function)
    D1s = D1s, # Smoothed estimates 
    eD1 = eD1, # Error of the drift coefficient
    D2 = D2, # Raw estimates (same output as from Langevin::Langevin1D function)
    D2s = D2s, # Smoothed estimates 
    eD2 = eD2, # Error of the diffusion coefficient
    D4 = D4, # Fourth moment coefficient
    mean_bin = mean_bin, # Mean of each bin
    density = dens, # Density of data points in each bin
    M1 = M1, # First moment
    eM1 = eM1, # Error of first moment
    M2 = M2, # Second moment
    eM2 = eM2, # Error of second moment
    M4 = M4, # Fourth moment
    U = U # Bin edges
  )
  
  # class(ret) <- "Langevin"
  return(ret)
}

# apply Langevin1D_adapted
apply_Langevin1D_adapted <- function(Ux, # Time series data
                                     bins, # Number of bins
                                     ntau, # Maximum time lag
                                     sf, # Sampling frequency
                                     bin_min = 50, # Minimum data points per bin
                                     bw_sd = 0.3 # Bandwidth for smoothing
) {
  
  DDLangout = Langevin1D_adapted(Ux, bins, ntau, sf, bw_sd, bin_min = 50) # Apply Langevin1D_adapted function
  
  # # Extract smoothed output
  D1s = DDLangout$D1s # Smoothed drift coefficients
  D2s = DDLangout$D2s # Smoothed diffusion coefficients
  bin.mid = DDLangout$mean_bin # Mean of each bin
  
  # Remove NAs and zeros
  find_idx <- function(x, tolerance = 1e-8) { # Function to find indices of non-NA and non-zero values
    na_idx <- which(!is.na(x))
    zero_idx <- which(abs(x) > tolerance)
    idx <- intersect(na_idx, zero_idx)
    return(idx)
  }
  D1sx_idx <- find_idx(D1s$x) # Find indices for D1s$x
  D1sy_idx <- find_idx(D1s$y) # Find indices for D1s$y
  D2sx_idx <- find_idx(D2s$x) # Find indices for D2s$x
  D2sy_idx <- find_idx(D2s$y) # Find indices for D2s$y
  idx = sort(dplyr::intersect(
    dplyr::intersect(D1sx_idx, D1sy_idx),
    dplyr::intersect(D2sx_idx, D2sy_idx)
  )) # Find common indices
  
  D1s = list(x = D1s$x[idx],
             y = D1s$y[idx]) # Subset D1s using common indices
  D2s = list(x = D2s$x[idx],
             y = D2s$y[idx]) # Subset D2s using common indices
  bin.mid = bin.mid[idx] # Subset bin.mid using common indices
  
  # Find roots - where D1 crosses x=0 (Fixed points)
  sdrift = sign(D1s$y) # Get the sign of D1s$y
  dsdrift = c(0,-diff(sdrift)) # Calculate the difference in sign, adding a 0 at the beginning
  ixeq = which(dsdrift != 0)  # Indices of the roots - where the sign changes
  xeq = bin.mid[ixeq] # Get the x-values (bin means) at the roots (fixed points)
  
  return(list(
    D1s = D1s, # Smoothed D1 estimates 
    D2s = D2s, # Smoothed D2 estimates
    DDLangout = DDLangout, # Both raw and smoothed estimates
    bin.mid = bin.mid,
    xeq = xeq, # Fixed points
    ixeq = ixeq # Indices of fixed points
  ))
}


# Interpolate with approxfun() for D1 and D2
D1fun = function(x, D1s) {
  yhat = stats::approx(
    x = D1s$x,
    y = D1s$y,
    xout = x,
    method = 'linear',
    rule = 2
  )$y
  return(yhat)
}

D2fun = function(x, D2s)  {
  yhat = stats::approx(
    x = D2s$x,
    y = D2s$y,
    xout = x,
    method = 'linear',
    rule = 2
  )$y
  return(yhat)
}


# Get effective potential
get_effective_potential = function(D1s,
                                   D2s,
                                   interpol_steps,
                                   xvec) {
  Dratio = D1s$y / D2s$y
  D1.over.D2 = function(x) {
    yhat = stats::approx(
      x = D1s$x,
      y = Dratio,
      xout = x,
      method = 'linear',
      rule = 2
    )$y
    return(yhat)
  }
  
  # Calculate and plot effective potential  ------------------------------------------
  eprange = range(xvec)
  xvec.ep = seq(eprange[1], eprange[2], length.out = interpol_steps)
  
  # Check D1/D2
  D1D2 = rep(0, interpol_steps)
  D1D2adj = rep(0, interpol_steps)
  for (i in 1:interpol_steps) {
    D1D2[i] = D1.over.D2(xvec.ep[i])
    D1D2adj[i] = log(0.5 * D2fun(xvec.ep[i], D2s) ^ 2) - D1D2[i]
  }
  
  # Plot effective potential function
  EPF = rep(0, interpol_steps - 1)
  EPF_err = EPF
  for (i in 1:(interpol_steps - 1)) {
    x0 = xvec.ep[1]
    x1 = xvec.ep[i + 1]
    xhalf = (x0 + x1) / 2
    # These integrate commands return 'bad behavior of integrand'
    #integral = integrate(D1.over.D2,lower=x0,upper=x1)$value
    #integral = integrate(Vectorize(D1.over.D2),lower=x0,upper=x1)$value
    # This function from the cubature package seems to work
    hcub = cubature::hcubature(f = D1.over.D2,
                               lowerLimit = x0,
                               upperLimit = x1)
    logdiff = log(0.5 * D2fun(xhalf, D2s) ^ 2)  # D2 = 0.5*sigma^2 for standard Langevin
    EPF[i] = -1 * hcub$integral + logdiff
    EPF_err[i] = hcub$error
  }
  
  return(
    list(
      EPF = EPF,
      EPF_err = EPF_err,
      Dratio = Dratio,
      D1D2 = D1D2,
      D1D2adj = D1D2adj
    )
  )
}


# Get potential function
get_potential = function(D1s, # Smoothed drift coefficients
                         D2s, # Smoothed diffusion coefficients
                         xeq, # Equilibria from D1 estimate
                         interpol_steps, # Number of interpolation steps
                         xvec # Vector of x values for interpolation
) {
  
  # Initialize vectors for drift and diffusion
  drift = rep(0, interpol_steps) 
  diff = rep(0, interpol_steps) 
  
  # Calculate drift and diffusion at each interpolation step
  for (i in 1:interpol_steps) {
    drift[i] = D1fun(xvec[i], D1s) # Calculate drift using D1fun
    diff[i] = D2fun(xvec[i], D2s) # Calculate diffusion using D2fun
  }
  
  # Calculate potential function (estimated)
  PF = rep(0, interpol_steps) 
  for (i in 2:interpol_steps) {
    x0 = xvec[1] 
    x1 = xvec[i] 
    # Integrate D1fun to get potential
    PF[i] = stats::integrate(
      D1fun,
      lower = x0,
      upper = x1,
      D1s = D1s,
      subdivisions = max(c(100, interpol_steps *6)) # To prevent error of maximum subdivisions reached
    )$value 
  }
  negPF = -1 * PF # Negate potential
  
  # Equilibria (fixed points) from potential estimate
  hills_idx = pracma::findpeaks(negPF)[, 2] # Unstable fixed points (peaks)
  valleys_idx = pracma::findpeaks(-negPF)[, 2] # Stable fixed points (valleys)
  idx_xeq_PF = c(hills_idx, valleys_idx) # Combine indices of peaks and valleys
  pot_xeq = sort(xvec[c(idx_xeq_PF)]) # Get x-values of equilibria from potential estimate
  
  # Stability of equilibria (from potential estimate)
  names_xeq_PF = c(rep("unstable", length(hills_idx)), rep("stable", length(valleys_idx))) # Assign stability labels
  
  # Find closest index for pot_xeq (equilibria from potential estimate)
  stab_pot_xeq = plyr::laply(pot_xeq, function(eq) {
    idx = which.min(abs(xvec - eq)) # Find index of closest value in xvec
    idx_match = which.min(abs(idx_xeq_PF - idx)) # Find closest index in idx_xeq_PF
    return(names_xeq_PF[idx_match]) # Return stability label
  })
  
  # Find closest index for xeq (equilibria from D1 estimate)
  stab_xeq = plyr::laply(xeq, function(eq) {
    idx = which.min(abs(xvec - eq)) # Find index of closest value in xvec
    idx_match = which.min(abs(idx_xeq_PF - idx)) # Find closest index in idx_xeq_PF
    return(names_xeq_PF[idx_match]) # Return stability label
  })
  
  print('Equilibria from Potential estimate', quote = F)
  print(pot_xeq, quote = F)
  print('Stabilities of equilibria from Potential estimate', quote = F)
  print(stab_pot_xeq, quote = F)
  
  print('Equilibria from D1 estimate', quote = F)
  print(xeq[order(xeq)], quote = F)
  print('Stabilities of equilibria from D1 estimate', quote = F)
  print(stab_xeq, quote = F) # this is the stability of equilibria from D1 estimate!!!
  
  # Calculate and plot effective potential  ------------------------------------------
  EP = get_effective_potential(D1s,
                               D2s,
                               interpol_steps,
                               xvec) # Calculate effective potential
  
  return(
    utils::modifyList(list(
      drift = drift,
      diff = diff,
      xvec = xvec,
      stab_xeq = stab_xeq,
      stab_pot_xeq = stab_pot_xeq, 
      pot_xeq = pot_xeq, 
      negPF = negPF
    ), EP) # Return a list of results
  )
}


# Get exit time
get_exit_time <- function(D1s, # Smoothed drift coefficients
                          D2s, # Smoothed diffusion coefficients
                          xeq, # Equilibria from D1 estimate
                          xvec, # Vector of x values for interpolation
                          atol = 1e-5,
                          interpol_steps # Number of interpolation steps
                          ) {
  # Calculate Exit Times =======================================================
  
  # function for solving the boundary value problem as two differential equations
  #  for T (col 1) and dT/dx (col 2)
  feval2 = function(x, y, plist) {
    out1 = y[2]
    out2 = -(D1fun(x, plist$D1s) * y[2] + 1) / D2fun(x, plist$D2s)
    return(list(c(out1, out2)))
  }
  
  
  for (basin in c("left", "right")) {
    if (basin == "left") {
      # set up for bvpSolve
      yini = c(NA, 0) # NA for an initial value which is not specified
      names(yini) <- c("T", "dT/dx")
      yend = c(0, NA) # NA for a final value which is not specified
      
      # solve the left basin from x = 0 (reflecting) to x=xeq[2] (absorbing)
      # x = seq(xeq[1]-1,xeq[2],length.out=30)  # x vector, original
      x = seq(min(xvec), xeq[2], length.out = ceiling(interpol_steps / 2)) # wider interval/discretized domain
      
    } else if (basin == "right") {
      # set up for bvpSolve
      yini = c(0, NA) # reversed from left basin; NA for an "final" value which is not specified
      names(yini) <- c("T", "dT/dx")
      yend = c(NA, 0) # reversed from left basin; # NA for an "initial" value which is not specified
      
      # solve the right basin from x=xeq[2] (absorbing) to x > xeq[3] (reflecting)
      # x = seq(xeq[2],xeq[3]+1,length.out=30)  # x vector original
      x = seq(xeq[2], max(xvec), length.out = ceiling(interpol_steps / 2))  # wider interval/discretized domain
      
    }
    
    # solve with bvpcol or bvptwp
    # KCE: bvptwp() can give an error that probably has something to do with how steep the function is. Adjust the tolerance parameter atol to prevent this.
    atol_ = atol
    trycol <- NA
    class(trycol) = "try-error"
    while ("try-error" %in% class(trycol) & atol_ <= 1e-2) {
      atol_ = 10 ** (log10(atol_) + 1) # Increase tolerance
      trycol <- try({
        bvpSolve::bvptwp(
          yini = yini,
          x = x,
          func = feval2,
          yend = yend,
          parm = list(D1s = D1s, D2s = D2s),
          atol = atol_
        )
      })
    }
    
    # Diagnostics
    diagn = list(df = attr(trycol, "istate"),
                 message = utils::capture.output(bvpSolve::diagnostics.bvpSolve(trycol, type = "message")))
    
    if (basin == "left") {
      atolL = atol_
      ETL = trycol # left exit time
      diagn_ETL = diagn
    } else if (basin == "right") {
      atolR = atol_
      ETR = trycol # right exit time
      diagn_ETR = diagn
    }
  }
  
  # CALCULATE WEIGHTS FOR AVERAGE EXIT TIME -----------------------------------------------------
  # Weight of x is p(x) from stationary distribution of Langevin equation
  # Based on appendix of Carpenter & Brock, Ecology Letters 2006 based on the book
  # 'Noise-Induced Transitions' by Horsthemke and Lefever 1984
  
  # Calculate weighted average ET for both basins ===================================
  # ETL[nrow(ETL),1] is the same as ETR[1,1]; they connect at xeq[2]
  nL = length(ETL[, 1])
  nR = length(ETR[, 1])
  x = c(ETL[1, 1] - 0.01, ETL[2:nL, 1], ETR[2:(nR), 1], ETR[nR, 1] + 0.01) # x has gaps that match ETL+ETR
  dx = diff(x)
  ETboth = c(ETL[1, 2] - 0.01, ETL[2:nL, 2], ETR[2:(nR), 2], ETR[nR, 2] +
               0.01) # Concatenated exit time
  nx = length(x)
  
  # Check behavior of finner over range of x
  yfinner = finner(x, D1s, D2s)
  
  # Weights by method in Arani et al. (2021)
  # See 2020-08-20 version of this code for the H&L version
  print("Computing weights for exit time")
  
  WTS = try(expr = get_weights(x, nx, nL, nR, D1s, D2s, ETL, ETR),
            silent = F)
  
  if ("try-error" %in% class(WTS)) {
    wts = NA
    wtraw_err = NA
    meanETl = sum(ETL[, 2])
    meanETr = sum(ETR[, 2])
  } else {
    wts = WTS$wts
    wtraw_err = WTS$wtraw_err
    meanETl = WTS$meanETl
    meanETr = WTS$meanETr
  }
  
  # Format results
  ET_df = rbind(ETL[1:nrow(ETL),],
                ETR[c(2:nrow(ETR),
                      nrow(ETR)),]) %>%
    as.data.frame() %>%
    dplyr::rename("ET" = "T", "dET/dx" = "dT/dx") %>%
    dplyr::mutate(
      yfinner = yfinner,
      wts = wts,
      wtraw_err = wtraw_err,
      LR = c(rep("left", nL), rep("right", nR))
    )  %>% tidyr::gather(variable, value, -c(x, wtraw_err)) %>%
    dplyr::rename(err = wtraw_err) %>% dplyr::mutate(err = ifelse(variable == "wts", NA, err))
  
  ET_df[ET_df$x == min(ET_df$x), "x"] = ET_df[ET_df$x == min(ET_df$x), "x"] - 0.01
  ET_df[ET_df$x == max(ET_df$x), "x"] = ET_df[ET_df$x == max(ET_df$x), "x"] + 0.01
  
  return(
    list(ET_df = ET_df, # Dataframe of exit times and weights from stationary probability distribution
         diagn_ETL = diagn_ETL,
         diagn_ETR = diagn_ETR,
         atolL = atolL,
         atolR = atolR,
         meanETl = meanETl, # Mean exit time left basin
         meanETr = meanETr # Mean exit time right basin
    )
  )
}


# Get resilience metrics from potential function (basin width and potential depth)
get_resilience_potential <- function(fp_df, # Data frame of fixed points
                                     compl_df # Data frame with potential function estimates
) {
  # Get distance to basin threshold (theoretical)
  theoretical_fps = fp_df %>% 
    filter(source == "Theoretical") %>% # Filter for theoretical fixed points
    pull(xintercept) # Extract x-intercepts (fixed point locations)
  
  # Dist_thr = |x_{att} - x_{saddle}|
  # Calculate basin width as the absolute difference between stable and unstable fixed points on x-axis
  basin_width_theo_right = abs(theoretical_fps[3] - theoretical_fps[2]) # Right basin
  basin_width_theo_left = abs(theoretical_fps[1] - theoretical_fps[2]) # Left basin
  
  # Get distance to basin threshold (estimated)
  # Check if the fixed points have the expected stability pattern
  if (identical(fp_df$color,
                c(
                  "stable",
                  "unstable",
                  "stable",
                  "stable",
                  "unstable",
                  "stable"
                ))) {
    estimated_fps = fp_df %>% 
      filter(source == "Estimated") %>% # Filter for estimated fixed points
      pull(xintercept) # Extract x-intercepts (fixed point locations)
    
    # Dist_thr = |x_{att} - x_{saddle}|
    # Calculate basin width as the absolute difference between stable and unstable fixed points on x-axis
    basin_width_est_right = abs(estimated_fps[3] - estimated_fps[2]) # Right basin
    basin_width_est_left = abs(estimated_fps[1] - estimated_fps[2]) # Left basin
  } else{
    basin_width_est_right = NULL # Set to NULL if the pattern doesn't match
    basin_width_est_left = NULL
  }
  
  # Get potential depth (theoretical)
  theoretical_potential = compl_df %>% 
    filter(variable == "potential", source == "Theoretical") # Filter for theoretical potential values
  
  # Create an interpolation function for the theoretical potential
  theo_potentialfun = approxfun(
    x = theoretical_potential$x,
    y = theoretical_potential$value,
    method = "linear",
    rule = 2
  )
  theo_fps_Ux = theo_potentialfun(theoretical_fps) # Get potential values at theoretical fixed points
  
  # Pot_dist = |U(x_{att})-U(x_{saddle})|
  # Calculate potential depth as the absolute difference in potential between stable and unstable fixed points (y-axis)
  potential_depth_theo_right = abs(theo_fps_Ux[3] - theo_fps_Ux[2]) # Right basin
  potential_depth_theo_left = abs(theo_fps_Ux[1] - theo_fps_Ux[2]) # Left basin
  
  # get potential depth (estimated)
  # Check if the fixed points have the expected stability pattern
  if (identical(fp_df$color,
                c(
                  "stable",
                  "unstable",
                  "stable",
                  "stable",
                  "unstable",
                  "stable"
                ))) {
    estimated_potential = compl_df %>% 
      filter(variable == "potential", source == "Estimated") # Filter for estimated potential values
    
    # Create an interpolation function for the estimated potential
    est_potentialfun = approxfun(
      x = estimated_potential$x,
      y = estimated_potential$value,
      method = "linear",
      rule = 2
    )
    est_fps_Ux = est_potentialfun(estimated_fps) # Get potential values at estimated fixed points
    
    # Pot_dist = |U(x_{att})-U(x_{saddle})|
    # Calculate potential depth as the absolute difference in potential between stable and unstable fixed points (y-axis)
    potential_depth_est_right = abs(est_fps_Ux[3] - est_fps_Ux[2]) # Right basin
    potential_depth_est_left = abs(est_fps_Ux[1] - est_fps_Ux[2]) # Left basin
  } else{
    potential_depth_est_right = NULL # Set to NULL if the pattern doesn't match
    potential_depth_est_left = NULL
  }
  
  # Return a list of calculated resilience metrics
  return(
    list(
      basin_width_theo_left = basin_width_theo_left,
      basin_width_theo_right = basin_width_theo_right,
      basin_width_est_left = basin_width_est_left,
      basin_width_est_right = basin_width_est_right,
      potential_depth_theo_left = potential_depth_theo_left,
      potential_depth_theo_right = potential_depth_theo_right,
      potential_depth_est_left = potential_depth_est_left,
      potential_depth_est_right = potential_depth_est_right
    )
  )
}


# Get variance and skewness around dominant modes of probability distribution of system states ("realizations")
get_resilience_prob_dist <- function(Ux, bw_sd) {
  
  # Estimate the density
  density_estimate <- density(Ux)
  
  # Extract x and y values
  x_values <- density_estimate$x
  y_values <- density_estimate$y
  
  # Identify the modes
  find_modes <- function(x, y) {
    # Find local maxima
    local_maxima <- which(diff(sign(diff(y))) == -2) + 1
    modes <- x[local_maxima]
    return(modes)
  }
  
  modes <- find_modes(x_values, y_values)
  
  # Calculate the variance and skewness around each mode
  calculate_moments <- function(x, y, mode) {
    variance <- sum((x - mode) ^ 2 * y) / sum(y)
    skewness <- sum((x - mode) ^ 3 * y) / sum(y)
    skewness <- skewness / (variance ^ (3 / 2))
    return(list(variance = variance, skewness = skewness))
  }
  
  moments <- lapply(modes, function(mode)
    calculate_moments(x_values, y_values, mode))
  
  return(
    list(
      est_var_mode_left = moments[[1]]$variance,
      est_var_mode_right = moments[[1]]$skewness,
      est_skewness_mode_left = moments[[2]]$variance,
      est_skewness_mode_right = moments[[2]]$skewness
    ))
}


# Function for inner integral
finner = function (z, D1s, D2s) {
  fx = D1fun(z, D1s) / (D2fun(z, D2s))
  return(fx)
}

# Function for g^2 weights
gg = function(z, D2s) {
  fx = 1 / (D2fun(z, D2s))
  return(fx)
}


# Get stationary probability density function to weight mean exit times of all initial conditions
get_weights <- function(x, nx, nL, nR, D1s, D2s, ETL, ETR) {
  wtraw = rep(0, nx)
  wtraw_err = wtraw
  for (i in 2:(nx)) {
    hcub = cubature::hcubature(
      f = finner,
      lowerLimit = x[1],
      upperLimit = x[i],
      D1s = D1s,
      D2s = D2s
    )
    epart = exp(hcub$integral)
    wtraw[i] = epart * gg(x[i], D2s)
    wtraw_err[i] = hcub$error
  }
  
  # normalize weights
  wts = wtraw / sum(wtraw)
  
  # artistic convenience on first weight
  wts[1] = 0.98 * wts[2]
  
  # Calculate mean exit time for left basin ---------------------------------------------
  # save axes
  xL = x[1:nL]
  wtsL = wts[1:nL]
  
  meanETl = sum(ETL[, 2] * wtsL) / sum(wtsL)
  print('', quote = F)
  print('Mean ET for left basin', quote = F)
  print(meanETl)
  print('-----------------------------------------------', quote = F)
  
  # Calculate weighted average ET for right basin ===========================================================
  xR = x[(nL + 1):nx]
  wtsR = wts[(nL + 1):nx]
  
  meanETr = sum(ETR[, 2] * wtsR) / sum(wtsR)
  print('', quote = F)
  print('Mean ET for right basin', quote = F)
  print(meanETr)
  print('-----------------------------------------------', quote = F)
  
  return(list(
    wts = wts,
    wtraw_err = wtraw_err,
    meanETl = meanETl,
    meanETr = meanETr
  ))
}


# Plot timeseries, probability distribution, and both theoretical drift, diffusion, potential function, effective potential, weights, and exit times
new_plot_overview <-
  function(out,
           filepath_image
  ) {
    with(out, {
      # Create a title for the plot using LaTeX formatting
      # The title includes the scenario, type_D2, strength_D2, sf, and N 
      title_plot = latex2exp::TeX(
        sprintf(
          "%s; %s $D_2$ (strength: %.2f); $f_s$ = %.2f, $N$ = %d",
          stringr::str_to_title(scenario),
          stringr::str_to_title(type_D2),
          strength_D2,
          sf,
          N
        ),
        bold = TRUE
      )
      
      # Set the linewidth for the fixed points in the plot
      lwd_FP = 0.8
      
      # Assign colors to vertical lines at theoretical stable and unstable fixed points
      col_FP_theoretical = dplyr::recode(
        stabs$stab_fps,
        "unstable" = scales::viridis_pal(option = "inferno")(50)[29],
        "stable" = scales::viridis_pal(option = "viridis")(50)[45]
      )
      
      # Assign colors to vertical lines at estimated stable and unstable fixed points
      col_FP_estimated = dplyr::recode(
        est_Carp$fp_df$color,
        "unstable" = scales::viridis_pal(option = "inferno")(50)[25],
        "stable" = scales::viridis_pal(option = "viridis")(50)[49]
      )
      
      # Check if mean exit time was estimated
      if(is.null(est_Carp$meanETl) == FALSE){
        # If estimated, create labels for the drift, diffusion, potential, EPF, ET, and wts variables
        variable.labs = c(
          "drift" = latex2exp::TeX("Drift", output = 'character'),
          "diffusion" = latex2exp::TeX("Diffusion", output =   'character'),
          "potential" = latex2exp::TeX("Potential", output = 'character'),
          "EPF" = latex2exp::TeX("Effective Potential", output = 'character'),
          "ET" = latex2exp::TeX("Exit Time", output = 'character'),
          "wts" = latex2exp::TeX("Weights", output = 'character')
        )
      } else { 
        
        # If not estimated, create labels only for drift, diffusion, potential, and EPF variables
        variable.labs = c(
          "drift" = latex2exp::TeX("Drift", output = 'character'),
          "diffusion" = latex2exp::TeX("Diffusion", output =   'character'),
          "potential" = latex2exp::TeX("Potential", output = 'character'),
          "EPF" = latex2exp::TeX("Effective Potential", output = 'character')
        ) # no "ET" and "wts" variables
      }
      
      # Determine the minimum length between the vector Ux and plot_t
      # plot_t =  min(length(as.vector(Ux)), plot_t)
      
      # Set the size of the points in the plot
      sz_point = 0.35
      
      # Plot timeseries (pl_ts) 
      pl_ts = 
        # Initialize the ggplot with a data frame containing time and Ux values
        ggplot2::ggplot(data.frame(t = 1:length(Ux) / sf / 60, x = as.vector(Ux))) +
        
        # Add a horizontal line at the y-intercept specified by stabs$fps
        # The line is dashed, with a specified color and linewidth
        ggplot2::geom_hline(
          yintercept = stabs$fps,
          color = col_FP_theoretical,
          linetype = 'dashed',
          linewidth = lwd_FP,
          alpha = 1
        ) +
        
        # Add points to the plot with specified aesthetics
        ggplot2::geom_point(
          ggplot2::aes(x = t, y = x),
          alpha = .75,
          col = "gray30",
          size = sz_point
        ) +
        
        # Add labels to the plot, including a title and subtitle
        # The subtitle includes the analytical and estimated fixed points
        ggplot2::labs(
          x = "Time (min)",
          y = latex2exp::TeX("$x$"),
          title = title_plot,
          subtitle =
            sprintf(
              "Timeseries with analytical fixed points: %s; Estimated fixed points: %s",
              paste0(sort(round(stabs$fps, 2)), collapse = ", "),
              paste0(sort(round(est_Carp$xeq, 2)), collapse = ", ")
            )
        ) +
        
        # Adjust the x-axis to start at 0
        ggplot2::scale_x_continuous(expand = c(0, 0))
      
      
      # Plot inset plot of timeseries (plot_ts_inset) with timeseries from first 1.67 min
      plot_t_inset = 100 * sf
      pl_ts_inset = ggplot2::ggplot(data.frame(t = 1:plot_t_inset / sf / 60, x = as.vector(Ux)[1:plot_t_inset])) +
        ggplot2::geom_hline(
          yintercept = stabs$fps,
          color = col_FP_theoretical,
          linetype = 'dashed',
          linewidth = lwd_FP,
          alpha = 1
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = t, y = x),
          size = .3,
          alpha = .75,
          col = "gray30"
        )  +
        ggplot2::labs(x = "", y = "", title = "") +
        ggplot2::scale_x_continuous(expand = c(.005, 0))
      
      # Add inset plot (pl_ts_inset) to timeseries plot (pl_ts) and style with style_plot function
      # The inset plot is positioned and sized within the main plot
      pl_ts_with_inset <-
        cowplot::ggdraw() +
        cowplot::draw_plot(style_plot(pl_ts)) +
        cowplot::draw_plot(
          style_plot(pl_ts_inset)  +
            # Style the inset plot: transparent background and no grid
            ggplot2::theme(
              plot.background = ggplot2::element_rect(fill = 'transparent', color =
                                                        NA),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank()
            ) +
            # Remove margin around the inset plot
            ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "pt")),
          x = 0.04,
          y = .12,
          width = .2,
          height = .5
        )
      
      # Create a density plot of timeseries states (pl_dens) 
      pl_dens <- ggplot2::ggplot(data.frame(x = as.vector(Ux))) +
        # Add a histogram with density on the y-axis
        ggplot2::geom_histogram(
          ggplot2::aes(x = x, y = ggplot2::after_stat(density)),
          bins = bins,
          colour = "lightgrey",
          fill = "white"
        ) +
        # Add a density line to the histogram
        ggplot2::geom_density(ggplot2::aes(x = x), linewidth = 0.8)  +
        # Add vertical lines at the theoretical fixed points 
        ggplot2::geom_vline(
          xintercept = stabs$fps,
          color = col_FP_theoretical,
          linetype = 'dashed',
          linewidth = lwd_FP,
          alpha = .75
        ) +
        # Adjust the x-axis to start at 0
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::labs(
          # Add labels to the plot, including number of modes of the distribution
          x = latex2exp::TeX("$x$"),
          y = "Density",
          title = "",
          subtitle = sprintf(
            "Density (%d mode%s)",
            length(LaplacesDemon::Modes(Ux)$modes),
            ifelse(length(LaplacesDemon::Modes(Ux)$modes) > 1, "s", "")
          )
        )
      
      # Plot theoretical drift, diffusion, and potential function (pl_Carp)
      # Set x-axis limit to maximum value of timeseries 
      xlim = max(abs(Ux))
      
      # pl_Carp with theoretical and estimated on the same facet plots
      pl_Carp =
        est_Carp$compl_df %>%
        # Filter the data to include only those variables present in variable.labs
        dplyr::filter(variable %in% names(variable.labs)) %>%
        # Convert the 'value' column to numeric and add a new column 'variable_name' by recoding 'variable' column based on variable.labs
        dplyr::mutate(
          value = as.numeric(as.character(value)),
          variable_name = dplyr::recode_factor(variable,!!!variable.labs, .ordered = T)
        ) %>%
        # Initialize a ggplot
        ggplot2::ggplot() +
        # Add a horizontal line at y=0
        ggplot2::geom_hline(
          yintercept = 0,
          color = 'grey10',
          linetype = "solid",
          linewidth = .2,
          alpha = 1
        ) +
        # Add vertical lines at both theoretical and estimated fixed points
        ggplot2::geom_vline(
          data = est_Carp$fp_df[4:6,],
          ggplot2::aes(xintercept = xintercept, color = interaction(color, source)),
          linetype = 'dashed',
          linewidth = lwd_FP,
          alpha = .75
        ) +
        # Define manual color scale for the vertical lines
        ggplot2::scale_color_manual(
          name = "",
          breaks = c(
            # "unstable.Theoretical",
            # "stable.Theoretical",
            "unstable.Estimated",
            "stable.Estimated"
          ), 
          values = c(
            # "unstable.Theoretical" = rgb(0, 114, 178, maxColorValue = 255),
            # "stable.Theoretical" = rgb(0, 114, 178, maxColorValue = 255),
            "unstable.Estimated" = scales::viridis_pal(option = "inferno")(50)[29],
            "stable.Estimated" = scales::viridis_pal(option = "viridis")(50)[45]
          ),
          guide = "none"
        ) +
        ggnewscale::new_scale_color() +
        # Add points representing estimated variables, color them based on 'source' (Theoretical or Estimated)
        ggplot2::geom_point(ggplot2::aes(x = x, y = value, color = source), size = .5) +
        # Define manual color scale for the points
        ggplot2::scale_color_manual(
          name = "",
          breaks = c("Theoretical", "Estimated"),
          values = c(
            "Theoretical" = 'grey10',
            "Estimated" = 'red3'),
          guide = "none"
        ) +
        # Adjust the x-axis to start at 0 
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        # Move the y-axis to the right
        ggplot2::scale_y_continuous(position = 'right') +
        # Create facets based on 'variable_name', with free y scales
        ggh4x::facet_grid2(. ~ variable_name,
                           scales = 'free_y', independent = 'y',
                           labeller = label_parsed, switch = "y", axes = "all") +
        # Remove the legend
        ggplot2::theme(legend.position = "none")
      
      
      # Check if exit time was estimated 
      if (is.null(est_Carp$meanETl) == FALSE) {
        # If estimated, add below labels to the plot 
        pl_Carp = pl_Carp +
          ggplot2::labs(
            y = "",
            # Y-axis label is empty
            x = latex2exp::TeX("$x$"),
            # X-axis label using LaTeX for formatting
            title = latex2exp::TeX(
              sprintf(
                # Title with theoretical coefficients and other parameters
                "Theoretical coefficients: $D_1(x) = %.2fx^3 + %.2fx^2 + %.2fx$; $D_2(x) = %.2fx^2 + %.2f$ (%d bins, $n_{tau}$ = %d, %d interpolation steps; smoothing = %.2f)",
                D$d13,
                D$d12,
                D$d11,
                D$d22,
                D$d20,
                bins,
                ntau,
                interpol_steps,
                bw_sd
              ),
              bold = T  # Make the title bold
            ),
            subtitle = latex2exp::TeX(
              sprintf(
                # Subtitle with theoretical and estimated exit times
                "Theoretical Exit Time $mu_{left} =$ %.2f, $mu_{right} =$ %.2f; Estimated Exit Time $mu_{left} =$ %.2f (tol = %.4f, %sconverged), $mu_{right} =$ %.2f  (tol = %.4f, %sconverged)",
                est_Carp$theo_meanETl,
                est_Carp$theo_meanETr,
                est_Carp$meanETl,
                est_Carp$atolL,
                ifelse(as.logical(est_Carp$diagn_ETL$df[["flag"]]), "not ", ""),
                est_Carp$meanETr,
                est_Carp$atolR,
                ifelse(as.logical(est_Carp$diagn_ETL$df[["flag"]]), "not ", "")
              )
            )
          )
      } else {
        # If mean exit time was not estimated, add labels stating mean exit time was not estimated
        pl_Carp = pl_Carp +
          ggplot2::labs(
            y = "",
            # Y-axis label is empty
            x = latex2exp::TeX("$x$"),
            # X-axis label using LaTeX for formatting
            title = latex2exp::TeX(
              sprintf(
                # Title with theoretical coefficients and other parameters
                "Theoretical coefficients: $D_1(x) = %.2fx^3 + %.2fx^2 + %.2fx$; $D_2(x) = %.2fx^2 + %.2f$ (%d bins, $n_{tau}$ = %d, %d interpolation steps; smoothing = %.2f)",
                D$d13,
                D$d12,
                D$d11,
                D$d22,
                D$d20,
                bins,
                ntau,
                interpol_steps,
                bw_sd
              ),
              bold = T  # Make the title bold
            ),
            subtitle = latex2exp::TeX(
              sprintf(
                # Subtitle with theoretical exit times and a message indicating that the mean exit time could not be estimated
                "Theoretical Exit Time $mu_{left} =$ %.2f, $mu_{right} =$ %.2f; Exit Time could not be estimated",
                est_Carp$theo_meanETl,
                est_Carp$theo_meanETr
              )
            )
          )
      }
      
      # Combine plots and style using the style_plot function
      pl_combo = cowplot::plot_grid(
        cowplot::plot_grid(
          pl_ts_with_inset,  # Timeseries plot with inset
          style_plot(pl_dens) + ggplot2::labs(x = ""),  # Styled density plot without x-axis label 
          rel_widths = c(1, .3),  
          nrow = 1  
        ),
        NULL,  # Add space between the top and bottom plots
        # Add the styled pl_Carp plot with a customized legend position
        style_plot(pl_Carp) +
          ggplot2::theme(legend.position = c(.87, 1)),  
        ncol = 1,  
        rel_heights = c(.75, .025, 1)  # Relative heights of the plots and space between them
      )
      
      # Close any open graphics devices
      grDevices::graphics.off()
      
      # Create a PDF file to save the combined plot
      filepath_image %>%
        grDevices::pdf(
          width = 16,  # Set the width of the PDF
          height = 8,  # Set the height of the PDF
          bg = "transparent"  # Set the background to transparent
        )
      
      # Draw the combined plot onto the grid
      grid::grid.draw(pl_combo)
      
      # Close the PDF device
      grDevices::dev.off()
      
      # Close any remaining open graphics devices
      grDevices::graphics.off()
      
      # Check if the file exists and print the file path if it does
      if (file.exists(filepath_image)) {
        print(filepath_image)
      }
      
    })
  }

# Style plot
style_plot <- function(pl) {
  return(
    pl + ggplot2::theme_bw() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
        axis.title = ggplot2::element_text(hjust = 0.5, size =
                                             16),
        axis.text = ggplot2::element_text(size =
                                            12),
        legend.title = ggplot2::element_text(hjust = 0.5, size =
                                               16),
        legend.text = ggplot2::element_text(hjust = 0.5, size =
                                              18),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 18),
        # Increase font size facet labels
        strip.text.y = ggplot2::element_text(size = 16,
                                             margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm")),
        strip.text.x = ggplot2::element_text(size = 16,
                                             margin = ggplot2::margin(0.1, 0.1, 0.1, 0.1, "cm")),
        plot.margin = unit(c(10, 10, 10, 10), "pt"),
        text = ggplot2::element_text(family = "serif")
      )
  )
}

setup_filepaths <- function(filepath_est,
                            filepath_figs,
                            scenario,
                            type_D2,
                            strength_D2,
                            N_high_res,
                            sf_high_res,
                            sf,
                            N,
                            noise_iter,
                            step_idx,
                            d13,
                            d12,
                            d11,
                            d10,
                            d22,
                            d21,
                            d20,
                            bins,
                            ntau,
                            interpol_steps,
                            bw_sd
) {
  # Set up filepath for estimates
  filepath_out = file.path(
    filepath_est,
    sprintf("%s", scenario),
    sprintf("D2strength_%s", strength_D2),
    sprintf("%s-D2", type_D2), 
    sprintf("N%s", N), 
    sprintf("sf%s", sf),
    sprintf("tau%s", ntau), 
    sprintf(
      "D2strength%.4f_sf%.2f_N%d_iter%04d_step%04d_pars%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_bins%d_ntau%d_interpol%d_bw%.2f_Downsampled_from_sf%d_N%d.RDS",
      strength_D2,
      sf,
      N, 
      noise_iter,
      step_idx,
      d13,
      d12,
      d11,
      d10,
      d22,
      d21,
      d20,
      bins, 
      ntau, 
      interpol_steps, 
      bw_sd,
      sf_high_res,
      N_high_res
    )
  )
  # Set up filepath for plot
  filepath_image = file.path(
    filepath_figs,
    sprintf("%s", scenario),
    sprintf("D2strength_%s", strength_D2),
    sprintf("%s-D2", type_D2), 
    sprintf("N%s", N), 
    sprintf("sf%s", sf), 
    sprintf("tau%s", ntau), 
    stringr::str_replace(basename(filepath_out), ".RDS", ".pdf")
  )
  
  if (!dir.exists(dirname(filepath_out))) {
    dir.create(dirname(filepath_out), recursive = T)
  }
  if (!dir.exists(dirname(filepath_image))) {
    dir.create(dirname(filepath_image), recursive = T)
  }
  return(list(filepath_out = filepath_out, filepath_image = filepath_image))
  
}

