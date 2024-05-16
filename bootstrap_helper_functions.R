# Created May 16th 

# Get "parent function"
## i.e. Fit polynomials to the estimated drift and diffusion coefficients
bootstrap_est_D <- function(D1s, D2s, eD1, eD2, mean_bin, bw_sd, Ux){
  # not smoothed
  # estD1 <- coef(lm(D1 ~ mean_bin + I(mean_bin^2) + I(mean_bin^3), weights = 1/eD1)) 
  # estD2 <- coef(lm(D2 ~ mean_bin+ I(mean_bin^2), weights = 1/eD2)) 
  
  bw <- bw_sd * sd(Ux)
  eD1s = ksmooth(x= mean_bin,y= eD1,kernel= 'normal',bandwidth= bw,
                 x.points=mean_bin)
  eD2s = ksmooth(x= mean_bin ,y= eD2,kernel= 'normal',bandwidth=bw,
                 x.points=mean_bin)
  
  # fit simple model
  # estD1_simple <- coef(lm(D1s$y ~ 0 + D1s$x + I(D1s$x^3), weights = 1/eD1s$y)) # specify only x and x^3?
  # estD2_simple <- coef(lm(D2s$y ~ 1, weights = 1/eD2s$y))
  
  # fit "complex" polynomial (smoothed)
  estD1 <- coef(lm(D1s$y ~ D1s$x + I(D1s$x^2) + I(D1s$x^3), weights = 1/eD1s$y))
  estD2 <- coef(lm(D2s$y ~ 1, weights = 1/eD2s$y))
  
  return(list(
    d13 = round(estD1[4], digits = 3),
    d12 = round(estD1[3], digits = 3),
    d11 = round(estD1[2], digits = 3),
    d10 = round(estD1[1], digits = 3),
    d22 = 0,
    d21 = 0,
    d20 = round(estD2[1], digits = 3)
    # d13 = estD1_simple[2],
    # d12 = 0,
    # d11 = estD1_simple[1],
    # d10 = 0,
    # d22 = 0,
    # d21 = 0,
    # d20 = estD2_simple[1]
  ))
}

# Bootstrap function: simulate time series from the reconstructed coefficients
## Written based on bootstrap steps from Arani supplementary materials
bootstrap_Lang <-
  function(bootstrap_n,
           bootstrapDDs,
           N,
           sf,
           ntau,
           bins,
           bw_sd,
           interpol_steps,
           compl_df) {
    execution_times <- list()
    full_boot_est_Carp <- list()
    full_boot_Ux <- list()
    
    estimated_prob_dist = compl_df %>%
      filter(variable == "wts", source == "Estimated")
    estimated_prob_dist$value <- as.numeric(estimated_prob_dist$value)
    
    random_initial_states <-
      sample(
        estimated_prob_dist$x,
        size = bootstrap_n,
        replace = TRUE,
        prob = estimated_prob_dist$value
      )
    # plot(random_initial_states)
    
    for (i in 1:bootstrap_n) {
      start_time <- Sys.time()
      
      boot_Ux <-
        do.call(timeseries1D, utils::modifyList(
          bootstrapDDs,
          list(
            N = N * sf,
            sf = sf,
            startpoint = random_initial_states[i]
          )
        ))
      
      if (any(is.na(boot_Ux) | is.infinite(boot_Ux))) {
        print("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
      } else {
        print("Timeseries simulated successfully")
      }
      
      full_boot_Ux[[i]] <- boot_Ux
      
      # Get stability of fixed points
      boot_stabs <- get_stability(bootstrapDDs)
      
      # Apply Exit Time Analysis (Carpenter's (2022) approach adapted )
      print("Estimate Langevin model & calculate exit time and alt. metrics")
      boot_est_Carp <- est_D_Carp(
        Ux = boot_Ux,
        sf = sf,
        D = bootstrapDDs,
        stabs = boot_stabs,
        Tstep = 1:length(as.vector(boot_Ux)) / sf,
        ntau = ntau,
        bins = bins,
        bw_sd = bw_sd,
        interpol_steps = interpol_steps
      )
      
      end_time <- Sys.time()
      execution_times[[i]] <- end_time - start_time
      
      full_boot_est_Carp[[i]] <- boot_est_Carp
    }
    
    return(
      list(
        full_boot_est_Carp = full_boot_est_Carp,
        full_boot_Ux =  full_boot_Ux,
        execution_times = execution_times
      )
    )
  }