### Bootstrap 
# load example data
example <-
  readRDS(
    "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/Rinn_est_scale_check/right-fp-gains-dominance/constant-D2/1000-N/D2strength0.3000_sf10_N1000_iter0005_step0001_pars-3.00_0.00_3.00_0.00_0.00_0.00_0.30_bins100_ntau3_interpol50_bw0.30.RDS"
  )
example_Ux <- example$Ux
plot(example_Ux)





# bootstrap: estimate metrics using different resamples of Ux
bootstrap_resample <-
  function(bootstrap_n,
           # bootstrapDDs,
           Ux,
           D,
           N,
           sf,
           ntau,
           bins,
           bw_sd,
           interpol_steps,
           compl_df) {
    execution_times <- list()
    full_boot_est_Carp <- list()
    
    
    
    # estimated_prob_dist = compl_df %>%
    #   filter(variable == "wts", source == "Estimated")
    # estimated_prob_dist$value <- as.numeric(estimated_prob_dist$value)
    
    # random_initial_states <-
    #   sample(
    #     estimated_prob_dist$x,
    #     size = bootstrap_n,
    #     replace = TRUE,
    #     prob = estimated_prob_dist$value
    #   )
    # plot(random_initial_states)
    
    for (i in 1:bootstrap_n) {
      start_time <- Sys.time()
      
      # Get stability of fixed points
      boot_stabs <- get_stability(D)
      
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
        # full_boot_Ux =  full_boot_Ux,
        execution_times = execution_times
      )
    )
    

  }