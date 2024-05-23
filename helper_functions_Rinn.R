### Modifications by Luiza 18.03.2024 ###
# (1) change in est_Carp_D function to calculate theoretical exit time when estimated fps are not "stable" "unstable" "stable"
# (2) specify function packages for parallel processing (removed 'pracma' for ceiling function)
# (3) removed if(datagen = "Langevin") in plot_overview_new

# dirname(rstudioapi::getActiveDocumentContext()$path)

# Generate timeseries using Langevin
generate_Langevin <- function(N,
                              sf,
                              noise_iter,
                              d13,
                              d12,
                              d11,
                              d10,
                              d22,
                              d21,
                              d20) {
  # Generate the time series
  set.seed(noise_iter)
  Ux <- Langevin::timeseries1D(
    N = N * sf,
    # startpoint = startpoint,
    d13 = d13,
    d12 = d12,
    d11 = d11,
    d10 = d10,
    d22 = d22,
    d21 = d21,
    d20 = d20,
    sf = sf,
    dt = 1 / sf
  )
  
  # Cut the timeseries if it went to infinity
  if (any(is.na(Ux) | is.infinite(Ux))) {
    message("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
    return()
  } else {
    return(Ux)
  }
}

# Polynomial (third-order)
poly3 <-
  function(x, a0, a1, a2, a3) {
    a0 + a1 * x + a2 * x ^ 2 + a3 * x ^ 3
  }

# Derive theoretical drift, diffusion, and potential analytically
get_theoretical_D <- function(D, min_x = -5, max_x = 5, interpol_steps) {
  with(D, {
    # Potential function from drift coefficients
    # drift = d10 + d11*x + d12*x^2 + d13*x^3 = negative derivative of potential
    # P = constant + a1*x + a2*x^2 + a3*x^3 + a4*x^4
    # d13 = a4 * 4
    # d12 = a3 * 3
    # d11 = a2 * 2
    # d10 = a1 * 1
    # constant falls away
    a4 = d13 / 4
    a3 = d12 / 3
    a2 = d11 / 2
    a1 = d10
    
    x = seq(min_x, max_x, length.out = interpol_steps)
    
    theoretical_df =
      data.frame(
        x = x,
        drift = d10 + d11 * x + d12 * x ^ 2 + d13 * x ^ 3,
        diffusion =  d20 + d21 * x + d22 * x ^ 2,
        potential = -(a1 * x + a2 * x ^ 2 + a3 * x ^ 3 + a4 *
                        x ^ 4)) %>%
      # Add effective potential
      dplyr::mutate(EPF = c(NA, get_effective_potential(list(x = x, y = drift),
                                                        list(x = x, y = diffusion),
                                                        interpol_steps = interpol_steps,
                                                        xvec = x)$EPF)) %>%
      tidyr::gather(variable, value,-x)
    
    return(theoretical_df)
  })
}

get_D <- function(nr_steps_bif,
                  scenario, 
                  # = c(
                  #   "2fps-balanced-deepening",
                  #   # "2fps-balanced-shallowing",
                  #   "left-fp-gains-dominance",
                  #   "right-fp-gains-dominance"
                  # )[1],
                  type_D2,
                  # = c("constant", "quadratic")[1],
                  strength_D2) {
  # if (scenario == "1fp") {
  #   # alphas = 0
  #   # betas = 1
  #   d13 = 0
  #   d12 = 0
  #   d11 = seq(3, 1, length.out = nr_steps_bif) # potential is plotted as negPF
  #   d10 = 0
  # } else  if (scenario == "2fps-balanced-shallowing") {
  #   # alphas = 0 #-.1
  #   # betas = seq(2, 0.1, length.out = nr_steps_bif)
  #   d13 = seq(-3,-1, length.out = nr_steps_bif)
  #   d12 = 0
  #   d11 = seq(3, 1, length.out = nr_steps_bif)
  #   d10 = 0
  # } else  
  if (scenario == "2fps-balanced-deepening") {
    # alphas = 0 #-.1
    # betas = rev(seq(2, 0.1, length.out = nr_steps_bif))
    d13 = rev(seq(-3,-1, length.out = nr_steps_bif))
    # d13 = rev(seq(-5,-1, length.out = nr_steps_bif))
    d12 = 0
    d11 = rev(seq(3, 1, length.out = nr_steps_bif))
    # d11 = rev(seq(5, 1, length.out = nr_steps_bif))
    d10 = 0
  } else if (scenario == "right-fp-gains-dominance") {
    # alphas = seq(-0.3, 0.3, length.out = nr_steps_bif)
    # betas = 1.5
    d13 = -3
    # d12 = seq(0.8, 0.5, length.out = nr_steps_bif)
    d12 = seq(0.05, 0.25, length.out = nr_steps_bif) #changed 25.04.2024
    d11 = 3
    d10 = 0
  } else if (scenario == "left-fp-gains-dominance") {
    # alphas = rev(seq(-0.3, 0.3, length.out = nr_steps_bif))
    # betas = 1.5
    d13 = -2.5
    d12 = seq(-0.8,-0.5, length.out = nr_steps_bif)
    d11 = 2.5
    d10 = 0
  }
  
  if (type_D2 == "constant") {
    d22 = 0
    d21 = 0
    d20 = 1
    # Cobb1981
    # cardan = 27 * alphas[1] ^ 2 - 4 * betas[1] ^ 3
    # modality = ifelse(cardan > 0, "unimodal", "bimodal")
    # alpha_meaning = ifelse(modality == "unimodal",
    #                        "skewness",
    #                        "the relative height of the two modes")
    # beta_meaning = ifelse(modality == "unimodal",
    #                       "kurtosis",
    #                       "the separation of the two modes")
    #
    # sprintf(
    #   "For the starting bifurcation parameter value, Cardan's discriminant is %.4f, meaning the PDF is %s if the diffusion is of noise type N (i.e. a constant). Alpha thus determines %s and beta determines %s.",
    #   cardan,
    #   modality,
    #   alpha_meaning,
    #   beta_meaning
    # )
    
    # } else if (type_D2 == "linear-positive") {
    #   d22 = 0
    #   d21 = 1
    #   d20 = 1
    # } else if (type_D2 == "linear-negative") {
    #   d22 = 0
    #   d21 = -1
    #   d20 = 1
  } else if (type_D2 == "quadratic") {
    d22 = 1
    d21 = 0
    d20 = 1
  }
  
  Ds = data.frame(
    d13 = d13,
    d12 = d12,
    d11 = d11,
    d10 = d10,
    # d13 = -1,
    # d12 = 0,
    # d11 = betas,
    # d10 = alphas,
    d22 = d22 * strength_D2,
    d21 = d21 * strength_D2,
    d20 = d20 * strength_D2
  )  %>% purrr::transpose() %>% unique() %>%
    # In case of no parameter change, repeat same sequence
    rep(ifelse(length(.) == 1, nr_steps_bif, 1))
  
  # return(list(d10=d10,d11=d11,d12=d12,d13=d13, d20=d20,d21=d21,d22=d22))
  return(Ds)
}


get_stability <- function(D) {
  with(D, {
    # Find fixed points
    fps = sort(unique(round(Re(polyroot(
      c(d10, d11, d12, d13)
    )), 10)))
    f_prime = D(expression(d10 + d11 * x + d12 * x ^ 2 + d13 * x ^ 3), "x")
    stab = plyr::laply(fps, function(x) {
      eval(f_prime)
    })
    stab_fps = dplyr::recode(as.numeric(stab < 0), "0" = "unstable", "1" = "stable")
    
    print('Theoretical Equilibria', quote = F)
    print(fps, quote = F)
    print(stab_fps, quote = F)
    return(list(fps = fps, stab_fps = stab_fps))
  })
}


# est_D_Rinn <- function(Ux,
#                        bins,
#                        steps,
#                        bin_min,
#                        kernel,
#                        d10,
#                        d11,
#                        d12,
#                        d13,
#                        d20,
#                        d21,
#                        d22) {
#   # Estimate drift and diffusion coefficients
#   ests <- Langevin::Langevin1D(Ux, bins, steps,
#                                bin_min = bin_min, # minimal number of events per bin, defaults to 100
#                                kernel = kernel)
#   # ". In a histogram based regression, the size of the bins located at xi is limited by the available amount of data. Therefore, a kernel based regression as described in [15] with use of the Nadaraya-Watson estimator is favorable, which results in a smooth curve." (Honisch)
#   est_df = data.frame(
#     bin_nr = 1:length(ests$mean_bin),
#     mean_bin = ests$mean_bin,
#     D1 = ests$D1,
#     eD1 = ests$eD1,
#     D2 = ests$D2,
#     eD2 = ests$eD2
#   )
#
#   estD1 <-
#     coef(stats::lm(
#       D1 ~ mean_bin + I(mean_bin ^ 2) + I(mean_bin ^ 3),
#       weights = 1 / eD1,
#       data = est_df %>% filter(eD1 != 0) %>% arrange(bin_nr)
#     ))
#   estD2 <-
#     coef(stats::lm(
#       D2 ~ mean_bin + I(mean_bin ^ 2) + I(mean_bin ^ 3),
#       weights = 1 / eD2,
#       data = est_df %>% filter(eD2 != 0) %>% arrange(bin_nr)
#     ))
#   mean_bin = ests$mean_bin
#   est_df_Lang = data.frame(
#     bin_nr = 1:length(mean_bin),
#     mean_bin = mean_bin,
#     D1XXestimated_raw = ests$D1,
#     D1XXerr = ests$eD1,
#     D2XXestimated_raw = ests$D2,
#     D2XXerr = ests$eD2,
#     D1XXestimated_fit = poly3(
#       mean_bin,
#       a0 = estD1[1],
#       a1 = estD1[2],
#       a2 = estD1[3],
#       a3 = estD1[4]
#     ),
#     D1XXtheoretical = poly3(
#       mean_bin,
#       a0 = d10,
#       a1 = d11,
#       a2 = d12,
#       a3 = d13
#     ),
#     D2XXestimated_fit = poly3(
#       mean_bin,
#       a0 = estD2[1],
#       a1 = estD2[2],
#       a2 = estD2[3],
#       a3 = 0
#     ),
#     D2XXtheoretical = poly3(
#       mean_bin,
#       a0 = d20,
#       a1 = d21,
#       a2 = d22,
#       a3 = 0
#     )
#   ) %>%
#     tidyr::pivot_longer(!c(bin_nr, mean_bin),
#                         names_to = "type",
#                         values_to = "value") %>%
#     tidyr::separate(type, into = c("type_D", "type_est"), sep = "XX") %>%
#     group_by(bin_nr, type_D) %>%
#     mutate(err = ifelse(type_est == "estimated_raw", value[type_est == "err"], NA)) %>%
#     ungroup() %>% filter(type_est != "err")
#   return(est_df_Lang)
# }


est_D_Carp <- function(Ux,
                       sf,
                       D,
                       stabs,
                       Tstep,
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
    Tstep = Tstep,
    bw_sd = bw_sd
  )
  
  # Create an interpolation vector specifying points on the x-axis. Exclude the outer edges of the data and make sure the interpolation vector falls within the range of the drift and diffusion function by a margin
  # nearest = .1
  # xvec = seq(-nearest * floor(abs(min(DD$D1s$x)) / nearest),
  #            nearest * floor(abs(max(DD$D1s$x)) / nearest),
  #            length.out = interpol_steps)
  xvec = seq(min(DD$D1s$x, na.rm = TRUE) - add_to_x, #added na.rm = TRUE
             max(DD$D1s$x, na.rm = TRUE) + add_to_x,
             length.out = interpol_steps)
  
  # Get (effective) potential
  Pot = get_potential(
    D1s = DD$D1s,
    D2s = DD$D2s,
    xeq = DD$xeq,
    interpol_steps = interpol_steps,
    xvec = xvec
  )
  
  Theoretical_df = get_theoretical_D(D, min_x = min(xvec), max_x = max(xvec), interpol_steps = interpol_steps) # KCE
  
  # Only if there are two stable equilibria on the outside and on unstable equilibrium on the inside, compute exit time
  # if (identical(stabs$stab_fps, c("stable", "unstable", "stable"))) {
  # Get exit time using theoretical drift and diffusion
  TheoreticalExitTime = get_exit_time(
    Theoretical_df %>% dplyr::filter(variable == "drift") %>% dplyr::rename(y = value) %>% dplyr::arrange(x) %>% dplyr::select(-variable),
    Theoretical_df %>% dplyr::filter(variable == "diffusion")  %>% dplyr::rename(y = value)  %>% dplyr::arrange(x) %>% dplyr::select(-variable),
    xeq = stabs$fps,
    # theoretical fixed points
    xvec = xvec,
    interpol_steps = interpol_steps
  ) # interpolation vector specifying points on the x-axis
  # } else{
  #   TheoreticalExitTime = list(Theoretical_df = data.frame())
  # }
  
  # Only if there are two stable equilibria on the outside and on unstable equilibrium on the inside, compute exit time
  if (identical(Pot$stab_xeq, c("stable", "unstable", "stable"))) {
    
    # Exit time calculated with estimated drift and diffusion
    EstimatedExitTime = get_exit_time(
      D1s = DD$D1s,
      D2s = DD$D2s,
      xeq = DD$xeq, # fixed points estimated from D1
      xvec = xvec, # interpolation vector specifying points on the x-axis
      interpol_steps = interpol_steps
    )
    
    fp_df = rbind(data.frame(xintercept = stabs$fps,
                             color = stabs$stab_fps, source = "Theoretical"),
                  data.frame(xintercept = DD$xeq,
                             color = Pot$stab_xeq, source = "Estimated")) %>%
      mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
    
  } else if(identical(Pot$stab_pot_xeq, c("stable", "unstable", "stable"))) {
    
    # Exit time calculated with estimated drift and diffusion
    EstimatedExitTime = get_exit_time(
      D1s = DD$D1s,
      D2s = DD$D2s,
      xeq = Pot$pot_xeq, # fixed points estimated from potential
      xvec = xvec, # interpolation vector specifying points on the x-axis
      interpol_steps = interpol_steps
    )
    
    fp_df = rbind(data.frame(xintercept = stabs$fps,
                             color = stabs$stab_fps, source = "Theoretical"),
                  data.frame(xintercept = Pot$pot_xeq,
                             color = Pot$stab_pot_xeq, source = "Estimated")) %>%
      mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
  } else {
    EstimatedExitTime = list(ET_df = data.frame())
    
    fp_df = rbind(data.frame(xintercept = stabs$fps,
                             color = stabs$stab_fps, source = "Theoretical"),
                  data.frame(xintercept = DD$xeq,
                             color = Pot$stab_pot_xeq, source = "Estimated")) %>%
      mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
  }
  
  # Format results (estimated D1 & D2, effective potential)
  est_df_Rinn = data.frame(
    x = DD$bin.mid,
    D1 = DD$D1s$y,
    D2 = DD$D2s$y
  ) %>% tidyr::gather(variable, value,-x) %>%
    dplyr::mutate(err = NA)
  Pot_df = data.frame(
    x = Pot$xvec,
    negPF = Pot$negPF,
    EPF = c(NA, Pot$EPF),
    EPF_err = c(NA, Pot$EPF_err),
    drift = Pot$drift,
    diffusion = Pot$diff
  ) %>% dplyr::rename(potential = negPF) %>% tidyr::gather(variable, value,-c(x, EPF_err)) %>%
    dplyr::rename(err = EPF_err) %>% dplyr::mutate(err = ifelse(variable == "EPF", NA, err))
  
  # Format results
  compl_df = rbind(Theoretical_df %>% mutate(source = "Theoretical", err = NA),
                   rbind(
                     # est_df_Carp,
                         Pot_df) %>% mutate(source = "Estimated"),
                   TheoreticalExitTime$ET_df %>% mutate(source = "Theoretical"),
                   EstimatedExitTime$ET_df %>% mutate(source = "Estimated")
  ) %>% mutate(source = factor(source, levels = c("Theoretical", "Estimated"), ordered = T))
  
  alt_metrics_potential = get_resilience_potential(fp_df = fp_df, compl_df = compl_df)
  alt_metrics_prob_dist = get_resilience_prob_dist(Ux = Ux, bw_sd = bw_sd)
  
  return(
    list(
      DD = DD, 
      est_df_Rinn = est_df_Rinn,
      compl_df = compl_df,
      fp_df = fp_df,
      xeq = DD$xeq,
      alt_metrics_potential = alt_metrics_potential,
      alt_metrics_prob_dist = alt_metrics_prob_dist,
      # est_df_Carp = est_df_Carp,
      # Pot_df = Pot_df,
      stab_xeq = Pot$stab_xeq,
      Dratio = Pot$Dratio,
      D1D2 = Pot$D1D2,
      D1D2adj = Pot$D1D2adj,
      # Theoretical_ET_df = Theoretical_ET_df,
      # ET_df = ET_df,
      theo_meanETl = TheoreticalExitTime$meanETl,
      theo_meanETr = TheoreticalExitTime$meanETr,
      theo_atolL = TheoreticalExitTime$atolL,
      theo_atolR = TheoreticalExitTime$atolR,
      theo_diagn_ETL = TheoreticalExitTime$diagn_ETL,
      theo_diagn_ETR = TheoreticalExitTime$diagn_ETR,
      meanETl = EstimatedExitTime$meanETl,
      meanETr = EstimatedExitTime$meanETr,
      # nL = EstimatedExitTime$nL,
      # nR = EstimatedExitTime$nR,
      atolL = EstimatedExitTime$atolL,
      atolR = EstimatedExitTime$atolR,
      diagn_ETL = EstimatedExitTime$diagn_ETL,
      diagn_ETR = EstimatedExitTime$diagn_ETR
    )
  )
}

# no theoretical (for bootstrapping)
# est_D_Carp_bootstrap <- function(Ux,
#                        sf,
#                        D,
#                        stabs,
#                        Tstep,
#                        ntau,
#                        bins,
#                        bw_sd,
#                        interpol_steps,
#                        add_to_x = .5) {
#   # Estimate drift and diffusion
#   DD = apply_Langevin1D_adapted(
#     Ux = Ux,
#     bins = bins,
#     ntau = ntau,
#     sf = sf,
#     bin_min = 50,
#     Tstep = Tstep,
#     bw_sd = bw_sd
#   )
# 
#   # Create an interpolation vector specifying points on the x-axis. Exclude the outer edges of the data and make sure the interpolation vector falls within the range of the drift and diffusion function by a margin
#   # nearest = .1
#   # xvec = seq(-nearest * floor(abs(min(DD$D1s$x)) / nearest),
#   #            nearest * floor(abs(max(DD$D1s$x)) / nearest),
#   #            length.out = interpol_steps)
#   xvec = seq(min(DD$D1s$x, na.rm = TRUE) - add_to_x, #added na.rm = TRUE
#              max(DD$D1s$x, na.rm = TRUE) + add_to_x,
#              length.out = interpol_steps)
# 
#   # Get (effective) potential
#   Pot = get_potential(
#     D1s = DD$D1s,
#     D2s = DD$D2s,
#     xeq = DD$xeq,
#     interpol_steps = interpol_steps,
#     xvec = xvec
#   )
# 
#   # Theoretical_df = get_theoretical_D(D, min_x = min(xvec), max_x = max(xvec), interpol_steps = interpol_steps) # KCE
# 
#   # Only if there are two stable equilibria on the outside and on unstable equilibrium on the inside, compute exit time
#   # if (identical(stabs$stab_fps, c("stable", "unstable", "stable"))) {
#   # # Get exit time using theoretical drift and diffusion
#   # TheoreticalExitTime = get_exit_time(
#   #   Theoretical_df %>% dplyr::filter(variable == "drift") %>% dplyr::rename(y = value) %>% dplyr::arrange(x) %>% dplyr::select(-variable),
#   #   Theoretical_df %>% dplyr::filter(variable == "diffusion")  %>% dplyr::rename(y = value)  %>% dplyr::arrange(x) %>% dplyr::select(-variable),
#   #   xeq = stabs$fps,
#   #   # theoretical fixed points
#   #   xvec = xvec,
#   #   interpol_steps = interpol_steps
#   # ) # interpolation vector specifying points on the x-axis
#   # } else{
#   #   TheoreticalExitTime = list(Theoretical_df = data.frame())
#   # }
# 
#   # Only if there are two stable equilibria on the outside and on unstable equilibrium on the inside, compute exit time
#   if (identical(Pot$stab_xeq, c("stable", "unstable", "stable"))) {
# 
#     # Exit time calculated with estimated drift and diffusion
#     EstimatedExitTime = get_exit_time(
#       D1s = DD$D1s,
#       D2s = DD$D2s,
#       xeq = DD$xeq, # fixed points estimated from D1
#       xvec = xvec, # interpolation vector specifying points on the x-axis
#       interpol_steps = interpol_steps
#     )
#   } else if(identical(Pot$stab_pot_xeq, c("stable", "unstable", "stable"))) {
# 
#     # Exit time calculated with estimated drift and diffusion
#     EstimatedExitTime = get_exit_time(
#       D1s = DD$D1s,
#       D2s = DD$D2s,
#       xeq = Pot$pot_xeq, # fixed points estimated from potential
#       xvec = xvec, # interpolation vector specifying points on the x-axis
#       interpol_steps = interpol_steps
#     )
#   } else {
#     EstimatedExitTime = list(ET_df = data.frame())
#   }
# 
#   # Format results (estimated D1 & D2, effective potential)
#   est_df_Rinn = data.frame(
#     x = DD$bin.mid,
#     D1 = DD$D1s$y,
#     D2 = DD$D2s$y
#   ) %>% tidyr::gather(variable, value,-x) %>%
#     dplyr::mutate(err = NA)
#   Pot_df = data.frame(
#     x = Pot$xvec,
#     negPF = Pot$negPF,
#     EPF = c(NA, Pot$EPF),
#     EPF_err = c(NA, Pot$EPF_err),
#     drift = Pot$drift,
#     diffusion = Pot$diff
#   ) %>% dplyr::rename(potential = negPF) %>% tidyr::gather(variable, value,-c(x, EPF_err)) %>%
#     dplyr::rename(err = EPF_err) %>% dplyr::mutate(err = ifelse(variable == "EPF", NA, err))
# 
#   # Format results
#   compl_df = rbind(Pot_df %>% mutate(source = "Estimated"),
#                    # TheoreticalExitTime$ET_df %>% mutate(source = "Theoretical"),
#                    EstimatedExitTime$ET_df %>% mutate(source = "Estimated")
#   ) %>% mutate(source = factor(source, levels = c("Estimated"), ordered = T))
#   fp_df = rbind(data.frame(xintercept = stabs$fps,
#                            color = stabs$stab_fps, source = "Theoretical"),
#                 data.frame(xintercept = DD$xeq,
#                            color = Pot$stab_xeq, source = "Estimated")) %>%
#     mutate(source = factor(source, levels = c("Estimated"), ordered = T))
# 
#   # alt_metrics_potential = get_resilience_potential(fp_df = fp_df, compl_df = compl_df)
#   # alt_metrics_prob_dist = get_resilience_prob_dist(fp_df = fp_df,compl_df = compl_df)
# 
#   return(
#     list(
#       DD = DD,
#       est_df_Rinn = est_df_Rinn,
#       compl_df = compl_df,
#       fp_df = fp_df,
#       xeq = DD$xeq,
#       # alt_metrics_potential = alt_metrics_potential,
#       # alt_metrics_prob_dist = alt_metrics_prob_dist,
#       # # est_df_Carp = est_df_Carp,
#       # Pot_df = Pot_df,
#       stab_xeq = Pot$stab_xeq,
#       Dratio = Pot$Dratio,
#       D1D2 = Pot$D1D2,
#       D1D2adj = Pot$D1D2adj,
#       # Theoretical_ET_df = Theoretical_ET_df,
#       # ET_df = ET_df,
#       # theo_meanETl = TheoreticalExitTime$meanETl,
#       # theo_meanETr = TheoreticalExitTime$meanETr,
#       # theo_atolL = TheoreticalExitTime$atolL,
#       # theo_atolR = TheoreticalExitTime$atolR,
#       # theo_diagn_ETL = TheoreticalExitTime$diagn_ETL,
#       # theo_diagn_ETR = TheoreticalExitTime$diagn_ETR,
#       meanETl = EstimatedExitTime$meanETl,
#       meanETr = EstimatedExitTime$meanETr,
#       # nL = EstimatedExitTime$nL,
#       # nR = EstimatedExitTime$nR,
#       # atolL = EstimatedExitTime$atolL,
#       # atolR = EstimatedExitTime$atolR,
#       # diagn_ETL = EstimatedExitTime$diagn_ETL,
#       # diagn_ETR = EstimatedExitTime$diagn_ETR
#     )
#   )
# }

# apply_DDbintau <- function(Ux,
#                            Tstep,
#                            # sf, # added to check scaling
#                            ntau,
#                            bins,
#                            bw_sd) {
#   # Step 3 Carpenter 2022
#   
#   # Set up for binning method
#   bw <-
#     bw_sd * sd(Ux)  # try bandwidth between 0.1*sd and 0.5*sd (according to Carpenter, 2022)
#   
#   # # Note that the bw also cuts of the beginning and end of the x vector. Make sure the bandwidth is not so large that there is no data left:
#   # while ((bins / 10 * bw) > abs(diff(range(Ux)))) {
#   #   bw_sd = bw_sd - .05
#   #   bw <- bw_sd * sd(Ux)  # try between 0.1*sd and 0.5*sd
#   # }
#   DDout = DDbins(Ux, bw, ntau, bins)
#   
#   # Extract smoothed output
#   D1s = DDout[[1]]
#   D2s = DDout[[2]]
#   # sigmas = DDout[[3]] # = sqrt(2*D2)
#   bin.mid = DDout[[4]]
#   
#   # Remove NAs - locations in either the drift or diffusion function where they couldn't be estimated
#   D1s_idx = purrr::map(D1s, function(x) {
#     which(!is.na(x))
#   })
#   D2s_idx = purrr::map(D2s, function(x) {
#     which(!is.na(x))
#   })
#   idx = sort(dplyr::intersect(
#     dplyr::intersect(D1s_idx$x, D1s_idx$y),
#     dplyr::intersect(D2s_idx$x, D2s_idx$y)
#   ))
#   D1s = list(x = D1s$x[idx],
#              # y = (D1s$y[idx])*sf, # changed: multiply by sf for scaling
#              y = D1s$y[idx])
#   D2s = list(x = D2s$x[idx],
#              # y = ((D2s$y[idx])*sf)/2, # changed: multiply by sf and divide by 2 for scaling
#              y = D2s$y[idx])
#   bin.mid = bin.mid[idx]
#   
#   # Find equilibria - where D1 crosses x=0
#   sdrift = sign(D1s$y)
#   dsdrift = c(0,-diff(sdrift))
#   ixeq = which(dsdrift != 0)  # indices of the equilibria
#   xeq = D1s$x[ixeq]
#   
#   print('Equilibria from D1 estimate', quote = F)
#   print(xeq[order(xeq)], quote = F)
#   return(list(
#     D1s = D1s,
#     D2s = D2s,
#     bin.mid = bin.mid,
#     xeq = xeq,
#     ixeq = ixeq
#   ))
# }

Langevin1D_adapted <- function(Ux, #changed from "data"
                               bins, 
                               ntau, #adapted (originally "steps")
                               sf, 
                               bw_sd, #added 
                               bin_min = 50
) {
  # if (!exists("omp_get_num_procs")) {
  #   stop("OpenMP is not available.")
  # }
  # 
  # if (!is.numeric(reqThreads)) {
  #   stop("reqThreads should be numeric")
  # }
  
  steps <- 1:ntau
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
    
    for (n in 1:(length(Ux) - max(steps))) {
      if (Ux[n] >= U[i] && Ux[n] < U[i + 1] && is.finite(Ux[n])) {
        for (s in 1:nsteps) {
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
      M1[i, ] <- sum_m1 / len_step
      M2[i, ] <- sum_m2 / len_step
      M4[i, ] <- sum_m4 / len_step
      
      eM1[i, ] <- sqrt((M2[i, ] - M1[i, ]^2) / len_step)
      eM2[i, ] <- sqrt((M4[i, ] - M2[i, ]^2) / len_step)
      
      coef <- lm(M1[i, ] ~ steps, weights = 1 / eM1[i, ])$coefficients
      D1[i] <- sf * coef[2]
      
      y <- M2[i, ] - (coef[2] * steps)^2
      coef <- lm(y ~ steps, weights = 1 / eM2[i, ])$coefficients
      D2[i] <- sf * coef[2] / 2
      
      coef <- lm(M4[i, ] ~ steps)$coefficients
      D4[i] <- sf * coef[2] / 24
      
      eD1[i] <- sqrt((2 * sf * D2[i] - D1[i]^2) / dens[i])
      eD2[i] <- sqrt((2 * sf * D4[i] - D2[i]^2) / dens[i])
    }
  }
  
  # added: fill values that couldn't be estimated with NA's instead of 0's
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
  
  # # added!!! (not in original Langevin1D)
  # # remove NA's (where D1 and D2 couldn't be estimated) before smoothing
  D1_idx = which(!is.na(D1))  #changed from is.na from the original apply_DDbintau
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
    D1 = D1,
    D1s = D1s, #smoothed
    eD1 = eD1,
    D2 = D2,
    D2s = D2s, #smoothed
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

# apply Langevin1D_adapted
apply_Langevin1D_adapted <- function(Ux,
                                     bins,
                                     ntau,
                                     sf, 
                                     bin_min = 50,
                                     Tstep,
                                     bw_sd = 0.3) {
  
  # Set up for binning method (as per Carpenter, 2022)
  # bw <-
  #   bw_sd * sd(Ux)  # try bandwidth between 0.1*sd and 0.5*sd (according to Carpenter, 2022)
  
  # # Note that the bw also cuts of the beginning and end of the x vector. Make sure the bandwidth is not so large that there is no data left:
  # while ((bins / 10 * bw) > abs(diff(range(Ux)))) {
  #   bw_sd = bw_sd - .05
  #   bw <- bw_sd * sd(Ux)  # try between 0.1*sd and 0.5*sd
  # }
  DDLangout = Langevin1D_adapted(Ux, bins, ntau, sf, bw_sd, bin_min = 50)
  
  # Extract raw output (instead of smoothed)
  # D1 = DDLangout$D1
  # D2 = DDLangout$D2
  # bin.mid = DDLangout$mean_bin
  # 
  # # Remove NAs - locations in either the drift or diffusion function where they couldn't be estimated
  # D1_idx = which(!is.na(D1)) #changed from is.na from the original apply_DDbintau
  # D2_idx = which(!is.na(D2))
  # idx = sort(
  #   dplyr::intersect(D1_idx,D2_idx)
  # )
  # D1 = list(x = bin.mid[idx],
  #           # y = (D1s$y[idx])*sf, # changed: multiply by sf for scaling
  #           y = D1[idx])
  # D2 = list(x = bin.mid[idx],
  #           # y = ((D2s$y[idx])*sf)/2, # changed: multiply by sf and divide by 2 for scaling
  #           y = D2[idx])
  # bin.mid = bin.mid[idx]
  
  # # Extract smoothed output
  D1s = DDLangout$D1s
  D2s = DDLangout$D2s
  # sigmas = DDout[[3]] # = sqrt(2*D2)
  bin.mid = DDLangout$mean_bin

  # Remove NAs and zeros
  find_idx <- function(x, tolerance = 1e-8) {
    na_idx <- which(!is.na(x))
    zero_idx <- which(abs(x) > tolerance)
    idx <- intersect(na_idx, zero_idx)
    return(idx)
  }
  D1sx_idx <- find_idx(D1s$x)
  D1sy_idx <- find_idx(D1s$y)
  D2sx_idx <- find_idx(D2s$x)
  D2sy_idx <- find_idx(D2s$y)
  idx = sort(dplyr::intersect(
    dplyr::intersect(D1sx_idx, D1sy_idx),
    dplyr::intersect(D2sx_idx, D2sy_idx)
  ))

  D1s = list(x = D1s$x[idx],
             # y = (D1s$y[idx])*sf, # changed: multiply by sf for scaling
             y = D1s$y[idx])
  D2s = list(x = D2s$x[idx],
             # y = ((D2s$y[idx])*sf)/2, # changed: multiply by sf and divide by 2 for scaling
             y = D2s$y[idx])
  bin.mid = bin.mid[idx]
  
  # Find equilibria - where D1 crosses x=0
  sdrift = sign(D1s$y)
  dsdrift = c(0,-diff(sdrift))
  ixeq = which(dsdrift != 0)  # indices of the equilibria
  xeq = bin.mid[ixeq]

  # print('Equilibria from D1 estimate', quote = F)
  # print(xeq[order(xeq)], quote = F)
  return(list(
    D1s = D1s,
    D2s = D2s,
    DDLangout = DDLangout,
    bin.mid = bin.mid,
    xeq = xeq,
    ixeq = ixeq
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


get_effective_potential = function(D1s,
                                   D2s,
                                   interpol_steps,
                                   xvec) {
  Dratio = D1s$y / D2s$y
  D1.over.D2 = function(x) {
    # D2 = 0.5*sigma^2(x) for standard Langevin
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

get_potential = function(D1s,
                         D2s,
                         xeq,
                         interpol_steps,
                         xvec) {
  # Step 4 Carpenter 2022
  # Set up for bvp
  
  # plot drift and diffusion (estimated)
  # xvec = seq(avec[1]-0.1,xeq[3]+0.2,length.out=100) # original range
  drift = rep(0, interpol_steps)
  diff = rep(0, interpol_steps)
  for (i in 1:interpol_steps) {
    drift[i] = D1fun(xvec[i], D1s)
    diff[i] = D2fun(xvec[i], D2s)
  }
  
  # Plot potential function (estimated)
  PF = rep(0, interpol_steps)
  for (i in 2:interpol_steps) {
    x0 = xvec[1]
    x1 = xvec[i]
    PF[i] = stats::integrate(
      D1fun,
      lower = x0,
      upper = x1,
      D1s = D1s,
      # rel.tol=.Machine$double.eps^.005,
      subdivisions = max(c(100, interpol_steps *6)) # To prevent error of maximum subdivisions reached
    )$value
  }
  negPF = -1 * PF
  
  # Equilibria from potential estimate
  hills_idx = pracma::findpeaks(negPF)[, 2] # Unstable FP
  valleys_idx = pracma::findpeaks(-negPF)[, 2] # Stable FP
  idx_xeq_PF = c(hills_idx, valleys_idx)
  pot_xeq = sort(xvec[c(idx_xeq_PF)])
  
  # Stability of equilibria (from potential estimate)
  names_xeq_PF = c(rep("unstable", length(hills_idx)), rep("stable", length(valleys_idx)))
  
  # Find closest index for pot_xeq (equilibria from potential estimate)
  stab_pot_xeq = plyr::laply(pot_xeq, function(eq) {
    idx = which.min(abs(xvec - eq))
    idx_match = which.min(abs(idx_xeq_PF - idx))
    return(names_xeq_PF[idx_match])
  })
  
  # Find closest index for xeq (equilibria from D1 estimate)
  stab_xeq = plyr::laply(xeq, function(eq) {
    idx = which.min(abs(xvec - eq))
    idx_match = which.min(abs(idx_xeq_PF - idx))
    return(names_xeq_PF[idx_match])
  })
  
  print('Equilibria from Potential estimate', quote = F)
  print(pot_xeq, quote = F)
  print('Stabilities of equilibria from Potential estimate', quote = F)
  print(stab_pot_xeq, quote = F)
  
  print('Equilibria from D1 estimate', quote = F)
  print(xeq[order(xeq)], quote = F)
  print('Stabilities of equilibria from D1 estimate', quote = F)
  print(stab_xeq, quote = F) # this is the stability of equilibria from D1 estimate!!!
  
  # quartz()
  # par(
  #   mfrow = c(3, 1),
  #   mar = c(4, 4.5, 2, 2) + 0.1,
  #   cex.axis = 2,
  #   cex.lab = 2
  # )
  # plot(
  #   xvec,
  #   drift,
  #   type = 'l',
  #   lwd = 2,
  #   col = 'black',
  #   xlab = 'Chl Level',
  #   ylab = 'Drift'
  # )
  # abline(h = 0, lty = 2)
  # grid()
  # plot(
  #   xvec,
  #   negPF,
  #   type = 'l',
  #   lwd = 2,
  #   col = 'black',
  #   xlab = 'Chl Level',
  #   ylab = '-Potential'
  # )
  # grid()
  # plot(
  #   xvec,
  #   diff,
  #   type = 'l',
  #   lwd = 2,
  #   col = 'black',
  #   xlab = 'Chl Level',
  #   ylab = 'Diffusion'
  # )
  # grid()
  #
  # Calculate and plot effective potential  ------------------------------------------
  EP = get_effective_potential(D1s,
                               D2s,
                               interpol_steps,
                               xvec)
  
  #  Dratio = D1s$y / D2s$y
  # D1.over.D2 = function(x) {
  #   # D2 = 0.5*sigma^2(x) for standard Langevin
  #   yhat = approx(
  #     x = D1s$x,
  #     y = Dratio,
  #     xout = x,
  #     method = 'linear',
  #     rule = 2
  #   )$y
  #   return(yhat)
  # }
  # eprange = range(xvec)
  # xvec.ep = seq(eprange[1], eprange[2], length.out = interpol_steps)
  #
  # # Check D1/D2
  # D1D2 = rep(0, interpol_steps)
  # D1D2adj = rep(0, interpol_steps)
  # for (i in 1:interpol_steps) {
  #   D1D2[i] = D1.over.D2(xvec.ep[i])
  #   D1D2adj[i] = log(0.5 * D2fun(xvec.ep[i], D2s) ^ 2) - D1D2[i]
  # }
  #
  # # quartz()
  # # par(
  # #   mfrow = c(1, 1),
  # #   mar = c(4, 4, 2, 2) + 0.1,
  # #   cex.lab = 1.5,
  # #   cex.axis = 1.5
  # # )
  # # plot(
  # #   xvec.ep,
  # #   D1D2,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'darkred',
  # #   xlab = 'x',
  # #   ylab = 'D1/D2'
  # # )
  # # abline(h = 0, lty = 3, lwd = 2)
  # # grid()
  # #
  # # quartz()
  # # par(
  # #   mfrow = c(1, 1),
  # #   mar = c(4, 4, 2, 2) + 0.1,
  # #   cex.lab = 1.5,
  # # #   cex.axis = 1.5
  # # )
  # # plot(
  # #   xvec.ep,
  # #   D1D2adj,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'darkred',
  # #   xlab = 'x',
  # #   ylab = '-(D1/D2)+log(D2)'
  # # )
  # # abline(h = 0, lty = 3, lwd = 2)
  # # grid()
  #
  # # Plot effective potential function
  # EPF = rep(0, interpol_steps - 1)
  # EPF_err = EPF
  # for (i in 1:(interpol_steps - 1)) {
  #   x0 = xvec.ep[1]
  #   x1 = xvec.ep[i + 1]
  #   xhalf = (x0 + x1) / 2
  #   # These integrate commands return 'bad behavior of integrand'
  #   #integral = integrate(D1.over.D2,lower=x0,upper=x1)$value
  #   #integral = integrate(Vectorize(D1.over.D2),lower=x0,upper=x1)$value
  #   # This function from the cubature package seems to work
  #   hcub = hcubature(f = D1.over.D2,
  #                        lowerLimit = x0,
  #                        upperLimit = x1)
  #   logdiff = log(0.5 * D2fun(xhalf, D2s) ^ 2)  # D2 = 0.5*sigma^2 for standard Langevin
  #   EPF[i] = -1 * hcub$integral + logdiff
  #   EPF_err[i] = hcub$error
  # }
  #
  # # quartz(width = 7, height = 10)
  # # par(
  # #   mfrow = c(2, 1),
  # #   mar = c(4, 4, 2, 2) + 0.1,
  # #   cex.lab = 1.5,
  # #   cex.axis = 1.5
  # # )
  # # plot(
  # #   xvec,
  # #   negPF,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'black',
  # #   xlab = 'Chl Level',
  # #   ylab = 'Potential'
  # # )
  # # plot(
  # #   xvec.ep[2:interpol_steps],
  # #   EPF,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'black',
  # #   xlab = 'Chl Level',
  # #   ylab = 'Effective Potential'
  # # )
  # #
  # # quartz(width = 8, height = 12)
  # # par(
  # #   mfrow = c(3, 1),
  # #   mar = c(4, 4.5, 2, 2) + 0.1,
  # #   cex.axis = 2,
  # #   cex.lab = 2
  # # )
  # # plot(
  # #   xvec,
  # #   drift,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'black',
  # #   xlab = 'Chl Level',
  # #   ylab = 'Drift'
  # # )
  # # abline(h = 0, lty = 3)
  # # abline(v = xeq, lty = 2)
  # # plot(
  # #   xvec,
  # #   diff,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'black',
  # #   xlab = 'Chl Level',
  # #   ylab = 'Diffusion'
  # # )
  # # abline(v = xeq, lty = 2)
  # # plot(
  # #   xvec.ep[2:interpol_steps],
  # #   EPF,
  # #   type = 'l',
  # #   lwd = 2,
  # #   col = 'black',
  # #   xlab = 'Chl Level',
  # #   ylab = 'Effective Potential'
  # # )
  # # abline(v = xeq, lty = 2)
  
  return(
    utils::modifyList(list(
      drift = drift,
      diff = diff,
      xvec = xvec,
      # xeq = xeq, 
      stab_xeq = stab_xeq,
      stab_pot_xeq = stab_pot_xeq, 
      pot_xeq = pot_xeq, 
      negPF = negPF
    ), EP)
  )
}

get_exit_time <- function(D1s,
                          D2s,
                          xeq,
                          xvec,
                          atol = 1e-5,
                          interpol_steps) {
  # Step 4 Carpenter (2022)
  # Calculate Exit Times =======================================================
  
  # function for solving the boundary value problem as two differential equations
  #  for T (col 1) and dT/dx (col 2)
  feval2 = function(x, y, plist) {
    out1 = y[2]
    out2 = -(D1fun(x, plist$D1s) * y[2] + 1) / D2fun(x, plist$D2s)
    return(list(c(out1, out2)))
  }
  
  # feval2new = function(x, y, plist) {
  #   # D2() is sigma not D2 = 0.5*sigma^2
  #   out1 = y[2]
  #   out2 = -(D1fun(x, plist$D1s) * y[2] + 1) / (0.5 * D2fun(x, plist$D2s) ^
  #                                                 2)
  #   return(list(c(out1, out2)))
  # }
  
  
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
          # epsini=.5,
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
  # quartz()
  # plot(
  #   x,
  #   yfinner,
  #   type = 'l',
  #   lwd = 2,
  #   col = 'darkred',
  #   xlab = 'x',
  #   ylab = 'finner'
  # )
  
  # Weights by method in Science paper
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
  
  # Format
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
    list(ET_df = ET_df,
         # ETL = ETL,
         # ETR = ETR,
         diagn_ETL = diagn_ETL,
         diagn_ETR = diagn_ETR,
         # yfinner = yfinner,
         # wts = wts,
         # wtraw_err = WTS$wtraw_err,
         # nL = nL,
         # nR = nR,
         atolL = atolL,
         atolR = atolR,
         meanETl = meanETl,
         meanETr = meanETr
    )
  )
}


get_resilience_potential <- function(fp_df, compl_df) {
  # get distance to basin threshold (theoretical)
  theoretical_fps = fp_df %>%
    filter(source == "Theoretical") %>%
    pull(xintercept)
  
  basin_width_theo_right = abs(theoretical_fps[3] - theoretical_fps[2]) # Dist_thr = |x_{att} - x_{saddle}|
  basin_width_theo_left = abs(theoretical_fps[1] - theoretical_fps[2])
  
  # get distance to basin threshold (estimated)
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
      filter(source == "Estimated") %>%
      pull(xintercept)
    
    basin_width_est_right = abs(estimated_fps[3] - estimated_fps[2]) # Dist_thr = |x_{att} - x_{saddle}|
    basin_width_est_left = abs(estimated_fps[1] - estimated_fps[2])
  } else{
    basin_width_est_right = NULL
    basin_width_est_left = NULL
  }
  
  # get potential depth (theoretical)
  theoretical_potential = compl_df %>%
    filter(variable == "potential", source == "Theoretical")
  # closest_indices = sapply(theoretical_fps, function(fps)
  #   which.min(abs(theoretical_potential$x - fps))) # find the indices of x values closest to fixed points
  
  theo_potentialfun = approxfun(
    x = theoretical_potential$x,
    y = theoretical_potential$value,
    method = "linear",
    rule = 2
  )
  theo_fps_Ux = theo_potentialfun(theoretical_fps)
  
  # theoretical_Ux = theoretical_potential %>%
  #   filter(row_number() %in% closest_indices) %>%
  #   pull(value) %>%
  #   as.numeric() # find U(x_{att}), U(x_{saddle}), and U(x_{att})
  
  potential_depth_theo_right = abs(theo_fps_Ux[3] - theo_fps_Ux[2]) # Pot_dist = |U(x_{att})-U(x_{saddle})|
  potential_depth_theo_left = abs(theo_fps_Ux[1] - theo_fps_Ux[2])
  
  # get potential depth (estimated)
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
      filter(variable == "potential", source == "Estimated")
    # closest_indices = sapply(estimated_fps, function(fps) which.min(abs(estimated_potential$x - fps))) # find the indices of x values closest to fixed points
    
    est_potentialfun = approxfun(
      x = estimated_potential$x,
      y = estimated_potential$value,
      method = "linear",
      rule = 2
    )
    est_fps_Ux = est_potentialfun(estimated_fps)
    
    # # find U(x_{att}), U(x_{saddle}), and U(x_{att})
    # estimated_Ux = estimated_potential %>%
    #   filter(row_number() %in% closest_indices) %>%
    #   pull(value) %>%
    #   as.numeric() # find U(x_{att}), U(x_{saddle}), and U(x_{att})
    
    potential_depth_est_right = abs(est_fps_Ux[3] - est_fps_Ux[2]) # Pot_dist = |U(x_{att})-U(x_{saddle})|
    potential_depth_est_left = abs(est_fps_Ux[1] - est_fps_Ux[2])
  } else{
    potential_depth_est_right = NULL
    potential_depth_est_left = NULL
  }
  
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
  # Get probability distribution of system states
  prob_dist <- density(Ux, n = 1000, bw = bw_sd * sd(Ux))
  
  # Extract the x and y values from the density estimate
  x_values <- prob_dist$x
  y_values <- prob_dist$y
  
  # Find the indice of the saddle point (local minima)
  saddle_point_index <- which(diff(sign(diff(y_values))) == 2) + 1
  
  # Divide probability distribution at saddle point to calculate variance and skewness around the dominant modes
  x_values_est_L <- x_values[1:saddle_point_index]
  x_values_est_R <- x_values[saddle_point_index:length(x_values)]
  
  prob_values_est_L <- y_values[1:saddle_point_index]
  prob_values_est_R <- y_values[saddle_point_index:length(y_values)]
  
  # Calculate variance around dominant modes (left and right)
  est_mode_left = x_values_est_L[which.max(prob_values_est_L)]
  est_var_mode_left =
    sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 2)
  
  est_mode_right = x_values_est_R[which.max(prob_values_est_R)]
  est_var_mode_right =
    sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 2)
  
  # Calculate the skewness around dominant modes (left and right)
  est_skewness_mode_right <-
    sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 3) / (sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 2)) ^ (3 / 2)
  est_skewness_mode_left <-
    sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 3) / (sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 2)) ^ (3 / 2)
  
  return(
    list(
      est_var_mode_left = est_var_mode_left,
      est_var_mode_right = est_var_mode_right,
      est_skewness_mode_left = est_skewness_mode_left,
      est_skewness_mode_right = est_skewness_mode_right
    )
  )
}


# Get variance and skewness around dominant modes of stationary probability distribution
get_resilience_stationary_prob_dist <- function(fp_df, compl_df) {
  
  # get stationary probability distribution (theoretical)
  theoretical_prob_dist = compl_df %>%
    filter(variable == "wts", source == "Theoretical")
  theoretical_saddle = fp_df %>%
    filter(source == "Theoretical") %>%
    slice(2) %>%
    pull(xintercept) # get saddle point
  theoretical_saddle_index = sapply(theoretical_saddle, function(fp)
    which.min(abs(theoretical_prob_dist$x - fp)))
  
  # divide stationary probability distribution at saddle point between basins for calculating variance and skewness around each mode (theoretical)
  theoretical_prob_dist_right =
    theoretical_prob_dist %>% slice(1:theoretical_saddle_index)
  theoretical_prob_dist_left =
    theoretical_prob_dist %>% slice((theoretical_saddle_index + 1):n())
  
  # calculate variance around modes (theoretical)
  x_values_theo_R = theoretical_prob_dist_right$x
  prob_values_theo_R =
    as.numeric(theoretical_prob_dist_right$value)
  
  # find the mode and calculate the variance around the mode (right)
  theo_mode_right =
    x_values_theo_R[which.max(prob_values_theo_R)] # mode of right basin
  theo_var_mode_right =
    sum(prob_values_theo_R * (x_values_theo_R - theo_mode_right) ^ 2) # variance around mode
  
  x_values_theo_L = theoretical_prob_dist_left$x
  prob_values_theo_L = as.numeric(theoretical_prob_dist_left$value)
  
  # find the mode and calculate the variance around the mode (left)
  theo_mode_left = x_values_theo_L[which.max(prob_values_theo_L)]
  theo_var_mode_left =
    sum(prob_values_theo_L * (x_values_theo_L - theo_mode_left) ^ 2)
  
  # calculate the skewness around the mode (theoretical)
  theo_skewness_mode_right <-
    sum(prob_values_theo_R * (x_values_theo_R - theo_mode_right) ^ 3) / (sum(prob_values_theo_R * (x_values_theo_R - theo_mode_right) ^ 2)) ^ (3 / 2)
  theo_skewness_mode_left <-
    sum(prob_values_theo_L * (x_values_theo_L - theo_mode_left) ^ 3) / (sum(prob_values_theo_L * (x_values_theo_L - theo_mode_left) ^ 2)) ^ (3 / 2)
  if (identical(fp_df$color,
                c(
                  "stable",
                  "unstable",
                  "stable",
                  "stable",
                  "unstable",
                  "stable"
                ))) {
    # get left vs. right basins of stationary probability distribution (estimated)
    estimated_prob_dist = compl_df %>%
      filter(variable == "wts", source == "Estimated")
    estimated_saddle = fp_df %>%
      filter(source == "Estimated") %>%
      slice(2) %>%
      pull(xintercept) # get saddle point
    estimated_saddle_index = sapply(estimated_saddle, function(fp)
      which.min(abs(estimated_prob_dist$x - fp)))
    
    estimated_prob_dist_right =
      estimated_prob_dist %>% slice(1:estimated_saddle_index)
    estimated_prob_dist_left =
      estimated_prob_dist %>% slice((estimated_saddle_index + 1):n())
    
    # calculate variance around modes (estimated)
    x_values_est_R = estimated_prob_dist_right$x
    prob_values_est_R = as.numeric(estimated_prob_dist_right$value)
    
    # find the mode and calculate the variance around the mode (right)
    est_mode_right = x_values_est_R[which.max(prob_values_est_R)]
    est_var_mode_right =
      sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 2)
    
    x_values_est_L = estimated_prob_dist_left$x
    prob_values_est_L = as.numeric(estimated_prob_dist_left$value)
    
    # find the mode and calculate the variance around the mode (right)
    est_mode_left = x_values_est_L[which.max(prob_values_est_L)]
    est_var_mode_left =
      sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 2)
    
    # calculate the skewness around the mode (estimated)
    est_skewness_mode_right <-
      sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 3) / (sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 2)) ^ (3 / 2)
    est_skewness_mode_left <-
      sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 3) / (sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 2)) ^ (3 / 2)
    
  }else{
    est_var_mode_left = NULL
    est_var_mode_right = NULL
    est_skewness_mode_left = NULL
    est_skewness_mode_right = NULL
  }
  
  return(
    list(
      theo_var_mode_left = theo_var_mode_left,
      theo_var_mode_right = theo_var_mode_right,
      est_var_mode_left = est_var_mode_left,
      est_var_mode_right = est_var_mode_right,
      theo_skewness_mode_left = theo_skewness_mode_left,
      theo_skewness_mode_right = theo_skewness_mode_right,
      est_skewness_mode_left = est_skewness_mode_left,
      est_skewness_mode_right = est_skewness_mode_right
    )
  )
}


# function for inner integral
finner = function (z, D1s, D2s) {
  fx = D1fun(z, D1s) / (D2fun(z, D2s))
  return(fx)
}

# function for g^2 weights
gg = function(z, D2s) {
  fx = 1 / (D2fun(z, D2s))
  return(fx)
}


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


# get_theoretical_exit_time <-
#   function(D, # d10, d11, d12, d13, d20, d21, d22 for getting theoretical drift and diffusion
#            D1s, # to get domain of x's to calculate analytically theoretical drift and diffusion
#            xeq, # theoretical fixed points
#            xvec # interpolation vector specifying points on the x-axis
#           ) {
#     # Derive theoretical drift and diffusion analytically (over same as range of x's for estimated D1 and D2)
#     # TheoreticalExitTime_df = get_theoretical_D(D, min_x = min(D1s$x), max_x = max(D1s$x), interpol_steps = length(xvec))
#     TheoreticalExitTime_df = get_theoretical_D(D, min_x = min(xvec), max_x = max(xvec), interpol_steps = length(xvec)) # KCE
#
#     # Format theoretical D1's
#     theoretical_D1 <-
#       TheoreticalExitTime_df %>% filter(variable == "drift_analytical") %>% rename(y = value) %>% arrange(x) %>% select(-variable)
#
#     # Format theoretical D2's
#     theoretical_D2 <-
#       TheoreticalExitTime_df %>% filter(variable == "diff_analytical")  %>% rename(y = value)  %>% arrange(x) %>% select(-variable)
#
#     # Get exit time using theoretical drift and diffusion
#     ExitTime = get_exit_time(theoretical_D1,
#                   theoretical_D2,
#                   xeq = xeq,   # theoretical fixed points
#                   xvec = xvec) # interpolation vector specifying points on the x-axis
#     return(ExitTime)
#   }


new_plot_overview <-
  function(out,
           filepath_image
           # plot_t = c(10000 * 100, Inf)[1]
  ) {
    with(out, {
      # Create a title for the plot using LaTeX formatting
      # The title includes the scenario, type_D2, strength_D2, sf, and N 
      title_plot = latex2exp::TeX(
        sprintf(
          "%s; %s $D_2$ (strength: %.2f); $f_s$ = %d, $N$ = %d",
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
      # xlim = max(abs(c(stabs$fps, Ux)))
      
      # pl_Carp with theoretical and estimated on different facet plots 
      # pl_Carp =
      #   est_Carp$compl_df %>%
      #   # Filter the data to include only those variables present in variable.labs
      #   # i.e. drift, diffusion, potential, effective potential, and if estimated, exit time and wts
      #   dplyr::filter(variable %in% names(variable.labs)) %>%
      #   # Convert the 'value' column to numeric and add a new column 'variable_name' by recoding 'variable' column based on variable.labs
      #   dplyr::mutate(
      #     value = as.numeric(as.character(value)),
      #     variable_name = dplyr::recode_factor(variable, !!!variable.labs, .ordered = T)
      #   ) %>%
      #   # Initialize ggplot
      #   ggplot2::ggplot() +
      #   # Add a horizontal line at y=0 as x-axis 
      #   ggplot2::geom_hline(
      #     yintercept = 0,
      #     color = 'grey10',
      #     linetype = "solid",
      #     linewidth = .2,
      #     alpha = 1
      #   ) +
      #   # Add vertical lines at the estimated fixed points 
      #   ggplot2::geom_vline(
      #     data = est_Carp$fp_df,
      #     ggplot2::aes(xintercept = xintercept, color = color),
      #     linetype = 'dashed',
      #     linewidth = lwd_FP,
      #     alpha = .75
      #   ) +
      #   # Define manual color scale for unstable and stable estimated fixed points
      #   ggplot2::scale_color_manual(
      #     name = "",
      #     breaks = c("unstable", "stable", unique(variable.labs)),
      #     values = c(
      #       "unstable" = scales::viridis_pal(option = "inferno")(50)[29],
      #       "stable" = scales::viridis_pal(option = "viridis")(50)[45]
      #     ),
      #     guide = "none"
      #   ) +
      #   ggnewscale::new_scale_color() +
      #   # Add points representing errors
      #   ggplot2::geom_point(
      #     ggplot2::aes(x = x, y = err, color = source),
      #     shape = 4,
      #     color = 'gray40',
      #     size = .5
      #   ) +
      #   # Add points representing estimated variables
      #   ggplot2::geom_point(ggplot2::aes(x = x, y = value, color = source), size = .5) +
      #   # Define manual color scale for the points
      #   ggplot2::scale_color_manual(
      #     name = "",
      #     breaks = c("Theoretical", "Estimated"),
      #     values = c(
      #       "Theoretical" = 'grey10',
      #       "Estimated" = 'red3'
      #     ),
      #     guide = "none"
      #   ) +
      #   # Adjust the x-axis to start at 0
      #   ggplot2::scale_x_continuous(expand = c(0, 0)) +
      #   # Move the y-axis to the right
      #   ggplot2::scale_y_continuous(position = 'right') +
      #   # Create facets based on 'source' and 'variable_name', with free y scales
      #   ggh4x::facet_grid2(
      #     source ~ variable_name,
      #     scales = 'free_y',
      #     independent = 'y',
      #     labeller = label_parsed,
      #     switch = "y",
      #     axes = "all"
      #   ) +
      #   # Remove the legend
      #   ggplot2::theme(legend.position = "none")
      
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
    # sprintf("%s-D2", type_D2), # removed since only constant noise
    sprintf("N%s", N), # added
    sprintf("sf%s", sf), # added
    sprintf("tau%s", ntau), # added
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
    # sprintf("%s-D2", type_D2), # removed since only constant noise
    sprintf("N%s", N), # added
    sprintf("sf%s", sf), # added
    sprintf("tau%s", ntau), # added
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

