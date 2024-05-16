### Assumption checks ###
# Set working directory
filepath_base = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(filepath_base)

# Load necessary libraries and helper functions
library(ggplot2)
library(dplyr)
library(plotrix)
source(file.path(filepath_base, "bootstrap_helper_functions.R"))

# Load example for demo purposes
example <-
  readRDS(
    "Rinn_est_scale_check/2fps-balanced-deepening/constant-D2/1000-N/D2strength0.3000_sf10_N1000_iter0077_step0003_pars-2.00_0.00_2.00_0.00_0.00_0.00_0.30_bins100_ntau3_interpol50_bw0.30.RDS"
  )

### N-step ahead 
bootstrapDDs <- bootstrap_est_D(
  example$est_Carp$DD$DDLangout$D1s,
  example$est_Carp$DD$DDLangout$D2s,
  example$est_Carp$DD$DDLangout$eD1,
  example$est_Carp$DD$DDLangout$eD2,
  example$est_Carp$DD$DDLangout$mean_bin,
  example$bw_sd,
  example$Ux
) # estimate parent model from estimated drift and diffusion coefficients 

#debug
Ux = example$Ux
N = example$N
sf = example$sf
ntau = example$ntau
bins = example$bins
bw_sd = example$bw_sd 
interpol_steps = example$interpol_steps
compl_df = example$est_Carp$compl_df
noise_iter = example$noise_iter
D = example$D
Nstep = 3
bootstrap_n = 2
Npredictions = 5
rm(list = c("Ux", "N", "sf", "ntau", "bins", "bw_sd", "interpol_steps", "compl_df", "noise_iter", "D", "Nstep", "boostrap_n", "Onestep_demo", "actual_timeseries", "i", "j", "last_state", "one_step_pred", "Onestep_predictions", "Onestep_Ux", "seeds", "MSE", "initial_states", "actual_Nstep_states", "actual_timeseries_states", "Nstep_demo", "Nstep_predictions", "residuals", "RMSE", "squared_residuals", "model_timeseries", "Npredictions", "start_states"))

# Nstep ahead function
Nstep_ahead <-
  function(Nstep,
           # bootstrap_n,
           Npredictions, 
           bootstrapDDs,
           Ux,
           N,
           sf,
           noise_iter,
           bins,
           ntau,
           bw_sd) {
    # Nstep_predictions <- matrix(NA, nrow = bootstrap_n + 1, ncol = Nstep)
    # last_state <- Ux[length(Ux)]
    # seeds <- c(seq(50, 149), noise_iter)
    
    # get initial starting states from original timeseries
    start_states <- Ux[round(seq(17, length(Ux)-Nstep, length.out = Npredictions))]
    
    # get actual nstep ahead points from original timeseries 
    actual_Nstep_states <- matrix(NA, nrow = Npredictions, ncol = Nstep)
    for(i in 1:Nstep){
      actual_Nstep_states[,i] <-
        Ux[round(seq(17, length(Ux) - Nstep, length.out = Npredictions)) + i]
    }
    
    # get nstep predictions
    Nstep_predictions <- matrix(NA, nrow = Npredictions, ncol = Nstep)
    for(i in 1:Nstep){
      for(j in 1:Npredictions){
        # print(c("start state is", start_states[j]))
        # set.seed(noise_iter)
        model_timeseries <-
          timeseries1D(
            N = 10,
            sf = 10,
            startpoint = start_states[j], 
            d13 = bootstrapDDs$d13,
            d12 = bootstrapDDs$d12,
            d11 = bootstrapDDs$d11,
            d10 = bootstrapDDs$d10,
            d22 = bootstrapDDs$d22,
            d21 = bootstrapDDs$d21,
            d20 = bootstrapDDs$d20
          )
        
        # print(c("prediction is", model_timeseries[i]))
        Nstep_predictions[j,i] <-
          model_timeseries[i]
      }
    }
    
    # calculate RMSE
    squared_residuals <- (actual_Nstep_states - Nstep_predictions)^2
    RMSE = matrix(NA, ncol = Nstep)
    for(i in 1: Nstep){
      RMSE[,i] <- sqrt(mean(squared_residuals[,i]))
    }
    
    # calculate MAD
    # median_Nstep_pred <- matrix(nrow = Npredictions, ncol = Nstep)
    # MAD = vector()
    # for(i in 1:Nstep){
    #   median_Nstep_pred[,i] <- rep(median(Nstep_predictions[,i]), times = Npredictions)
    #   abs_dev <- abs(Nstep_predictions - median_Nstep_pred)
    #   MAD[i] <- median(abs_dev[,i])
    # }
    
    # calculate median absolute error
    # abs_error <- abs(actual_Nstep_states - Nstep_predictions)
    # MAE = vector()
    # for(i in 1: Nstep){
    #   MAE[i] <- median(abs_error[,i])
    # }
    
    return(
      list(
        Npredictions = Npredictions,
        actual_Nstep_states = actual_Nstep_states,
        Nstep_predictions = Nstep_predictions,
        RMSE = RMSE
        # MAD = MAD
        # MAE = MAE
      )
    )
  }

Nstep_ahead(Nstep = 5,
            # bootstrap_n = 100,
            Npredictions = 15,
            bootstrapDDs = bootstrapDDs,
            Ux = example$Ux,
            N = example$N,
            sf = example$sf,
            noise_iter = example$noise_iter,
            bins = example$bins,
            ntau = example$ntau,
            bw_sd = example$bw_sd)

## Demo
Nstep_demo <- matrix(nrow = length(seq(2, 200, by = 10)), ncol = 3)
for (i in 1:3){
  for(j in seq(2, 200, by = 10)){
    n <- which(seq(2, 200, by = 10) == j)
    Nstep_demo[n,i] <- Nstep_ahead(Nstep = 3,
                                   # bootstrap_n = 100,
                                   Npredictions = j,
                                   bootstrapDDs = bootstrapDDs,
                                   Ux = example$Ux,
                                   N = example$N,
                                   sf = example$sf,
                                   noise_iter = example$noise_iter,
                                   bins = example$bins,
                                   ntau = example$ntau,
                                   bw_sd = example$bw_sd)$MAE[i]
  }
}

Nstep_demo_df <- data.frame(Npredictions = rep(seq(2, 200, by = 10), times = 3), MAE = c(Nstep_demo[,1], Nstep_demo[,2], Nstep_demo[,3]), step = c(rep("step1", times = 20), rep("step2", times = 20), rep("step3", times = 20)))

## Facet plot
ggplot(Nstep_demo_df, aes(x = Npredictions, y = MAE)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ step, scales = "fixed", ncol = 1) +
  labs(title = "N = 1000, sf = 10, tau = 3; no set.seed", x = "Npredictions", y = "MAE") +
  theme_bw()

## Check effect of setting different seeds 
Nstep_ahead_different_seeds <- function(Nstep, Npredictions, Nseeds, bootstrapDDs, Ux, N, sf, noise_iter, bins, ntau, bw_sd) {
  # Initialize a 3D array to store the predictions
  Nstep_predictions <- array(NA, dim = c(Nseeds, Npredictions, Nstep))
  
  # Get initial starting states from original timeseries
  start_states <- Ux[round(seq(17, length(Ux)-Nstep, length.out = Npredictions))]
  
  # Get actual nstep ahead points from original timeseries 
  actual_Nstep_states <- matrix(NA, nrow = Npredictions, ncol = Nstep)
  for(i in 1:Nstep){
    actual_Nstep_states[,i] <- Ux[round(seq(17, length(Ux) - Nstep, length.out = Npredictions)) + i]
  }
  
  # Iterate over the seeds
  for (s in 1:Nseeds) {
    # Set the seed
    set.seed(s)
    # Get nstep predictions
    for(i in 1:Nstep){
      for(j in 1:Npredictions){
        model_timeseries <- timeseries1D(
          N = 10,
          sf = 10,
          startpoint = start_states[j], 
          d13 = bootstrapDDs$d13,
          d12 = bootstrapDDs$d12,
          d11 = bootstrapDDs$d11,
          d10 = bootstrapDDs$d10,
          d22 = bootstrapDDs$d22,
          d21 = bootstrapDDs$d21,
          d20 = bootstrapDDs$d20
        )
        
        Nstep_predictions[s, j, i] <- model_timeseries[i]
      }
    }
  }
  
  # Calculate 95% confidence intervals
  CI_lower <- array(NA, dim = c(Npredictions, Nstep))
  CI_upper <- array(NA, dim = c(Npredictions, Nstep))
  for(i in 1:Nstep){
    for(j in 1:Npredictions){
      CI_lower[j, i] <- quantile(Nstep_predictions[, j, i], 0.025)  # 2.5 percentile
      CI_upper[j, i] <- quantile(Nstep_predictions[, j, i], 0.975)  # 97.5 percentile
    }
  }
  
  return(
    list(
      # Npredictions = Npredictions,
      actual_Nstep_states = actual_Nstep_states,
      Nstep_predictions = Nstep_predictions,
      CI_lower = CI_lower,
      CI_upper = CI_upper
      # RMSE = RMSE
    )
  )
}

## Demo
Nstep = 8
Npredictions = 15
Nseeds = 1000
example$noise_iter

demo <- Nstep_ahead_different_seeds(
  Nstep = Nstep,
  Npredictions = Npredictions,
  Nseeds = Nseeds,
  bootstrapDDs = bootstrapDDs,
  Ux = example$Ux[1:200],
  # Ux = example$Ux,
  N = example$N,
  sf = example$sf,
  noise_iter = example$noise_iter,
  bins = example$bins,
  ntau = example$ntau,
  bw_sd = example$bw_sd
)

## Plot
plot(
  1:length(example$Ux[1:200]),
  example$Ux[1:200],
  type = "l",
  # cex = 0.5,
  main = "Timeseries with Nstep prediction CI's",
  xlab = "Time",
  ylab = "Ux"
)

# Get the indices for the actual N-step states
start_state_indices <- round(seq(17, length(example$Ux[1:200]) - Nstep, length.out = Npredictions))

# Add the confidence intervals to the plot
for (i in 1:Nstep) {
  for (j in 1:Npredictions) {
    # Get the current index
    index <- start_state_indices[j] + i
    # print(index)
    
    # Add a line for the lower confidence interval
    lines(c(index, index), c(demo$CI_lower[j, i], demo$CI_upper[j, i]), col = "red")
    # print(demo$CI_lower[j, i])
    # print(demo$CI_upper[j, i])
    
    # Add a point for the actual N-step state
    # points(index, demo$actual_Nstep_states[j, i], col = "blue", cex = 0.5)
    
    # Add a point for Nstep predictions with same set.seed as original timeseries
    points(index, demo$Nstep_predictions[example$noise_iter, j, i], col = "blue", cex = 0.5)
  }
}

### Plot using ggplot
# Get the indices for the actual N-step states
Nstep = 8
Npredictions = 15
start_state_indices <- round(seq(17, length(example$Ux[1:200]) - Nstep, length.out = Npredictions))

# Create a data frame
df <- data.frame(
  Nstep_Index = rep(start_state_indices, each = Nstep) + rep(1:Nstep, times = Npredictions),
  Value = as.vector(t(demo$actual_Nstep_states)),
  CI_Lower = as.vector(t(demo$CI_lower)),
  CI_Upper = as.vector(t(demo$CI_upper))
)

# Plot using ggplot
ggplot(df, aes(x = Nstep_Index, y = Value)) +
  geom_line() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, color = "red") +
  labs(x = "Time", y = "Ux", title = "Time Series with Confidence Intervals") +
  theme_bw()

### Plot using plotCI
plotCI(df$Nstep_Index,as.vector(t(demo$actual_Nstep_states)), ui = as.vector(t(demo$CI_upper)), li = as.vector(t(demo$CI_lower)), main = "Timeseries with Nstep confidence intervals", xlab = "", ylab = "Ux", col = "red", scol = "red", pch = 1, cex = 0.5)

### Moving window
library(zoo)
#debug
Ux = example$Ux
nwindow = 100
length(rol_mean)

# function to check stationarity via moving window
stationary_check <- function(Ux, sf, nwindow, lags = NULL) {
  # calculate rolling mean and standard deviation
  # rol_mean <- rollmean(Ux, 50, fill = NA)
  # rol_std <- rollapply(Ux, nwindow, sd, fill = NA)
  rol_mean <- stats::filter(Ux, rep(1/nwindow, nwindow), sides=2)
  rol_std <- sqrt(stats::filter(Ux^2, rep(1/nwindow, nwindow), sides=2) - rol_mean^2)
  
  # plot timeseries with rolling window statistics
  data <- data.frame(Time = 1:length(Ux)/ sf / 60, 
                     # Timeseries = Ux, 
                     RollingMean = rol_mean, 
                     RollingStd = rol_std)
  
  rol_plot <- ggplot(data, aes(x = Time)) +
    geom_line(aes(y = RollingMean, colour = "rolling mean")) +
    geom_line(aes(y = RollingStd, colour = "rolling stdev")) +
    labs(title = "moving window", x = "time (min)", y = "", ) +
    theme_minimal() +
    scale_colour_manual(
      "",
      breaks = c("rolling mean", "rolling stdev"),
      values = c("deeppink3", "darkblue")
    )
  
  # Augmented Dickey-Fuller Test
  adf_test <- adf.test(Ux, alternative = "stationary")
  
  # autocorrelation structure/function for stationary time series 
  acf_stationarity <- acf(Ux, lag.max = length(Ux), plot = FALSE) 
  plot(acf_stationarity , main = 'autocorrelation function', cex.main = 1.2,  
       cex.lab = 1.1) # ACF should decay rapidly for stationary timeseries, indicating a lack of long-term dependencies
  
  print(adf_test)
  print(rol_plot)
}

# demo
stationary_check(example$Ux, example$sf, nwindow = 100, lags = 1000)

#### Check stationary (augmented Dickey-Fuller (ADF) test)
?adf.test
tseries::adf.test(example$Ux)

#### Check Markov property
spgs::markov.test(example$Ux, type = "lb.test")
spgs::markov.test(Ux, type = "rank.test")


#### "Recycle bin"
# # n step ahead (recursive; doesn't work anymore due to edits)
# Nstep_ahead <-
#   function(Nstep,
#            # bootstrap_n,
#            Npredictions, 
#            bootstrapDDs,
#            Ux,
#            N,
#            sf,
#            noise_iter,
#            bins,
#            ntau,
#            bw_sd) {
#     Nstep_predictions <- matrix(NA, nrow = Npredictions, ncol = Nstep)
#     # Nstep_predictions <- matrix(NA, nrow = bootstrap_n + 1, ncol = Nstep)
#     # last_state <- Ux[length(Ux)]
#     # seeds <- c(seq(50, 149), noise_iter)
#     
#     # set.seed(noise_iter)
#     initial_states <- Ux[round(seq(17, length(Ux)-3, length.out = 100))]
#     
#     actual_timeseries_states <- list()
#     for(i in 1:Nstep){
#       actual_timeseries_states[[i]] <-
#         Ux[round(seq(17, length(Ux) - 3, length.out = 100)) + i]
#       }
#     
#     model_timeseries <-
#       timeseries1D(
#         N = N*sf,
#         sf = sf,
#         startpoint = Ux[1], # note: this is the first point of Ux, not initial state of Ux
#         d13 = bootstrapDDs$d13,
#         d12 = bootstrapDDs$d12,
#         d11 = bootstrapDDs$d11,
#         d10 = bootstrapDDs$d10,
#         d22 = bootstrapDDs$d22,
#         d21 = bootstrapDDs$d21,
#         d20 = bootstrapDDs$d20
#       )
#     
#     # print(c("actual timeseries states are:", actual_timeseries_states))
#     
#     for (i in 1:Nstep) {
#       # NstepDDLangout <- apply_Langevin1D_adapted(
#       #   Ux = Ux,
#       #   bins = bins,
#       #   ntau = ntau,
#       #   sf = sf,
#       #   bin_min = 50,
#       #   Tstep = 1:length(as.vector(Ux)) / sf,
#       #   bw_sd = bw_sd
#       # )$DDLangout
#       # 
#       # NstepDDs <-
#       #   bootstrap_est_D(
#       #     NstepDDLangout$D1s,
#       #     NstepDDLangout$D2s,
#       #     NstepDDLangout$eD1,
#       #     NstepDDLangout$eD2,
#       #     NstepDDLangout$mean_bin,
#       #     bw_sd = bw_sd,
#       #     Ux = Ux)
#       # 
#       # print(unlist(NstepDDs))
#       
#       for (j in 1:(bootstrap_n+1)) {
#         print(c("last state is", last_state))
#         
#         set.seed(seeds[j])
#         # Nstep_Ux <-
#         #   do.call(timeseries1D, utils::modifyList(NstepDDs, list(
#         #     N = 10,
#         #     sf = 10,
#         #     startpoint = last_state
#         #   )))
#         Nstep_Ux <-
#           timeseries1D( 
#             N = 10,
#             sf = 10,
#             startpoint = last_state,
#             d13 = NstepDDs$d13,
#             d12 = NstepDDs$d12,
#             d11 = NstepDDs$d11,
#             d10 = NstepDDs$d10,
#             d22 = NstepDDs$d22,
#             d21 = NstepDDs$d21,
#             d20 = NstepDDs$d20
#           )
#         
#         next_prediction <- Nstep_Ux[i]
#         print(c("the next point estimate is", next_prediction))
#         Nstep_predictions[j, i] <- next_prediction
#         
#         if(j == bootstrap_n+1){
#           Ux <-
#             c(Ux, next_prediction) #update Ux with prediction
#           last_state <- Ux[length(Ux)]
#           print(c("length of new_Ux is", length(Ux)))
#           print(c("new last state is", last_state))
#         }
#       }
#     }
#     return(
#       list(
#         actual_timeseries_states = actual_timeseries_states,
#         Nstep_predictions = Nstep_predictions
#       )
#     )
#   }
# 
# Nstep_demo <- Nstep_ahead(Nstep = 3,
#                bootstrap_n = 100,
#                bootstrapDDs = bootstrapDDs,
#                Ux = example$Ux,
#                N = example$N,
#                sf = example$sf,
#                noise_iter = example$noise_iter,
#                bins = example$bins,
#                ntau = example$ntau,
#                bw_sd = example$bw_sd) 
# 
# # Calculate confidence intervals
# lower_ci <- apply(Nstep_demo$Nstep_predictions[-nrow(Nstep_demo$Nstep_predictions), ], 2, function(x) quantile(x, probs = 0.025))
# upper_ci <- apply(Nstep_demo$Nstep_predictions[-nrow(Nstep_demo$Nstep_predictions), ], 2, function(x) quantile(x, probs = 0.975))
# 
# # plot
# Nstep = 3
# for(i in 1:Nstep) {
#   hist(
#     Nstep_demo$Nstep_predictions[-nrow(Nstep_demo$Nstep_predictions), i],
#     breaks = 30,
#     prob = TRUE,
#     col = "lightblue",
#     xlab = paste(c(
#       i, "-step ahead predictions"
#     )),
#     ylab = "Density",
#     main = "Distribution of bootstrapped estimates"
#   )
#   lines(
#     density(Nstep_demo$Nstep_predictions[-nrow(Nstep_demo$Nstep_predictions), i]),
#     col = "darkblue",
#     lwd = 2
#   )
#   
#   # actual timeseries state
#   abline(
#     v = Nstep_demo$actual_timeseries_states[i],
#     col = "red",
#     lwd = 2
#   )
#   # point estimate
#   abline(
#     v = Nstep_demo$Nstep_predictions[nrow(Nstep_demo$Nstep_predictions), i],
#     col = "purple",
#     lwd = 2
#   )
#   # vertical lines for confidence interval
#   abline(v = lower_ci[i],
#          col = "lightblue",
#          lty = 2)
#   abline(v = upper_ci[i],
#          col = "lightblue",
#          lty = 2)
#   legend(
#     "topright",
#     legend = c("Actual state", "Point Estimate", "95% CI"),
#     col = c("red", "purple", "lightblue"),
#     lty = c(1, 1, 2),
#     lwd = c(2, 2, 1),
#     x.intersp = 0.8,
#     y.intersp = 0.8,
#     xjust = 2,
#     yjust = 2
#   )
# }

# # One-step ahead
# Onestep_ahead <- function(bootstrapDDs, bootstrap_n, Ux, N, sf, noise_iter){
#   Onestep_predictions <- vector()
#   last_state <- Ux[length(Ux)]
#   seeds <- c(noise_iter, seq(50,149))
#   
#   set.seed(noise_iter)
#   # actual_timeseries_state <- do.call(timeseries1D, utils::modifyList(D, list(
#   #   N = 10,
#   #   sf = 10,
#   #   startpoint = last_state
#   # )))[1]
#   actual_timeseries_state <-
#     timeseries1D(
#       N = 10,
#       sf = 10,
#       startpoint = last_state,
#       d13 = bootstrapDDs$d13,
#       d12 = bootstrapDDs$d12,
#       d11 = bootstrapDDs$d11,
#       d10 = bootstrapDDs$d10,
#       d22 = bootstrapDDs$d22,
#       d21 = bootstrapDDs$d21,
#       d20 = bootstrapDDs$d20
#     )[1]
#   
#   actual_timeseries_state <-
#     timeseries1D(
#       N = 10,
#       sf = 10,
#       startpoint = example$Ux[218],
#       d13 = bootstrapDDs$d13,
#       d12 = bootstrapDDs$d12,
#       d11 = bootstrapDDs$d11,
#       d10 = bootstrapDDs$d10,
#       d22 = bootstrapDDs$d22,
#       d21 = bootstrapDDs$d21,
#       d20 = bootstrapDDs$d20
#     )[1]
#   
#   for(i in 1:(bootstrap_n+1)){
#     
#     set.seed(seeds[i])
#     Onestep_Ux <- do.call(timeseries1D, utils::modifyList(bootstrapDDs, list(
#       N = 10,
#       sf = 10,
#       startpoint = last_state
#     )))
#     
#     if (any(is.na(Onestep_Ux) | is.infinite(Onestep_Ux))) {
#       print("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
#     } else {
#       print("Timeseries simulated successfully")
#     }
#     
#     one_step_pred <- Onestep_Ux[1]
#     Onestep_predictions[i] <- one_step_pred
#   }
#   
#   return(list(actual_timeseries_state = actual_timeseries_state, Onestep_predictions = Onestep_predictions))
# }
# 
# Onestep_demo <-
#   Onestep_ahead(
#     bootstrapDDs = bootstrapDDs,
#     bootstrap_n = 100,
#     Ux = example$Ux,
#     N = example$N,
#     sf = example$sf,
#     noise_iter = example$noise_iter
#   )
# 
# # plot CI
# confidence_interval <- quantile(Onestep_demo$Onestep_predictions, c(0.025, 0.975))
# hist(Onestep_demo$Onestep_predictions[-1], breaks = 30, prob = TRUE, col = "lightblue", 
#      xlab = "Bootstrapped 1-step ahead predictions", ylab = "Density", main = "Distribution of bootstrapped estimates")
# lines(density(Onestep_demo$Onestep_predictions[-1]), col = "darkblue", lwd = 2)
# 
# # actual timeseries state
# abline(v = Onestep_demo$actual_timeseries_state, col = "red", lwd = 2)
# # point estimate
# abline(v = Onestep_demo$Onestep_predictions[1], col = "purple", lwd = 2)
# # vertical lines for confidence interval
# abline(v = confidence_interval[1], col = "lightblue", lty = 2)
# abline(v = confidence_interval[2], col = "lightblue", lty = 2)
# legend(
#   "topright",
#   legend = c("Actual state", "Point Estimate", "95% CI"),
#   col = c("red", "purple", "lightblue"),
#   lty = c(1, 1, 2),
#   lwd = c(2, 2, 1),
#   x.intersp = 0.8,
#   y.intersp = 0.8,
#   xjust = 2,
#   yjust = 2
# )


