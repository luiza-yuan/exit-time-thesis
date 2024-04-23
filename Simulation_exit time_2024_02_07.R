# Application of Carpenter's (2022) exit time estimation to simulated data using Rinn's (2016) Langevin package
# Kyra Evers, University of Amsterdam, 17.01.2024

######### Modifications by Luiza Yuan, 29.01.2024 #########
# (1) modified get_exit_time function ('yini' and 'yend' set up for bvpSolve for left vs. right basins)
# (2) modified get_D function(from specifying alpha/beta to specifying d10,d11,d12,d13)
# (3) modified plotting function (new_plot_overview function)
# (4) added get_theoretical_exit_time function (included in est_D_Carp function)
# (5) notes at the bottom regarding scaling and theoretical exit time estimation issues (still in progress)
#################################

######### Modifications by Kyra Evers, 07.02.2024 #########
# (1) helper_function.R
# (2) plot formatting changed
#################################

######### Modifications by Luiza Yuan, 09.02.2024 #########
# (1) changes to forloop (e.g. added datagen, nr_steps_bif, modified some arguments)
#################################

######### Modifications by Luiza Yuan, 05.03.2024 & 15.03.2024 #########
# (1) alternative resilience metrics calculation: (1) basin width (2) potential depth (3) variance around modes of probability distribution and (4) skewness around modes of probability distribution
#################################

######### Modifications by Luiza Yuan, 10.03.2024 #########
# (1) Setup github 
# (2) These are comments added from RStudio
#################################

######### Modifications by Luiza Yuan, 20.02.2024 #########
# (1) parallel processing 
#################################

rm(list = ls())
graphics.off()

# Create necessary directories
filepath_base = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(filepath_base)
filepath_figs = file.path(filepath_base, "Rinn_figs_scale_check") #figs_theoreticalET_scalecheck")
filepath_est = file.path(filepath_base, "Rinn_est_scale_check")# "est_theoreticalET_scalecheck")
if (!dir.exists(filepath_figs)) {
  dir.create(filepath_figs, recursive = T)
}

library(bvpSolve)
library(cubature)
library(stats)
library(Langevin)
library(boot)
library(tseries)
library(spgs)
library(dplyr)
library(ggplot2)
library(parallel) # added 05.03.2024
library(doParallel) # added 05.03.2024
library(foreach)
# library(tuneR)
library(cowplot)
library(ggnewscale)
library(latex2exp)

# source("ExitTime_BinMethod_PeterLakeExample-main/DDbintau.R")
# source(file.path(dirname(filepath_base), "ExitTime_BinMethod_PeterLakeExample-main/DDbintau.R"))
# source(file.path(filepath_base, "helper_functions.R")) # attention: make sure that it is the version updated 20.03.2024
source(file.path(filepath_base, "helper_functions_Rinn.R"))

# Choose parameters to loop through
forloop = tidyr::expand_grid(
  # datagen = "Langevin",
  nr_steps_bif = 5, #length.out for deepening, asymmetry, etc.
  type_D2 = c("constant"),
  #"quadratic",
  scenario = c("2fps-balanced-deepening", "right-fp-gains-dominance"), #"left-fp-gains-dominance"
  strength_D2 = c(.3),
  sf = c(10), #c(10, 100),
  N = c(600), #c(500, 100000),
  bins = c(60), #c(30, 40, 100),
  interpol_steps = 50,# c(50, 100, 500),
  ntau = c(3), # c(3, 5, 10),
  bw_sd = .3,
  #10000
  noise_iter = c(2), #c(1:5)
  # noise_iter_concat_times = 5
) %>% purrr::transpose() %>% unique()

# Check each item in the forloop list and keep only the ones where N*sf = bins*100
filtered_forloop <- list()
for (i in seq_along(forloop)) {
  if (forloop[[i]]$N * forloop[[i]]$sf == forloop[[i]]$bins * 100) {
    filtered_forloop[[i]] <- forloop[[i]]
  }
}
# Remove NULL elements from the list
filtered_forloop <- filtered_forloop[!sapply(filtered_forloop, is.null)]

# Debug
#datagen = "Langevin"
nr_steps_bif = 5

step_idx = 3
for_par=forloop[[3]]
type_D2 = for_par$type_D2
scenario = for_par$scenario
bins = for_par$bins
strength_D2 = for_par$strength_D2
bw_sd = for_par$bw_sd
ntau= for_par$ntau
interpol_steps = for_par$interpol_steps
sf = for_par$sf
N = for_par$N
noise_iter = for_par$noise_iter
add_to_x = .5
rm(list = c("nr_steps_bif", "step_idx", "type_D2", "scenario", "mean_bin", "bw_sd", "ntau", "interpol_steps", "sf", "N", "noise_iter", "add_to_x", "bins", "i", "strength_D2" ,"basin", "D", "DD", "Ds", "est_Carp", "Pot", "pl_Carp", "stabs", "Theoretical_df", "yend", "yini", "xvec", "Ux"))

# noise_iter_concat_times = for_par$noise_iter_concat_times

# # Set up cluster
# cl <- parallel::makeCluster(detectCores() - 1, type = "PSOCK")
# doParallel::registerDoParallel(cl)
# 
# getDoParWorkers() #check how many workers 'foreach' is going to use

# Set up cluster (forking)
cl <- parallel::makeForkCluster(detectCores() - 1) # using forking
doParallel::registerDoParallel(cl)

getDoParWorkers() #check how many workers 'foreach' is going to use

# Packages, functions, variables 
# packages <-
#   c(
#     "bvpSolve",
#     "cubature",
#     "stats",
#     "Langevin",
#     "dplyr",
#     "ggplot2",
#     "parallel",
#     "doParallel",
#     "foreach",
#     "cowplot",
#     "ggnewscale",
#     "latex2exp"
#   )
# variables <-
#   c(
#     "datagen",
#     "nr_steps_bif",
#     "type_D2",
#     "scenario",
#     "strength_D2",
#     "sf",
#     "N",
#     "bins",
#     "interpol_steps",
#     "ntau",
#     "bw_sd",
#     "noise_iter"
#   )
functions <-
  c(
    # "apply_DDbintau",
    "Langevin1D_adapted",
    "apply_Langevin1D_adapted",
    "D1fun",
    "D2fun",
    # "DDbins",
    "est_D_Carp",
    "finner",
    "generate_Langevin",
    "get_D",
    "get_effective_potential",
    "get_exit_time",
    "get_potential",
    "get_stability",
    "get_theoretical_D",
    "get_weights",
    "gg",
    "new_plot_overview",
    "poly3",
    "setup_filepaths",
    "style_plot"
  )

for (fn in functions) {
  rm(list = fn, envir = .GlobalEnv)
}

# Loop through scenarios
execution_times <- list()
foreach(for_par = forloop) %do% {
  with(for_par, {
    print(as.data.frame(for_par))
    #print("step 1")
    
    # Get polynomial coefficients for forloop
    Ds = get_D(nr_steps_bif,
               scenario,
               type_D2,
               strength_D2)

    # Loop through steps in bifurcation parameter
    foreach(step_idx = 1:nr_steps_bif,
            D = Ds,
            # .packages = c("ggplot2"),
            # .export = c("interpol_steps"),
            .export = c(functions)) %dopar% {
              #print("step 2")
              # start_time <- Sys.time()
              # Setup and create filepaths
              paths = do.call(setup_filepaths, utils::modifyList(
                D,
                list(
                  filepath_est = filepath_est,
                  filepath_figs = filepath_figs,
                  scenario = scenario,
                  type_D2 = type_D2,
                  strength_D2 = strength_D2,
                  sf = sf,
                  N = N,
                  noise_iter = noise_iter,
                  step_idx = step_idx,
                  bins = bins,
                  ntau = ntau,
                  interpol_steps = interpol_steps,
                  bw_sd = bw_sd
                )
              ))
              
              # print(paths)

              if (!file.exists(paths$filepath_out)) {
                
                # print("step 3")
                
                # Generate timeseries
                Ux = do.call(generate_Langevin, utils::modifyList(D, list(
                  N = N,
                  sf = sf,
                  noise_iter = noise_iter
                )))
                
                # Generate concatenated timeseries
                # Ux <- NULL
                # for (i in 0:(noise_iter_concat_times-1)){
                #   Ux_iter <- do.call(generate_Langevin, utils::modifyList(D, list(
                #     N = N/noise_iter_concat_times,
                #     sf = sf,
                #     noise_iter = noise_iter + i
                #   )))
                #   # print(noise_iter + i)
                #   Ux <- c(Ux, Ux_iter)
                # }

                # Get stability of fixed points
                print("Get theoretical fixed points and their stability")
                stabs = get_stability(D)

                # Apply Exit Time Analysis (Carpenter's (2022) approach adapted )
                print("Estimate Langevin model & calculate exit time and alt. metrics")
                est_Carp = est_D_Carp(
                  Ux,
                  sf, 
                  D,
                  stabs,
                  # Tstep = 1:length(as.vector(Ux)),
                  Tstep = 1:length(as.vector(Ux)) / sf,
                  ntau = ntau,
                  bins = bins,
                  bw_sd = bw_sd,
                  interpol_steps = interpol_steps
                )

                # Save results
                out = c(as.list(environment()))  # Gather environment
                out[unlist(lapply(out, class)) == "function"] = NULL # Remove functions
                # saveRDS(out, paths$filepath_out)
                
                saveRDS(out, paste0(paths$filepath_out))
                # saveRDS(out, paste0())
                
                # end_time <- Sys.time()  
                # execution_times[[i]] <- end_time - start_time 
              }
              out = readRDS(paths$filepath_out)
              # out = readRDS()
              print("Plot results")
              new_plot_overview(out, paths$filepath_image, plot_t = ifelse(N*sf < 100000, Inf, 100000))
              graphics.off()
              # end_time <- Sys.time()
              # execution_times[[i]] <- end_time - start_time
            }
  })
}
parallel::stopCluster(cl) # End cluster  

### Checking helper_functions_Rinn
# load example data
example <-
  readRDS(
    "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/Rinn_est_scale_check/2fps-balanced-deepening/constant-D2/D2strength0.3000_sf10_N1000_iter0003_step0002_pars-1.50_0.00_1.50_0.00_0.00_0.00_0.30_bins100_ntau3_interpol50_bw0.30.RDS"
  )

# D1 and D2 ("raw") estimated from Langevin1D_adapted
out$est_Carp$est_df_Rinn[out$est_Carp$est_df_Rinn$variable == "D1",]

out$est_Carp$compl_df[out$est_Carp$compl_df$variable == "drift",]
out_D1 <- as.numeric(out$est_Carp$compl_df[out$est_Carp$compl_df$variable == "drift" & out$est_Carp$compl_df$source == "Estimated",]$value)
out_D1x <- out$est_Carp$compl_df[out$est_Carp$compl_df$variable == "drift" & out$est_Carp$compl_df$source == "Estimated",]$x

out_D1_theo <- as.numeric(out$est_Carp$compl_df[out$est_Carp$compl_df$variable == "drift" & out$est_Carp$compl_df$source == "Theoretical",]$value)
out_D1x_theo <- out$est_Carp$compl_df[out$est_Carp$compl_df$variable == "drift" & out$est_Carp$compl_df$source == "Theoretical",]$x

plot(out_D1x, out_D1)
points(out_D1x_theo,out_D1_theo)

plot(xvec, drift)
points(DD$D1s$x, DD$D1s$y, col = "red")
points(DD$D1s$x, DD$D1s$x - DD$D1s$x^3, col = "blue")

plot(xvec, negPF)
plot(Theoretical_df$x[Theoretical_df$variable == "potential"], Theoretical_df$value[Theoretical_df$variable == "potential"])

### Bootstrap 
# load example data
example <-
  readRDS(
    "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/Rinn_est_scale_check/2fps-balanced-deepening/constant-D2/D2strength0.3000_sf10_N1000_iter0003_step0002_pars-1.50_0.00_1.50_0.00_0.00_0.00_0.30_bins100_ntau3_interpol50_bw0.30.RDS"
  )

# get random states
example$est_Carp$
random_initial_states <- sample(estimated_prob_dist$x, size = 100, replace = TRUE, prob = estimated_prob_dist$value)
plot(random_initial_states)

plot(estimated_prob_dist$x, estimated_prob_dist$value)

#debug 
D1s = example$est_Carp$DD$DDLangout$D1s
D2s = example$est_Carp$DD$DDLangout$D2s
eD1 = example$est_Carp$DD$DDLangout$eD1
eD2 = example$est_Carp$DD$DDLangout$eD2
mean_bin = example$est_Carp$DD$DDLangout$mean_bin
bw_sd = example$bw_sd
Ux = example$Ux

rm(list = c("D1s", "D2s", "eD1", "eD2", "mean_bin", "bw_sd", "Ux", "eD1s", "eD2s"))

# Fit polynomials to the estimated drift and diffusion coefficients
bootstrap_est_D <- function(D1s, D2s, eD1, eD2, mean_bin, bw_sd, Ux){
  # estD1 <- coef(lm(D1 ~ mean_bin + I(mean_bin^2) + I(mean_bin^3), weights = 1/eD1)) # specify only x and x^3?
  # estD2 <- coef(lm(D2 ~ mean_bin+ I(mean_bin^2), weights = 1/eD2)) # depends on type of noise? 
  bw <- bw_sd * sd(Ux)
  eD1s = ksmooth(x= mean_bin,y= eD1,kernel= 'normal',bandwidth= bw,
                x.points=mean_bin)
  eD2s = ksmooth(x= mean_bin ,y= eD2,kernel= 'normal',bandwidth=bw,
                x.points=mean_bin)
  
  estD1_simple <- coef(lm(D1s$y ~ 0 + D1s$x + I(D1s$x^3), weights = 1/eD1s$y)) # specify only x and x^3?
  estD2_simple <- coef(lm(D2s$y ~ 1, weights = 1/eD2s$y))
  
  return(list(
    # d13 = estD1[4],
    # d12 = estD1[3],
    # d11 = estD1[2],
    # d10 = estD1[1],
    # d22 = estD2[3],
    # d21 = estD2[2],
    # d20 = estD2[1]
    d13 = estD1_simple[2],
    d12 = 0,
    d11 = estD1_simple[1],
    d10 = 0,
    d22 = 0,
    d21 = 0,
    d20 = estD2_simple[1]
  ))
}

bootstrapDDs <- bootstrap_est_D(example$est_Carp$DD$DDLangout$D1s, example$est_Carp$DD$DDLangout$D2s, example$est_Carp$DD$DDLangout$eD1, example$est_Carp$DD$DDLangout$eD2, example$est_Carp$DD$DDLangout$mean_bin, example$bw_sd, example$Ux)

# check: plot estD1 and estD2
x_values <- seq(min(example$est_Carp$DD$DDLangout$D1s$x, na.rm = TRUE), max(example$est_Carp$DD$DDLangout$D1s$x, na.rm = TRUE), length.out = 100)

# predicted values for D1
y_pred_D1_simple <- estD1_simple[1] * x_values + estD1_simple[2] * x_values^3
y_pred_D1 <- estD1[1] + estD1[2]*x_values + estD1[3]*x_values^2 + estD1[4]*x_values^3
# predicted values for D2
y_pred_D2_simple <- estD2_simple[1]
y_pred_D2 <- estD2[1] + estD2[2]*x_values + estD2[3]*x_values^2

#predicted values bootstrap
y_pred_D1_boot <- bootstrapDDs$d11 * x_values + bootstrapDDs$d13 * x_values^3
y_pred_D2_boot <- bootstrapDDs$d20

plot(example$est_Carp$DD$DDLangout$D1s$x, example$est_Carp$DD$DDLangout$D1s$y, xlab = "mean_bin", ylab = "D1", main = "Polynomial Regression") #actual estimated
lines(x_values, y_pred_D1_boot, col = "red")
# lines(x_values, y_pred_D1, col = "red")
# lines(x_values, y_pred_D1_simple, col = "pink")
points(example$est_Carp$DD$DDLangout$D2s$x, example$est_Carp$DD$DDLangout$D2s$y, col = "blue") #actual estimated
abline(h = y_pred_D2_boot, col = "purple")
# lines(x_values, y_pred_D2, col = "blue")
# abline(h = y_pred_D2_simple, col = "purple")

# bootstrap: simulate time series from the reconstructed coefficients
bootstrap_Lang <- function(bootstrap_n, bootstrapDDs, N, sf, ntau, bins, bw_sd, interpol_steps, compl_df, noise_iter){
  full_boot_est_Carp <- list()
  execution_times <- vector("list", length = bootstrap_n)
  
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
  plot(random_initial_states)
  
  for(i in 1:bootstrap_n){
    start_time <- Sys.time()  # Record start time
    
    set.seed(noise_iter)
    boot_Ux <- do.call(timeseries1D, utils::modifyList(bootstrapDDs, list(
      N = N*sf,
      sf = sf,
      startpoint = random_initial_states[i]
    )))
    
    if (any(is.na(boot_Ux) | is.infinite(boot_Ux))) {
      print("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
    } else {
      print("Timeseries simulated successfully")
    }
    
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
    
    end_time <- Sys.time()  # Record end time
    execution_times[[i]] <- end_time - start_time  # Calculate and store execution time
    
    full_boot_est_Carp[[i]] <- boot_est_Carp
  }
  
  return(list(full_boot_est_Carp = full_boot_est_Carp, execution_times = execution_times))
}

#debug
bootstrap_n = 100
bootstrapDDs = bootstrapDDs
N = example$N
sf = example$sf
ntau = example$ntau
bins = example$bins
bw_sd = example$bw_sd 
interpol_steps = example$interpol_steps
noise_iter = example$noise_iter
compl_df = example$est_Carp$compl_df

rm(list = c("bootstrap_n", "N", "sf", "ntau", "bins", "bw_sd", "interpol_steps", "noise_iter", "compl_df"))
rm(list = c("bw", "D1", "D2", "i", "D", "boot_D1", "reconstructed_D1s_per_bin"))

estimated_prob_dist = example$est_Carp$compl_df %>%
  filter(variable == "wts", source == "Estimated")
estimated_prob_dist$value <- as.numeric(estimated_prob_dist$value)

# Set up cluster (forking)
# cl <- parallel::makeForkCluster(detectCores() - 1) # using forking
# doParallel::registerDoParallel(cl)
# 
# getDoParWorkers() 

# Initialize a vector to store execution times
estDDs <- list(
  d13 = estD1[4],
  d12 = estD1[3],
  d11 = estD1[2],
  d10 = estD1[1],
  d22 = estD2[3],
  d21 = estD2[2],
  d20 = estD2[1]
)

estDDs2 <- list(
  d13 = estD1[4] -1.,
  d12 = estD1[3] + 1,
  d11 = estD1[2] + 1,
  d10 = estD1[1] -1,
  d22 = estD2[3] -1,
  d21 = estD2[2] - 1,
  d20 = estD2[1] + 1
)

# check: plot estD1 and estD2
x_values <- seq(min(example$est_Carp$DD$DDLangout$D1s$x, na.rm = TRUE), max(example$est_Carp$DD$DDLangout$D1s$x, na.rm = TRUE), length.out = 100)

# predicted values for D1
y_pred_D1_test<- estDDs2[4]$d10 + estDDs2[3]$d11*x_values + estDDs2[2]$d12*x_values^2 + estDDs2[1]$d13*x_values^3
# predicted values for D2
y_pred_D2_test <- estDDs2[7]$d20 + estDDs2[6]$d21*x_values + estDDs2[5]$d22*x_values^2

plot(example$est_Carp$DD$DDLangout$D1s$x, example$est_Carp$DD$DDLangout$D1s$y, xlab = "mean_bin", ylab = "D1", main = "Polynomial Regression") #actual estimated
lines(x_values, y_pred_D1_test, col = "red")
# lines(x_values, y_pred_D1, col = "red")
# lines(x_values, y_pred_D1_simple, col = "pink")
points(example$est_Carp$DD$DDLangout$D2s$x, example$est_Carp$DD$DDLangout$D2s$y, col = "blue") #actual estimated
abline(h =y_pred_D2_test, col = "purple")

# generating timeseries from this goes to infinity
bootstrap_complex_model <-
  bootstrap_Lang(bootstrap_n = 20,
                 bootstrapDDs = estDDs,
                 example$N,
                 example$sf, 
                 example$ntau, 
                 example$bins, 
                 example$bw_sd, 
                 example$interpol_steps,
                 example$est_Carp$compl_df,
                 example$noise_iter)

# length(bootstrap_demo$execution_times)
# length(bootstrap_demo$full_boot_est_Carp[[12]]$meanETl)

saveRDS(bootstrap_demo_smoothed, paste0(filepath_base, "bootstrap_demo_smoothed.RDS"))

# Stop the cluster
# stopCluster(cl)

### checking how to adjust for bias 
bootstrap_demo <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesisbootstrap_demo.RDS") # this I ran 100 times, the parent model was reconstructed based on NOT smoothed D1's and D2's 
bootstrap_demo2 <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesisbootstrap_demo2.RDS") # 20 times, the parent model is the theoretical model
bootstrap_demo_smoothed<- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesisbootstrap_demo_smoothed.RDS") # 20 times, the parent model was reconstructed based on smoothed D1's and D2's 

bootstrap_demo$full_boot_est_Carp[[1]]$fp_df
bootstrap_demo2$full_boot_est_Carp[[1]]$fp_df
bootstrap_demo_smoothed$full_boot_est_Carp[[1]]$fp_df

# bootstrap mean exit times from parent model reconstructed from NOT smoothed D1's and D2's
boot_meanETl <- vector()
for(i in 1:100){
  boot_meanETl[i] <- bootstrap_demo$full_boot_est_Carp[[i]]$meanETl
}
confidence_interval <- quantile(boot_meanETl, c(0.025, 0.975))
confidence_interval # the confidence interval doesn't include the estimated mean exit time
example$est_Carp$meanETl

ggplot(data = data.frame(x = boot_meanETl), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  geom_density(color = "darkblue", size = 0.5) +
  geom_vline(aes(xintercept = example$est_Carp$meanETl), linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Distribution mean exit time (left) estimated from reconstructed models (not smoothed)", x = "Values", y = "Density")

boot_meanETr <- vector()
for(i in 1:100){
  boot_meanETr[i] <- bootstrap_demo$full_boot_est_Carp[[i]]$meanETr
}
confidence_interval <- quantile(boot_meanETr, c(0.025, 0.975))
example$est_Carp$meanETr

ggplot(data = data.frame(x = boot_meanETr), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  geom_density(color = "darkblue", size = 0.5) +
  geom_vline(aes(xintercept = example$est_Carp$meanETr), linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Distribution mean exit time (right) estimated from reconstructed models (not smoothed)", x = "Values", y = "Density")

# bootstrap mean exit times from theoretical (i.e. if parent model WAS theoretical model)
boot_meanETl <- vector()
for(i in 1:20){
  boot_meanETl[i] <- bootstrap_demo2$full_boot_est_Carp[[i]]$meanETl
}
confidence_interval <- quantile(boot_meanETl, c(0.025, 0.975))
example$est_Carp$meanETl

ggplot(data = data.frame(x = boot_meanETl), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  geom_density(color = "darkblue", size = 0.5) +
  geom_vline(
    aes(xintercept = example$est_Carp$meanETl),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  labs(title = "Distribution mean exit time (left) estimated with theoretical model as parent model", x = "Values", y = "Density")

boot_meanETr <- vector()
for(i in 1:20){
  boot_meanETr[i] <- bootstrap_demo2$full_boot_est_Carp[[i]]$meanETr
}
confidence_interval <- quantile(boot_meanETr, c(0.025, 0.975))
example$est_Carp$meanETr

ggplot(data = data.frame(x = boot_meanETr), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  geom_density(color = "darkblue", size = 0.5) +
  geom_vline(
    aes(xintercept = example$est_Carp$meanETr),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  labs(title = "Distribution mean exit time (right) estimated with theoretical model as parent model", x = "Values", y = "Density")

# bootstrap mean exit times from smoothed reconstructed model (the smoothed D1s and D2s are used to estimate the parent model)
boot_meanETl <- vector()
for(i in 1:20){
  boot_meanETl[i] <- bootstrap_demo_smoothed$full_boot_est_Carp[[i]]$meanETl
}
confidence_interval <- quantile(boot_meanETl, c(0.025, 0.975))
example$est_Carp$meanETl

ggplot(data = data.frame(x = boot_meanETl), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  # geom_density(color = "darkblue", size = 0.5) +
  geom_vline(aes(xintercept = example$est_Carp$meanETl), linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Distribution mean exit time (left) estimated from reconstructed models (smoothed)", x = "Values", y = "Density")

boot_meanETr <- vector()
for(i in 1:20){
  boot_meanETr[i] <- bootstrap_demo_smoothed$full_boot_est_Carp[[i]]$meanETr
}
confidence_interval <- quantile(boot_meanETr, c(0.025, 0.975))
example$est_Carp$meanETr

ggplot(data = data.frame(x = boot_meanETr), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  # geom_density(color = "darkblue", size = 0.5) +
  geom_vline(aes(xintercept = example$est_Carp$meanETr), linetype = "dashed", color = "red", alpha = 0.5) +
  labs(title = "Distribution mean exit time (right) estimated from reconstructed models (smoothed)", x = "Values", y = "Density")

##### how to adjust bias for D1?
boot_D1 <- list()
for(i in 1:100){
  boot_D1[[i]] <- bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$D1
}

transposed_D1 <- list()
median_D1_reconstructed <- vector()
for(i in seq_along(boot_D1[[1]])) {
  transposed_D1[[i]] <- sapply(boot_D1, function(x) x[[i]])
  median_D1_reconstructed[i] <- median(transposed_D1[[i]])
}
plot(example$est_Carp$DD$DDLangout$D1, main = "bootstrapped D1s (not smoothed)")
points(median_D1_reconstructed, col = "red", pch = 4)

reconstructed_D1s_perbin <- data.frame(group = rep(1:100, each = 100), value = unlist(transposed_D1))

# Combine reconstructed D1's per bin with parent D1 estimates
reconstructed_D1s_perbin <- data.frame(group = rep(1:100, each = 100), value = unlist(transposed_D1), parent_D1 = rep(example$est_Carp$DD$DDLangout$D1, each = 100))

# Density plot of D1 estimates per bin from reconstructed models (reconstructed on parent model, not smoothed) 
plot_notsmoothed_parent <- ggplot(reconstructed_D1s_perbin, aes(x = value)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  geom_density(color = "darkblue", size = 0.5) +
  facet_wrap( ~ group, scales = "free") +
  geom_vline(
    data = subset(reconstructed_D1s_perbin, !is.na(parent_D1)),
    aes(xintercept = parent_D1),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    axis.title = element_blank(),
    # axis.text = element_blank(),
    axis.ticks = element_blank()
  )
# ggsave("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/bootstrapped_D1s_notsmoothed.pdf", plot = plot_notsmoothed_parent)

# Density plot of D1 estimates per bin from reconstructed models (reconstructed on parent model, smoothed) 
boot_D1 <- list()
for(i in 1:20){
  boot_D1[[i]] <- bootstrap_demo_smoothed$full_boot_est_Carp[[i]]$DD$DDLangout$D1
}

transposed_D1 <- list()
median_D1_reconstructed <- vector()
for(i in seq_along(boot_D1[[1]])) {
  transposed_D1[[i]] <- sapply(boot_D1, function(x) x[[i]])
  median_D1_reconstructed[i] <- median(transposed_D1[[i]])
}
plot(example$est_Carp$DD$DDLangout$D1, main = "bootstrapped D1s (smoothed)")
points(median_D1_reconstructed, col = "red", pch = 4)

reconstructed_D1s_perbin <- data.frame(group = rep(1:100, each = 100), value = unlist(transposed_D1))

# Combine reconstructed D1's per bin with parent D1 estimates
reconstructed_D1s_perbin <- data.frame(group = rep(1:100, each = 100), value = unlist(transposed_D1), parent_D1 = rep(example$est_Carp$DD$DDLangout$D1, each = 100))

ggplot(reconstructed_D1s_perbin, aes(x = value)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  geom_density(color = "darkblue", size = 0.5) +
  facet_wrap( ~ group, scales = "free") +
  geom_vline(
    data = subset(reconstructed_D1s_perbin, !is.na(parent_D1)),
    aes(xintercept = parent_D1),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  theme_minimal() +
  theme(
    panel.spacing = unit(0, "points"),
    strip.background = element_blank(),
    axis.title = element_blank(),
    # axis.text = element_blank(),
    axis.ticks = element_blank()
  )

# checks
boot_D1_bin28 <- c()
for(i in 1:100){
  boot_D1_bin28[i] <- bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$D1[28]
}

transposed_D1[[28]] == boot_D1_bin28 #yes 
median(boot_D1_bin28)
mean(sort(boot_D1_bin28)[50], sort(boot_D1_bin28)[51])

example$est_Carp$DD$DDLangout$mean_bin
example$est_Carp$DD$DDLangout$D1
example$est_Carp$DD$DDLangout$D2

parent_EPF$y <- as.numeric(bootstrap_demo$full_boot_est_Carp[[1]]$compl_df$value[bootstrap_demo$full_boot_est_Carp[[1]]$compl_df$variable == "EPF" & bootstrap_demo$full_boot_est_Carp[[1]]$compl_df$source == "Theoretical"])

parent_EPF$x <- bootstrap_demo$full_boot_est_Carp[[1]]$compl_df$x[bootstrap_demo$full_boot_est_Carp[[1]]$compl_df$variable == "EPF" & bootstrap_demo$full_boot_est_Carp[[1]]$compl_df$source == "Theoretical"]

plot(parent_EPF$x, parent_EPF$y)

# Call the function in parallel
foreach(i = 1:bootstrap_n, .combine = c) %dopar% {
  execution_times[[i]] <- system.time({
    bootstrap_Lang(bootstrap_n = bootstrap_n,
                   bootstrapDDs = bootstrapDDs,
                   N = example$N,
                   sf = example$sf, 
                   ntau = example$ntau, 
                   bins = example$bins, 
                   bw_sd = example$bw_sd, 
                   interpol_steps = example$interpol_steps,
                   compl_df = example$est_Carp$compl_df)
  })[["elapsed"]]
  
  # Print the execution times
  print(execution_times)
}

# Plot bootstrap time series and probability density function (PDF)
try <- timeseries1D(1000, startpoint = random_initial_states[1], d13 = estD1_simple[2], d12 = 0, d11 = estD1_simple[1], d22 = 0, d21 = 0, d20 = estD2_simple, sf = 10)

op <- par(no.readonly = TRUE)
par(mar = c(4.2, 5.4, 0.5, 0.5))
layout(matrix(c(1, 1, 2), 1, 3))
plot((1:length(try)/sf), try, xlab = "t [s]", ylab = "x [a.u.]", t = 'l')
plot(density(try), xlab = "x [a.u.]", ylab = "density", t = 'l',main = "")
par(op)

# Plot Langevin bootstrap CI of D1's and D2's 
library("plotrix")
for(i in 1:5){
  par(mfrow = c(1, 2))
  plotCI(
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin,
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$D1,
    uiw = bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$eD1,
    xlab = "x [a.u.]",
    ylab = expression(paste("Drift coefficient",
                            +Dˆ(1), "(x) [a.u.]")),
    cex = 2,
    pch = 20
  )
  lines(
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin,
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin - bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin ^
      3,
    col = "red",
    lwd = 3,
    lty = 2
  )
  plotCI(
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin,
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$D2,
    uiw = bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$eD2,
    xlab = "x [a.u.]",
    ylab = expression(paste(
      "Diffusion coefficient", Dˆ(2), "(x) [a.u.]"
    )),
    cex = 2,
    pch = 20
  )
  lines(
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin,
    bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$mean_bin ^ 2 + 1,
    col = "red",
    lwd = 3,
    lty = 2
  )
}

### Assumption checks ###

### N-step
bootstrapDDs <- bootstrap_est_D(example$est_Carp$DD$DDLangout$D1s, example$est_Carp$DD$DDLangout$D2s, example$est_Carp$DD$DDLangout$eD1, example$est_Carp$DD$DDLangout$eD2, example$est_Carp$DD$DDLangout$mean_bin, example$bw_sd, example$Ux)

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

# One step
Onestep_ahead <- function(bootstrapDDs, bootstrap_n, Ux, N, sf, noise_iter){
  Onestep_predictions <- vector()
  last_state <- Ux[length(Ux)]
  seeds <- c(noise_iter, seq(50,149))
  
  set.seed(noise_iter)
  # actual_timeseries_state <- do.call(timeseries1D, utils::modifyList(D, list(
  #   N = 10,
  #   sf = 10,
  #   startpoint = last_state
  # )))[1]
  actual_timeseries_state <-
    timeseries1D(
      N = 10,
      sf = 10,
      startpoint = last_state,
      d13 = bootstrapDDs$d13,
      d12 = bootstrapDDs$d12,
      d11 = bootstrapDDs$d11,
      d10 = bootstrapDDs$d10,
      d22 = bootstrapDDs$d22,
      d21 = bootstrapDDs$d21,
      d20 = bootstrapDDs$d20
    )[1]
  
  actual_timeseries_state <-
    timeseries1D(
      N = 10,
      sf = 10,
      startpoint = example$Ux[218],
      d13 = bootstrapDDs$d13,
      d12 = bootstrapDDs$d12,
      d11 = bootstrapDDs$d11,
      d10 = bootstrapDDs$d10,
      d22 = bootstrapDDs$d22,
      d21 = bootstrapDDs$d21,
      d20 = bootstrapDDs$d20
    )[1]
  
  for(i in 1:(bootstrap_n+1)){
    
    set.seed(seeds[i])
    Onestep_Ux <- do.call(timeseries1D, utils::modifyList(bootstrapDDs, list(
      N = 10,
      sf = 10,
      startpoint = last_state
    )))
    
    if (any(is.na(Onestep_Ux) | is.infinite(Onestep_Ux))) {
      print("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
    } else {
      print("Timeseries simulated successfully")
    }
    
    one_step_pred <- Onestep_Ux[1]
    Onestep_predictions[i] <- one_step_pred
  }
  
  return(list(actual_timeseries_state = actual_timeseries_state, Onestep_predictions = Onestep_predictions))
}

Onestep_demo <-
  Onestep_ahead(
    bootstrapDDs = bootstrapDDs,
    bootstrap_n = 100,
    Ux = example$Ux,
    N = example$N,
    sf = example$sf,
    noise_iter = example$noise_iter
  )

# plot CI
confidence_interval <- quantile(Onestep_demo$Onestep_predictions, c(0.025, 0.975))
hist(Onestep_demo$Onestep_predictions[-1], breaks = 30, prob = TRUE, col = "lightblue", 
     xlab = "Bootstrapped 1-step ahead predictions", ylab = "Density", main = "Distribution of bootstrapped estimates")
lines(density(Onestep_demo$Onestep_predictions[-1]), col = "darkblue", lwd = 2)

# actual timeseries state
abline(v = Onestep_demo$actual_timeseries_state, col = "red", lwd = 2)
# point estimate
abline(v = Onestep_demo$Onestep_predictions[1], col = "purple", lwd = 2)
# vertical lines for confidence interval
abline(v = confidence_interval[1], col = "lightblue", lty = 2)
abline(v = confidence_interval[2], col = "lightblue", lty = 2)
legend(
  "topright",
  legend = c("Actual state", "Point Estimate", "95% CI"),
  col = c("red", "purple", "lightblue"),
  lty = c(1, 1, 2),
  lwd = c(2, 2, 1),
  x.intersp = 0.8,
  y.intersp = 0.8,
  xjust = 2,
  yjust = 2
)

# n step (stationarity)
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
    
    start_states <- Ux[round(seq(17, length(Ux)-3, length.out = Npredictions))]
    
    actual_Nstep_states <- matrix(NA, nrow = Npredictions, ncol = Nstep)
    for(i in 1:Nstep){
      actual_Nstep_states[,i] <-
        Ux[round(seq(17, length(Ux) - 3, length.out = Npredictions)) + i]
    }
    
    Nstep_predictions <- matrix(NA, nrow = Npredictions, ncol = Nstep)
    for(i in 1:Nstep){
      for(j in 1:Npredictions){
        # print(c("start state is", start_states[j]))
        set.seed(noise_iter)
        model_timeseries <-
          timeseries1D(
            N = N*sf,
            sf = sf,
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
    
    return(
      list(
        actual_Nstep_states = actual_Nstep_states,
        Nstep_predictions = Nstep_predictions,
        RMSE = RMSE
      )
    )
  }

Nstep_demo <- Nstep_ahead(Nstep = 3,
                          # bootstrap_n = 100,
                          Npredictions = 100,
                          bootstrapDDs = bootstrapDDs,
                          Ux = example$Ux,
                          N = example$N,
                          sf = example$sf,
                          noise_iter = example$noise_iter,
                          bins = example$bins,
                          ntau = example$ntau,
                          bw_sd = example$bw_sd) 

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

#### Check stationary
?adf.test
tseries::adf.test(example$Ux)

#### Check Markov property
spgs::markov.test(example$Ux, type = "lb.test")

# Notes from Luiza Yuan, 03.2024:
### Check parallel processing output ###
attempt<-
  readRDS(
    "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/est_parallel_check/2fps-balanced-deepening/constant-D2/D2strength0.5000_sf10_N1000_iter0001_step0001_pars-1.00_0.00_1.00_0.00_0.00_0.00_0.50_bins100_ntau10_interpol100_bw0.30.RDS"
  )
attempt$est_Carp$fp_df$xintercept[4]

### Check new_plot_overview function when no mean exit times estimated ###
out <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Modified_for_Luiza_2024_02_07/est_parallel_check/2fps-balanced-deepening/constant-D2/N1000sf10interpol500/D2strength0.5000_sf10_N1000_iter0001_step0005_pars-3.00_0.00_3.00_0.00_0.00_0.00_0.50_bins100_ntau10_interpol500_bw0.30.RDS")
filepath_image <- "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Modified_for_Luiza_2024_02_07/figs_parallel_check/2fps-balanced-deepening/constant-D2/N1000sf10interpol500/D2strength0.5000_sf10_N1000_iter0001_step0005_pars-3.00_0.00_3.00_0.00_0.00_0.00_0.50_bins100_ntau10_interpol500_bw0.30.pdf"
new_plot_overview(out, paths$filepath_image, plot_t = ifelse(N*sf < 100000, Inf, 100000))
# debug
est_Carp = out$est_Carp
stabs = out$stabs
Ux = out$Ux

# Notes from Luiza Yuan, 29.01.2024:

#################################################################################################################
#I just wanted to check whether differences in the plots for theoretical drift, diffusion, and potential were simply due to differences in the min and max x's used in sf10 and sf100
#################################################################################################################

# Get two examples of sf10 vs sf100 (in this case iter = 1, step = 4 i.e. drift = -2.5x^3 + 2.5x)
N400_sf10 <-
  readRDS(
    "/Users/luizashen58/Documents/Documents - Luiza’s MacBook Air/UvA Master Yr 2/Thesis/theoretical_ET_2024_01_29/est_theoreticalET/2fps-balanced-deepening/constant-D2/D2strength_0.3000_sf10_N400_iter0001_step0004_pars-2.50_0.00_2.50_0.00_0.00_0.00_0.30.RDS"
  )

N400_sf100 <-
  readRDS(
    "/Users/luizashen58/Documents/Documents - Luiza’s MacBook Air/UvA Master Yr 2/Thesis/est_theoreticalET_scalecheck/2fps-balanced-deepening/constant-D2/D2strength_0.3000_sf100_N400_iter0001_step0004_pars-2.50_0.00_2.50_0.00_0.00_0.00_0.30.RDS"
  )

## The theoretical drift, diffusion, and potentials were not saved (I need to fix this), so I recreated them here in order to check some things

# original code for plotting theoretical drift/diffusion/potential:
# xlim = max(abs(c(stabs$fps, Ux)))
# theoretical_df = get_theoretical_D(D, min_x = -xlim, max_x = xlim)

get_theoretical_df <- function(D, data) {
  xlim = max(abs(c(data$stabs$fps, data$Ux)))
  get_theoretical_D(D, min_x = -xlim, max_x = xlim)
}

N400_sf10_theoretical_df <-
  get_theoretical_df(N400_sf10$D, N400_sf10)
N400_sf100_theoretical_df <-
  get_theoretical_df(N400_sf100$D, N400_sf100)

# check
N400_sf10_min_x <-
  min(N400_sf10_theoretical_df$x[N400_sf10_theoretical_df$variable == "drift_analytical"])
N400_sf10_max_x <-
  max(N400_sf10_theoretical_df$x[N400_sf10_theoretical_df$variable == "drift_analytical"])

N400_sf100_min_x <-
  min(N400_sf100_theoretical_df$x[N400_sf100_theoretical_df$variable == "drift_analytical"])
N400_sf100_max_x <-
  max(N400_sf100_theoretical_df$x[N400_sf100_theoretical_df$variable == "drift_analytical"])

- 2.5 * (N400_sf10_min_x ^ 3) + 2.5 * (N400_sf10_min_x) == max(N400_sf10_theoretical_df$value[N400_sf10_theoretical_df$variable == "drift_analytical"])
- 2.5 * (N400_sf10_max_x ^ 3) + 2.5 * (N400_sf10_max_x) == min(N400_sf10_theoretical_df$value[N400_sf10_theoretical_df$variable == "drift_analytical"])

- 2.5 * (N400_sf100_min_x ^ 3) + 2.5 * (N400_sf100_min_x) == max(N400_sf100_theoretical_df$value[N400_sf100_theoretical_df$variable == "drift_analytical"])
- 2.5 * (N400_sf100_max_x ^ 3) + 2.5 * (N400_sf100_max_x) == min(N400_sf100_theoretical_df$value[N400_sf100_theoretical_df$variable == "drift_analytical"])


#################################################################################################################
#I wanted to check how the 'bins' argument affects theoretical exit time estimation because 'bins' should affect only the estimation of drift and diffusion and not anything to do with the theoretical drift/diffusion/exit times
################################################################################################################

get_theoretical_df <- function(D, data) {
  xlim = max(abs(c(data$stabs$fps, data$Ux)))
  get_theoretical_D(D, min_x = -xlim, max_x = xlim)
} # copy/pasted here just for clarity

# read-in example for N = 400, sf = 10, bins = 100
N400_sf10 <-
  readRDS(
    "/Users/luizashen58/Documents/Documents - Luiza’s MacBook Air/UvA Master Yr 2/Thesis/theoretical_ET_2024_01_29/est_theoreticalET/2fps-balanced-deepening/constant-D2/D2strength_0.3000_sf10_N400_iter0001_step0004_pars-2.50_0.00_2.50_0.00_0.00_0.00_0.30.RDS"
  ) # copy/pasted here just for clarity
N400_sf10_theoretical_df <-
  get_theoretical_df(N400_sf10$D, N400_sf10)

# read-in example for N = 400, sf = 10, bins = 40
N400_sf10_bin40 <-
  readRDS(
    "/Users/luizashen58/Documents/Documents - Luiza’s MacBook Air/UvA Master Yr 2/Thesis/est_theoreticalET_scalecheck/2fps-balanced-deepening/constant-D2/N400_sf10_bin40/D2strength_0.3000_sf10_N400_iter0001_step0004_pars-2.50_0.00_2.50_0.00_0.00_0.00_0.30.RDS"
  )

N400_sf10_bin40_theoretical_df <-
  get_theoretical_df(N400_sf10_bin40$D, N400_sf10_bin40)

# the values for theoretical drift, diffusion, and potential are the same regardless of bin size (as expected)
N400_sf10_bin40_theoretical_df$value[N400_sf10_bin40_theoretical_df$variable == "drift_analytical"] == N400_sf10_theoretical_df$value[N400_sf10_theoretical_df$variable == "drift_analytical"]

N400_sf10_bin40_theoretical_df$value[N400_sf10_bin40_theoretical_df$variable == "diff_analytical"] == N400_sf10_theoretical_df$value[N400_sf10_theoretical_df$variable == "diff_analytical"]

N400_sf10_bin40_theoretical_df$value[N400_sf10_bin40_theoretical_df$variable == "potential_analytical"] == N400_sf10_theoretical_df$value[N400_sf10_theoretical_df$variable == "potential_analytical"]

## checking 'xvec' (interpolation vector specifying points on the x-axis)
# (this is the only thing I found in the theoretical exit time estimation/get_theoretical_exit_time function that is related to bins/estimated drift/diffusion)

# xvec specification in get_exit_time function (and hence in get_theoretical_exit_time function):

# nearest = .1
# xvec = seq(-nearest * floor(abs(min(DD$D1s$x)) / nearest),
#            nearest * floor(abs(max(DD$D1s$x)) / nearest),
#            length.out = interpol_steps)

## recreating DD for N400_sf10 (b/c not saved, need to fix)
sf = N400_sf10$sf
steps = c(1:3)
bins = 100 #note
bw_sd = .3
interpol_steps = 100
Tstep = 1:length(as.vector(N400_sf10$Ux)) / sf
ntau = max(steps)

N400_sf10_DD = apply_DDbintau(
  Ux = N400_sf10$Ux,
  Tstep = Tstep,
  ntau = ntau,
  bins = bins,
  bw_sd = bw_sd
)

# xvec for N400_sf10
nearest = .1
N400_sf10_xvec <-
  seq(-nearest * floor(abs(min(N400_sf10_DD$D1s$x)) / nearest),
      nearest * floor(abs(max(N400_sf10_DD$D1s$x)) / nearest),
      length.out = interpol_steps)

## recreating DD for N400_sf10_bin40 (b/c not saved)
sf = N400_sf10_bin40$sf
steps = c(1:3)
bins = 40 # note
bw_sd = .3
interpol_steps = 100
Tstep = 1:length(as.vector(N400_sf10_bin40$Ux)) / sf
ntau = max(steps)

N400_sf10_bin40_DD = apply_DDbintau(
  Ux = N400_sf10_bin40$Ux,
  Tstep = Tstep,
  ntau = ntau,
  bins = bins,
  bw_sd = bw_sd
)

# xvec for N400_sf10_bin40
nearest = .1
N400_sf10_bin40_xvec <-
  seq(-nearest * floor(abs(min(
    N400_sf10_bin40_DD$D1s$x
  )) / nearest),
  nearest * floor(abs(max(
    N400_sf10_bin40_DD$D1s$x
  )) / nearest),
  length.out = interpol_steps)

## Now, inspect setups for bvpsolve!!!

# original specification for discretized domain over which the values of the two-point boundary value problem solution variables are calculated:

# solve the left basin from x = 0 (reflecting) to x=xeq[2] (absorbing)
# x = seq(xeq[1]-1,xeq[2],length.out=30)  # x vector, original (Carpenter)
# x = seq(min(xvec), xeq[2], length.out = ceiling(interpol_steps / 2)) # wider interval/discretized domain

# N400_sf10_x_original = seq(N400_sf10$est_Carp$xeq[1]-1,N400_sf10$est_Carp$xeq[2],length.out=30)
# N400_sf10_bin40_x_original = seq(N400_sf10$est_Carp$xeq[1]-1,N400_sf10$est_Carp$xeq[2],length.out=30)
# N400_sf10_x_original == N400_sf10_bin40_x_original # these are the same (ofc)

## this is the actual setup used (wider interval/discretized domain)
N400_sf10_left_x = seq(min(N400_sf10_xvec),
                       N400_sf10$est_Carp$xeq[2],
                       length.out = ceiling(interpol_steps / 2))
N400_sf10_bin40_left_x = seq(
  min(N400_sf10_bin40_xvec),
  N400_sf10_bin40$est_Carp$xeq[2],
  length.out = ceiling(interpol_steps / 2)
)

c(min(N400_sf10_left_x), max(N400_sf10_left_x))
c(min(N400_sf10_bin40_left_x), max(N400_sf10_bin40_left_x)) # this interval is wider!

# right basin from x=xeq[2] (absorbing) to x > xeq[3] (reflecting)
# x = seq(xeq[2],xeq[3]+1,length.out=30)  # x vector original (Carpenter)
# x = seq(xeq[2], max(xvec), length.out = ceiling(interpol_steps / 2))  # wider interval/discretized domain

# N400_sf10_x_original = seq(N400_sf10$est_Carp$xeq[2],N400_sf10$est_Carp$xeq[3]+1,length.out=30)
# N400_sf10_bin40_x_original = seq(N400_sf10$est_Carp$xeq[2],N400_sf10$est_Carp$xeq[3]+1, length.out=30)
# N400_sf10_x_original == N400_sf10_bin40_x_original # again, these are the same (ofc)

## this is the actual setup used (wider interval/ discretized domain)
N400_sf10_right_x = seq(N400_sf10$est_Carp$xeq[2],
                        max(N400_sf10_xvec),
                        length.out = ceiling(interpol_steps / 2))
N400_sf10_bin40_right_x = seq(
  N400_sf10_bin40$est_Carp$xeq[2],
  max(N400_sf10_bin40_xvec),
  length.out = ceiling(interpol_steps / 2)
)

c(min(N400_sf10_right_x), max(N400_sf10_right_x))
c(min(N400_sf10_bin40_right_x), max(N400_sf10_bin40_right_x)) # this interval is wider!
