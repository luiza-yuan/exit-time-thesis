# Timeseries simulation and estimation of mean exit time and alternative resilience metrics using the Langevin approach
# When simulating the timeseries, we first generate a high resolution timeseries close to the conditions of Arani et al. (2021), with 12000 timeunits with timestep of 0.0001. This forms the "complete" timeseries (e.g. one can imagine this as the complete timeseries of an individual). From here, we downsample before conducting the analysis. This is similar to in real life, when we take only a sample of the complete timeseries when collecting data. 
# As an example below, after burning the first 10,000 datapoints (as best practice to ensure convergence for stationarity and Markovian process), we take every 200th datapoint. In other words, for each timeunit, we now have "timestep" of 0.2. 

# Set working directory
rm(list = ls())
graphics.off()

# Get/set working directory
if (interactive() &&
    requireNamespace("rstudioapi", quietly = TRUE)) {
  # RStudio-specific code
  filepath_base <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(filepath_base)
} else {
  # Alternative code for non-interactive environments (e.g. Snellius)
  filepath_base <- getwd()  # Use the current working directory
}

# Create necessary directories for estimates and figures
filepath_figs = file.path(filepath_base, "figs") 
filepath_est = file.path(filepath_base, "est")
if (!dir.exists(filepath_figs)) {
  dir.create(filepath_figs, recursive = T)
}

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
source(file.path(filepath_base, "helper_functions_2024_05_22.R"))

# Choose parameters to loop through
forloop = tidyr::expand_grid(
  # number of bifurcation steps 
  nr_steps_bif = 3, 
  # diffusion type: (1) constant (2) quadratic 
  type_D2 = c("constant"), 
  # bifurcation scenario: (1) "2fps-balanced-deepening" (balanced double basin), (2) "right-fp-gains-dominance" (asymmetrical basin depth)
  scenario = c("2fps-balanced-deepening"), 
  # strength of diffusion (noise intensity)
  strength_D2 = c(.3), 
  # number of time units for simulating timeseries
  N_high_res = c(12000), 
  # sampling frequency (i.e. 1 / number of timesteps) for simulating timeseries
  sf_high_res = c(1000), 
  # length of timeseries final sample
  N = c(11990),
  # new "sampling frequency" for timeseries final sample
  sf = c(5), 
  # number of bins for range of states for estimating drift and diffusion functions
  bins = c(50), 
  # number of interpolation steps 
  interpol_steps = 50,
  # maximum timelag for estimating conditional moments
  ntau = c(2), 
  # length to add to discretized domain over which to solve boundary value problem for exit time
  add_to_x = 0.5,
  # bin width for smoothing D1 and D2 estimates
  bw_sd = 0.3,
  # number of noise iterations
  noise_iter = c(1:5), 
) %>% purrr::transpose() %>% unique()

# Set up cluster (forking)
cl <- parallel::makeForkCluster(detectCores() - 1) # using forking
doParallel::registerDoParallel(cl)

getDoParWorkers() #check how many workers 'foreach' is going to use

# Names of functions for parallel processing 
functions <-
  c(
    "apply_Langevin1D_adapted",
    "D1fun",
    "D2fun",
    "est_D_Carp",
    "finner",
    "generate_Langevin",
    "downsample_Langevin",
    "get_D",
    "get_effective_potential",
    "get_exit_time",
    "get_potential",
    "get_resilience_potential",
    "get_resilience_prob_dist",
    "get_stability",
    "get_theoretical_D",
    "get_weights",
    "gg",
    "Langevin1D_adapted",
    "new_plot_overview",
    "poly3",
    "setup_filepaths",
    "style_plot"
  )

# Loop through scenarios
start_time <- Sys.time()
foreach(for_par = forloop) %do% {
  with(for_par, {
    print(as.data.frame(for_par))
    
    # Get polynomial coefficients for generating timeseries from Langevin model
    Ds = get_D(nr_steps_bif,
               scenario,
               type_D2,
               strength_D2)
    
    # Loop through steps in bifurcation parameter
    foreach(step_idx = 1:nr_steps_bif,
            D = Ds,
            .export = c(functions)) %dopar% {
              
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
                  bw_sd = bw_sd,
                  sf_high_res = sf_high_res,
                  N_high_res = N_high_res
                )
              ))
              
              if (!file.exists(paths$filepath_out)) {
              
              # Generate high resolution timeseries with Langevin model 
              # Note: Arani et al. (2021) 12,000 time units with time step of 0.001
              print("Generate high resolution timeseries")
              Ux_original = do.call(generate_Langevin, utils::modifyList(D, list(
                N = N_high_res,
                sf = sf_high_res,
                noise_iter = noise_iter
              )))
              
              # Burn first 10000 datapoints (burn)
              Ux_original_burned = window(Ux_original, start = c(11, 1), end = c(12000, 1000))
              
              # Downsample from original high resolution timeseries to get final Ux
              Ux <- downsample_Langevin(Ux_original_burned, n_points = N*sf)
              
              # Get stability of fixed points
              print("Get theoretical fixed points and their stability")
              stabs = get_stability(D)
              
              # Apply Exit Time Analysis (Carpenter's (2022) approach adapted )
              print("Estimate Langevin model & calculate exit times and alt. metrics")
              est_Carp = est_D_Carp(
                Ux,
                frequency(Ux),
                D,
                stabs,
                ntau = ntau,
                bins = bins,
                bw_sd = bw_sd,
                interpol_steps = interpol_steps
              )
              
              # Save results
              out = c(as.list(environment()))  # Gather environment
              
              # Filter out original high resolution timeseries (due to large size)
              # Filter out "out" variable as to not save duplicate of output
              variables_to_exclude = c("Ux_original", "Ux_original_burned", "out", "filtered_out", "cl", "for_par", "forloop", "start_time", "functions")
              filtered_out = out[!names(out) %in% variables_to_exclude]
              
              filtered_out[unlist(lapply(filtered_out, class)) == "function"] = NULL # Remove functions
              
              # Save the filtered variables to an RDS file
              saveRDS(filtered_out, paths$filepath_out)
              # saveRDS(out, paths$filepath_out)
              
              }
              out = readRDS(paths$filepath_out)
              
              new_plot_overview(out, paths$filepath_image)
              
              graphics.off()
            }
  })
}

end_time <- Sys.time()
execution_time <- end_time - start_time
print(execution_time)

# Define log file path and append the execution time to the log file
log_file <- "execution_time_log.txt"
write(
  paste(
    "Execution Time for condition: N =",
    forloop[[1]]$N,
    ", sf =",
    forloop[[1]]$sf,
    "for",
    length(forloop),
    "iterations is",
    round(difftime(end_time, start_time, units = "hours"), 3),
    "hours or",
    round(difftime(end_time, start_time, units = "mins"), 3),
    "minutes"
  ),
  file = log_file,
  append = TRUE
)

# End cluster
parallel::stopCluster(cl)  
