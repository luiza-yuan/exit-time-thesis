# Time series simulation and estimation of mean exit time, potential depth, basin width, and variance around maxima using the Langevin approach
# When simulating the time series, we first generate a high resolution time series close to the conditions of Arani et al. (2021), with 12000 time units with time step of 0.0001. This forms the "complete" timeseries (e.g. one can imagine this as the complete time series of an individual). From here, we downsample before conducting the analysis. This is similar to in real life, when we take only a sample of the complete time series when collecting data.

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
filepath_figs = file.path(filepath_base, "27_conditions_figs")
filepath_est = file.path(filepath_base, "27_conditions_est")
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
source(file.path(filepath_base, "helper_functions_July.R"))

# Choose parameters to loop through
forloop = tidyr::expand_grid(
  # number of bifurcation steps (i.e. number of times basins are deepened)
  nr_steps_bif = 3,
  # constant or quadratic diffusion function
  type_D2 = c("constant"),
  # scenario options: "2fps-balanced-deepening" (balanced double basin), "right-fp-gains-dominance" (asymmetrical basin depth)
  scenario = c("2fps-balanced-deepening", "right-fp-gains-dominance"),
  # noise intensity
  strength_D2 = c(.2),
  # number of time units for simulating timeseries
  N_high_res = c(12000),
  # sampling frequency (i.e. 1 / number of timesteps) for simulating timeseries
  sf_high_res = c(1000),
  # burn first 10 timeunits (i.e. 10*1000 = 10000 datapoints)
  N_burned = c(10),
  # sampling frequency 4 (interval = 250); sampling frequency 8 (interval = 125); sampling frequency 10 (interval = 100);
  interval = c(250, 125, 100),
  # c(250, 125, 100)
  timeseries_length = c(7000, 12000, 20000),
  # c(7000, 10000, 60000)
  bins = c(10),
  #c(10, 50),
  interpol_steps = 100,
  # c(50, 100, 500),
  ntau = c(2),
  # c(2, 3, 4),
  add_to_x = 0.5,
  # 0.5
  bw_sd = 0.3,
  noise_iter = c(401:500),
) %>% purrr::transpose() %>% unique()

# Functions
functions <-
  c(
    "D1fun",
    "D2fun",
    "downsample_Langevin",
    "est_D_Carp",
    "finner",
    "generate_Langevin",
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

# Set up cluster (forking)
cl <- parallel::makeForkCluster(detectCores() - 1) # using forking
doParallel::registerDoParallel(cl)

getDoParWorkers() #check how many workers 'foreach' is going to use

# Get start time
start_time <- Sys.time()

# Loop through scenarios
foreach(for_par = forloop) %do% {
  with(for_par, {
    print(as.data.frame(for_par))
    
    # Get polynomial coefficients for generating timeseries from Langevin model
    Ds = get_D(nr_steps_bif, scenario, type_D2, strength_D2)
    
    # Loop through steps in bifurcation parameter
    foreach(
      step_idx = 1:nr_steps_bif,
      D = Ds,
      .export = c(functions)
    ) %dopar% {
      # Setup and create filepaths
      paths = do.call(setup_filepaths, utils::modifyList(
        D,
        list(
          filepath_est = filepath_est,
          filepath_figs = filepath_figs,
          scenario = scenario,
          type_D2 = type_D2,
          strength_D2 = strength_D2,
          timeseries_length = timeseries_length,
          interval = interval,
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
        # Arani et al. (2021) 12,000 time units with time step of 0.001
        print("Generate true timeseries (high resolution)")
        Ux_original = do.call(generate_Langevin, utils::modifyList(
          D,
          list(
            N = N_high_res,
            sf = sf_high_res,
            noise_iter = noise_iter
          )
        ))
        
        # Burn first 10000 datapoints (burn)
        Ux_original_burned = window(Ux_original,
                                    start = c(N_burned + 1, 1),
                                    end = c(12000, 1000))
        
        # Downsample from original high resolution timeseries to get final Ux
        Ux <- downsample_Langevin(
          Ux = Ux_original_burned,
          interval = interval,
          timeseries_length = timeseries_length
        )
        
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
          interpol_steps = interpol_steps,
          add_to_x = add_to_x
        )
        
        # Save results
        out = c(as.list(environment()))  # Gather environment
        
        # Filter out original high resolution timeseries (due to large size)
        # Filter out "out" variable as to not save duplicate of output
        variables_to_exclude = c(
          "Ux_original",
          "Ux_original_burned",
          "out",
          "filtered_out",
          "cl",
          "for_par",
          "forloop",
          "start_time",
          "functions"
        )
        filtered_out = out[!names(out) %in% variables_to_exclude]
        
        filtered_out[unlist(lapply(filtered_out, class)) == "function"] = NULL # Remove functions
        
        # Save the filtered variables to an RDS file
        saveRDS(filtered_out, paths$filepath_out)
        
      }
    }
  })
}

end_time <- Sys.time() # Get end time
execution_time <- end_time - start_time # Get execution time

# Define log file path and append the execution time to the log file
#log_file <- "execution_time_log.txt"
#write(
# paste(
#  "Execution Time for condition: timeseries_length =",
# forloop[[1]]$timeseries_length,
# ", interval =",
# forloop[[1]]$interval,
# "for",
# length(forloop),
# "iterations is",
# round(difftime(end_time, start_time, units = "hours"), 3),
# "hours or",
# round(difftime(end_time, start_time, units = "mins"), 3),
# "minutes"
# ),
# file = log_file,
# append = TRUE
#)

# End cluster
parallel::stopCluster(cl)
