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
filepath_figs = file.path(filepath_base, "SLURM_test_figs") # Snellius
filepath_est = file.path(filepath_base, "SLURM_test_est") # Snellius
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
  nr_steps_bif = 3, #length.out for deepening, asymmetry, etc.
  type_D2 = c("constant"),
  #"quadratic",
  scenario = c("2fps-balanced-deepening"), #"left-fp-gains-dominance" "2fps-balanced-deepening", "right-fp-gains-dominance" 
  strength_D2 = c(.3),
  sf = c(10), #c(10, 100),
  N = c(500), #c(500, 100000),
  bins = c(50), #c(30, 40, 100),
  interpol_steps = 50,# c(50, 100, 500),
  ntau = c(2), # c(3, 5, 10),
  bw_sd = c(0.3),
  #10000
  noise_iter = c(1:3), #c(1:5)
  # noise_iter_concat_times = 5
) %>% purrr::transpose() %>% unique()

# Check each item in the forloop list and keep only the ones where N*sf = bins*100
# filtered_forloop <- list()
# for (i in seq_along(forloop)) {
#   if (forloop[[i]]$N * forloop[[i]]$sf == forloop[[i]]$bins * 100) {
#     filtered_forloop[[i]] <- forloop[[i]]
#   }
# }
# # Remove NULL elements from the list
# filtered_forloop <- filtered_forloop[!sapply(filtered_forloop, is.null)]

# Debug
#datagen = "Langevin"
nr_steps_bif = 3

step_idx = 1
for_par=forloop[[1]]
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
# noise_iter_concat_times = for_par$noise_iter_concat_times
rm(list = c("nr_steps_bif", "step_idx", "for_par", "type_D2", "scenario", "bins", "strength_D2", "bw_sd", "ntau", "interpol_steps", "sf", "N", "noise_iter", "add_to_x", "noise_iter_concat_times"))

# Functions
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
    "get_resilience_potential",
    "get_resilience_prob_dist",
    "get_stability",
    "get_theoretical_D",
    "get_weights",
    "gg",
    "new_plot_overview",
    "poly3",
    "setup_filepaths",
    "style_plot"
  )

packages <-
  c(
    "bvpSolve",
    "cubature",
    "stats",
    "Langevin",
    "dplyr",
    "ggplot2",
    "parallel",
    "doParallel",
    "foreach",
    "cowplot",
    "ggnewscale",
    "latex2exp"
  )

variables <-
  c(
    "nr_steps_bif",
    "type_D2",
    "scenario",
    "strength_D2",
    "sf_high_res",
    "N_high_res",
    "sf",
    "N",
    "bins",
    "interpol_steps",
    "ntau",
    "bw_sd",
    "noise_iter"
  )


# # Set up cluster (forking)
cl <- parallel::makeForkCluster(detectCores() - 1) # using forking

# # Set up cluster (no forking)
# cl <- parallel::makeCluster(detectCores() - 1) # no forking
# cl <- makeCluster(no_cores, outfile = "cluster_log.txt")

doParallel::registerDoParallel(cl)
getDoParWorkers() #check how many workers 'foreach' is going to use

# Loop through scenarios
start_time <- Sys.time()
foreach(for_par = forloop) %do% {
  with(for_par, {
    print(as.data.frame(for_par))
    
    # Get polynomial coefficients for forloop
    Ds = get_D(nr_steps_bif,
               scenario,
               type_D2,
               strength_D2)
    
    # Loop through steps in bifurcation parameter
    foreach(
      step_idx = 1:nr_steps_bif,
      D = Ds,
      .packages = c(packages),
      .export = c(functions, "filepath_est", "filepath_figs")
      # .export = c(functions)
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
          sf_high_res = 0,
          N_high_res = 0,
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
      
      if (!file.exists(paths$filepath_out)) {
      
      # Generate timeseries
      print("Generate timeseries")
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
      print("Estimate Langevin model & calculate exit times and alt. metrics")
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
      variables_to_exclude = c(
        "cl",
        "for_par",
        "forloop",
        "execution_time",
        "filtered_out",
        "start_time",
        "end_time",
        "functions",
        "packages",
        "variables"
      )
      
      filtered_out = out[!names(out) %in% variables_to_exclude]
      
      filtered_out[unlist(lapply(filtered_out, class)) == "function"] = NULL # Remove functions
      # out[unlist(lapply(out, class)) == "function"] = NULL # Remove functions
      saveRDS(filtered_out, paths$filepath_out)
      
      }
      out = readRDS(paths$filepath_out)
      # # print("Plot results")
      # # # new_plot_overview(out, paths$filepath_image, plot_t = ifelse(N*sf < 100000, Inf, 100000))
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
write(paste("Execution Time for condition: N =", forloop[[1]]$N, ", sf =", forloop[[1]]$sf, "for", length(forloop), "iterations is", round(difftime(end_time, start_time, units = "hours"), 3), "hours or", round(difftime(end_time, start_time, units = "mins"), 3), "minutes"), file = log_file, append = TRUE)

parallel::stopCluster(cl) # End cluster  
