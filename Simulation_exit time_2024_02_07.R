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
# (1) parallel processing (still needs edits to plot when mean exit time is not estimated)
#################################

rm(list = ls())
graphics.off()

# Create necessary directories
filepath_base = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(filepath_base)
filepath_figs = file.path(filepath_base, "figs_sim_sf_change") #figs_theoreticalET_scalecheck")
filepath_est = file.path(filepath_base, "est_sim_sf_change")# "est_theoreticalET_scalecheck")
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
source(file.path(dirname(filepath_base), "ExitTime_BinMethod_PeterLakeExample-main/DDbintau.R"))
source(file.path(filepath_base, "helper_functions.R")) # attention: make sure that it is the version updated 20.03.2024

# Choose parameters to loop through
forloop = tidyr::expand_grid(
  # datagen = "Langevin",
  nr_steps_bif = 5, #length.out for deepening, asymmetry, etc.
  type_D2 = c("constant"),
  #"quadratic",
  scenario = c("2fps-balanced-deepening"),
  # "left-fp-gains-dominance",
  # "right-fp-gains-dominance"),
  strength_D2 = c(.3),
  sf = 40, #c(10, 100),
  N = c(200), #c(500, 100000),
  bins = c(80), #c(30, 40, 100),
  interpol_steps = 100,# c(50, 100, 500),
  ntau = 10, # c(3, 5, 10),
  bw_sd = .3,
  #10000
  noise_iter = c(1) #c(1:5)
) %>% purrr::transpose() %>% unique()

# Debug
#datagen = "Langevin"
nr_steps_bif = 5

step_idx = 2
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
    "apply_DDbintau",
    "D1fun",
    "D2fun",
    "DDbins",
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

# Loop through scenarios
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

                # Get stability of fixed points
                stabs = get_stability(D)

                # Apply Exit Time Analysis (Carpenter's (2022) approach)
                print("Estimate D via Carpenter's (2022) method")
                est_Carp = est_D_Carp(
                  Ux,
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
                saveRDS(out, paths$filepath_out)
              }
              out = readRDS(paths$filepath_out)
              print("Plot results")
              new_plot_overview(out, paths$filepath_image, plot_t = ifelse(N*sf < 100000, Inf, 100000))
              graphics.off()
            }
  })
}
parallel::stopCluster(cl) # End cluster

### attempt Bootstrap 
boot::tsboot(attempt$Ux)
?tsboot()

### Assumption checks ###
# Check stationary
?adf.test
tseries::adf.test(attempt$Ux)

# Check Markov property
spgs::markov.test(attempt$Ux, type = "lb.test")

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
