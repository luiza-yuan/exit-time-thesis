# Application of Carpenter's (2022) exit time estimation to simulated data using Rinn's (2016) Langevin package
# Kyra Evers, University of Amsterdam, 17.01.2024

######### Modifications by Luiza Yuan, 29.01.2024 #########
# (1) modified get_exit_time function ('yini' and 'yend' set up for bvpSolve for left vs. right basins)
# (2) modified get_D function(from specifying alpha/beta to specifying d10,d11,d12,d13)
# (3) modified plotting function (new_plot_overview function)
# (4) added get_theoretical_exit_time function (included in est_D_Carp function)
# (5) notes at the bottom regarding scaling and theoretical exit time estimation issues (still in progress)
#################################

rm(list = ls())
graphics.off()

# Create necessary directories
filepath_base = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(filepath_base)
filepath_figs = file.path(filepath_base, "figs_domain_check") #figs_theoreticalET_scalecheck")
filepath_est = file.path(filepath_base, "est_domain_check")# "est_theoreticalET_scalecheck")
if (!dir.exists(filepath_figs)) {
  dir.create(filepath_figs, recursive = T)
}

library(bvpSolve)
library(cubature)
library(stats)
library(Langevin)
library(dplyr)
library(ggplot2)
library(foreach)
# library(tuneR)
library(cowplot)
library(ggnewscale)
library(latex2exp)

# source("ExitTime_BinMethod_PeterLakeExample-main/DDbintau.R")
source(file.path(dirname(filepath_base), "ExitTime_BinMethod_PeterLakeExample-main/DDbintau.R"))
source("helper_functions.R")

# Choose parameters to loop through
forloop = tidyr::expand_grid(
  type_D2 = c("constant"),
  #"quadratic",
  scenario = c("2fps-balanced-deepening"),
  # "2fps-balanced-shallowing",
  # "left-fp-gains-dominance",
  # "right-fp-gains-dominance"),
  strength_D2 = c(.3),
  sf = 2, #c(10, 100),
  N = 100, #c(500, 100000),
  bins = 100, #c(30,40,100),
  interpol_steps = 20,# c(50,100),
  ntau = 10, # c(3,5, 10),
  bw_sd = .3,
  #10000
  noise_iter = c(1) #c(1:5)
) %>% purrr::transpose() %>% unique()



# Parameters
datagen = "Langevin"
nr_steps_bif = 5

# Debug
step_idx = 1
for_par=forloop[[1]]
type_D2 = for_par$type_D2
scenario = for_par$scenario
bins = for_par$bins
strength_D2 = for_par$strength_D2
bw_sd = for_par$bw_sd
ntau= for_par$ntau
interpol_steps = for_par$interpol_steps#100 #30, 60
sf = for_par$sf
N = for_par$N
noise_iter = for_par$noise_iter
add_to_x = .5

# Set up cluster
# cl <- parallel::makeCluster(4)
# doParallel::registerDoParallel(cl)

# Loop through scenarios
foreach(for_par = forloop) %do% {
  with(for_par, {
    print(as.data.frame(for_par))

    # Get polynomial coefficients for forloop
    Ds = get_D(nr_steps_bif,
               scenario,
               type_D2,
               strength_D2)

    # Loop through steps in bifurcation parameter
    foreach(step_idx = 1:nr_steps_bif,
            D = Ds) %do% {
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
                  step_idx = step_idx, bins=bins, ntau=ntau, interpol_steps=interpol_steps, bw_sd=bw_sd
                )
              ))

              if (!file.exists(paths$filepath_out)) {
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
# parallel::stopCluster(cl) # End cluster



# Notes:

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
