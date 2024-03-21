##### Alternative resilience metrics #####
filepath_base_demo <-
  "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Modified_for_Luiza_2024_02_07/Alt_metrics_demo_11_03_2024"
setwd(filepath_base_demo)

file_names <-
  c(
    "Mult_noise_D2strength0.3000_sf10_N1000_iter0001_step0001_pars-1.00_0.00_1.00_0.00_0.00_0.00_0.30_bins100_ntau10_interpol100_bw0.30.RDS",
    "Mult_noise_D2strength0.3000_sf10_N1000_iter0001_step0002_pars-1.50_0.00_1.50_0.00_0.00_0.00_0.30_bins100_ntau10_interpol100_bw0.30.RDS",
    "D2strength0.3000_sf10_N1000_iter0001_step0002_pars-1.50_0.00_1.50_0.00_0.00_0.00_0.30_bins100_ntau10_interpol100_bw0.30.RDS"
  )

loaded_data <- list()
counter = 1

# Loop through each file and load the data
for (file in file_names) {
  # Use readRDS to load the RDS file
  data <- readRDS(file)
  
  # Add the loaded data to the list
  loaded_data[[counter]] <- data
  
  counter = counter + 1
}

# Distance to basin threshold, potential depth
get_resilience_potential <- function(fp_df, compl_df) {
  # get distance to basin threshold (theoretical)
  theoretical_fps = fp_df %>%
    filter(source == "Theoretical") %>%
    pull(xintercept)
  
  basin_width_theo_right = abs(theoretical_fps[3] - theoretical_fps[2]) # Dist_thr = |x_{att} - x_{saddle}|
  basin_width_theo_left = abs(theoretical_fps[1] - theoretical_fps[2])
  
  # get distance to basin threshold (estimated)
  estimated_fps = fp_df %>%
    filter(source == "Estimated") %>%
    pull(xintercept)
  
  basin_width_est_right = abs(estimated_fps[3] - estimated_fps[2]) # Dist_thr = |x_{att} - x_{saddle}|
  basin_width_est_left = abs(estimated_fps[1] - estimated_fps[2])
  
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

# demonstrate use
get_resilience_potential(fp_df = loaded_data[[2]]$est_Carp$fp_df,
                         compl_df = loaded_data[[2]]$est_Carp$compl_df) # concat multiple noise

get_resilience_potential(fp_df = loaded_data[[3]]$est_Carp$fp_df,
                         compl_df = loaded_data[[3]]$est_Carp$compl_df) # noise = 0.3

# Get variance and skewness around dominant modes of probability distribution
get_resilience_prob_dist <- function(fp_df, compl_df) {
  # get left vs. right basins of stationary probability distribution (theoretical)
  theoretical_prob_dist = compl_df %>%
    filter(variable == "wts", source == "Theoretical")
  theoretical_saddle = fp_df %>%
    filter(source == "Theoretical") %>%
    slice(2) %>%
    pull(xintercept) # get saddle point
  theoretical_saddle_index = sapply(theoretical_saddle, function(fp)
    which.min(abs(theoretical_prob_dist$x - fp)))
  
  theoretical_prob_dist_right =
    theoretical_prob_dist %>% slice(1:theoretical_saddle_index)
  theoretical_prob_dist_left =
    theoretical_prob_dist %>% slice((theoretical_saddle_index + 1):n())
  
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
  
  # calculate the skewness around the mode (theoretical)
  theo_skewness_mode_right <-
    sum(prob_values_theo_R * (x_values_theo_R - theo_mode_right) ^ 3) / (sum(prob_values_theo_R * (x_values_theo_R - theo_mode_right) ^
                                                                               2)) ^ (3 / 2)
  theo_skewness_mode_left <-
    sum(prob_values_theo_L * (x_values_theo_L - theo_mode_left) ^ 3) / (sum(prob_values_theo_L * (x_values_theo_L - theo_mode_left) ^
                                                                              2)) ^ (3 / 2)
  
  # calculate the skewness around the mode (estimated)
  est_skewness_mode_right <-
    sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^ 3) / (sum(prob_values_est_R * (x_values_est_R - est_mode_right) ^
                                                                            2)) ^ (3 / 2)
  est_skewness_mode_left <-
    sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^ 3) / (sum(prob_values_est_L * (x_values_est_L - est_mode_left) ^
                                                                           2)) ^ (3 / 2)
  
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

# demonstrate use
get_resilience_prob_dist(fp_df = loaded_data[[2]]$est_Carp$fp_df,
                         compl_df = loaded_data[[2]]$est_Carp$compl_df) # concat multiple noise

get_resilience_prob_dist(fp_df = loaded_data[[3]]$est_Carp$fp_df,
                         compl_df = loaded_data[[3]]$est_Carp$compl_df) # noise = 0.3
