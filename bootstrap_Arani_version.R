### Bootstrap 
# load example data
example <-
  readRDS(
    "/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/Rinn_est_scale_check/right-fp-gains-dominance/constant-D2/1000-N/D2strength0.3000_sf10_N1000_iter0005_step0001_pars-3.00_0.00_3.00_0.00_0.00_0.00_0.30_bins100_ntau3_interpol50_bw0.30.RDS"
  )

# get random states

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
  # not smoothed
  # estD1 <- coef(lm(D1 ~ mean_bin + I(mean_bin^2) + I(mean_bin^3), weights = 1/eD1)) 
  # estD2 <- coef(lm(D2 ~ mean_bin+ I(mean_bin^2), weights = 1/eD2)) 
  
  bw <- bw_sd * sd(Ux)
  eD1s = ksmooth(x= mean_bin,y= eD1,kernel= 'normal',bandwidth= bw,
                 x.points=mean_bin)
  eD2s = ksmooth(x= mean_bin ,y= eD2,kernel= 'normal',bandwidth=bw,
                 x.points=mean_bin)
  
  # fit simple model
  # estD1_simple <- coef(lm(D1s$y ~ 0 + D1s$x + I(D1s$x^3), weights = 1/eD1s$y)) # specify only x and x^3?
  # estD2_simple <- coef(lm(D2s$y ~ 1, weights = 1/eD2s$y))
  
  # fit "complex" polynomial (smoothed)
  estD1 <- coef(lm(D1s$y ~ D1s$x + I(D1s$x^2) + I(D1s$x^3), weights = 1/eD1s$y))
  estD2 <- coef(lm(D2s$y ~ 1, weights = 1/eD2s$y))
  
  return(list(
    d13 = round(estD1[4], digits = 3),
    d12 = round(estD1[3], digits = 3),
    d11 = round(estD1[2], digits = 3),
    d10 = round(estD1[1], digits = 3),
    d22 = 0,
    d21 = 0,
    d20 = round(estD2[1], digits = 3)
    # d13 = estD1_simple[2],
    # d12 = 0,
    # d11 = estD1_simple[1],
    # d10 = 0,
    # d22 = 0,
    # d21 = 0,
    # d20 = estD2_simple[1]
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
y_pred_D2 <- estD2[1]

#predicted values bootstrap
y_pred_D1_boot <- bootstrapDDs$d10 + bootstrapDDs$d11 * x_values + bootstrapDDs$d12 * x_values^2 + bootstrapDDs$d13 * x_values^3
y_pred_D2_boot <- bootstrapDDs$d20

plot(example$est_Carp$DD$DDLangout$D1s$x, example$est_Carp$DD$DDLangout$D1s$y, xlab = "mean_bin", ylab = "D1", main = "Polynomial Regression") #actual estimated
points(x_values, y_pred_D1_boot, col = "red")
# lines(x_values, y_pred_D1, col = "red")
# lines(x_values, y_pred_D1_simple, col = "pink")
points(example$est_Carp$DD$DDLangout$D2s$x, example$est_Carp$DD$DDLangout$D2s$y, col = "blue") #actual estimated
abline(h = y_pred_D2_boot, col = "purple")
# lines(x_values, y_pred_D2, col = "blue")
abline(h = y_pred_D2, col = "purple")
# abline(h = y_pred_D2_simple, col = "purple")

# bootstrap: simulate time series from the reconstructed coefficients
bootstrap_Lang <-
  function(bootstrap_n,
           bootstrapDDs,
           N,
           sf,
           ntau,
           bins,
           bw_sd,
           interpol_steps,
           compl_df) {
    execution_times <- list()
    full_boot_est_Carp <- list()
    full_boot_Ux <- list()
    
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
    # plot(random_initial_states)
    
    for (i in 1:bootstrap_n) {
      start_time <- Sys.time()
      
      # set.seed(noise_iter)
      boot_Ux <-
        do.call(timeseries1D, utils::modifyList(
          bootstrapDDs,
          list(
            N = N * sf,
            sf = sf,
            startpoint = random_initial_states[i]
          )
        ))
      
      if (any(is.na(boot_Ux) | is.infinite(boot_Ux))) {
        print("Timeseries went to infinity or contains NaN values. Try other parameter settings.")
      } else {
        print("Timeseries simulated successfully")
      }
      
      full_boot_Ux[[i]] <- boot_Ux
      
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
      
      end_time <- Sys.time()
      execution_times[[i]] <- end_time - start_time
      
      full_boot_est_Carp[[i]] <- boot_est_Carp
    }
    
    return(
      list(
        full_boot_est_Carp = full_boot_est_Carp,
        full_boot_Ux =  full_boot_Ux,
        execution_times = execution_times
      )
    )
  }

#debug
bootstrap_n = 20
bootstrapDDs = estDDs
N = example$N
sf = example$sf
ntau = example$ntau
bins = example$bins
bw_sd = example$bw_sd 
interpol_steps = example$interpol_steps
noise_iter = example$noise_iter
compl_df = example$est_Carp$compl_df
estimated_prob_dist = example$est_Carp$compl_df %>%
  filter(variable == "wts", source == "Estimated")
estimated_prob_dist$value <- as.numeric(estimated_prob_dist$value)

rm(list = c("bootstrap_n", "N", "sf", "ntau", "bins", "bw_sd", "interpol_steps", "noise_iter", "compl_df", "bw", "D1", "D2", "i", "D", "boot_D1", "D1s","D2s","eD1s", "eD2s","reconstructed_D1s_per_bin", "bootstrapDDs", "estDDs"))

estDDs <- list(
  d13 = round(estD1[4], digits = 2),
  d12 = round(estD1[3], digits = 2),
  d11 = round(estD1[2], digits = 2),
  d10 = round(estD1[1], digits = 2),
  d22 = 0,
  d21 = 0,
  d20 = round(estD2[1], digits = 2)
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

# D1 all coefficients, constant D2
bootstrap_complex <-
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

saveRDS(bootstrap_largeN_with_Ux, paste0(filepath_base, "_bootstrap_largeN_with_Ux.RDS"))

bootstrap_largeN_with_Ux <- 
  bootstrap_Lang(bootstrap_n = 20,
                 bootstrapDDs = bootstrapDDs,
                 example$N,
                 example$sf, 
                 example$ntau, 
                 example$bins, 
                 example$bw_sd, 
                 example$interpol_steps,
                 example$est_Carp$compl_df,
                 example$noise_iter)

### checking how to adjust for bias 
bootstrap_demo <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_demo.RDS") # this I ran 100 times, the parent model was reconstructed based on NOT smoothed D1's and D2's 
bootstrap_demo2 <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_demo2.RDS") # 20 times, the parent model is the theoretical model
bootstrap_demo_smoothed<- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_demo_smoothed.RDS") # 20 times, the parent model was reconstructed based on smoothed D1's and D2's 
bootstrap_largeN <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_largeN.RDS") # N = 5000; simple model
bootstrap_complex <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_complex_N1000.RDS") #N = 1000 sf = 10
bootstrap_largeN_complex <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_largeN_with_Ux.RDS") # N = 5000; complex model
bootstrap_complex_N1000 <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_complex_N1000.RDS") # N1000, complex parent, with Ux

bootstrap_complex_N1000_random_noise <- readRDS("/Users/luizashen58/Library/CloudStorage/OneDrive-UvA/Exit Time Project Master's Thesis/exit-time-thesis/exit-time-thesis_bootstrap_complex_N1000_random_noise.RDS")

bootstrap_demo$full_boot_est_Carp[[1]]$fp_df
bootstrap_demo2$full_boot_est_Carp[[1]]$fp_df
bootstrap_demo_smoothed$full_boot_est_Carp[[1]]$fp_df
bootstrap_complex$full_boot_est_Carp[[1]]$fp_df
bootstrap_complex_N1000$full_boot_est_Carp[[1]]$fp_df
bootstrap_largeN$full_boot_est_Carp[[1]]$fp_df
bootstrap_largeN_with_Ux$full_boot_est_Carp[[1]]$fp_df
bootstrap_largeN_complex$full_boot_est_Carp[[1]]$fp_df

# bootstrap mean exit times from parent model reconstructed from NOT smoothed D1's and D2's
boot_meanETl <- vector()
for(i in 1:100){
  boot_meanETl[i] <- bootstrap_demo$full_boot_est_Carp[[i]]$meanETl
}
confidence_interval <- quantile(boot_meanETl, c(0.025, 0.975))
confidence_interval # the confidence interval doesn't include the estimated mean exit time
example$est_Carp$meanETl

ggplot(data = data.frame(x = boot_meanETl[1:20]), aes(x = x)) +
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

ggplot(data = data.frame(x = boot_meanETr[1:20]), aes(x = x)) +
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
  boot_meanETl[i] <- bootstrap_complex_N1000_random_noise$full_boot_est_Carp[[i]]$meanETl
}
confidence_interval <- quantile(boot_meanETl, c(0.025, 0.975))
example$est_Carp$meanETl

ggplot(data = data.frame(x = boot_meanETl), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  # geom_density(color = "darkblue", size = 0.5) +
  geom_vline(
    aes(xintercept = example$est_Carp$meanETl),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  labs(title = "Distribution bootstrapped mean exit time (left)", x = "Values", y = "Density")

boot_meanETr <- vector()
for(i in 1:20){
  boot_meanETr[i] <- bootstrap_complex_N1000_random_noise$full_boot_est_Carp[[i]]$meanETr
}
confidence_interval <- quantile(boot_meanETr, c(0.025, 0.975))
example$est_Carp$meanETr

ggplot(data = data.frame(x = boot_meanETr), aes(x = x)) +
  # geom_density(fill = "skyblue", color = "black") +
  geom_histogram(fill = "skyblue",
                 color = "black",
                 bins = 20) +
  # geom_density(color = "darkblue", size = 0.5) +
  geom_vline(
    aes(xintercept = example$est_Carp$meanETr),
    linetype = "dashed",
    color = "red",
    alpha = 0.5
  ) +
  labs(title = "Distribution bootstrapped mean exit time (right)", x = "Values", y = "Density")

##### how to adjust bias for D1?
boot_D1s <- list()
for(i in 1:100){
  boot_D1s[[i]] <- bootstrap_demo$full_boot_est_Carp[[i]]$DD$DDLangout$D1s$y
}

transposed_D1s <- list()
median_D1s_reconstructed <- vector()
for(i in seq_along(boot_D1s[[1]])) {
  transposed_D1s[[i]] <- sapply(boot_D1s, function(x) x[[i]])
  median_D1s_reconstructed[i] <- median(transposed_D1s[[i]])
}

plot(example$est_Carp$DD$DDLangout$D1s$y, main = "bootstrapped D1s (not smoothed)")
points(median_D1s_reconstructed, col = "red", pch = 4)

reconstructed_D1s_perbin <- data.frame(group = rep(1:example$bins, each = 100), value = unlist(transposed_D1))

# Combine reconstructed D1's per bin with parent D1 estimates
reconstructed_D1s_perbin <- data.frame(group = rep(1:example$bins, each = 100), value = unlist(transposed_D1), parent_D1 = rep(example$est_Carp$DD$DDLangout$D1, each = 100))

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

# bootstrapped D1 and D2 estimates per bin from reconstructed models (reconstructed on smoothed parent)
library(plotrix)
boot_D2s_y <- list()
for(i in 1:20){
  boot_D2s_y[[i]] <- bootstrap_complex_N1000_random_noise$full_boot_est_Carp[[i]]$DD$DDLangout$D2s$y
}

max_length <- max(sapply(boot_D2s_y, length))
boot_D2s_y <- lapply(boot_D2s_y, function(vec) {
  if (length(vec) < max_length) {
    c(vec, rep(NA, max_length - length(vec)))
  } else {
    vec
  }
})

transposed_D2s_y <- list()
median_D2s_y_reconstructed <- vector()
lower_D2s_y_reconstructed <- vector()
upper_D2s_y_reconstructed <- vector()
for(i in seq_along(boot_D2s_y[[1]])) {
  # for(i in seq(22,76)) {
  transposed_D2s_y[[i]] <- sapply(boot_D2s_y, function(x) x[[i]])
  median_D2s_y_reconstructed[i] <- median(transposed_D2s_y[[i]])
  lower_D2s_y_reconstructed[i] <- quantile(transposed_D2s_y[[i]], 0.025, na.rm = TRUE)
  upper_D2s_y_reconstructed[i] <- quantile(transposed_D2s_y[[i]], 0.975, na.rm = TRUE)
}

plot_D2s <-
  data.frame(
    bins = c(1:100),
    median_D2s_y_reconstructed = median_D2s_y_reconstructed,
    upper_D2s_y_reconstructed = upper_D2s_y_reconstructed,
    lower_D2s_y_reconstructed = lower_D2s_y_reconstructed,
    original_est_D2s = c(example$est_Carp$DD$DDLangout$D2s$y, NA, NA, NA)
  )

boot_D1s_y <- list()
for(i in 1:20){
  boot_D1s_y[[i]] <- bootstrap_complex_N1000_random_noise$full_boot_est_Carp[[i]]$DD$DDLangout$D1s$y
}

max_length <- max(sapply(boot_D1s_y, length))
boot_D1s_y <- lapply(boot_D1s_y, function(vec) {
  if (length(vec) < max_length) {
    c(vec, rep(NA, max_length - length(vec)))
  } else {
    vec
  }
})

transposed_D1s_y <- list()
median_D1s_y_reconstructed <- vector()
lower_D1s_y_reconstructed <- vector()
upper_D1s_y_reconstructed <- vector()
for(i in seq_along(boot_D1s_y[[1]])) {
  # for(i in seq(22,76)) {
  transposed_D1s_y[[i]] <- sapply(boot_D1s_y, function(x) x[[i]])
  median_D1s_y_reconstructed[i] <- median(transposed_D1s_y[[i]])
  lower_D1s_y_reconstructed[i] <- quantile(transposed_D1s_y[[i]], 0.025, na.rm = TRUE)
  upper_D1s_y_reconstructed[i] <- quantile(transposed_D1s_y[[i]], 0.975, na.rm = TRUE)
}

plot_D1s <-
  data.frame(
    bins = c(1:100),
    median_D1s_y_reconstructed = median_D1s_y_reconstructed,
    upper_D1s_y_reconstructed = upper_D1s_y_reconstructed,
    lower_D1s_y_reconstructed = lower_D1s_y_reconstructed,
    original_est_D1s = c(example$est_Carp$DD$DDLangout$D1s$y, NA, NA, NA)
  )

plotCI(plot_D2s$bins, median_D2s_y_reconstructed, ui = upper_D2s_y_reconstructed, li = lower_D2s_y_reconstructed, main = "median D2s per bin", xlab = "bins", ylab = "estimated D2's", col = "red", scol = "red", pch = 1, cex = 0.5)
points(example$est_Carp$DD$DDLangout$D2s$y, col = "black", pch = 4)

plotCI(plot_D1s$bins, median_D1s_y_reconstructed, ui = upper_D1s_y_reconstructed, li = lower_D1s_y_reconstructed, main = "median D1s per bin", xlab = "bins", ylab = "estimated D1's", col = "red", scol = "red", pch = 1, cex = 0.5)
points(example$est_Carp$DD$DDLangout$D1s$y, col = "black", pch = 4, cex = 0.5)

ggplot(plot_D2s,
       aes(bins, median_D2s_y_reconstructed)) + geom_point(color = "red") + geom_point(aes(bins, original_est_D2s), color = "black") +
  geom_errorbar(aes(ymin = lower_D2s_y_reconstructed, ymax = upper_D2s_y_reconstructed, color = "red"))+theme_minimal() + labs(color = "median estimated D2's per bin") + xlab("bins") + ylab("estimated D2's")

ggplot(plot_D1s,
       aes(bins, median_D1s_y_reconstructed)) + geom_point(color = "red", shape = 1, show.legend = FALSE) + geom_point(aes(bins, original_est_D1s), color = "black") +
  geom_errorbar(aes(ymin = lower_D1s_y_reconstructed, ymax = upper_D1s_y_reconstructed, color = "red"))+theme_minimal() + xlab("bins") + ylab("estimated D1's")

plot(example$est_Carp$DD$DDLangout$D2s$y, main = "median D2s per bin", cex = 0.8, xlab = "bins", ylab = "estimated D2's", pch = 4)
points(median_D2s_y_reconstructed, col = "red", pch = 1, cex =0.8)

reconstructed_D1s_perbin <- data.frame(group = rep(1:example$bins, each = 100), value = unlist(transposed_D1))

# Combine reconstructed D1's per bin with parent D1 estimates
reconstructed_D1s_perbin <- data.frame(group = rep(1:example$bins, each = 100), value = unlist(transposed_D1), parent_D1 = rep(example$est_Carp$DD$DDLangout$D1, each = 100))

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

# Density plot of D1 estimates per bin from large N reconstructed models (reconstructed on parent model, smoothed) 
boot_D1 <- list()
for(i in 1:100){
  boot_D1[[i]] <- bootstrap_largeN$full_boot_est_Carp[[i]]$DD$DDLangout$D1[72:171]
}

transposed_D1 <- list()
median_D1_reconstructed <- vector()
for(i in seq_along(boot_D1[[1]])) {
  transposed_D1[[i]] <- sapply(boot_D1, function(x) x[[i]])
  median_D1_reconstructed[i] <- median(transposed_D1[[i]])
}
plot(example$est_Carp$DD$DDLangout$D1, main = "bootstrapped D1s")
points(median_D1_reconstructed, col = "red", pch = 4)

reconstructed_D1s_perbin <- data.frame(group = rep(1:100, each = 100), value = unlist(transposed_D1))

# Combine reconstructed D1's per bin with parent D1 estimates
reconstructed_D1s_perbin <- data.frame(group = rep(1:100, each = 100), value = unlist(transposed_D1), parent_D1 = rep(example$est_Carp$DD$DDLangout$D1[72:171], each = 100))

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

library(tidyr)
# plot bootstrapped Ux's 
full_boot_Ux <- data.frame(
  time = 1:length(as.vector(example$Ux))/example$sf/60,
  bootUx1 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[1]],
  bootUx2 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[2]],
  bootUx3 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[6]],
  bootUx4 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[4]],
  bootUx5 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[5]],
  bootUx6 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[6]],
  bootUx6 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[7]],
  bootUx6 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[8]],
  bootUx6 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[9]],
  bootUx6 = bootstrap_complex_N1000_random_noise$full_boot_Ux[[10]])

full_boot_Ux_long <- full_boot_Ux %>% 
  pivot_longer(
    cols = -time, 
    names_to = "series", 
    values_to = "x"
  )

ggplot2::ggplot(data.frame(
  time = 1:length(example$Ux) / example$sf / 60,
  x = as.vector(example$Ux)[1:length(Ux)]
)) +
  ggplot2::geom_point(
    ggplot2::aes(x = time, y = x),
    alpha = .75,
    col = "gray30",
    size = sz_point
  )

# Example plotting code using facet_grid
ggplot(full_boot_Ux_long, aes(x = time, y = x)) +
  geom_point(
    ggplot2::aes(x = time, y = x), size = 0.02, color = "gray30") +
  facet_grid(rows = vars(series)) +  # Arrange plots in a grid by series
  labs(x = "time", y = "x", title = "Faceted Time Series Plot") +
  theme_minimal()

# Plot density
pl_dens <- ggplot2::ggplot(data.frame(x = as.vector(Ux))) +
  ggplot2::geom_histogram(
    ggplot2::aes(x = x, y = ggplot2::after_stat(density)),
    bins = bins,
    colour = "lightgrey",
    fill = "white"
  ) +
  ggplot2::geom_density(ggplot2::aes(x = x), linewidth = 0.8)  +
  ggplot2::geom_vline(
    xintercept = stabs$fps,
    color = col_FP,
    linetype = 'dashed',
    linewidth = lwd_FP,
    alpha = .75
  ) +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::labs(
    x = latex2exp::TeX("$x$"),
    y = "Density",
    title = "",
    subtitle = sprintf(
      "Density (%d mode%s)",
      length(LaplacesDemon::Modes(Ux)$modes),
      ifelse(length(LaplacesDemon::Modes(Ux)$modes) > 1, "s", "")
    )
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