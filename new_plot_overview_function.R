### new_plot_overview function
# This document contains (1) the older new_plot_overview function, with original comments and changes from the very beginning of the project (2) updated, cleaned up version of the function with unnecessary code removed, commented, with both versions where the theoretica and estimated are plotted on the same facet plots and on different facet plots

example <-
  readRDS(
    "Rinn_est_scale_check/2fps-balanced-deepening/constant-D2/1000-N/D2strength0.3000_sf10_N1000_iter0077_step0003_pars-2.00_0.00_2.00_0.00_0.00_0.00_0.30_bins100_ntau3_interpol50_bw0.30.RDS"
  )

# This is the new_plot_overview function (with all comments) plotting theoretical and estimated separately
new_plot_overview <-
  function(out,
           filepath_image,
           plot_t = c(10000 * 100, Inf)[1]) {
    with(out, {
      #if (datagen == "Langevin") {
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
      # } else if (grepl("noise", datagen, fixed = T)) {
      #   title_plot = latex2exp::TeX(sprintf(
      #     "%s noise; $f_s$ = %d, $N$ = %d",
      #     stringr::str_to_title(sub("_noise", "", datagen)),
      #     sf,
      #     N
      #   ))
      # }
      
      # Plotting parameters
      lwd_FP = 0.8
      col_FP = dplyr::recode(
        stabs$stab_fps,
        "unstable" = scales::viridis_pal(option = "inferno")(50)[29],
        "stable" = scales::viridis_pal(option = "viridis")(50)[45]
      )
      
      # variables if exit time is estimated vs. not:
      if(is.null(est_Carp$meanETl) == FALSE){
        variable.labs = c(
          "drift" = latex2exp::TeX("Drift", output = 'character'),
          "diffusion" = latex2exp::TeX("Diffusion", output =   'character'),
          "potential" = latex2exp::TeX("Potential", output = 'character'),
          "EPF" = latex2exp::TeX("Effective Potential", output = 'character'),
          "ET" = latex2exp::TeX("Exit Time", output = 'character'),
          "wts" = latex2exp::TeX("Weights", output = 'character')
        )
      } else { #else? (not sure if other scenarios are possible, no standardized error messages yet)
        variable.labs = c(
          "drift" = latex2exp::TeX("Drift", output = 'character'),
          "diffusion" = latex2exp::TeX("Diffusion", output =   'character'),
          "potential" = latex2exp::TeX("Potential", output = 'character'),
          "EPF" = latex2exp::TeX("Effective Potential", output = 'character')
        ) # no "ET" and "wts" variables
      }
      
      # variable.labs = c(
      #   "potential_analytical" = latex2exp::TeX("Potential (analytical)", output =
      #                                             'character'),
      #   "negPF" = latex2exp::TeX("- Potential (estimated)", output = 'character'),
      #   "EPF" = latex2exp::TeX("Effective Potential (estimated)", output = 'character'),
      #   "drift_analytical" = latex2exp::TeX("Drift (theoretical)", output = 'character'),
      #   "D1" = latex2exp::TeX("$D_1$ (estimated)", output = 'character'),
      #   "ET" = latex2exp::TeX("Exit Time (estimated)", output = 'character'),
      #   "theo_ET" = latex2exp::TeX("Exit Time (theoretical)", output =
      #                                'character'),
      #   "diff_analytical" = latex2exp::TeX("Diffusion (theoretical)", output =
      #                                        'character'),
      #   "D2" = latex2exp::TeX("$D_2$ (estimated)", output = 'character'),
      #   "wts" = latex2exp::TeX("Weights (estimated)", output = 'character'),
      #   "theo_wts" = latex2exp::TeX("Weights (theoretical)", output = 'character')
      #   # "drift" = latex2exp::TeX("Drift (estimated)", output = 'character'),
      #   # "diff" = latex2exp::TeX("Diffusion (estimated)", output = 'character')
      # )
      
      # Plot generated timeseries 
      plot_t =  min(length(as.vector(Ux)), plot_t)
      #sz_point = ifelse(N * sf <= 10000, 1, ifelse(N * sf <= 1000000, .1, .01))
      sz_point = 0.35
      pl_ts <-
        ggplot2::ggplot(data.frame(t = 1:plot_t / sf / 60, x = as.vector(Ux)[1:plot_t])) +
        ggplot2::geom_hline(
          yintercept = stabs$fps,
          color = col_FP,
          linetype = 'dashed',
          linewidth = lwd_FP,
          alpha = 1
        ) +
        ggplot2::geom_point(
          ggplot2::aes(x = t, y = x),
          alpha = .75,
          col = "gray30",
          size = sz_point
        )    +
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
        ggplot2::scale_x_continuous(expand = c(0, 0))
      
      plot_t_inset = 100 * sf
      pl_ts_inset = ggplot2::ggplot(data.frame(t = 1:plot_t_inset / sf / 60, x = as.vector(Ux)[1:plot_t_inset])) +
        ggplot2::geom_hline(
          yintercept = stabs$fps,
          color = col_FP,
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
      
      pl_ts_with_inset <-
        cowplot::ggdraw() +
        cowplot::draw_plot(style_plot(pl_ts)) +
        cowplot::draw_plot(
          style_plot(pl_ts_inset)  +
            ggplot2::theme(
              # panel.background = element_rect(fill='transparent'),
              plot.background = ggplot2::element_rect(fill = 'transparent', color =
                                                        NA),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank()
            ) +
            ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "pt")),
          x = 0.04,
          y = .12,
          width = .2,
          height = .5
        )
      
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
      
      # Plot theoretical drift, diffusion, and potential function
      xlim = max(abs(c(stabs$fps, Ux)))
      
      # Carpenter (2022)
      # design <-
      #   matrix(c(1:3, NA, 4:11), 3, 4, byrow = TRUE) # Layout facets
      
      pl_Carp =
        est_Carp$compl_df %>%
        # rbind(
        # theoretical_df %>% mutate(err = NA), # %>% mutate(value = ifelse(
        #   # variable != "potential_analytical", value, value
        # # )),
        # est_Carp$est_df_Carp,
        # est_Carp$Pot_df,
        # est_Carp$ET_df,
        # est_Carp$Theoretical_ET_df
        # ) %>%
        dplyr::filter(variable %in% names(variable.labs)) %>%
        dplyr::mutate(
          value = as.numeric(as.character(value)),
          variable_name = dplyr::recode_factor(variable,!!!variable.labs, .ordered = T)
        ) %>%
        ggplot2::ggplot() +
        # Add zero line
        ggplot2::geom_hline(
          yintercept = 0,
          color = 'grey10',
          linetype = "solid",
          linewidth = .2,
          alpha = 1
        ) +
        # Add fixed points
        ggplot2::geom_vline(
          data = est_Carp$fp_df,
          ggplot2::aes(xintercept = xintercept,
                       color = color),
          linetype = 'dashed',
          linewidth = lwd_FP,
          alpha = .75
        ) +
        ggplot2::scale_color_manual(
          name = "",
          breaks = c("unstable", "stable", unique(variable.labs)),
          values = c(
            "unstable" = scales::viridis_pal(option = "inferno")(50)[29],
            "stable" = scales::viridis_pal(option = "viridis")(50)[45]
            # "'Potential (analytical)'" = "black",
            # "'- Potential (estimated)'" = "#CC0000",
            # "'Effective Potential'" = "#CC0000",
            # "'Drift (theoretical)'" = "black",
            # "D[1]*' (estimated)'" = "#CC0000",
            # "'Exit Time (estimated)'" = "#CC0000",
            # "'Exit Time (theoretical)'" = "black",
            # "'Diffusion (theoretical)'" = "black",
            # "D[2]*' (estimated)'" = "#CC0000",
            # "'Weights'" = "#CC0000",
            # "'theo_wts'" = "black"
          ),
          guide = "none"
        ) +
        ggnewscale::new_scale_color() +
        # Add error
        ggplot2::geom_point(ggplot2::aes(x = x, y = err, color = source), shape = 4, color = 'gray40', size = .5) +
        # Add estimated variables
        ggplot2::geom_point(ggplot2::aes(x = x, y = value, color = source), size = .5) +
        ggplot2::scale_color_manual(
          name = "",
          breaks = c("Theoretical", "Estimated"),
          values = c(
            "Theoretical" = 'grey10',
            "Estimated" = 'red3'),
          guide = "none"
        ) +
        ggplot2::scale_x_continuous(expand = c(0, 0)) +
        ggplot2::scale_y_continuous(position = 'right') +
        # ggh4x::facet_manual(
        #   vars(variable_name),
        #   scales = "free_y",
        #   design = design,
        #   labeller = label_parsed
        # ) +
        ggh4x::facet_grid2(source ~ variable_name,
                           scales = 'free_y', independent = 'y',
                           labeller = label_parsed, switch = "y", axes = "all") +
        # scale_fill_manual(values = "black") +  # Use a single color for fill
        ggplot2::theme(legend.position = "none")
      
      if(is.null(est_Carp$meanETl) == FALSE){
        pl_Carp = pl_Carp +
          ggplot2::labs(
            y = "",
            x = latex2exp::TeX("$x$"),
            title = latex2exp::TeX(
              sprintf(
                "Theoretical coefficients: $D_1(x) = %.2fx^3 + %.2fx^2 + %.2fx$; $D_2(x) = %.2fx^2 + %.2f$ (%d bins, $n_{tau}$ = %d, %d interpolation steps; smoothing = %.2f)",
                D$d13,
                D$d12, 
                D$d11,
                D$d22,
                D$d20,  bins, ntau, interpol_steps, bw_sd
              ),
              bold = T
            ),
            subtitle =
              latex2exp::TeX(
                sprintf(
                  "Theoretical Exit Time $mu_{left} =$ %.2f, $mu_{right} =$ %.2f; Estimated Exit Time $mu_{left} =$ %.2f (tol = %.4f, %sconverged), $mu_{right} =$ %.2f  (tol = %.4f, %sconverged)",
                  est_Carp$theo_meanETl,
                  est_Carp$theo_meanETr,
                  est_Carp$meanETl, est_Carp$atolL, ifelse(as.logical(est_Carp$diagn_ETL$df[["flag"]]),"not ", ""),
                  est_Carp$meanETr, est_Carp$atolR, ifelse(as.logical(est_Carp$diagn_ETL$df[["flag"]]),"not ", "")
                )
              )
          ) 
      } else {
        pl_Carp = pl_Carp +
          ggplot2::labs(
            y = "",
            x = latex2exp::TeX("$x$"),
            title = latex2exp::TeX(
              sprintf(
                "Theoretical coefficients: $D_1(x) = %.2fx^3 + %.2fx^2 + %.2fx$; $D_2(x) = %.2fx^2 + %.2f$ (%d bins, $n_{tau}$ = %d, %d interpolation steps; smoothing = %.2f)",
                D$d13,
                D$d12, 
                D$d11,
                D$d22,
                D$d20,  bins, ntau, interpol_steps, bw_sd
              ),
              bold = T
            ),
            subtitle =
              latex2exp::TeX(
                sprintf(
                  "Theoretical Exit Time $mu_{left} =$ %.2f, $mu_{right} =$ %.2f; Exit Time could not be estimated",
                  est_Carp$theo_meanETl,
                  est_Carp$theo_meanETr
                )
              )
          )
      }
      
      
      # Combine plots
      pl_combo = cowplot::plot_grid(
        cowplot::plot_grid(
          pl_ts_with_inset,
          style_plot(pl_dens) + ggplot2::labs(x = ""),
          rel_widths = c(1, .3),
          nrow = 1
        ),
        NULL,
        # Add space between plots
        style_plot(pl_Carp) +
          ggplot2::theme(legend.position = c(.87, 1)),
        ncol = 1,
        rel_heights = c(.75, .025, 1)
      )
      
      grDevices::graphics.off()
      filepath_image %>%
        grDevices::pdf(width = 16,
                       height = 8,
                       bg = "transparent")
      grid::grid.draw(pl_combo) # Make plot
      grDevices::dev.off()
      grDevices::graphics.off()
      if (file.exists(filepath_image)) {
        print(filepath_image)
      }
      # print(pl_combo)
      
    })
  }

# Demo
## source style_plot first
new_plot_overview(example, example$paths$filepath_image, plot_t = ifelse(N*sf < 100000, Inf, 100000))


# New plot function cleaned up and commented (May 16th) 

# debug
out <- example
filepath_image <- example$paths$filepath_image
plot_t <- ifelse(N*sf < 100000, Inf, 100000)

scenario = example$scenario
type_D2 = example$type_D2
strength_D2 = example$strength_D2
stabs = example$stabs
est_Carp = example$est_Carp
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

# Demo 
new_plot_overview(example, example$paths$filepath_image)

