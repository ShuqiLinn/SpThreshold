library(ggplot2)
library(dplyr)

results          <- readRDS("simulation_results.rds")

m_levels  <- c(1, 2, 5, 10, 20, 50, 80, 100, 200)
m_pos_map <- setNames(1:length(m_levels), m_levels)

m_to_pos <- function(m_val){
   if(is.na(m_val)) return(NA_real_)
   if(is.infinite(m_val)) return(NA_real_)
   if(m_val > max(m_levels)) return(NA_real_)
   if(m_val <= min(m_levels)) return(1)
   approx(m_levels, 1:length(m_levels), xout = m_val)$y
   }

color_palette <- c("var0"   = "#E69F00",
                   "var1"   = "#56B4E9",
                   "varInf" = "#009E73")

shape_palette <- c("var0"   = 16,
                   "var1"   = 17,
                   "varInf" = 15)

x_var_legend_labels <- c("var0"   = "C1",
                         "var1"   = "C2",
                         "varInf" = "C3")

gamma_threshold <- 0.05


add_facet_labels <- function(df){
   df %>%
      mutate(
         m_pos = m_pos_map[as.character(m)],
         rho_label = factor(
            paste0("rho == ", rho_true),
            levels = c("rho == 0.05", "rho == 0.5", "rho == 0.95")),
         k_label = factor(
            case_when(
               abs(k - 0.0526) < 0.01 ~ "tau^2/sigma^2 == 0.05/0.95",
               abs(k - 1)      < 0.01 ~ "tau^2/sigma^2 == 0.50/0.50",
               abs(k - 19)     < 0.1  ~ "tau^2/sigma^2 == 0.95/0.05",
               TRUE ~ NA_character_),
            levels = c("tau^2/sigma^2 == 0.05/0.95",
                       "tau^2/sigma^2 == 0.50/0.50",
                       "tau^2/sigma^2 == 0.95/0.05")),
         x_var_label = factor(
            case_when(
               x_loc_var == 0 ~ "var0",
               x_loc_var == 1 ~ "var1",
               is.infinite(x_loc_var) ~ "varInf"),
            levels = c("var0", "var1", "varInf")))
   }


## ----------------------------------------------------------------------------
## Plot 1: relative variance difference, with m* line under C2
## ----------------------------------------------------------------------------

rel_var_summary <- results %>%
   group_by(n_locs, m, rho_true, k, x_loc_var) %>%
   summarise(mean_rel_var = mean(rel_var_diff, na.rm = TRUE),
             sd_rel_var   = sd(rel_var_diff,   na.rm = TRUE),
             n_obs        = n(),
             .groups      = "drop") %>%
   mutate(se_rel_var = sd_rel_var / sqrt(n_obs),
          ci_lower   = mean_rel_var - 1.96 * se_rel_var,
          ci_upper   = mean_rel_var + 1.96 * se_rel_var)

plot_data_1 <- add_facet_labels(rel_var_summary)


## m* vertical line metadata (only computed under C2, which is x_loc_var == 1)
m_star_vline <- results %>%
   filter(x_loc_var == 1) %>%
   distinct(n_locs, rho_true, k, x_loc_var, m_star_approx_asymp) %>%
   mutate(
      m_star_approx_asymp = pmax(m_star_approx_asymp, 2),
      out_of_range        = m_star_approx_asymp > max(m_levels),
      m_star_pos          = ifelse(out_of_range, NA, sapply(m_star_approx_asymp, m_to_pos)),
      near_right_edge     = !out_of_range & !is.na(m_star_pos) & m_star_pos > 7.5,
      m_star_label = case_when(
         out_of_range ~ paste0("m* = ", round(m_star_approx_asymp), " \u2192"),
         TRUE         ~ paste0("m* = ", round(m_star_approx_asymp))),
      label_hjust = case_when(
         out_of_range    ~ 1,
         near_right_edge ~ 1.1,
         TRUE            ~ -0.1),
      label_x = case_when(
         out_of_range ~ length(m_levels) + 0.5,
         TRUE         ~ m_star_pos),
      rho_label = factor(
         paste0("rho == ", rho_true),
         levels = c("rho == 0.05", "rho == 0.5", "rho == 0.95")),
      k_label = factor(
         case_when(
            abs(k - 0.0526) < 0.01 ~ "tau^2/sigma^2 == 0.05/0.95",
            abs(k - 1)      < 0.01 ~ "tau^2/sigma^2 == 0.50/0.50",
            abs(k - 19)     < 0.1  ~ "tau^2/sigma^2 == 0.95/0.05"),
         levels = c("tau^2/sigma^2 == 0.05/0.95",
                    "tau^2/sigma^2 == 0.50/0.50",
                    "tau^2/sigma^2 == 0.95/0.05")))


for(n_val in c(25, 100, 400)){

   plot_data_n <- plot_data_1 %>% filter(n_locs == n_val)
   y_max       <- max(plot_data_n$ci_upper, na.rm = TRUE) * 1.2

   m_star_vline_n        <- m_star_vline %>% filter(n_locs == n_val)
   m_star_vline_in_range <- m_star_vline_n %>% filter(!out_of_range)

   p1 <- ggplot(plot_data_n,
                aes(x = m_pos, y = mean_rel_var,
                    group = x_var_label, color = x_var_label,
                    fill = x_var_label)) +
      geom_hline(aes(yintercept = gamma_threshold,
                     linetype = "gamma == 0.05"),
                 color = "red", linewidth = 0.5) +
      geom_vline(data = m_star_vline_in_range,
                 aes(xintercept = m_star_pos, linetype = "m* (C2)"),
                 color = "#56B4E9", linewidth = 0.6) +
      geom_text(data = m_star_vline_n,
                aes(x = label_x, y = y_max,
                    label = m_star_label, hjust = label_hjust),
                color = "#56B4E9", vjust = 1, size = 5,
                inherit.aes = FALSE) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.15, linewidth = 0.5) +
      geom_line(linewidth = 0.6) +
      geom_point(aes(shape = x_var_label), size = 2) +
      facet_grid(rho_label ~ k_label, labeller = label_parsed) +
      scale_x_continuous(breaks = 1:length(m_levels),
                         labels = m_levels,
                         expand = expansion(mult = c(0.02, 0.05))) +
      scale_y_continuous(limits = c(NA, y_max)) +
      scale_color_manual(name = NULL,
                         values = color_palette,
                         labels = x_var_legend_labels) +
      scale_fill_manual(name = NULL,
                        values = color_palette,
                        labels = x_var_legend_labels) +
      scale_shape_manual(name = NULL,
                         values = shape_palette,
                         labels = x_var_legend_labels) +
      scale_linetype_manual(
         name = NULL,
         values = c("gamma == 0.05" = "dashed",
                    "m* (C2)"       = "dashed"),
         labels = c("gamma == 0.05" = expression(gamma == 0.05),
                    "m* (C2)"       = "m* (C2)"),
         guide = guide_legend(override.aes = list(color = c("red", "#56B4E9")))) +
      coord_cartesian(xlim = c(0.5, 9.5)) +
      theme_bw(base_size = 18, base_family = "serif") +
      theme(
         legend.position  = "bottom",
         legend.text      = element_text(size = 14),
         legend.box       = "horizontal",
         strip.background = element_rect(fill = "white", color = "grey80"),
         strip.text       = element_text(face = "bold", size = 16),
         axis.title       = element_text(size = 16),
         axis.text        = element_text(size = 14, color = "grey20"),
         panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
         panel.grid.minor = element_blank()) +
      guides(color    = guide_legend(order = 1),
             shape    = guide_legend(order = 1),
             fill     = guide_legend(order = 1),
             linetype = guide_legend(order = 2)) +
      labs(x = "Sample Size per Location (m)",
           y = expression(frac(group("|", Var(beta[1]*"|"*bold(Y)*","*rho) -
                                        Var(beta[1]*"|"*bold(Y)*","*0), "|"),
                               Var(beta[1]*"|"*bold(Y)*","*0))))

   print(p1)

   ggsave(filename = paste0("plot1_rel_var_diff_J", n_val, ".png"),
          plot     = p1,
          width    = 12,
          height   = 8,
          dpi      = 300)
   cat(sprintf("Saved plot1_rel_var_diff_J%d.png\n", n_val))
   }


## ----------------------------------------------------------------------------
## Plot 2: absolute mean difference
## ----------------------------------------------------------------------------

mean_diff_summary <- results %>%
   group_by(n_locs, m, rho_true, k, x_loc_var) %>%
   summarise(mean_abs_mean_diff = mean(abs(bias_mean), na.rm = TRUE),
             sd_abs_mean_diff   = sd(abs(bias_mean),   na.rm = TRUE),
             n_obs              = n(),
             .groups            = "drop") %>%
   mutate(se_abs_mean_diff = sd_abs_mean_diff / sqrt(n_obs),
          ci_lower         = mean_abs_mean_diff - 1.96 * se_abs_mean_diff,
          ci_upper         = mean_abs_mean_diff + 1.96 * se_abs_mean_diff)

plot_data_2 <- add_facet_labels(mean_diff_summary)


for(n_val in c(25, 100, 400)){

   plot_data_n <- plot_data_2 %>% filter(n_locs == n_val)

   p2 <- ggplot(plot_data_n,
                aes(x = m_pos, y = mean_abs_mean_diff,
                    group = x_var_label, color = x_var_label,
                    fill = x_var_label)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = "grey50", linewidth = 0.5) +
      geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                    width = 0.15, linewidth = 0.5) +
      geom_line(linewidth = 0.6) +
      geom_point(aes(shape = x_var_label), size = 2) +
      facet_grid(rho_label ~ k_label, labeller = label_parsed) +
      scale_x_continuous(breaks = 1:length(m_levels),
                         labels = m_levels,
                         expand = expansion(mult = c(0.02, 0.05))) +
      scale_color_manual(name = NULL,
                         values = color_palette,
                         labels = x_var_legend_labels) +
      scale_fill_manual(name = NULL,
                        values = color_palette,
                        labels = x_var_legend_labels) +
      scale_shape_manual(name = NULL,
                         values = shape_palette,
                         labels = x_var_legend_labels) +
      coord_cartesian(xlim = c(0.5, 9.5)) +
      theme_bw(base_size = 18, base_family = "serif") +
      theme(
         legend.position  = "bottom",
         legend.text      = element_text(size = 14),
         legend.box       = "horizontal",
         strip.background = element_rect(fill = "white", color = "grey80"),
         strip.text       = element_text(face = "bold", size = 16),
         axis.title       = element_text(size = 16),
         axis.text        = element_text(size = 14, color = "grey20"),
         panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
         panel.grid.minor = element_blank()) +
      labs(x = "Sample Size per Location (m)",
           y = expression(group("|", E(beta[1]*"|"*bold(Y)*","*rho) -
                                   E(beta[1]*"|"*bold(Y)*","*0),"|")))

   print(p2)

   ggsave(filename = paste0("plot2_mean_diff_J", n_val, ".png"),
          plot     = p2,
          width    = 12,
          height   = 8,
          dpi      = 300)
   cat(sprintf("Saved plot2_mean_diff_J%d.png\n", n_val))
   }
