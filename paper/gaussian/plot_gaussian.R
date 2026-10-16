library(ggplot2)
library(dplyr)

m_values <- c(2, 5, 10, 20, 50, 80, 100, 200)
rho_values <- c(0.05, 0.50, 0.95)
tau_values <- c(0.05, 0.50, 0.95)
colours <- c(C1 = "#E69F00", C2 = "#56B4E9", C3 = "#009E73")

variance_label <- expression(frac(
  group("|", Var(beta[1]*"|"*bold(Y)*","*rho) -
          Var(beta[1]*"|"*bold(Y)*","*0), "|"),
  Var(beta[1]*"|"*bold(Y)*","*0)))

mean_label <- expression(group("|",
                               E(beta[1]*"|"*bold(Y)*","*rho) - E(beta[1]*"|"*bold(Y)*","*0), "|"))

plot_comparison <- function(s, metric, y_label, reference, reference_colour) {
  s$estimate <- s[[metric]]
  s$se <- s[[paste0(metric, "_se")]]
  observed <- filter(s, method == "MCMC")
  
  ggplot(s, aes(m_pos, estimate, colour = covariate,
                group = interaction(covariate, method))) +
    geom_hline(
      yintercept = reference, colour = reference_colour,
      linetype = "dashed", linewidth = 0.5
    ) +
    geom_errorbar(
      data = observed,
      aes(ymin = pmax(0, estimate - 1.96 * se),
          ymax = estimate + 1.96 * se),
      width = 0.15, linewidth = 0.4, show.legend = FALSE
    ) +
    geom_line(
      aes(linetype = method), linewidth = 0.7,
      show.legend = c(colour = FALSE, linetype = TRUE)
    ) +
    geom_point(data = observed, aes(shape = covariate), size = 2) +
    facet_grid(rho_label ~ k_label, labeller = label_parsed) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = c(C1 = 16, C2 = 17, C3 = 15)) +
    scale_linetype_manual(
      values = c(MCMC = "solid", Exact = "dashed"),
      labels = c(MCMC = "Marginal MCMC", Exact = "Exact conditional")
    ) +
    scale_x_continuous(
      breaks = seq_along(m_values), labels = m_values
    ) +
    expand_limits(y = 0) +
    labs(x = "Sample Size per Location (m)", y = y_label) +
    guides(
      colour = guide_legend(
        order = 1,
        override.aes = list(
          shape = c(16, 17, 15), linetype = 0, size = 3
        )
      ),
      shape = "none",
      linetype = guide_legend(
        order = 2,
        override.aes = list(colour = "grey20", linewidth = 0.8)
      )
    ) +
    theme_bw(base_size = 20, base_family = "serif") +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.key.width = grid::unit(1.8, "cm"),
      strip.background = element_rect(fill = "white", colour = "grey80"),
      strip.text = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 14, colour = "grey20"),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}

plots <- list()

for (n_locs in c(25, 100, 400)) {
  d <- readRDS(paste0("n", n_locs, "_gaussian.rds")) |>
    filter(
      m %in% m_values,
      rho_true %in% rho_values,
      tau0_2_true %in% tau_values
    ) |>
    mutate(
      covariate = factor(covariate, levels = names(colours)),
      m_pos = match(m, m_values),
      rho_label = factor(
        rho_true, levels = rho_values,
        labels = paste0("rho == ", rho_values)
      ),
      k_label = factor(
        tau0_2_true, levels = tau_values,
        labels = sprintf(
          "tau[0]^2/sigma^2 == %.2f/%.2f",
          tau_values, 1 - tau_values
        )
      )
    )
  
  s <- bind_rows(
    d |> transmute(
      rho_label, k_label, covariate, m_pos, method = "MCMC",
      variance_difference = beta1_marginal_variance_relative_abs_difference,
      mean_difference = beta1_marginal_mean_abs_difference
    ),
    d |> transmute(
      rho_label, k_label, covariate, m_pos, method = "Exact",
      variance_difference = beta1_conditional_variance_relative_abs_difference,
      mean_difference = beta1_conditional_mean_abs_difference
    )
  ) |>
    group_by(rho_label, k_label, covariate, m_pos, method) |>
    summarise(
      variance = mean(variance_difference),
      variance_se = sd(variance_difference) / sqrt(n()),
      mean = mean(mean_difference),
      mean_se = sd(mean_difference) / sqrt(n()),
      .groups = "drop"
    ) |>
    mutate(method = factor(method, levels = c("MCMC", "Exact")))
  
  thresholds <- d |>
    filter(covariate == "C2", m == 200) |>
    group_by(rho_label, k_label) |>
    summarise(
      m_approx = ceiling(mean(m_approx)),
      m_exact = ceiling(mean(m_exact)),
      .groups = "drop"
    ) |>
    mutate(
      approx_pos = approx(m_values, seq_along(m_values), m_approx)$y,
      exact_pos = approx(m_values, seq_along(m_values), m_exact)$y,
      approx_label = paste0(
        "bar(m)[approx]^'*' == ", m_approx,
        ifelse(m_approx > 200, " ~ '(>200)'", "")
      ),
      exact_label = paste0(
        "bar(m)[exact]^'*' == ", m_exact,
        ifelse(m_exact > 200, " ~ '(>200)'", "")
      )
    )
  
  observed <- filter(s, method == "MCMC")
  curve_top <- max(
    s$variance,
    observed$variance + 1.96 * observed$variance_se
  )
  y_max <- max(1, ceiling(5 * curve_top / 0.75) / 5)
  
  p_variance <- plot_comparison(
    s, "variance", variance_label, 0.05, "red"
  ) +
    geom_vline(
      data = filter(thresholds, !is.na(approx_pos)),
      aes(xintercept = approx_pos),
      colour = colours["C2"], linetype = "longdash", linewidth = 0.65
    ) +
    geom_vline(
      data = filter(thresholds, !is.na(exact_pos)),
      aes(xintercept = exact_pos),
      colour = "grey35", linetype = "dotdash", linewidth = 0.65
    ) +
    geom_text(
      data = thresholds,
      aes(x = Inf, label = paste0("bold(", approx_label, ")")),
      y = 0.92 * y_max, inherit.aes = FALSE, parse = TRUE,
      hjust = 1.05, family = "serif", colour = colours["C2"],
      size = 4.0
    ) +
    geom_text(
      data = thresholds,
      aes(x = Inf, label = paste0("bold(", exact_label, ")")),
      y = 0.80 * y_max, inherit.aes = FALSE, parse = TRUE,
      hjust = 1.05, family = "serif", colour = "grey35",
      size = 4.0
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      breaks = seq(0, y_max, 0.2),
      labels = function(x) sprintf("%.1f", x),
      expand = expansion(mult = c(0.02, 0))
    )
  
  p_mean <- plot_comparison(
    s, "mean", mean_label, 0, "grey50"
  ) +
    scale_y_continuous(
      labels = function(x) format(x, trim = TRUE, scientific = FALSE),
      expand = expansion(mult = c(0.02, 0.15))
    )
  
  plots[[paste0("variance_n", n_locs)]] <- p_variance
  plots[[paste0("mean_n", n_locs)]] <- p_mean
  
  print(p_variance)
  print(p_mean)
  
  ggsave(
    paste0("gaussian_variance_n", n_locs, ".pdf"),
    p_variance, width = 12, height = 8,
    device = grDevices::pdf, family = "serif", useDingbats = FALSE
  )
  ggsave(
    paste0("gaussian_mean_n", n_locs, ".pdf"),
    p_mean, width = 12, height = 8,
    device = grDevices::pdf, family = "serif", useDingbats = FALSE
  )
}