library(ggplot2)
library(dplyr)

d <- readRDS("n100_gaussian.rds")

m_values <- c(2, 5, 10, 20, 50, 80, 100, 200)
rho_values <- c(0.05, 0.50, 0.95)
tau_values <- c(0.05, 0.50, 0.95)
colours <- c(C1 = "#E69F00", C2 = "#56B4E9", C3 = "#009E73")
z <- qnorm(0.975)

d <- d |>
  filter(m %in% m_values, rho_true %in% rho_values,
         tau0_2_true %in% tau_values)

s <- bind_rows(
  d |> transmute(m, rho_true, tau0_2_true, covariate,
                 model = "Spatial", error = beta1_spatial_error,
                 covered = beta1_spatial_covered),
  d |> transmute(m, rho_true, tau0_2_true, covariate,
                 model = "Nonspatial", error = beta1_iid_error,
                 covered = beta1_iid_covered)
) |>
  group_by(rho_true, tau0_2_true, covariate, m, model) |>
  summarise(
    R = n(),
    bias = mean(error),
    bias_mcse = sd(error) / sqrt(R),
    coverage = mean(covered),
    coverage_mcse = sqrt(coverage * (1 - coverage) / R),
    .groups = "drop"
  ) |>
  mutate(
    bias_lower = bias - z * bias_mcse,
    bias_upper = bias + z * bias_mcse,
    coverage_lower = pmax(0, coverage - z * coverage_mcse),
    coverage_upper = pmin(1, coverage + z * coverage_mcse),
    covariate = factor(covariate, levels = names(colours)),
    model = factor(model, levels = c("Spatial", "Nonspatial")),
    rho_label = factor(rho_true, levels = rho_values,
                       labels = paste0("rho == ", rho_values)),
    k_label = factor(tau0_2_true, levels = tau_values,
                     labels = sprintf("tau[0]^2/sigma^2 == %.2f/%.2f",
                                      tau_values, 1 - tau_values)),
    m_pos = match(m, m_values) + ifelse(model == "Spatial", -0.07, 0.07)
  )

make_plot <- function(metric, y_label, reference, reference_colour, limits) {
  ggplot(s, aes(m_pos, .data[[metric]], colour = covariate,
                group = interaction(covariate, model))) +
    geom_hline(yintercept = reference, colour = reference_colour,
               linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(
      aes(ymin = .data[[paste0(metric, "_lower")]],
          ymax = .data[[paste0(metric, "_upper")]]),
      width = 0.12, linewidth = 0.4, show.legend = FALSE
    ) +
    geom_line(aes(linetype = model), linewidth = 0.7,
              show.legend = c(colour = FALSE, linetype = TRUE)) +
    geom_point(aes(shape = covariate), size = 2) +
    facet_grid(rho_label ~ k_label, labeller = label_parsed) +
    scale_colour_manual(values = colours) +
    scale_shape_manual(values = c(C1 = 16, C2 = 17, C3 = 15)) +
    scale_linetype_manual(values = c(Spatial = "solid", Nonspatial = "dashed")) +
    scale_x_continuous(breaks = seq_along(m_values), labels = m_values) +
    scale_y_continuous(
      limits = limits, breaks = pretty(limits, n = 5),
      labels = function(x) format(x, trim = TRUE, scientific = FALSE),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    labs(x = "Sample Size per Location (m)", y = y_label) +
    guides(
      colour = guide_legend(order = 1,
                            override.aes = list(shape = c(16, 17, 15), linetype = 0, size = 3)),
      shape = "none",
      linetype = guide_legend(order = 2,
                              override.aes = list(colour = "grey20", linewidth = 0.8))
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

bias_limit <- max(abs(c(s$bias_lower, s$bias_upper)))
coverage_min <- max(0, floor(20 * min(s$coverage_lower, 0.90)) / 20)

p_bias <- make_plot(
  "bias", expression("Bias of " * hat(beta)[1]),
  0, "grey50", c(-bias_limit, bias_limit)
)

p_coverage <- make_plot(
  "coverage", "Coverage", 0.95, "red", c(coverage_min, 1)
)

print(p_bias)
print(p_coverage)

ggsave("gaussian_bias_n100.pdf", p_bias, width = 12, height = 8,
       device = grDevices::pdf, family = "serif", useDingbats = FALSE)
ggsave("gaussian_coverage_n100.pdf", p_coverage, width = 12, height = 8,
       device = grDevices::pdf, family = "serif", useDingbats = FALSE)