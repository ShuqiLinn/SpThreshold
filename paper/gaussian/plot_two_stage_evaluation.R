library(ggplot2)
library(dplyr)
library(ggh4x)

d <- readRDS("n100_gaussian.rds")

m_values <- c(2, 5, 10, 20, 50, 80, 100, 200)
rho_values <- c(0.05, 0.50, 0.95)
tau_values <- c(0.05, 0.50, 0.95)
methods <- c("Generating parameters", "Two-stage estimates")

d <- d |>
  filter(covariate == "C2", m %in% m_values,
         rho_true %in% rho_values, tau0_2_true %in% tau_values)

s <- bind_rows(
  d |> transmute(rho_true, tau0_2_true, m, method = methods[1],
                 threshold = two_stage_m_exact_true),
  d |> transmute(rho_true, tau0_2_true, m, method = methods[2],
                 threshold = two_stage_m_exact_estimated)
) |>
  group_by(rho_true, tau0_2_true, m, method) |>
  summarise(
    estimate = mean(threshold),
    lower = quantile(threshold, 0.025, names = FALSE, type = 7),
    upper = quantile(threshold, 0.975, names = FALSE, type = 7),
    .groups = "drop"
  ) |>
  mutate(
    method = factor(method, levels = methods),
    m_pos = match(m, m_values) + ifelse(method == methods[1], -0.06, 0.06),
    rho_label = factor(rho_true, levels = rho_values,
                       labels = sprintf("rho == %.2f", rho_values)),
    k_label = factor(tau0_2_true, levels = tau_values,
                     labels = c("kappa == 0.05/0.95", "kappa == 0.5/0.5",
                                "kappa == 0.95/0.05"))
  )

bounds <- s |>
  group_by(rho_label, k_label) |>
  summarise(top = max(4, estimate, upper), .groups = "drop")

p <- ggplot(s, aes(m_pos, estimate, colour = method,
                   linetype = method, group = method)) +
  geom_blank(data = bounds, aes(x = 1, y = top), inherit.aes = FALSE) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.14,
                linewidth = 0.45, linetype = "solid", show.legend = FALSE) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  facet_grid2(rho_label ~ k_label, scales = "free_y", independent = "y",
              axes = "y", labeller = label_parsed) +
  scale_colour_manual(values = setNames(c("#444444", "#0072B2"), methods)) +
  scale_linetype_manual(values = setNames(c("solid", "dashed"), methods)) +
  scale_x_continuous(breaks = seq_along(m_values), labels = m_values,
                     limits = c(0.65, 8.35), expand = expansion(mult = 0)) +
  scale_y_continuous(
    limits = function(x) c(0, x[2]),
    breaks = function(x) {
      b <- pretty(c(0, max(x)), n = 4)
      b[b >= 0 & b == floor(b)]
    },
    labels = function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE),
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  labs(x = expression("Sample Size per Location (" * m[fit] * ")"),
       y = expression(m[exact]^"*"), colour = NULL, linetype = NULL) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(1.8, "cm"),
    strip.background = element_rect(fill = "white", colour = "grey80"),
    strip.text = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, colour = "grey20"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

print(p)

ggsave("two_stage_evaluation.pdf", p, width = 12, height = 8,
       device = grDevices::pdf, family = "serif", useDingbats = FALSE)