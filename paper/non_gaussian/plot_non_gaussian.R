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
    ggplot(s, aes(m_pos, estimate, colour = covariate, group = covariate)) +
        geom_hline(yintercept = reference, colour = reference_colour,
                   linetype = "dashed", linewidth = 0.5) +
        geom_errorbar(aes(ymin = pmax(0, estimate - 1.96 * se),
                          ymax = estimate + 1.96 * se),
                      width = 0.15, linewidth = 0.4, show.legend = FALSE) +
        geom_line(linewidth = 0.7, show.legend = FALSE) +
        geom_point(aes(shape = covariate), size = 2) +
        facet_grid(rho_label ~ tau_label, labeller = label_parsed) +
        scale_colour_manual(values = colours) +
        scale_shape_manual(values = c(C1 = 16, C2 = 17, C3 = 15)) +
        scale_x_continuous(breaks = seq_along(m_values), labels = m_values) +
        expand_limits(y = 0) +
        labs(x = "Sample Size per Location (m)", y = y_label) +
      guides(
        colour = guide_legend(
          override.aes = list(shape = c(16, 17, 15), size = 3)
        ),
        shape = "none"
      ) +
        theme_bw(base_size = 20, base_family = "serif") +
        theme(legend.position = "bottom", legend.title = element_blank(),
              legend.key.width = grid::unit(1.8, "cm"),
              strip.background = element_rect(fill = "white", colour = "grey80"),
              strip.text = element_text(face = "bold", size = 16),
              axis.text = element_text(size = 14, colour = "grey20"),
              panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
              panel.grid.minor = element_blank())
}

plots <- list()

for (family in c("bernoulli", "nb")) {
    headers <- sprintf("tau[0]^2 == %.2f", tau_values)
    if (family == "nb")
        headers <- sprintf("tau[0]^2 == %.2f ~ ',' ~ r == %d", tau_values, c(1, 10, 100))

    d <- readRDS(paste0("n100_", family, ".rds")) |>
        filter(m %in% m_values, rho_true %in% rho_values,
               tau0_2_true %in% tau_values) |>
        mutate(covariate = factor(covariate, levels = names(colours)),
               m_pos = match(m, m_values),
               rho_label = factor(rho_true, levels = rho_values,
                                  labels = sprintf("rho == %.2f", rho_values)),
               tau_label = factor(tau0_2_true, levels = tau_values, labels = headers))

    stopifnot(all(is.finite(d$beta1_marginal_mean_abs_difference)),
              all(is.finite(d$beta1_marginal_variance_relative_abs_difference)))

    s <- d |>
        group_by(rho_label, tau_label, covariate, m_pos) |>
        summarise(datasets = n(),
            variance = mean(beta1_marginal_variance_relative_abs_difference),
            variance_se = sd(beta1_marginal_variance_relative_abs_difference) / sqrt(n()),
            mean = mean(beta1_marginal_mean_abs_difference),
            mean_se = sd(beta1_marginal_mean_abs_difference) / sqrt(n()), .groups = "drop")
    stopifnot(nrow(s) == 216L, all(s$datasets == 100L))

    y_max <- max(1, 1.08 * max(s$variance + 1.96 * s$variance_se))
    p_variance <- plot_comparison(s, "variance", variance_label, 0.05, "red") +
        scale_y_continuous(limits = c(0, y_max),
            labels = function(x) format(x, trim = TRUE, scientific = FALSE),
            expand = expansion(mult = c(0.02, 0)))
    p_mean <- plot_comparison(s, "mean", mean_label, 0, "grey50") +
        scale_y_continuous(labels = function(x) format(x, trim = TRUE, scientific = FALSE),
                           expand = expansion(mult = c(0.02, 0.15)))

    plots[[paste0(family, "_variance")]] <- p_variance
    plots[[paste0(family, "_mean")]] <- p_mean
    print(p_variance)
    print(p_mean)

    ggsave(paste0(family, "_variance_n100.pdf"), p_variance,
           width = 12, height = 8, device = grDevices::pdf,
           family = "serif", useDingbats = FALSE)
    ggsave(paste0(family, "_mean_n100.pdf"), p_mean,
           width = 12, height = 8, device = grDevices::pdf,
           family = "serif", useDingbats = FALSE)
}
