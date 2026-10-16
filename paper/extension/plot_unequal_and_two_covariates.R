library(ggplot2)
library(dplyr)

unequal <- readRDS("n100_unequal.rds")
two_covariates <- readRDS("n100_two_covariates.rds")

plot_extensions <- function(d, type) {
    unequal <- type == "unequal"
    m_values <- c(2, 5, 10, 20, 50, 80, 100, 200)
    rho_values <- c(0.50, 0.95)
    panel_levels <- if (unequal) c("C1", "C2", "C3") else
        c("C1-C1", "C1-C2", "C1-C3", "C2-C2", "C2-C3", "C3-C3")
    series_levels <- if (unequal) c("Balanced", "Unequal") else c("beta1", "beta2")
    series_labels <- if (unequal) series_levels else expression(beta[1], beta[2])
    colours <- setNames(
        if (unequal) c("#444444", "#0072B2") else c("#0072B2", "#D55E00"),
        series_levels
    )

    d <- d |> filter(m %in% m_values, rho_true %in% rho_values)
    coefficients <- if (unequal) "beta1" else c("beta1", "beta2")

    b <- bind_rows(lapply(coefficients, function(coefficient) {
        data.frame(
            rho_true = d$rho_true,
            m = d$m,
            panel = if (unequal) d$covariate else
                paste(d$covariate1, d$covariate2, sep = "-"),
            series = if (unequal)
                ifelse(d$analysis == "balanced_control", "Balanced", "Unequal") else
                coefficient,
            variance = d[[paste0(coefficient, "_marginal_variance_relative_abs_difference")]],
            mean = d[[paste0(coefficient, "_marginal_mean_abs_difference")]]
        )
    }))

    panel_codes <- unlist(lapply(rho_values, function(rho)
        paste(panel_levels, rho, sep = "_")))
    panel_labels <- unlist(lapply(rho_values, function(rho)
        sprintf('atop("%s", rho == %.2f ~ "," ~ kappa == 0.5/0.5)', panel_levels, rho)))

    index <- if (unequal) 1 else quote(j)
    y_labels <- list(
        variance = bquote(frac(
            group("|", Var(beta[.(index)]*"|"*bold(Y)*","*rho) -
                         Var(beta[.(index)]*"|"*bold(Y)*","*0), "|"),
            Var(beta[.(index)]*"|"*bold(Y)*","*0))),
        mean = bquote(group("|",
            E(beta[.(index)]*"|"*bold(Y)*","*rho) -
            E(beta[.(index)]*"|"*bold(Y)*","*0), "|"))
    )
    references <- c(variance = 0.05, mean = 0)
    plots <- summaries <- list()

    for (metric in names(y_labels)) {
        s <- b |>
            group_by(rho_true, panel, series, m) |>
            summarise(
                R = n(),
                estimate = mean(.data[[metric]]),
                mcse = sd(.data[[metric]]) / sqrt(R),
                .groups = "drop"
            ) |>
            mutate(
                lower = pmax(0, estimate - qnorm(0.975) * mcse),
                upper = estimate + qnorm(0.975) * mcse,
                panel = factor(panel, levels = panel_levels),
                series = factor(series, levels = series_levels),
                rho_label = factor(rho_true, levels = rho_values,
                    labels = sprintf("atop(rho == %.2f, kappa == 0.5/0.5)", rho_values)),
                panel_label = factor(paste(panel, rho_true, sep = "_"),
                                     levels = panel_codes, labels = panel_labels),
                m_pos = match(m, m_values) + 0.10 * (as.integer(series) - 1.5)
            )

        limits <- c(0, 1.05 * max(s$upper, references[[metric]]))

        p <- ggplot(s, aes(m_pos, estimate, colour = series, group = series)) +
            geom_hline(yintercept = references[[metric]], colour = "grey50",
                       linetype = "dotted", linewidth = 0.4) +
            geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.09,
                          linewidth = 0.3, show.legend = FALSE) +
            geom_line(aes(linetype = series), linewidth = 0.5) +
            geom_point(aes(shape = series), size = 1.4) +
            scale_colour_manual(values = colours, labels = series_labels) +
            scale_shape_manual(values = setNames(c(16, 17), series_levels),
                               labels = series_labels) +
            scale_linetype_manual(values = setNames(c("solid", "dashed"), series_levels)) +
            scale_x_continuous(breaks = seq_along(m_values), labels = m_values) +
            scale_y_continuous(
                limits = limits, breaks = pretty(limits, n = 4),
                labels = function(x) format(x, trim = TRUE, scientific = FALSE),
                expand = expansion(mult = c(0.025, 0.04))
            ) +
            guides(
                colour = guide_legend(order = 1, override.aes = list(
                    shape = c(16, 17),
                    linetype = c("solid", "dashed"))),
                shape = "none", linetype = "none"
            ) +
            labs(x = if (unequal) "Expected sample size per location (m)" else
                     "Sample size per location (m)", y = y_labels[[metric]]) +
            theme_bw(base_size = 11, base_family = "serif") +
            theme(
                legend.position = "bottom",
                legend.title = element_blank(),
                legend.key.width = grid::unit(1, "cm"),
                legend.text = element_text(size = 9),
                axis.text = element_text(size = 8, colour = "grey20"),
                strip.text = element_text(size = 9),
                strip.background = element_rect(fill = "white", colour = "grey75", linewidth = 0.4),
                panel.border = element_rect(colour = "grey65", linewidth = 0.4),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25),
                panel.grid.minor = element_blank(),
                panel.spacing = grid::unit(0.8, "lines")
            )

        if (unequal) {
            p <- p + facet_grid(rho_label ~ panel,
                                labeller = labeller(rho_label = label_parsed))
        } else {
            p <- p + facet_wrap(~panel_label, ncol = 3, labeller = label_parsed)
        }

        ggsave(paste0(type, "_", metric, "_n100.pdf"), p,
               width = 7.2, height = if (unequal) 4.8 else 8.2,
               device = grDevices::pdf, family = "serif", useDingbats = FALSE)
        plots[[metric]] <- p
        summaries[[metric]] <- mutate(s, metric = metric)
    }

    list(plots = plots, summary = bind_rows(summaries))
}

unequal_figures <- plot_extensions(unequal, "unequal")
two_covariate_figures <- plot_extensions(two_covariates, "two_covariates")

print(unequal_figures$plots$variance)
print(unequal_figures$plots$mean)
print(two_covariate_figures$plots$variance)
print(two_covariate_figures$plots$mean)
