library(ggplot2)

design <- readRDS("n100_nonmonotone.rds")
n <- design$n_locs
lambda <- design$eigenvalues
d <- design$d
within <- n - sum(d)
gamma <- 0.05
m_max <- 1500
c2_colour <- "#56B4E9"
point_m <- c(2, 5, 10, 20, 50, 80, 100, 200, seq(300, m_max, by = 100))

settings <- expand.grid(
  tau0_2 = c(0.05, 0.50, 0.95),
  rho = c(0.05, 0.50, 0.95)
)
settings$sigma2 <- 1 - settings$tau0_2
settings$kappa <- settings$tau0_2 / settings$sigma2
settings$rho_label <- factor(
  sprintf("rho == %.2f", settings$rho),
  levels = sprintf("rho == %.2f", c(0.05, 0.50, 0.95))
)
settings$kappa_label <- factor(
  sprintf("kappa == %.2f/%.2f", settings$tau0_2, settings$sigma2),
  levels = c("kappa == 0.05/0.95", "kappa == 0.50/0.50", "kappa == 0.95/0.05")
)

curves <- vector("list", nrow(settings))
points <- vector("list", nrow(settings))
thresholds <- vector("list", nrow(settings))

for (j in seq_len(nrow(settings))) {
  z <- settings[j, , drop = FALSE]
  c_i <- 1 - z$rho + z$rho * lambda
  r <- 1 / (mean(1 / c_i) * c_i)
  
  exact_gap <- function(m) {
    vapply(m, function(value) {
      t <- value * z$kappa
      numerator <- sum(d * t * (r - 1) / ((1 + t) * (1 + t * r)))
      denominator <- within + sum(d / (1 + t * r))
      abs(numerator) / denominator
    }, numeric(1))
  }
  
  leading <- abs(sum(d * (1 - 1 / r))) / (z$kappa * within)
  tail_bound <- sum(d * abs(1 - 1 / r)) / (z$kappa * within)
  search_end <- max(2L, ceiling(tail_bound / gamma))
  integer_m <- seq.int(2L, search_end)
  violations <- integer_m[exact_gap(integer_m) > gamma]
  m_exact <- max(c(1L, violations)) + 1L
  m_approx <- max(2L, ceiling(leading / gamma))
  
  grid <- sort(unique(c(
    seq(2, 20, length.out = 600),
    seq(20, 200, length.out = 900),
    seq(200, m_max, length.out = 1301),
    m_exact, m_approx
  )))
  
  curves[[j]] <- data.frame(
    z, m = grid, difference = exact_gap(grid), row.names = NULL
  )
  points[[j]] <- data.frame(
    z, m = point_m, difference = exact_gap(point_m), row.names = NULL
  )
  thresholds[[j]] <- data.frame(
    z, m_approx = m_approx, m_exact = m_exact
  )
}

curves <- do.call(rbind, curves)
points <- do.call(rbind, points)
thresholds <- do.call(rbind, thresholds)
rownames(curves) <- rownames(points) <- rownames(thresholds) <- NULL

panel_top <- tapply(curves$difference, curves$rho_label, max)
thresholds$panel_top <- pmax(
  ifelse(thresholds$rho > 0.05, gamma, 0),
  unname(panel_top[as.character(thresholds$rho_label)])
)
thresholds$approx_y <- 1.24 * thresholds$panel_top
thresholds$exact_y <- 1.10 * thresholds$panel_top
thresholds$approx_label <- paste0("m[approx]^'*' == ", thresholds$m_approx)
thresholds$exact_label <- paste0("m[exact]^'*' == ", thresholds$m_exact)

p <- ggplot(curves, aes(m, difference)) +
  geom_hline(
    data = settings[settings$rho > 0.05, ],
    aes(yintercept = gamma),
    colour = "red", linetype = "dashed", linewidth = 0.5
  ) +
  geom_vline(
    data = thresholds, aes(xintercept = m_approx),
    colour = c2_colour, linetype = "longdash", linewidth = 0.65
  ) +
  geom_vline(
    data = thresholds, aes(xintercept = m_exact),
    colour = "grey35", linetype = "dotdash", linewidth = 0.65
  ) +
  geom_line(colour = c2_colour, linewidth = 0.8) +
  geom_point(data = points, colour = c2_colour, shape = 17, size = 1.5) +
  geom_text(
    data = thresholds,
    aes(x = Inf, y = approx_y, label = paste0("bold(", approx_label, ")")),
    inherit.aes = FALSE, parse = TRUE, hjust = 1.08,
    family = "serif", colour = c2_colour, size = 4.0
  ) +
  geom_text(
    data = thresholds,
    aes(x = Inf, y = exact_y, label = paste0("bold(", exact_label, ")")),
    inherit.aes = FALSE, parse = TRUE, hjust = 1.08,
    family = "serif", colour = "grey35", size = 4.0
  ) +
  facet_grid(
    rho_label ~ kappa_label, scales = "free_y", labeller = label_parsed
  ) +
  scale_x_continuous(
    breaks = c(2, 200, 500, 1000, 1500),
    labels = c("2", "200", "500", "1,000", "1,500"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    labels = function(x) format(x, trim = TRUE, scientific = FALSE),
    expand = expansion(mult = c(0.025, 0.10))
  ) +
  expand_limits(y = 0) +
  labs(
    x = "Sample Size per Location (m)",
    y = expression(frac(
      group("|", Var(beta[1]*"|"*bold(Y)*","*rho) -
              Var(beta[1]*"|"*bold(Y)*","*0), "|"),
      Var(beta[1]*"|"*bold(Y)*","*0)))
  ) +
  theme_bw(base_size = 20, base_family = "serif") +
  theme(
    strip.background = element_rect(fill = "white", colour = "grey80"),
    strip.text = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14, colour = "grey20"),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(8, 10, 8, 8)
  )

print(p)

ggsave("nonmonotone.pdf", p, width = 12, height = 8,
       device = grDevices::pdf, family = "serif", useDingbats = FALSE)

print(thresholds[, c("rho", "kappa", "m_approx", "m_exact")], row.names = FALSE)