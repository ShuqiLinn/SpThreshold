library(ggplot2)

results <- readRDS("radon_results.rds")

draws <- list(
  Spatial = do.call(rbind, results$draws$spatial),
  Nonspatial = do.call(rbind, results$draws$iid)
)

parameters <- c("floor", "log_uranium")

densities <- do.call(rbind, lapply(parameters, function(parameter) {
  bandwidth <- max(vapply(draws, function(x) {
    bw.nrd0(x[[parameter]])
  }, numeric(1)))
  
  limits <- range(unlist(lapply(draws, function(x) {
    x[[parameter]]
  }))) + c(-3, 3) * bandwidth
  
  do.call(rbind, lapply(names(draws), function(model) {
    d <- density(
      draws[[model]][[parameter]],
      bw = bandwidth,
      from = limits[1],
      to = limits[2],
      n = 1024
    )
    data.frame(
      parameter = parameter,
      model = model,
      value = d$x,
      density = d$y
    )
  }))
}))

densities$parameter <- factor(densities$parameter, levels = parameters)
densities$model <- factor(densities$model, levels = names(draws))

labels <- c(
  "floor" = "'Floor'~(beta[1])",
  "log_uranium" = "'Log uranium'~(beta[2])"
)

p <- ggplot(
  densities,
  aes(value, density, colour = model, linetype = model)
) +
  geom_line(linewidth = 0.85) +
  facet_wrap(
    ~parameter,
    nrow = 1,
    scales = "free",
    labeller = as_labeller(labels, label_parsed)
  ) +
  scale_colour_manual(
    values = c(Spatial = "#0072B2", Nonspatial = "#D55E00")
  ) +
  scale_linetype_manual(
    values = c(Spatial = "solid", Nonspatial = "dashed")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Regression coefficient",
    y = "Posterior density",
    colour = NULL,
    linetype = NULL
  ) +
  theme_bw(base_size = 12, base_family = "serif") +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "white", colour = "grey70"),
    strip.text = element_text(size = 12),
    axis.text = element_text(size = 10, colour = "black"),
    legend.position = "bottom",
    legend.key.width = grid::unit(1.4, "cm"),
    panel.spacing = grid::unit(0.8, "lines")
  )

print(p)

ggsave(
  "radon_posterior_densities.pdf",
  p,
  width = 8,
  height = 3.3,
  device = grDevices::pdf,
  family = "serif",
  useDingbats = FALSE,
  bg = "white"
)