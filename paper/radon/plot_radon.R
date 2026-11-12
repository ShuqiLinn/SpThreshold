library(ggplot2)
library(patchwork)
library(sf)

homes <- read.csv("srrs2.dat", strip.white = TRUE)
county <- read.csv("cty.dat", strip.white = TRUE)
names(homes) <- tolower(trimws(names(homes)))
names(county) <- tolower(trimws(names(county)))

homes <- homes[trimws(homes$state) == "MN" & homes$activity > 0, ]
homes$fips <- 1000L * homes$stfips + homes$cntyfips
homes$log_radon <- log(homes$activity)

county <- county[county$stfips == 27L, ]
county$fips <- 1000L * county$stfips + county$ctfips
county$log_uranium <- log(county$uppm)
county <- unique(county[c("fips", "log_uranium")])

# Average the log concentrations within each county.
radon_summary <- aggregate(log_radon ~ fips, homes, mean)
names(radon_summary)[2] <- "mean_log_radon"
county <- merge(county, radon_summary, by = "fips", all.x = TRUE)

floor_summary <- aggregate(floor ~ fips, homes, mean)
names(floor_summary)[2] <- "floor_prop"
county <- merge(county, floor_summary, by = "fips", all.x = TRUE)
county$no_floor_variation <- county$floor_prop %in% c(0, 1)

map <- st_read("mn_counties_2023.gpkg", quiet = TRUE)
ix <- match(as.integer(map$GEOID), county$fips)
map$mean_log_radon <- county$mean_log_radon[ix]
map$floor_prop <- county$floor_prop[ix]
map$log_uranium <- county$log_uranium[ix]
map$no_floor_variation <- county$no_floor_variation[ix]

map_theme <- theme_void(base_size = 18, base_family = "serif") +
  theme(
    plot.title = element_text(size = 17, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    plot.background = element_rect(fill = "white", colour = NA)
  )

compact_guide <- guide_colourbar(
  direction = "vertical",
  theme = theme(
    legend.key.height = grid::unit(30, "mm"),
    legend.key.width = grid::unit(4.5, "mm")
  )
)

radon_plot <- ggplot(map) +
  geom_sf(aes(fill = mean_log_radon), colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Mean log\nradon", option = "D", na.value = "grey85"
  ) +
  guides(fill = compact_guide) +
  coord_sf(datum = NA) +
  labs(title = "Radon concentration") +
  map_theme

floor_plot <- ggplot(map) +
  geom_sf(aes(fill = floor_prop), colour = "white", linewidth = 0.2) +
  geom_sf(
    data = map[which(map$no_floor_variation), ],
    fill = NA, colour = "black", linewidth = 0.6, show.legend = FALSE
  ) +
  scale_fill_viridis_c(
    name = "Proportion\nfirst floor", option = "D", na.value = "grey85",
    breaks = seq(0, 0.6, by = 0.2)
  ) +
  guides(fill = compact_guide) +
  coord_sf(datum = NA) +
  labs(title = "Floor of measurement") +
  map_theme

uranium_plot <- ggplot(map) +
  geom_sf(aes(fill = log_uranium), colour = "white", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "Log uranium", option = "D", na.value = "grey85",
    breaks = c(-0.5, 0, 0.5)
  ) +
  guides(fill = compact_guide) +
  coord_sf(datum = NA) +
  labs(title = "County uranium") +
  map_theme

p <- wrap_plots(radon_plot, floor_plot, uranium_plot, nrow = 1) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    plot.tag.position = "bottom",
    plot.tag = element_text(size = 20, family = "serif"),
    plot.margin = margin(5.5, 5.5, 10, 5.5)
  )
print(p)

ggsave("radon_outcome_floor_uranium.pdf", p, width = 12, height = 4.6,
       device = grDevices::pdf, family = "serif", useDingbats = FALSE)
