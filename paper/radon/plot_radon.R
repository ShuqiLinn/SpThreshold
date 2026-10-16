library(ggplot2)
library(patchwork)
library(sf)

homes <- read.csv("srrs2.dat", strip.white = TRUE)
county <- read.csv("cty.dat", strip.white = TRUE)
names(homes) <- tolower(trimws(names(homes)))
names(county) <- tolower(trimws(names(county)))

homes <- homes[trimws(homes$state) == "MN" & homes$activity > 0, ]
homes$fips <- 1000L * homes$stfips + homes$cntyfips

county <- county[county$stfips == 27L, ]
county$fips <- 1000L * county$stfips + county$ctfips
county$log_uranium <- log(county$uppm)
county <- unique(county[c("fips", "log_uranium")])

floor_summary <- aggregate(floor ~ fips, homes, mean)
names(floor_summary)[2] <- "floor_prop"
county <- merge(county, floor_summary, by = "fips", all.x = TRUE)
county$no_floor_variation <- county$floor_prop %in% c(0, 1)

download.file(
  "https://www2.census.gov/geo/tiger/GENZ2023/shp/cb_2023_us_county_500k.zip",
  "cb_2023_us_county_500k.zip", mode = "wb"
)

boundaries <- st_read("/vsizip/cb_2023_us_county_500k.zip", quiet = TRUE)
map <- boundaries[boundaries$STATEFP == "27", ]
ix <- match(as.integer(map$GEOID), county$fips)
map$floor_prop <- county$floor_prop[ix]
map$log_uranium <- county$log_uranium[ix]
map$no_floor_variation <- county$no_floor_variation[ix]

map_theme <- theme_void(base_size = 16, base_family = "serif") +
  theme(
    plot.title = element_text(size = 13, hjust = 0, margin = margin(b = 8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    plot.background = element_rect(fill = "white", colour = NA)
  )

compact_guide <- guide_colourbar(
  direction = "vertical",
  theme = theme(
    legend.key.height = grid::unit(28, "mm"),
    legend.key.width = grid::unit(4, "mm")
  )
)

floor_plot <- ggplot(map) +
  geom_sf(aes(fill = floor_prop), colour = "white", linewidth = 0.2) +
  geom_sf(
    data = map[map$no_floor_variation, ],
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

p <- wrap_plots(floor_plot, uranium_plot, nrow = 1)
print(p)

ggsave("radon_covariates.pdf", p, width = 9, height = 5,
       device = grDevices::pdf, family = "serif", useDingbats = FALSE)