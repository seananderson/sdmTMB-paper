library(sdmTMB)
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

if (!require("rnaturalearthhires", quietly = TRUE)) {
  stop("Please install 'rnaturalearthhires'.\n",
  "`remotes::install_github('ropensci/rnaturalearthhires')`")
}
if (!require("ggsidekick", quietly = TRUE)) {
  stop("Please install 'ggsidekick'.\n",
  "`remotes::install_github('seananderson/ggsidekick')`")
}

theme_set(ggsidekick::theme_sleek())

mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)

fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on", silent = FALSE
)
print(fit)
p <- predict(fit, newdata = subset(qcs_grid, year == 2003))

nd <- data.frame(
  depth =
    seq(min(pcod$depth), max(pcod$depth), length.out = 300)
)
p_smooth <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA)

fit_spatiotemporal <- sdmTMB(
  density ~ s(depth),
  family = tweedie(link = "log"), data = pcod, mesh = mesh,
  time = "year", spatial = "on", spatiotemporal = "ar1", silent = FALSE
)
p_tv <- predict(fit_spatiotemporal,
  newdata = qcs_grid, return_tmb_object = TRUE,
  area = 4
)
ind <- get_index(p_tv)
cog <- get_cog(p_tv, format = "wide")

# just for plotting in UTMs
pcod$X1000 <- pcod$X * 1000
pcod$Y1000 <- pcod$Y * 1000
mesh_1000 <- make_mesh(pcod, xy_cols = c("X1000", "Y1000"), cutoff = 10 * 1000)

# plot -------------------

lims_x <- c(230957.7 + 105000, 1157991 - 570000) + c(-10000, 10000)
lims_y <- c(5366427 + 270000, 6353456 - 513000) + c(-10000, 10000)
land <- "grey86"
land_border <- "grey86"

map_data <- rnaturalearth::ne_countries(
  scale = "large",
  returnclass = "sf", country = "canada"
)
# Crop the polygon for plotting and efficiency:
# st_bbox(map_data) # find the rough coordinates
bc_coast <- suppressWarnings(suppressMessages(
  st_crop(
    map_data,
    c(xmin = -134, ymin = 46, xmax = -120, ymax = 57)
  )
))
utm_zone9 <- 3156
bc_coast_proj <- sf::st_transform(bc_coast, crs = utm_zone9)

.lab <- expression("Density" ~ (kg / km^2))

# base plot:
g_pred <- ggplot(bc_coast_proj) +
  geom_sf(colour = land_border, lwd = 0.3, fill = land) +
  geom_raster(data = p, aes(x = X * 1000, y = Y * 1000, fill = exp(est)))

greys <- RColorBrewer::brewer.pal(9, name = "Greys")
.n <- length(unique(cog$year))
cross_cols <- rev(colorRampPalette(c("#FFFFFF", greys[length(greys)]))(.n+3))
cross_cols <- cross_cols[-c(1:3)]
names(cross_cols) <- unique(cog$year)

# build in year order to crosses overlap in right order
for (i in unique(cog$year)) {
  .d <- filter(cog, year == i)
  .d$year <- as.character(.d$year)
  g_pred <- g_pred +
    # # white outline:
    # geom_segment(aes(y = lwr_y * 1000, yend = upr_y * 1000, x = est_x * 1000, xend = est_x * 1000), .d, lwd = 1.1, colour = "grey90") +
    # geom_segment(aes(y = est_y * 1000, yend = est_y * 1000, x = lwr_x * 1000, xend = upr_x * 1000), .d, lwd = 1.1, colour = "grey90") +
    # coloured lines:
    geom_segment(aes(y = lwr_y * 1000, yend = upr_y * 1000, x = est_x * 1000, xend = est_x * 1000, colour = year), .d, lwd = 0.8) +
    geom_segment(aes(y = est_y * 1000, yend = est_y * 1000, x = lwr_x * 1000, xend = upr_x * 1000, colour = year), .d, lwd = 0.8) +
    # scale_colour_viridis_c(option = "D")
    scale_colour_manual(values = cross_cols)
}

# theme it etc.
g_pred <- g_pred + theme(
  panel.grid.major = element_line(colour = "grey90"),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
) +
  scale_x_continuous(
    limits = lims_x
  ) +
  scale_y_continuous(
    limits = lims_y
  ) +
  scale_fill_viridis_c(trans = "sqrt", option = "C", breaks = c(50, 200, 400)) +
  labs(fill = .lab, x = "UTMs E (km)", y = "UTMs N (km)", colour = "Year") +
  guides(colour = "none") +
  theme(
    legend.position = c(0.2, 0.2), legend.key.height = unit(10, "pt")
  )
g_pred

blues <- RColorBrewer::brewer.pal(5, "Blues")
reds <- RColorBrewer::brewer.pal(5, "Reds")

g_mesh <- ggplot(bc_coast_proj) +
  geom_sf(colour = land_border, lwd = 0.3, fill = land) +
  inlabru::gg(mesh_1000$mesh, edge.color = "grey61") +
  coord_sf(xlim = lims_x, ylim = lims_y) +
  geom_point(data = pcod, alpha = 0.8, mapping = aes(X * 1000, Y * 1000, size = density), colour = blues[4]) +
  scale_size_area(max_size = 10) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(legend.position = "none") +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
g_mesh

g_smooth <- ggplot(
  p_smooth,
  aes(depth, exp(est),
    ymin = exp(est - 1.96 * est_se),
    ymax = exp(est + 1.96 * est_se)
  )
) +
  geom_ribbon(fill = blues[2]) +
  geom_line(colour = blues[5], lwd = 1) +
  labs(x = "Depth (m)", y = expression("Density" ~ (kg / km^2))) +
  xlim(min(pcod$depth), 300) +
  ylim(0, 179) +
  coord_cartesian(expand = FALSE)
g_smooth

# ggplot(bc_coast_proj) +
#   geom_sf(colour = land_border, lwd = 0.3, fill = land) +
#   coord_sf(xlim = lims_x, ylim = lims_y) +
#   geom_point(aes(est_x, est_y, colour = year), cog) +
#   geom_segment(aes(y = lwr_y * 1000, yend = upr_y * 1000, x = est_x * 1000, xend = est_x * 1000, colour = year), cog) +
#   geom_segment(aes(y = est_y * 1000, yend = est_y * 1000, x = lwr_x * 1000, xend = upr_x * 1000, colour = year), cog) +
#   scale_colour_viridis_c(option = "D") +
#   theme(panel.grid.major = element_line(colour = "grey90")) +
#   theme(legend.position = "none") +
#   theme(
#     panel.grid.major = element_line(colour = "grey90"),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank()
#   )


g_ind <- ggplot(ind, aes(year, est / 1000)) +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  ylab("Biomass (t)") +
  xlab("Year") +
  ylim(0, 1700) +
  coord_cartesian(expand = FALSE)
g_ind

g_omega <- ggplot(bc_coast_proj) +
  geom_sf(colour = land_border, lwd = 0.3, fill = land) +
  geom_raster(data = p, aes(x = X * 1000, y = Y * 1000, fill = omega_s)) +
  coord_sf() +
  scale_x_continuous(
    limits = lims_x
  ) +
  scale_y_continuous(
    limits = lims_y
  ) +
  scale_fill_gradient2(high = reds[4], low = blues[4]) +
  labs(fill = "Spatial random\nfield deviations", x = "UTMs E (km)", y = "UTMs N (km)") +
  theme(
    legend.position = c(0.2, 0.2), legend.key.height = unit(10, "pt")
  ) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
g_omega

g3 <- g_smooth / g_ind
layout <- "
      12
      34
"
g <- wrap_plots(g_mesh, g_omega, g_pred, g3) + plot_layout(design = layout) +
  plot_annotation(tag_levels = c('A'))
# g
ggsave("figs/pcod-fig.pdf", width = 9, height = 7)
