library(sdmTMB)
library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

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

# p3 <- predict(fit, newdata = subset(qcs_grid, year == 2003), sims = 500L)
# se <- apply(p3, 1, sd)
# p$se <- se
# g_pred_se <- ggplot(p, aes(X, Y, fill = se)) +
#   geom_raster() +
#   scale_fill_viridis_c(trans = "log10", option = "C") +
#   coord_fixed() +
#   labs(fill = "Predicted\ndensity (units)", x = "UTMs E (km)", y = "UTMs N (km)") +
#   theme(legend.position = c(0.2, 0.2), legend.key.height = unit(10, "pt"))

nd <- data.frame(
  depth =
    seq(min(pcod$depth), max(pcod$depth), length.out = 300)
)
p_smooth <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA)

fit_spatiotemporal <- sdmTMB(
  density ~ s(depth),
  family = tweedie(link = "log"), data = pcod, mesh = mesh,
  time = "year", spatial = "off", spatiotemporal = "ar1", silent = FALSE
)
p_tv <- predict(fit_spatiotemporal,
  newdata = qcs_grid, return_tmb_object = TRUE,
  area = 4
)
ind <- get_index(p_tv)

# just for plotting in UTMs
pcod$X1000 <- pcod$X * 1000
pcod$Y1000 <- pcod$Y * 1000
mesh_1000 <- make_mesh(pcod, xy_cols = c("X1000", "Y1000"), cutoff = 10 * 1000)

# plot -------------------

lims_x <- c(230957.7 + 105000, 1157991 - 570000) + c(-10000, 10000)
lims_y <- c(5366427 + 270000, 6353456 - 513000) + c(-10000, 10000)
land <- "grey86"
land_border <- "grey86"

# devtools::install_github("ropensci/rnaturalearthhires") # install if needed
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

# .lab <- expression(atop("Predicted density", (km/m^2)))

g_pred <- ggplot(bc_coast_proj) +
  geom_sf(colour = land_border, lwd = 0.3, fill = land) +
  geom_raster(data = p, aes(x = X * 1000, y = Y * 1000, fill = exp(est))) +
  # coord_sf(datum = st_crs(utm_zone9)) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  ) +
  scale_x_continuous(
    limits = lims_x
    # labels = function(x) paste0(x / 1000)
  ) +
  scale_y_continuous(
    limits = lims_y
    # labels = function(x) paste0(x / 1000)
  ) +
  scale_fill_viridis_c(trans = "sqrt", option = "C", breaks = c(50, 200, 400)) +
  # scale_fill_viridis_c(trans = "sqrt", option = "B", begin = 0.3) +
  labs(fill = .lab, x = "UTMs E (km)", y = "UTMs N (km)") +
  theme(
    # panel.background = element_rect(fill = "slategray1", colour= "white"),
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
  # labs(x = "UTMs E (km)", y = "UTMs N (km)") +
  # scale_x_continuous(
  #   limits = lims_x
  #   # labels = function(x) paste0(x / 1000)
  # ) +
  # scale_y_continuous(
  #   limits = lims_y
  #   # labels = function(x) paste0(x / 1000)
  # ) +
  scale_size_area(max_size = 10) +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(legend.position = "none") +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )
# g_mesh

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
  # coord_fixed(expand = FALSE, ratio = 1.28) +
  coord_cartesian(expand = FALSE)
# theme(panel.grid.major = element_line(colour = "grey90"))
# g_smooth

g_ind <- ggplot(ind, aes(year, est / 1000)) +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  # theme(panel.grid.major = element_line(colour = "grey90")) +
  # theme(axis.title.x = element_blank()) +
  ylab("Biomass (t)") +
  xlab("Year") +
  ylim(0, 1700) +
  coord_cartesian(expand = FALSE)

g_omega <- ggplot(bc_coast_proj) +
  geom_sf(colour = land_border, lwd = 0.3, fill = land) +
  geom_raster(data = p, aes(x = X * 1000, y = Y * 1000, fill = omega_s)) +
  coord_sf() +
  # coord_sf(datum = st_crs(utm_zone9)) +
  scale_x_continuous(
    limits = lims_x
    # labels = function(x) paste0(x / 1000)
  ) +
  scale_y_continuous(
    limits = lims_y
    # labels = function(x) paste0(x / 1000)
  ) +
  # scale_fill_gradient2(high = reds[4], low = scales::muted("purple")) +
  scale_fill_gradient2(high = reds[4], low = blues[4]) +
  labs(fill = "Spatial random\nfield deviations", x = "UTMs E (km)", y = "UTMs N (km)") +
  theme(
    # panel.background = element_rect(fill = "slategray1", colour= "white"),
    legend.position = c(0.2, 0.2), legend.key.height = unit(10, "pt")
  ) +
  theme(
    panel.grid.major = element_line(colour = "grey90"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

g3 <- g_smooth / g_ind
layout <- "
      12
      34
"
g <- wrap_plots(g_mesh, g_omega, g_pred, g3) + plot_layout(design = layout)
ggsave("figs/pcod-fig.pdf", width = 9, height = 7)
