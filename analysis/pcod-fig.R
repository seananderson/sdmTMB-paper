library("sdmTMB")
library("dplyr")
library("ggplot2")
library("sf")
library("patchwork")
dir.create("figs", showWarnings = FALSE)

if (!require("rnaturalearthhires", quietly = TRUE)) {
  message(
    "Please install 'rnaturalearthhires' to use a high-resolution map.\n",
    "`remotes::install_github('ropensci/rnaturalearthhires')`"
  )
  rnaturalearthhires_installed <- FALSE
} else {
  rnaturalearthhires_installed <- TRUE
}

# from remotes::install_github('seananderson/ggsidekick')
theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size / 2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      axis.ticks.length = grid::unit(half_line / 2.2, "pt"), strip.background = element_rect(
        fill = NA,
        colour = NA
      ), strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"), axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"), legend.title = element_text(
        colour = "grey30",
        size = rel(0.9)
      ), panel.border = element_rect(
        fill = NA,
        colour = "grey70", linewidth = 1
      ), legend.key.size = grid::unit(
        0.9,
        "lines"
      ), legend.text = element_text(
        size = rel(0.7),
        colour = "grey30"
      ), legend.key = element_rect(
        colour = NA,
        fill = NA
      ), legend.background = element_rect(
        colour = NA,
        fill = NA
      ), plot.title = element_text(
        colour = "grey30",
        size = rel(1)
      ), plot.subtitle = element_text(
        colour = "grey30",
        size = rel(0.85)
      )
    )
}
theme_set(theme_sleek())

# Pacific cod example (Figure 2) ------------------------------------------

mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on", silent = FALSE
)
print(fit)
p <- predict(fit, newdata = qcs_grid)

nd <- data.frame(
  depth =
    seq(min(pcod$depth), max(pcod$depth), length.out = 300)
)
p_smooth <- predict(fit, newdata = nd, se_fit = TRUE, re_form = NA)

fit_spatiotemporal <- sdmTMB(
  density ~ 0 + s(depth),
  time_varying = ~1,
  family = tweedie(link = "log"), data = pcod, mesh = mesh,
  time = "year", spatial = "on", spatiotemporal = "iid", silent = FALSE
)
survey_grid <- replicate_df(qcs_grid, "year", unique(pcod$year))
p_tv <- predict(fit_spatiotemporal,
  newdata = survey_grid, return_tmb_object = TRUE
)
ind <- get_index(p_tv, area = 4, bias_correct = TRUE)
cog <- get_cog(p_tv, format = "wide", area = 4)

# Plot --------------------------------------------------------------------

# for plotting in UTMs
pcod$X1000 <- pcod$X * 1000
pcod$Y1000 <- pcod$Y * 1000
mesh_1000 <- make_mesh(pcod, xy_cols = c("X1000", "Y1000"), cutoff = 10 * 1000)

lims_x <- c(230957.7 + 105000, 1157991 - 570000) + c(-10000, 10000)
lims_y <- c(5366427 + 270000, 6353456 - 513000) + c(-10000, 10000)
land <- "grey86"
land_border <- "grey86"

map_data <- rnaturalearth::ne_countries(
  scale = if (rnaturalearthhires_installed) "large" else "medium",
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

# base plot
g_pred <- ggplot(bc_coast_proj) +
  geom_sf(colour = land_border, lwd = 0.3, fill = land) +
  geom_raster(data = p, aes(x = X * 1000, y = Y * 1000, fill = exp(est)))

greys <- RColorBrewer::brewer.pal(9, name = "Greys")
.n <- length(unique(cog$year))

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
    legend.position = c(0.2, 0.2), legend.key.height = grid::unit(10, "pt")
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

g_ind <- ggplot(ind, aes(year, est / 1000)) +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000), fill = "grey90") +
  geom_line(lwd = 1, colour = "grey30") +
  ylab("Biomass (t)") +
  xlab("Year") +
  ylim(0, max(ind$upr) / 1000 * 1.05) +
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
  guides(fill = "none") +
  theme(
    legend.position = c(0.2, 0.2), legend.key.height = grid::unit(10, "pt")
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
  plot_annotation(tag_levels = c("A"))

print(g)
ggsave("figs/pcod-fig.pdf", width = 9, height = 7)

sessionInfo()
