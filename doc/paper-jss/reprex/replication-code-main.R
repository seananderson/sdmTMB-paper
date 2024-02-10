## ----preliminaries, echo=FALSE, results='hide', include=FALSE, cache=FALSE----
library(knitr)
opts_chunk$set(engine = 'R', tidy = FALSE)
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


## ----main-setup, include=FALSE, cache=FALSE---------------------------
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = "#>",
 # cache = TRUE,
 # autodep = TRUE,
  fig.width = 7,
  fig.asp = 0.618,
  fig.pos = "ht",
  fig.align = "center",
  cache.comments = TRUE,
  dev = "png",
  dpi = 140,
  # optipng = "-strip all"
  highlight = FALSE
)
# knitr::knit_hooks$set(optipng = knitr::hook_optipng)
options(prompt = "R> ", continue = "+  ", width = 72, useFancyQuotes = FALSE)
opts_chunk$set(prompt = TRUE)
options(replace.assign = TRUE, width = 72, prompt = "R> ")
# knitr::render_sweave()


## ----libraries, echo=FALSE, warning=FALSE, message=FALSE--------------
library(ggplot2)
library(dplyr)
library(sdmTMB)


## ----matern-range, warning=FALSE, message=FALSE, fig.width=9, out.width="\\textwidth", fig.asp=0.33, fig.align='center', fig.cap="Example Gaussian random fields for two range values. The range describes the distance at which spatial correlation decays to $\\approx 0.13$ in coordinate units (i.e., the distance at which two points are effectively independent). Panel (a) shows a shorter range than panel (b), which results in a ``wigglier'' surface. Panel (c) shows the Mat\\'ern function for these two range values. The dashed horizontal line shows a correlation of 0.13."----
predictor_dat <- expand.grid(
  x = seq(0, 1, length.out = 100),
  y = seq(0, 1, length.out = 100),
  year = seq_len(6)
)
sim_mesh <- make_mesh(predictor_dat, xy_cols = c("x", "y"), cutoff = 0.01)
s1 <- sdmTMB_simulate(
  formula = ~1,
  data = predictor_dat,
  mesh = sim_mesh,
  range = 0.2,
  phi = 0.1,
  sigma_O = 0.2,
  seed = 1,
  B = 0
)
sim_g1 <- ggplot(s1, aes(x, y, fill = mu)) +
  geom_raster(show.legend = FALSE) +
  scale_fill_viridis_c(option = "C") +
  coord_fixed(expand = FALSE) +
  theme_light() +
  # theme(axis.title = element_blank()) +
  ggtitle("(a) Range = 0.2") +
  labs(x = "X", y = "Y")

s2 <- sdmTMB_simulate(
  formula = ~1,
  data = predictor_dat,
  mesh = sim_mesh,
  range = 0.6,
  phi = 0.1,
  sigma_O = 0.2,
  seed = 1,
  B = 0
)
sim_g2 <- ggplot(s2, aes(x, y, fill = mu)) +
  geom_raster(show.legend = FALSE) +
  scale_fill_viridis_c(option = "C") +
  coord_fixed(expand = FALSE) +
  theme_light() +
  # theme(axis.title = element_blank()) +
  ggtitle("(b) Range = 0.6") +
  labs(x = "X", y = "Y")

x <- seq(0, 1, length.out = 200)
r <- seq(0.2, 1, 0.2)
r <- c(0.2, 0.6)
df <- data.frame(
  x = rep(x, length(r)),
  range = rep(r, each = length(x))
)
matern <- function(h, sigma = 1, kappa, nu = 1) {
  ret <- (sigma^2/(2^(nu - 1) * gamma(nu))) *
    ((kappa * abs(h))^nu) *
    besselK(kappa * abs(h), nu)
  ret
  ret[x == 0] <- sigma^2
  ret
}
blues <- RColorBrewer::brewer.pal(length(r) + 1, "Blues")[-1]
df$cor <- matern(df$x, kappa = sqrt(8) / df$range)
sim_g3 <- ggplot(df, aes(x, cor, col = as.factor(range), group = as.factor(range))) +
  geom_line() +
  theme_light() +
  xlab("Distance") +
  ylab("Correlation") +
  labs(colour = "Range") +
  coord_cartesian(expand = FALSE, ylim = c(0, 1)) +
  scale_colour_manual(values = blues) +
  geom_hline(yintercept = 0.13, col = "grey50", lty = 2) +
  scale_x_continuous(breaks = r) +
  theme(legend.position = c(0.7, 0.7)) +
  ggtitle("(c) Matérn correlation function")

cowplot::plot_grid(sim_g1, sim_g2, sim_g3, align = "h", ncol = 3L)


## ----sdmTMB-install, eval=FALSE, echo=TRUE----------------------------
## install.packages("sdmTMB")


## ----sdmTMB-install-suggest, eval=FALSE, echo=TRUE--------------------
## install.packages("sdmTMB", dependencies = TRUE)


## ----remotes, eval=FALSE, echo=TRUE-----------------------------------
## install.packages("remotes")
## remotes::install_github("pbs-assess/sdmTMB", dependencies = TRUE)


## ----libs, warning=FALSE, message=FALSE, echo=TRUE, cache=FALSE-------
library(sdmTMB)
library(dplyr)
library(ggplot2)


## ----setoptions, echo=FALSE-------------------------------------------
options(
  pillar.print_max = 3,
  pillar.print_min = 3,
  pillar.advice = FALSE,
  pillar.width = 80
)
options(width = 80)


## ----libs-extras, cache=FALSE, echo=FALSE-----------------------------
theme_set(theme_light())
options(ggplot2.continuous.fill = "viridis")


## ----pcod-head, echo=TRUE---------------------------------------------
select(pcod, lat, lon, X, Y, depth, present)


## ----pcod-utms-eval, echo=TRUE, eval=FALSE----------------------------
## pcod <- add_utm_columns(pcod, c("lon", "lat"), units = "km")


## ----dog-binomial-mesh, results='hide', message=FALSE, warning=FALSE, echo=TRUE----
mesh_pcod <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 8)


## ----dog-binomial-mesh2, results='hide', message=FALSE, warning=FALSE, echo=TRUE, fig.cap="SPDE mesh (lines) combined with the trawl survey observations (points). The locations where lines intersect are referred to as ``vertices'' or ``knots''. Finer meshes will be slower to fit but generally increase the accuracy of the SPDE approximation, to a point. A greater degree of control over the mesh construction can be achieved by using \\pkg{fmesher} or \\proglang{R}-\\pkg{INLA} directly and supplying the object to \\code{make\\_mesh()}.", fig.width=4.5, fig.asp=1, out.width="3in"----
mesh_pcod2 <- make_mesh(
  pcod,
  xy_cols = c("X", "Y"),
  fmesher_func = fmesher::fm_mesh_2d_inla,
  cutoff = 8,
  max.edge = c(10, 40),
  offset = c(10, 40)
)
plot(mesh_pcod2)


## ----pcod-fit, echo=TRUE----------------------------------------------
fit_bin_rf <- sdmTMB(
  present ~ poly(log(depth), 2),
  data = pcod,
  mesh = mesh_pcod2,
  spatial = "on",
  family = binomial(link = "logit")
)


## ----pcod-fit-off, echo=TRUE------------------------------------------
fit_bin <- update(fit_bin_rf, spatial = "off")


## ----pcod-eg1-sanity, eval=TRUE, echo=TRUE, message=FALSE-------------
sanity(fit_bin_rf)


## ----pcod-bin-summary, echo=TRUE--------------------------------------
summary(fit_bin_rf)


## ----pcod-tidy, eval=TRUE, echo=TRUE, results='markup'----------------
tidy(fit_bin_rf, conf.int = TRUE)
tidy(fit_bin, conf.int = TRUE)


## ----pcod-eg1-tidy-re, eval=TRUE, echo=TRUE---------------------------
tidy(fit_bin_rf, effects = "ran_pars", conf.int = TRUE)


## ----morans, echo=FALSE, results='hide'-------------------------------
test_autocor <- function(obj) {
  set.seed(1)
  s <- simulate(obj, nsim = 500)
  pr <- predict(obj, type = "response")$est
  r <- DHARMa::createDHARMa(
    simulatedResponse = s,
    observedResponse = pcod$present,
    fittedPredictedResponse = pr
  )
  DHARMa::testSpatialAutocorrelation(r, x = pcod$X, y = pcod$Y, plot = FALSE)
}
(t_no_rf <- test_autocor(fit_bin))
(t_rf <- test_autocor(fit_bin_rf))
p_rf <- round(t_rf$p.value, 2)


## ----pcod-aic, eval=TRUE, echo=TRUE, results='markup'-----------------
AIC(fit_bin_rf, fit_bin)


## ----pcod-cv-future, eval=FALSE, echo=TRUE----------------------------
## library(future)
## plan(multisession)


## ----pcod-cv, eval=TRUE, echo=TRUE------------------------------------
set.seed(12928)
cv_bin_rf <- sdmTMB_cv(present ~ poly(log(depth), 2),
  data = pcod, mesh = mesh_pcod, spatial = "on",
  family = binomial(), k_folds = 10
)
set.seed(12928)
cv_bin <- sdmTMB_cv(present ~ poly(log(depth), 2),
  data = pcod, mesh = mesh_pcod, spatial = "off",
  family = binomial(), k_folds = 10
)


## ----pcod-cv-out, eval=TRUE, echo=TRUE, results='markup'--------------
cv_bin_rf$sum_loglik
cv_bin$sum_loglik


## ----pcod-predict, echo=TRUE, results='markup'------------------------
p <- predict(fit_bin_rf, newdata = qcs_grid)
select(p, X, Y, depth, est, est_non_rf, omega_s) |>
  as_tibble() |> head(n = 2)


## ----pcod-predict-maps, fig.width=10, fig.asp=0.4, out.width="6.1in", fig.cap="Prediction components from the binomial species distribution model of Pacific Cod. Shown are (a) the quadratic effect of bottom depth, (b) the spatial random field in link (logit) space, and (c) the overall prediction, which here is the combination of panels a and b. The spatial random field represents spatially correlated latent effects not accounted for by the fixed effects. Note the difference between predictions from depth alone (a) and predictions including a spatial random field (c)."----
plot_spatial_map <- function(dat, column, title) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    coord_fixed() +
    theme(legend.position= "bottom") +
    ggtitle(title) +
    theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
}
g1 <- plot_spatial_map(p, plogis(est_non_rf), "(a) Fixed effects")
g2 <- plot_spatial_map(p, omega_s, "(b) Spatial random field") +  scale_fill_gradient2()
g3 <- plot_spatial_map(p, plogis(est), "(c) Combined prediction")
cowplot::plot_grid(g1, g2, g3, ncol = 3)


## ----dog-head-dat, echo=TRUE------------------------------------------
dat <- select(dogfish, lon = longitude, lat = latitude, year,
  catch_weight, area_swept, depth)
dat


## ----dog-utms, echo=TRUE----------------------------------------------
dat <- add_utm_columns(dat, c("lon", "lat"),
  units = "km", utm_crs = 32609)
dat$log_depth <- log(dat$depth)
mesh <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 8)


## ----dog-mesh2, eval=FALSE--------------------------------------------
## plot(mesh)
## mesh$mesh$n


## ----dog-tw, results='hide', message=FALSE, warning=FALSE, echo=TRUE----
fit_tw <- sdmTMB(
  catch_weight ~ s(log_depth),
  data = dat,
  mesh = mesh,
  family = tweedie(),
  offset = log(dat$area_swept),
  time = "year",
  time_varying = ~ 1,
  time_varying_type = "ar1",
  spatial = "on",
  spatiotemporal = "iid",
  anisotropy = TRUE,
  silent = FALSE
)


## ----check-version, results='hide', echo=FALSE, include=FALSE---------
# delta_poisson_link_gamma() deprecated in favor of delta_gamma(type = "poisson-link")
new_deltas <- packageVersion("sdmTMB") > '0.4.2.9000'


## ----dog-update, results='hide', message=FALSE, warning=FALSE, echo=new_deltas, eval=new_deltas, include=new_deltas----
fit_dg <- update(fit_tw, family = delta_gamma(),
  spatiotemporal = list("off", "iid"))
fit_dl <- update(fit_dg, family = delta_lognormal())
fit_dpg <- update(fit_dg, family = delta_gamma(type = "poisson-link"))
fit_dpl <- update(fit_dg, family = delta_lognormal(type = "poisson-link"))


## ----dog-update-old, results='hide', message=FALSE, warning=FALSE, echo=!new_deltas, eval=!new_deltas, include=!new_deltas----
## fit_dg <- update(fit_tw, family = delta_gamma(),
##   spatiotemporal = list("off", "iid"))
## fit_dl <- update(fit_dg, family = delta_lognormal())
## fit_dpg <- update(fit_dg, family = delta_poisson_link_gamma())
## fit_dpl <- update(fit_dg, family = delta_poisson_link_lognormal())


## ----dog-aic, echo=TRUE-----------------------------------------------
AIC(fit_tw, fit_dg, fit_dl, fit_dpg, fit_dpl) |>
  mutate(delta_AIC = AIC - min(AIC)) |>
  arrange(delta_AIC)


## ----dog-ar1, results='hide', echo=TRUE-------------------------------
fit_dpl_iso <- update(fit_dpl, anisotropy = FALSE)
fit_dpl_ar1 <- update(fit_dpl, spatiotemporal = list("off", "ar1"))


## ----dog-aci2, echo=TRUE----------------------------------------------
AIC(fit_dpl_ar1, fit_dpl, fit_dpl_iso)


## ----dog-aniso, echo=TRUE, fig.cap= "A visualization of anisotropy from the function \\code{plot\\_anisotropy()}. Ellipses are centered at coordinates of zero in the units that the X-Y coordinates are modeled. The ellipses show the spatial and spatiotemporal range (distance at which correlation is effectively independent) in any direction from the center (zero).", out.width="4in"----
plot_anisotropy(fit_dpl)


## ----dog-ar1-2, echo=TRUE, results="hide"-----------------------------
fit_dpl_ar1_only <- update(fit_dpl_ar1, spatial = list("on", "off"))


## ----dog-aic3, echo=TRUE----------------------------------------------
AIC(fit_dpl_ar1_only, fit_dpl_ar1, fit_dpl)


## ----dog-print, echo=TRUE---------------------------------------------
fit <- fit_dpl_ar1_only
sanity(fit)
summary(fit)


## ----dog-grid, echo=TRUE----------------------------------------------
grid <- replicate_df(wcvi_grid, "year", time_values = unique(dat$year))
grid$log_depth <- log(grid$depth)
head(grid, n = 2)


## ----dog-pred1, echo=TRUE---------------------------------------------
pred <- predict(fit, newdata = grid, type = "response")


## ----dog-pred2, echo=TRUE---------------------------------------------
names(pred)


## ----plot-map, echo=FALSE---------------------------------------------
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(vars(year)) +
    coord_fixed() +
    theme(legend.position= "bottom") +
    theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
}


## ----dog-plot1, echo=FALSE--------------------------------------------
g_nonrf <- pred |> filter(year %in% c(2004, 2022)) |>
  plot_map(est_non_rf1) +
  scale_fill_viridis_c(trans = "log10") +
  ggtitle("(a) Non-random-field components; first delta model")


## ----dog-plot2, echo=FALSE--------------------------------------------
g_omega <- pred |> filter(year %in% c(2004, 2022)) |>
  plot_map(omega_s1) +
  scale_fill_gradient2() +
  ggtitle("(b) Spatial random field; first delta model")


## ----dog-plot3, echo=FALSE--------------------------------------------
g_eps <- pred |> filter(year %in% c(2004, 2022)) |>
  plot_map(epsilon_st2) +
  scale_fill_gradient2() +
  ggtitle("(c) Spatiotemporal random field; second delta model")


## ----dog-plot4, echo=FALSE--------------------------------------------
g_est <- pred |> filter(year %in% c(2004, 2022)) |>
  plot_map(est) +
  # labs(fill = "Density") +
  scale_fill_viridis_c(trans = "log10") +
  ggtitle("(d) Overall prediction")


## ----dog-wcvi-pred, fig.asp = 0.8, fig.cap="Example prediction elements from the spatiotemporal model of Pacific Dogfish biomass density. Throughout, two example years are shown. (a) \\code{est\\_non\\_rf1} refers to the prediction from all non-random-field elements (here, a smoother for bottom depth and the time-varying year effect) from the first delta model component, (b) \\code{omega\\_s1} refers to the spatial random field from the first delta model component, (c) \\code{epsilon\\_st2} refers to spatiotemporal random fields from the second delta model component, and (d) \\code{est} refers to the overall prediction estimate combining all effects. The spatial random field is constant through time (i.e., the two panels in b are identical) and represents static biotic or abiotic features not included as covariates (e.g., habitat). The spatiotemporal random fields are different each time step and here are constrained to follow an AR(1) process. They represent temporal variability in the spatial patterning of Pacific Spiny Dogfish (e.g., resulting from movement or local changes in population density).", fig.width=9, out.width="6in"----
cowplot::plot_grid(
  g_nonrf,
  g_omega,
  g_eps,
  g_est,
  ncol = 2L
)


## ----dog-depth, echo=TRUE, out.width="5in"----------------------------
nd <- data.frame(
  log_depth = seq(min(dat$log_depth), max(dat$log_depth), length.out = 100),
  year = max(dat$year)
)
pred_depth <- predict(
  fit, newdata = nd,
  model = NA, re_form = NA, se_fit = TRUE
)


## ----dog-depth-plot, echo=FALSE, fig.cap="The conditional effect of ocean bottom depth on Pacific Spiny Dogfish population density. The line and shaded ribbon represent the mean and 95\\% confidence interval, respectively. Other fixed effects are held at constant values and the random fields are set to their expected value (zero).", out.width="4in"----
ggplot(pred_depth, aes(
  exp(log_depth), exp(est),
  ymin = exp(est - 2 * est_se),
  ymax = exp(est + 2 * est_se))) +
  geom_ribbon(fill = "grey90") +
  geom_line() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
  labs(x = "Depth (m)", y = "Density")


## ----dog-index, echo=TRUE---------------------------------------------
pred2 <- predict(fit, newdata = grid, return_tmb_object = TRUE)
ind <- get_index(pred2, bias_correct = TRUE, area = rep(4, nrow(grid)))


## ----dog-index-plot, echo=FALSE, fig.cap="Area-weighted index of relative biomass over time for Pacific Spiny Dogfish. Dots and line segments represent means and 95\\% confidence intervals.", out.width="4in"----
ggplot(ind, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_pointrange() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
  labs(y = "Biomass", x = "Year")


## ----owl-data---------------------------------------------------------
snow <- readRDS("snow-data.rds")


## ----owl-data-head, echo=TRUE-----------------------------------------
select(snow, X, Y, year, year_f, nao, count) |> head()


## ----owl-fit, echo=TRUE, results = "hide"-----------------------------
mesh_snow <- make_mesh(snow, xy_cols = c("X", "Y"), cutoff = 1.5)
fit_owl <- sdmTMB(
  count ~ 1 + nao + (1 | year_f),
  spatial_varying = ~ nao,
  time = "year",
  data = snow,
  mesh = mesh_snow,
  family = nbinom2(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",
  reml = TRUE,
  silent = FALSE
)


## ----owl-sanity, eval=FALSE-------------------------------------------
## sanity(fit_owl)


## ----owl-print, echo=TRUE---------------------------------------------
summary(fit_owl)


## ----owl-tidy, echo=TRUE----------------------------------------------
tidy(fit_owl, conf.int = TRUE)


## ----owl-p, message=FALSE, echo=FALSE, eval=FALSE, cache=TRUE---------
## # not currently including, but might be helpful
## snow <- predict(fit_owl, newdata = snow)
## snow |> select(X, Y, year, nao, count, est, omega_s, epsilon_st, zeta_s_nao) |>
##   head(n = 2)


## ----zeta-effect, message=FALSE, echo=TRUE, eval=TRUE, cache=TRUE-----
zeta_s <- predict(fit_owl, newdata = snow, nsim = 200, sims_var = "zeta_s")
dim(zeta_s)
sims <- spread_sims(fit_owl, nsim = 200)
dim(sims)
combined <- sims$nao + t(zeta_s)
snow$nao_effect <- exp(apply(combined, 2, median))
snow$nao_effect_lwr <- exp(apply(combined, 2, quantile, probs = 0.1))
snow$nao_effect_upr <- exp(apply(combined, 2, quantile, probs = 0.9))


## ----owl-plot-basic, echo=TRUE, eval=FALSE----------------------------
## ggplot(snow, aes(X, Y)) + geom_point(aes(colour = nao_effect))


## ----shapes, echo=FALSE-----------------------------------------------
if (!file.exists("ne_10m_lakes")) {
  zip_file <- paste0("https://www.naturalearthdata.com/http//www.naturalearthdata.com/",
    "download/10m/physical/ne_10m_lakes.zip")
  download.file(zip_file, destfile = "ne_10m_lakes.zip")
  unzip("ne_10m_lakes.zip", exdir = "ne_10m_lakes")
}


## ----shapes-read, echo=FALSE------------------------------------------
Albers <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"
coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(crs = Albers)
lakes <- sf::st_read("ne_10m_lakes", quiet = TRUE)
lakes <- lakes[lakes$scalerank == 0, ] %>% sf::st_transform(crs = Albers)
land <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
land <- land %>% sf::st_transform(crs = Albers)


## ----owl-proj2--------------------------------------------------------
# project to get nice axis limits
snow2 <- snow |> mutate(X = X * 100000, Y = Y * 100000)
snow2 <- snow2 |> mutate(x = X, y = Y) |>
  sf::st_as_sf(coords = c("x", "y"), crs = Albers)


## ----owl-plot-fancy, echo=FALSE, eval=TRUE, fig.cap="Spatially varying effect of mean annual NAO (North Atlantic Oscillation) on counts of Snowy Owls observed on annual Christmas Bird Counts from 1979--2020 in Canada and the US. The effect is multiplicative on owl count per NAO unit. In the west, the lower bound of values overlaps 1 implying no effect, whereas in the southeast the effect becomes positive. Point size is scaled to the mean counts in each location.", fig.width=5, fig.asp=1, fig.pos="htbp", fig.align='center'----
nao_effect_df <- select(snow2, X, Y, count, nao_effect, nao_effect_lwr, nao_effect_upr) |>
  group_by(X, Y) |>
  summarise_all(mean)

snow_g1 <- ggplot(data = nao_effect_df) +
  geom_sf(data = land, fill = "white", colour = "white", lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray23", fill = "grey90", lwd = 0.35) +
  geom_point(aes(X, Y, colour = nao_effect, size = count), alpha = 0.5) +
  geom_sf(data = coast, colour = "gray50", fill = NA, lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray50", fill = NA, lwd = 0.35) +
  coord_sf(
    xlim = c(min(snow2$X), max(snow2$X)),
    ylim = c(min(snow2$Y), max(snow2$Y))
  ) +
  scale_colour_viridis_c(
    limit = c(min(snow$nao_effect_lwr), max(snow$nao_effect_upr)),
    guide = guide_colourbar(direction = "horizontal", title.vjust = 1, title.position = "top", label.position = "bottom")
  ) +
  guides(size = "none") +
  ggtitle(paste0("(a) Median multiplicative effect: ", round(min(snow$nao_effect), 2), " to ", round(max(snow$nao_effect), 2))) +
  labs(colour = "NAO effect\non Snowy Owl counts") +
  theme_bw() +
  theme(
    legend.position = c(0.25, 0.17),
    legend.title = element_text(size = 9, hjust = 0),
    plot.margin = unit(c(0, 0.5, 0, 0), "cm"),
    panel.background = element_rect(fill = "grey90", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.5),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(size = 10, hjust = -0.15),
    axis.title = element_blank()
  )


snow_g2 <- ggplot(data = nao_effect_df) +
  geom_sf(data = land, fill = "white", colour = "white", lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray23", fill = "grey90", lwd = 0.35) +
  geom_point(aes(X, Y, colour = nao_effect_lwr, size = count), alpha = 0.5) +
  geom_sf(data = coast, colour = "gray50", lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray50", fill = NA, lwd = 0.35) +
  coord_sf(
    xlim = c(min(snow2$X), max(snow2$X)),
    ylim = c(min(snow2$Y), max(snow2$Y))
  ) +
  scale_colour_viridis_c(
    limit = c(min(snow$nao_effect_lwr), max(snow$nao_effect_upr))
  ) +
  ggtitle(paste0("(b) Lower 80% CI: ", round(min(snow$nao_effect_lwr), 2), " to ", round(max(snow$nao_effect_lwr), 2))) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "grey90", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.5),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(size = 10), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank()
  )

snow_g3 <- ggplot(data = nao_effect_df) +
  geom_sf(data = land, fill = "white", colour = "white", lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray23", fill = "grey90", lwd = 0.35) +
  geom_point(aes(X, Y, colour = nao_effect_upr, size = count), alpha = 0.5) +
  geom_sf(data = coast, colour = "gray50", lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray50", fill = NA, lwd = 0.35) +
  coord_sf(
    xlim = c(min(snow2$X), max(snow2$X)),
    ylim = c(min(snow2$Y), max(snow2$Y))
  ) +
  scale_colour_viridis_c(
    limit = c(min(snow$nao_effect_lwr), max(snow$nao_effect_upr))
  ) +
  ggtitle(paste0("(c) Upper 80% CI: ", round(min(snow$nao_effect_upr), 2), " to ", round(max(snow$nao_effect_upr), 2))) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "grey90", colour = NA),
    panel.grid.major = element_line(colour = "white", size = 0.5),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    plot.title = element_text(size = 10), axis.ticks = element_blank(),
    axis.title = element_blank(), axis.text = element_blank()
  )

bottom_row <- cowplot::plot_grid(snow_g2, snow_g3, label_size = 12)
cowplot::plot_grid(snow_g1, bottom_row, ncol = 1, rel_heights = c(2.18, 1))


## ----pkg-versions, echo=FALSE, eval=FALSE-----------------------------
## packageVersion("INLA")
## packageVersion("inlabru")
## packageVersion("mgcv")
## packageVersion("spaMM")
## packageVersion("sdmTMB")


## ----pkgs, warning=FALSE, message=FALSE, cache=FALSE, echo=FALSE------
library("INLA")
library("inlabru")
library("ggplot2")
library("sdmTMB")
library("mgcv")
library("spaMM")


## ----mesh-timing, echo=TRUE-------------------------------------------
max_edge <- 0.06
loc_bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
segm_bnd <- INLA::inla.mesh.segment(loc_bnd)
mesh <- INLA::inla.mesh.2d(
  boundary = segm_bnd,
  max.edge = c(max_edge, 0.2),
  offset = c(0.1, 0.05)
)


## ----get-vertices, echo=FALSE-----------------------------------------
vertices <- mesh$n


## ----mesh-vis-timing, fig.width=7, fig.asp=1.1, out.width="4.3in", echo=FALSE, fig.cap="The meshes used in simulations from least to most vertices.", fig.pos='ht', fig.align='center'----
max_edge_vec <- rev(c(0.06, 0.075, 0.1, 0.15, 0.2))
g <- list()
for (i in seq_along(max_edge_vec)) {
  .mesh <- INLA::inla.mesh.2d(
    boundary = segm_bnd,
    max.edge = c(max_edge_vec[i], 0.2),
    offset = c(0.1, 0.05)
  )
  g[[i]] <- ggplot() +
    gg(.mesh) +
    theme_light() +
    coord_equal(expand = FALSE) +
    ggtitle(paste("Mesh vertices =", .mesh$n))
}
cowplot::plot_grid(plotlist = g, ncol = 2)


## ----sim-timing, echo=TRUE--------------------------------------------
set.seed(123)
n_obs <- 1000
predictor_dat <- data.frame(
  X = runif(n_obs), Y = runif(n_obs),
  a1 = rnorm(n_obs)
)
mesh_sdmTMB <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), mesh = mesh)

sim_dat <- sdmTMB_simulate(
  formula = ~ 1 + a1,
  data = predictor_dat,
  mesh = mesh_sdmTMB,
  family = gaussian(),
  range = 0.5,
  phi = 0.3,
  sigma_O = 0.2,
  B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
)


## ----sim-dat-plot-timing, echo=FALSE, fig.cap="Example simulated dataset with a spatial random field.", fig.pos='ht', fig.align='center'----
ggplot(sim_dat, aes(X, Y, colour = observed, size = abs(observed))) +
  geom_point() +
  scale_color_gradient2() +
  coord_fixed() +
  theme_light()


## ----sdmTMBfit-timing, echo=TRUE--------------------------------------
fit_sdmTMB <- sdmTMB(
  observed ~ a1,
  data = sim_dat,
  mesh = mesh_sdmTMB,
  family = gaussian(),
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 0.05, sigma_lt = 2)
  )
)


## ----spaMMfit-timing, echo=TRUE, cache=TRUE---------------------------
spde <- INLA::inla.spde2.pcmatern(
  mesh = mesh,
  prior.range = c(0.05, 0.05),
  prior.sigma = c(2, 0.05)
)

fit_spaMM <- fitme(
  observed ~ a1 + IMRF(1 | X + Y, model = spde),
  family = gaussian(),
  data = sim_dat
)


## ----inlabrufit-timing, echo=TRUE, cache=TRUE, warning=FALSE, eval=FALSE----
## dat_sp <- sp::SpatialPointsDataFrame(
##   cbind(sim_dat$X, sim_dat$Y),
##   proj4string = sp::CRS(
##     "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000
## + +y_0=0 +datum=NAD83 +units=km +no_defs"
##   ), data = sim_dat
## )
## components <- observed ~ -1 + Intercept(1) + a1 +
##   spatrf(main = coordinates, model = spde)
## like <- like(observed ~ Intercept + a1 + spatrf,
##   family = "gaussian", data = dat_sp
## )
## fit_bru <- bru(
##   like,
##   components = components,
##   options = bru_options(
##     control.inla = list(int.strategy = "eb", strategy = "gaussian"),
##     bru_max_iter = 1, num.threads = "1:1"
##   )
## )


## ----mgcv-spde-funcs-timing, echo=FALSE, eval=TRUE--------------------
# -------------------------------------------------------------------------
# Code in this chunk is from:
#
# Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic Partial
# Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020).
# https://doi.org/10.1007/s13253-019-00377-z
#
# Re-used here under a Creative Commons Attribution 4.0 International License
# http://creativecommons.org/licenses/by/4.0/
smooth.construct.spde.smooth.spec <- function(object, data, knots) {
  dim <- length(object$term)
  if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
  if (dim == 1) {
    x <- data[[object$term]]
  } else {
    x <- matrix(0, nr = length(data[[1]]), nc = 2)
    x[, 1] <- data[[object$term[1]]]
    x[, 2] <- data[[object$term[2]]]
  }
  if (is.null(object$xt)) {
    if (dim == 1) {
      t <- seq(min(x), max(x), len = object$bs.dim)
      mesh <- INLA::inla.mesh.1d(loc = t, degree = 2, boundary = "free")
    } else {
      stop("For 2D, mesh must be supplied as argument xt$mesh in s(...,xt = )")
    }
  } else {
    if (class(object$xt$mesh) != "inla.mesh") stop("xt must be NULL or an inla.mesh object")
    mesh <- object$xt$mesh
  }
  object$X <- as.matrix(INLA::inla.spde.make.A(mesh, x))
  inlamats <- INLA::inla.mesh.fem(mesh)
  object$S <- list()
  object$S[[1]] <- as.matrix(inlamats$c1)
  object$S[[2]] <- 2 * as.matrix(inlamats$g1)
  object$S[[3]] <- as.matrix(inlamats$g2)
  object$L <- matrix(c(2, 2, 2, 4, 2, 0), ncol = 2)
  object$rank <- rep(object$bs.dim, 3)
  object$null.space.dim <- 0
  object$mesh <- mesh
  object$df <- ncol(object$X)
  class(object) <- "spde.smooth"
  return(object)
}
Predict.matrix.spde.smooth <- function(object, data) {
  dim <- length(object$term)
  if (dim > 2 | dim < 1) stop("SPDE Matern can only be fit in 1D or 2D.")
  if (dim == 1) {
    x <- data[[object$term]]
  } else {
    x <- matrix(0, nr = length(data[[1]]), nc = 2)
    x[, 1] <- data[[object$term[1]]]
    x[, 2] <- data[[object$term[2]]]
  }
  Xp <- INLA::inla.spde.make.A(object$mesh, x)
  return(as.matrix(Xp))
}
# End of code from Miller et al. 2020.
# -------------------------------------------------------------------------


## ----mgcvfit-timing, echo=TRUE, eval=TRUE-----------------------------
class(mesh) <- "inla.mesh"
fit_bam <- bam(
  observed ~ a1 + s(X, Y,
    bs = "spde", k = mesh$n,
    xt = list(mesh = mesh)
  ),
  data = sim_dat,
  family = gaussian(),
  method = "fREML",
  control = gam.control(scalePenalty = FALSE),
  discrete = TRUE
)

sessionInfo()
