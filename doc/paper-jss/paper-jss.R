## ----main-setup, include=FALSE----------------------------------------------------------------------------------------------
options(prompt = 'R> ', continue = '+ ')
knitr::opts_chunk$set(
  echo = FALSE,
  message = FALSE,
  warning = FALSE,
  collapse = TRUE,
  comment = "#>",
  prompt = FALSE,
  cache = TRUE,
  autodep = TRUE,
  fig.width = 7,
  fig.asp = 0.618,
  fig.pos = "ht",
  cache.comments = TRUE,
  dev = "pdf",
  dpi = 140
  # optipng = "-strip all"
  # R.options = list(prompt = "R> ", continue = "> ")
)
# knitr::knit_hooks$set(optipng = knitr::hook_optipng)


## ----sdmTMB-lib0, eval=TRUE, echo=FALSE-------------------------------------------------------------------------------------
library(sdmTMB)
# simplify for paper:
pcod_original <- pcod # save it for later
pcod <- dplyr::select(pcod, year, density, depth, lon, lat) |>
  tibble::as_tibble()


## ----sdmTMB-lib, eval=TRUE, echo=TRUE---------------------------------------------------------------------------------------
library(sdmTMB)
head(pcod, n = 3)


## ----pcod-utms-eval, echo=FALSE, eval=TRUE----------------------------------------------------------------------------------
pcod <- add_utm_columns(pcod, c("lon", "lat"), units = "km")


## ----pcod-utms-echo, echo=TRUE, eval=FALSE----------------------------------------------------------------------------------
#> pcod <- add_utm_columns(pcod, c("lon", "lat"), units = "km")
#> #> Detected UTM zone 9N; CRS = 32609.
#> #> Visit https://epsg.io/32609 to verify.


## ----pcod-head2, echo=TRUE--------------------------------------------------------------------------------------------------
head(pcod, n = 3)


## ---- echo=TRUE, eval=TRUE, cache=TRUE--------------------------------------------------------------------------------------
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)


## ----pcod-eg1-fit, warning=FALSE, message=FALSE, echo=TRUE, eval=TRUE, cache=TRUE-------------------------------------------
fit <- sdmTMB(
  density ~ s(depth),
  data = pcod,
  family = tweedie(link = "log"),
  mesh = mesh,
  spatial = "on"
)


## ----pcod-eg1-summary,eval=TRUE, echo=TRUE----------------------------------------------------------------------------------
summary(fit)


## ----pcod-eg1-tidy-fe, eval=TRUE, echo=TRUE---------------------------------------------------------------------------------
tidy(fit, conf.int = TRUE)


## ----pcod-eg1-tidy-re, eval=TRUE, echo=TRUE---------------------------------------------------------------------------------
tidy(fit, effects = "ran_pars", conf.int = TRUE)


## ----pcod-eg1-sanity, eval=TRUE, echo=TRUE----------------------------------------------------------------------------------
sanity(fit)


## ----pcod-eg1-suppress-depth2, eval=TRUE, echo=FALSE------------------------------------------------------------------------
qcs_grid_original <- qcs_grid
qcs_grid$depth_scaled <- NULL # for simplicity in paper
qcs_grid$depth_scaled2 <- NULL # for simplicity in paper


## ----qcs-head, echo=TRUE----------------------------------------------------------------------------------------------------
head(qcs_grid, n = 3)


## ----pcod-eg1-predict, echo=TRUE, eval=TRUE, cache=TRUE---------------------------------------------------------------------
p <- predict(fit, newdata = qcs_grid)


## ----pcod-p-tibble----------------------------------------------------------------------------------------------------------
p <- tibble::as_tibble(p)


## ----pcod-head-predict, echo=TRUE-------------------------------------------------------------------------------------------
head(p, n = 3)


## ----time1, cache=TRUE------------------------------------------------------------------------------------------------------
t1 <- Sys.time()


## ----fit-tv, eval=TRUE, echo=TRUE, cache=TRUE, results='hide'---------------------------------------------------------------
fit_spatiotemporal <- sdmTMB(
  density ~ 0 + s(depth, k = 5),
  time_varying = ~ 1,
  family = tweedie(link = "log"),
  data = pcod,
  mesh = mesh,
  time = "year",
  spatial = "on",
  spatiotemporal = "iid",
  time_varying_type = "rw",
  anisotropy = TRUE,
  silent = FALSE
)


## ----time2, cache=TRUE------------------------------------------------------------------------------------------------------
t2 <- Sys.time()


## ----fit-tv-print, eval=TRUE, echo=TRUE-------------------------------------------------------------------------------------
summary(fit_spatiotemporal)


## ----pcod-restore, echo=FALSE-----------------------------------------------------------------------------------------------
pcod <- pcod_original
qcs_grid <- qcs_grid_original


## ----owl-fit, warning=FALSE, message=FALSE, cache=TRUE, echo=TRUE, eval=FALSE-----------------------------------------------
#> mesh <- make_mesh(snow, xy_cols = c("X", "Y"), cutoff = 1.5)
#> fit_owls <- sdmTMB(
#>   count ~ nao + (1 | year_factor),
#>   spatial_varying = ~ nao,
#>   family = nbinom2(link = "log"),
#>   data = snow,
#>   mesh = mesh,
#>   time = "year",
#>   spatial = "on",
#>   spatiotemporal = "iid"
#> )


## ----owl-nao, fig.cap="Spatially varying coefficient for effect of mean annual NAO (North Atlantic Oscillation) on counts of Snowy Owls observed on annual Christmas Bird Counts 1979--2020 in Canada and the US. Points represent all count locations and circle area is scaled to the mean number of owls observed per year (range: 0 to 8). The effect is multiplicative on owl count per NAO unit.", out.width="4.1in", fig.align='center'----
owl <- here::here("figs", "owl-nao-effect.png")
knitr::include_graphics(owl)


## ----compare-table, results='asis', eval=FALSE------------------------------------------------------------------------------
#> library(kableExtra)
#> suppressMessages(suppressWarnings(d <- readr::read_csv(here::here("doc/software-comparison.csv"), show_col_types = FALSE, comment = "#")))
#> d <- as.matrix(d)
#> # d[d == "1"] <- "\\checkmark"
#> d <- gsub("^1", "\\\\checkmark", d)
#> # d[d == "1*"] <- "\\checkmark\\footnote{a}"
#> # d[d == "GAMs*"] <- "GAMs\\footnote{a}"
#> # d[d == "0"] <- "--"
#> d <- gsub("^0", "--", d)
#> # d[d == "0*"] <- "--\\footnote{a}"
#> d <- as.data.frame(d)
#> names(d)[1] <- ""
#> row.names(d) <- NULL
#> d <- d[,c(1, 2, 3, 5, 6, 4, 7, 8, 9)]
#> d$Hmsc <- NULL
#>
#> d$INLA <- NULL
#> names(d)[names(d) == "inlabru"] <- "INLA/inlabru"
#>
#> knitr::kable(d, format = "latex", caption = "Comparison of functionality between several R packages that can fit geostatistical GLMMs.",
#>   booktabs = TRUE,
#>   linesep = c(
#>     rep('', 4), # coefs
#>     '\\addlinespace',
#>     rep('', 6), # correlation stuff
#>     '\\addlinespace',
#>     rep('', 12), # likelihood
#>     '\\addlinespace',
#>     rep('', 2), # bayes
#>     '\\addlinespace'
#>   ),
#>   escape = FALSE, row.names = FALSE) %>%
#>   kableExtra::kable_styling(font_size = 8) %>%
#>   kableExtra::footnote(
#>     general = paste(
#>       "$^1$Technically possible but non-trivial.",
#>       "$^2$Penalized smoother GAMs that determine `wiggliness'.",
#>       "$^3$inlabru but not INLA.",
#>       "$^4$Spatiotemporal fields as random walk only.",
#>       "$^5$SPDE approach as in Miller et al. (2019).",
#>       "$^6$Does have predictive process knots.",
#>       "$^7$Zero-inflated NB2 only.",
#>       "$^8$Tweedie power parameter fixed for \\\\texttt{mgcv::gamm()}.",
#>       "$^9$Possible as log-linked Poisson GLMM with aggregated data.",
#>       "$^{10}$Hurdle models possible by fitting components separately."
#>     ),
#>     escape = FALSE, threeparttable = TRUE
#>   )




## ----packages, message=FALSE, warning=FALSE, cache=FALSE--------------------------------------------------------------------
library("ggplot2")
library("dplyr")
library("sdmTMB")
theme_set(theme_minimal())


## ----plot-mesh, fig.align='center', fig.asp=1, out.width="3in", echo=TRUE, fig.cap="Delaunay triangulation mesh.", fig.pos='ht', fig.align='center'----
mesh <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
plot(mesh)


## ----fit-pcod, warning=FALSE, message=FALSE, cache=TRUE, echo=TRUE, results='hide'------------------------------------------
fit <- sdmTMB(
  density ~ 0 + as.factor(year),
  data = pcod,
  time = "year",
  mesh = mesh,
  family = tweedie(link = "log"),
  silent = FALSE
)


## ----print-fit, echo=TRUE---------------------------------------------------------------------------------------------------
fit


## ----pcod-tidy, echo=TRUE---------------------------------------------------------------------------------------------------
tidy(fit, conf.int = TRUE)
tidy(fit, effects = "ran_pars", conf.int = TRUE)


## ----pcod-predict, cache=TRUE, echo=TRUE------------------------------------------------------------------------------------
pred <- predict(fit)


## ----sim-resids, out.width="4in", warning=FALSE, results='hide', cache=TRUE, echo=TRUE, fig.cap="Simulated residuals from the DHARMa package.", fig.pos='ht', fig.align='center'----
set.seed(123)
simulate(fit, nsim = 500) %>%
  dharma_residuals(fit)


## ----mcmc-resids, warning=FALSE, results='hide', cache=TRUE-----------------------------------------------------------------
set.seed(123)
r <- residuals(fit,
  type = "mle-mcmc",
  mcmc_iter = 200, mcmc_warmup = 100
)


## ----mcmc-resids-vis, out.width="4in", fig.cap="MCMC residuals.", fig.pos='ht', fig.align='center'--------------------------
qqnorm(r)
qqline(r)


## ----predict-newdata--------------------------------------------------------------------------------------------------------
survey_grid <- replicate_df(
  qcs_grid,
  time_name = "year",
  time_values = unique(pcod$year)
)

pred_qcs <- predict(fit, newdata = survey_grid)


## ----plot-map---------------------------------------------------------------------------------------------------------------
plot_map <- function(dat, column) {
  ggplot(dat, aes_string("X", "Y", fill = column)) +
    geom_raster() +
    facet_wrap(~year) +
    coord_fixed()
}


## ----plot-all-effects, fig.cap="Model predictions based on all fixed effects and random effects.", fig.pos='ht', fig.align='center'----
plot_map(pred_qcs, "exp(est)") +
  # ggtitle("Prediction (fixed effects + all random effects)") +
  scale_fill_viridis_c(trans = "log10")


## ----plot-fix-defects, fig.cap="Predictions based on only the non-random field elements.", fig.pos='ht', fig.align='center'----
plot_map(pred_qcs, "exp(est_non_rf)") +
  # ggtitle("Prediction (fixed effects only)") +
  scale_fill_viridis_c(trans = "sqrt")


## ----plot-spatial-effects, fig.cap="Spatial random effects only.", fig.pos='ht', fig.align='center'-------------------------
plot_map(pred_qcs, "omega_s") +
  # ggtitle("Spatial random effects only") +
  scale_fill_gradient2()


## ----plot-spatiotemporal-effects, fig.cap="Spatiotemporal random effects only.", fig.pos='ht', fig.align='center'-----------
plot_map(pred_qcs, "epsilon_st") +
  # ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()


## ----pcod-sims, cache=TRUE--------------------------------------------------------------------------------------------------
pred_sims <- predict(fit, newdata = survey_grid, nsim = 200)
dim(pred_sims)


## ----sims-cv-plot, cache=TRUE, fig.cap="Spatiotemporal coefficient of variation.", fig.pos='ht', fig.align='center'---------
survey_grid$cv <- apply(pred_sims, 1, function(x) sd(exp(x)) / mean(exp(x)))
ggplot(survey_grid, aes(X, Y, fill = cv)) +
  geom_raster() +
  facet_wrap(~year) +
  coord_fixed() +
  scale_fill_viridis_c(trans = "log10", option = "D")


## ----pcod-app-index, fig.asp=0.45, out.width="5in", cache=TRUE, fig.cap="Geostatistical population index.", fig.pos='ht', fig.align='center'----
survey_grid$area <- 4 # all 2 x 2km
pred2 <- predict(
  fit,
  newdata = survey_grid,
  return_tmb_object = TRUE
)
ind <- get_index(pred2, area = survey_grid$area, bias_correct = TRUE)

ggplot(ind, aes(year, est / 1000)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr / 1000, ymax = upr / 1000), alpha = 0.4) +
  ylab("Biomass (t)")




## ----setup-owls, include = FALSE, cache=FALSE-------------------------------------------------------------------------------
theme_set(theme_bw() + theme(
  plot.title = element_text(size = 12),
  legend.position = "top", legend.justification = c(0, 0), # move legend to top left
  legend.title = element_text(size = 10, hjust = 0), legend.key.height = unit(0.3, "cm"),
  strip.background = element_rect(fill = NA, colour = NA), # de-emphasize facet labels
  panel.grid.major = element_line(colour = "grey90"), # for lat/lons
  axis.title = element_blank() # lat/lon is clear from axis text
))


## ----packages-owls, message=FALSE, warning=FALSE, cache=FALSE---------------------------------------------------------------
library("ggplot2")
library("patchwork")
library("dplyr")
library("sf")
library("sdmTMB")


## ----rnaturalearthdata, echo=FALSE------------------------------------------------------------------------------------------
if (!require("rnaturalearthdata", quietly = TRUE)) {
  stop(
    "Please install 'rnaturalearthdata'.\n",
    "`remotes::install_github('ropensci/rnaturalearthdata')`"
  )
}


## ----load-data, cache=FALSE-------------------------------------------------------------------------------------------------
snow <- readRDS(here::here("data/SNOW_data.rds"))


## ----head-snow, echo=TRUE---------------------------------------------------------------------------------------------------
head(snow, n = 3)


## ----proj, echo=TRUE--------------------------------------------------------------------------------------------------------
Albers <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0
+datum=NAD83 +units=m +no_defs"


## ----mesh-owls, echo=TRUE, out.width="4in", fig.cap="This mesh uses an Albers projection divided by 100000 to give units of 100 km and a cutoff distance of 1.5 units (or 150 km).", fig.pos='ht', fig.align='center', dev='png', dpi=100----
mesh <- make_mesh(snow, xy_cols = c("X", "Y"), cutoff = 1.5)

# built-in plot() method:
# plot(mesh)

# or with ggplot:
ggplot() + inlabru::gg(mesh$mesh) +
  geom_point(aes(X, Y), data = snow) +
  theme_light()


## ----fit-negb2, echo=TRUE, warning=FALSE------------------------------------------------------------------------------------
fit_owl <- sdmTMB(count ~ 1 + nao + (1 | year_f),
  spatial_varying = ~ nao,
  time = "year",
  data = snow,
  family = nbinom2(link = "log"),
  spatial = "on",
  spatiotemporal = "iid",
  mesh = mesh,
  reml = TRUE,
  silent = FALSE
)


## ----sanity, echo=TRUE------------------------------------------------------------------------------------------------------
sanity(fit_owl)


## ---------------------------------------------------------------------------------------------------------------------------
fit_owl$gradients


## ----owl-newton, results='hide'---------------------------------------------------------------------------------------------
fit_owl2 <- run_extra_optimization(fit_owl, newton_loops = 1)


## ----sanity2, echo=TRUE-----------------------------------------------------------------------------------------------------
sanity(fit_owl2)


## ----print-fit-owls, warning=FALSE, echo=TRUE, message=FALSE----------------------------------------------------------------
print(fit_owl2)


## ----tidy-owls--------------------------------------------------------------------------------------------------------------
tidy(fit_owl2, conf.int = TRUE)
tidy(fit_owl2, effects = "ran_pars", conf.int = TRUE)


## ----predict-owls, message=FALSE--------------------------------------------------------------------------------------------
pred <- predict(fit_owl2)
set.seed(19208)
pred$resid <- residuals(fit_owl2)


## ----qqnorm-owls, out.width="60%", fig.cap="Randomized quantile residuals.", fig.pos='ht', fig.align='center', dpi=90, dev='png'----
qqnorm(pred$resid)
qqline(pred$resid)


## ----stan-resids-owls, eval=FALSE-------------------------------------------------------------------------------------------
#> fit_ml <- update(fit_owl, reml = FALSE)
#> mcmc_resids <- residuals(fit_owl,
#>   type = "mle-mcmc",
#>   mcmc_iter = 101,
#>   mcmc_warmup = 100
#> )


## ----zero-test-owls, echo=TRUE, message=FALSE-------------------------------------------------------------------------------
s_nb2 <- simulate(fit_owl2, nsim = 400)


## ----zero-test-owls2, echo=TRUE---------------------------------------------------------------------------------------------
mean(s_nb2 == 0)
mean(snow$count == 0)


## ----zero-test-owls3, echo=TRUE---------------------------------------------------------------------------------------------
pred_fixed <- fit_owl2$family$linkinv(pred$est_non_rf)
r_nb2 <- DHARMa::createDHARMa(
  simulatedResponse = s_nb2,
  observedResponse = fit_owl2$data$count,
  fittedPredictedResponse = pred_fixed
)
DHARMa::testZeroInflation(r_nb2, plot = FALSE)


## ----proj-p-----------------------------------------------------------------------------------------------------------------
p <- pred %>% mutate(X = X * 100000, Y = Y * 100000)
p_proj <- p %>% mutate(x = X, y = Y) %>%
  sf::st_as_sf(coords = c("x", "y"), crs = Albers)


## ----shapes, echo=FALSE-----------------------------------------------------------------------------------------------------
if (!file.exists(here::here("data/ne_10m_lakes"))) {
  zip_file <- paste0("https://www.naturalearthdata.com/http//www.naturalearthdata.com/",
    "download/10m/physical/ne_10m_lakes.zip")
  download.file(zip_file, destfile = here::here("data/ne_10m_lakes.zip"))
  unzip(here::here("data/ne_10m_lakes.zip"), exdir = here::here("data/ne_10m_lakes"))
}


## ----shapes-read, echo=TRUE-------------------------------------------------------------------------------------------------
coast <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf") %>%
  sf::st_transform(crs = Albers)
lakes <- sf::st_read(here::here("data/ne_10m_lakes"), quiet = TRUE)
lakes <- lakes[lakes$scalerank == 0, ] %>% sf::st_transform(crs = Albers)


## ----p-mean, message=FALSE, echo=TRUE, eval=TRUE, cache=TRUE----------------------------------------------------------------
b <- tidy(fit_owl2, conf.int = TRUE)

nsim <- 200
zeta_s <- predict(fit_owl2, nsim = nsim, sims_var = "zeta_s")
sims <- spread_sims(fit_owl2, nsim = nsim)
beta <- sims$nao
combined <- beta + t(zeta_s)

p$zeta_s <- as.numeric(apply(t(zeta_s), 2, median))
p$nao_effect_sim <- as.numeric(apply(combined, 2, median))
p$nao_effect_lwr <- as.numeric(apply(combined, 2, quantile, probs = 0.10))
p$nao_effect_upr <- as.numeric(apply(combined, 2, quantile, probs = 0.90))

p_mean <- p %>%
  group_by(CBCID) %>% summarise(
    year = max(year),
    X = mean(X), Y = mean(Y),
    mean_est_count = exp(mean(est)),
    nao_effect = b$estimate[b$term == "nao"] + mean(zeta_s),
    nao_sim = mean(nao_effect_sim),
    nao_lwr = mean(nao_effect_lwr),
    nao_upr = mean(nao_effect_upr)
  )


## ---- echo=FALSE------------------------------------------------------------------------------------------------------------
nao_effect_lwr95 <- as.numeric(apply(combined, 2, quantile, probs = 0.025))


## ----plot-map-owls, echo=FALSE, eval=TRUE, out.width="65%", dpi=140, dev="png"----------------------------------------------
ggplot(data = p_proj) +
  geom_point(
    data = p_mean, aes(X, Y, colour = nao_effect, size = mean_est_count), alpha = 0.5
  ) +
  geom_sf(data = coast, colour = "gray50") +
  geom_sf(data = lakes, colour = "gray50", fill = NA) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000), # adjusts space on sides
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_viridis_c(
    limit = c(min(p_mean$nao_lwr), max(p_mean$nao_upr)),
    guide = guide_colourbar(direction = "horizontal", title.position = "top")
  ) +
  guides(size = "none") +
  labs(colour = "Estimated effect in log space") +
  # ggtitle("Map of combined main and spatially varying effects of NAO on Snowy Owl count") +
  theme(legend.position = c(0.1, 0.07), axis.title = element_blank())


## ----plot-map-ci, echo=FALSE, eval=TRUE, out.width="80%",fig.asp=0.3, fig.cap="Map of combined main and spatially varying effects of NAO on Snowy Owl count with confidence intervals. In the west, the lower values span zero implying a no evidence of an effect at this confidence level, while in the southeast the lowest estimates are still positive. Note that the values are in log space and thus additive. Points represent all count locations and circle area is scaled to the mean number of owls observed per year (range: 0 to 8).", fig.pos='ht', fig.align='center', dpi=140, dev="png"----
p1 <- ggplot(data = p_proj) +
  geom_point(
    data = p_mean,
    aes(X, Y, colour = nao_lwr, size = mean_est_count), alpha = 0.5
  ) +
  geom_sf(data = coast, colour = "gray50") +
  geom_sf(data = lakes, colour = "gray50", fill = NA) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000), # adjusts space on sides
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_viridis_c(
    limit = c(min(p_mean$nao_lwr), max(p_mean$nao_upr)),
    guide = guide_colourbar(direction = "horizontal", title.position = "top")
  ) +
  guides(size = "none", colour = "none") +
  # labs(colour = "Lower 90% CI") +
  ggtitle("Lower 80% CI") +
  theme(legend.position = c(0.1, 0.1), axis.title = element_blank())

p2 <- ggplot(data = p_proj) +
  geom_point(
    data = p_mean,
    aes(X, Y, colour = nao_upr, size = mean_est_count), alpha = 0.5
  ) +
  geom_sf(data = coast, colour = "gray50") +
  geom_sf(data = lakes, colour = "gray50", fill = NA) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000), # adjusts space on sides
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_viridis_c(
    limit = c(min(p_mean$nao_lwr), max(p_mean$nao_upr)),
    guide = guide_colourbar(direction = "horizontal", title.position = "top")
  ) +
  guides(size = "none", colour = "none") +
  # labs(colour = "Upper 90% CI") +
  ggtitle("Upper 80% CI") +
  theme(legend.position = c(0.1, 0.1), axis.title = element_blank())

p1 + p2 + patchwork::plot_layout()


## ----plot-all-effects-owls, out.width="90%", fig.asp=0.8, echo=FALSE, eval=TRUE, dpi=140, fig.cap="Predicted winter owl counts that incorporate all fixed effects and random effects.", fig.pos='ht', fig.align='center', dev="png"----
ggplot(data = p_proj) +
  geom_point(data = p, aes(X, Y, colour = exp(est)), size = 0.4, alpha = 0.5) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000),
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_viridis_c(
    trans = "sqrt", # makes variation in low numbers visible
    # can also trim the extreme high values
    limits = c(0, quantile(exp(p$est), 0.995)), na.value = "yellow",
    guide = guide_colourbar(title.position = "top")
  ) +
  labs(colour = "Predicted counts") +
  theme(axis.text = element_blank(),
  legend.title = element_text(size = 8, hjust = 0), legend.key.height = unit(0.25, "cm"),
  axis.ticks = element_blank()) +
  facet_wrap(~year_f)


## ----plot-spatial-effects-owls, echo=FALSE, eval=TRUE, out.width="3in", dpi=130, fig.cap="Spatial random effects that apply across all years.", fig.pos='ht', fig.align='center', dev = "png"----
ggplot(data = p_proj) +
  geom_point(data = p, aes(X, Y, colour = omega_s), alpha = 0.5) +
  geom_sf(data = coast, colour = "gray50") +
  geom_sf(data = lakes, colour = "gray50", fill = NA) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000),
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_gradient2(guide = guide_colourbar(
    direction = "horizontal", title.position = "top"
  )) +
  theme(legend.position = c(0.1, 0.1), axis.title = element_blank())


## ----plot-spatiotemporal-effects-owls, out.width="75%", fig.asp=0.8, echo=FALSE, eval=TRUE, dpi=120, fig.cap="Spatiotemporal random effects.", fig.pos='ht', fig.align='center', dev="png"----
ggplot(data = filter(p_proj, year > 1978)) +
  geom_point(
    data = filter(p, year > 1978), aes(X, Y, colour = epsilon_st),
    size = 0.4, alpha = 0.5
  ) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000),
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_gradient2(guide = guide_colourbar(title.position = "top")) +
  theme(axis.text = element_blank(),
  legend.title = element_text(size = 8, hjust = 0), legend.key.height = unit(0.25, "cm"),
  axis.ticks = element_blank()) +
  facet_wrap(~year_f)


## ----plot-fixed-effects-owls, out.width="75%", fig.asp=0.8, echo=FALSE, eval=TRUE, dpi=120, fig.cap="Predictions for each year based on only the non-random field elements.", fig.pos='ht', fig.align='center', dev="png"----
ggplot(data = filter(p_proj, year > 1978)) +
  geom_point(
    data = filter(p, year > 1978), aes(X, Y, colour = exp(est_non_rf)),
    size = 0.4, alpha = 0.5
  ) +
  coord_sf(
    xlim = c(min(p_proj$X) - 50000, max(p_proj$X) - 50000),
    ylim = c(min(p_proj$Y), max(p_proj$Y))
  ) +
  scale_colour_viridis_c(
    breaks = c(0.03, 0.05, 0.07),
    guide = guide_colourbar(title.position = "top")
  ) +
  theme(axis.text = element_blank(),
  legend.title = element_text(size = 8, hjust = 0), legend.key.height = unit(0.25, "cm"),
  axis.ticks = element_blank()) +
  facet_wrap(~year_f)

## ----packages-inla-comparison, message=FALSE, warning=FALSE, cache=FALSE, echo=TRUE-----------------------------------------
library("ggplot2")
library("dplyr")
library("sdmTMB")
library("INLA")
library("inlabru")
theme_set(theme_light())


## ----inla-threads, echo=TRUE------------------------------------------------------------------------------------------------
INLA::inla.setOption(num.threads = "1:1")


## ----pcod-inla-comparison-data, echo=TRUE-----------------------------------------------------------------------------------
d <- pcod %>% filter(year >= 2007)
yr_lu <- data.frame(year = unique(d$year), i_year = seq_along(unique(d$year)))
d <- dplyr::left_join(d, yr_lu, by = "year")


## ----set-priors, echo = TRUE------------------------------------------------------------------------------------------------
range_min <- 5
sigma_max <- 5
prior_prob <- 0.05


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------
h_spec <- list(rho = list(prior = "pccor1", param = c(0, 0.9)))


## ----sdm-mesh, echo = TRUE--------------------------------------------------------------------------------------------------
sdmTMB_mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 10)


## ----inla-mesh, echo = TRUE-------------------------------------------------------------------------------------------------
loc_xy <- as.matrix(d[, c("X", "Y"), drop = FALSE])
inla_mesh <- INLA::inla.mesh.create(loc_xy, refine = TRUE, cutoff = 10)
spde <- inla.spde2.pcmatern(
  mesh = inla_mesh,
  prior.range = c(range_min, prior_prob),
  prior.sigma = c(sigma_max, prior_prob)
)


## ----plot-inla-mesh, fig.asp=0.9, out.height="2.5in",echo=FALSE, eval=FALSE-------------------------------------------------
#> plot(sdmTMB_mesh)
#> title(main = "sdmTMB mesh")
#> plot(spde$mesh, main = "", asp = 1)
#> title(main = "INLA mesh")


## ----inla, echo=TRUE--------------------------------------------------------------------------------------------------------
tictoc::tic("INLA") # start timing
nyear <- length(unique(d$year))
iset <- inla.spde.make.index("i", n.spde = spde$n.spde, n.group = nyear)
A <- inla.spde.make.A(mesh = inla_mesh, loc = loc_xy, group = d$i_year)
sdat <- inla.stack(
  data = list(y = d$present),
  A = list(A, 1),
  effects = list(iset, year_factor = as.factor(d$year)),
  tag = "stdata"
)

formula <- y ~ 0 + year_factor + f(i,
  model = spde, group = i.group,
  control.group = list(model = "ar1", hyper = h_spec)
)

m_inla <- inla(formula,
  data = inla.stack.data(sdat),
  control.predictor = list(compute = TRUE, A = inla.stack.A(sdat)),
  family = "binomial",
  control.family = list(link = "logit"),
  control.fixed = list(
    expand.factor.strategy = "inla",
    mean = 0, prec = 1 / (100 * 100)
  ),
  control.inla = list(int.strategy = "eb", strategy = "gaussian")
)
tictoc::toc() # stop timing


## ---- echo=FALSE, eval=FALSE------------------------------------------------------------------------------------------------
#> summary(m_inla)


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------
dat <- as.data.frame(d)
sp::coordinates(dat) <- c("X", "Y")
dat <- as(dat, "SpatialPointsDataFrame")


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------
tictoc::tic("bru")
components <- present ~ 0 +
  f(
    main = coordinates, model = spde, group = i_year,
    control.group = list(model = "ar1", hyper = h_spec)
  ) +
  fac(main = year, model = "factor_full")

m_bru <- bru(
  components = components,
  family = "binomial",
  data = dat,
  options =
    list(
      bru_verbose = TRUE,
      control.fixed = list(
        expand.factor.strategy = "inla",
        mean = 0, prec = 1 / (100 * 100)
      ),
      control.family = list(link = "logit"),
      control.inla = list(int.strategy = "eb", strategy = "gaussian")
    )
)
tictoc::toc()


## ---- echo=TRUE, eval=FALSE-------------------------------------------------------------------------------------------------
#> summary(m_bru)


## ---- echo=TRUE-------------------------------------------------------------------------------------------------------------
m_bru$summary.random$fac[, c("ID", "mean", "sd", "0.025quant", "0.975quant")]


## ----sdmtmb, warning=FALSE, message=FALSE, echo=TRUE------------------------------------------------------------------------
tictoc::tic("sdmTMB")
m <- sdmTMB(present ~ 0 + as.factor(year),
  data = d,
  mesh = sdmTMB_mesh,
  time = "year",
  spatial = "off",
  spatiotemporal = "ar1",
  family = binomial(link = "logit"),
  priors = sdmTMBpriors(
    matern_st = pc_matern(
      range_gt = range_min, range_prob = prior_prob,
      sigma_lt = sigma_max, sigma_prob = prior_prob
    ),
    b = normal(0, 100),
    ar1_rho = normal(0, 2)
  ),
  reml = TRUE # to match INLA "eb" `int.strategy`
)
tictoc::toc()


## ----sdmtmb-m, echo=TRUE----------------------------------------------------------------------------------------------------
print(m)


## ----sdmtmb-coefs, echo=TRUE------------------------------------------------------------------------------------------------
sdmTMB_coefs <- tidy(m, effects = "fixed", conf.int = TRUE)
sdmTMB_coefs


## ----comparison, echo=TRUE--------------------------------------------------------------------------------------------------
# INLA
m_inla$summary.hyperpar[, c("mean", "sd", "0.025quant", "0.975quant")]

# inlabru
m_bru$summary.hyperpar[, c("mean", "sd", "0.025quant", "0.975quant")]

# sdmTMB
tidy(m, "ran_pars", conf.int = TRUE)


## ----coef-plot, fig.asp=0.4, out.width="85%", echo=FALSE, fig.cap="Fixed effect estimates from all three models with 95 percent CI.", fig.pos='ht', fig.align='center'----
inla_coefs <- m_inla$summary.fixed[, c("0.025quant", "mean", "0.975quant")] %>%
  rename(estimate = mean, conf.low = `0.025quant`, conf.high = `0.975quant`) %>%
  mutate(term = unique(sdmTMB_coefs$term), model = "INLA")

bru_coefs <- m_bru$summary.random$fac[, c("0.025quant", "mean", "0.975quant")] %>%
  rename(estimate = mean, conf.low = `0.025quant`, conf.high = `0.975quant`) %>%
  mutate(term = unique(sdmTMB_coefs$term), model = "inlabru")

dplyr::bind_rows(mutate(sdmTMB_coefs, model = "sdmTMB"), inla_coefs, bru_coefs) %>%
  ggplot(aes(x = estimate, y = term, xmin = conf.low, xmax = conf.high, colour = model)) +
  geom_pointrange(position = position_dodge(width = 0.2)) +
  theme(axis.title.y = element_blank())


## ----sdmTMB-sims, echo=FALSE------------------------------------------------------------------------------------------------
set.seed(1)
post_sdmTMB <- gather_sims(m, nsim = 2000)
post_sdmTMB <- post_sdmTMB %>%
  mutate(
    variable =
      factor(`.variable`,
        levels = unique(post_sdmTMB$`.variable`),
        labels = c(seq(2007, 2017, 2), "Range", "AR1 rho", "Random Field SD")
      )
  )


## ----inla-densities, echo=FALSE---------------------------------------------------------------------------------------------
old_inla <- any(grepl(
  "year_factoryear_factor",
  names(m_inla$marginals.fixed)
))
yrs <- seq(2007, 2017, 2)
if (old_inla) {
  yf <- "year_factoryear_factor"
} else {
  yf <- "year_factor"
}
x <- list()
for (i in seq_along(yrs)) {
  x[[i]] <- inla.smarginal(m_inla$marginals.fixed[[paste0(yf, yrs[i])]]) %>%
    as_tibble() %>%
    mutate(variable = paste0("as.factor.year.", yrs[i]))
}
x[[7]] <- inla.smarginal(m_inla$marginals.hyperpar$`Range for i`) %>%
  as_tibble() %>%
  mutate(variable = "range")
x[[8]] <- inla.smarginal(m_inla$marginals.hyperpar$`Stdev for i`) %>%
  as_tibble() %>%
  mutate(variable = "sigma_E")
x[[9]] <- inla.smarginal(m_inla$marginals.hyperpar$`GroupRho for i`) %>%
  as_tibble() %>%
  mutate(variable = "ar1_rho")
post_inla1 <- dplyr::bind_rows(x)
post_inla2 <- post_inla1 %>% mutate(variable = factor(variable,
  levels = unique(post_inla1$variable),
  labels = c(yrs, "Range", "Random Field SD", "AR1 rho")
))


## ----compare-plot, out.width="6in", fig.cap="Blue histograms represent simulated draws from the joint precision matrix of the sdmTMB model compared with the marginal posterior densities from the INLA model.", fig.align='center'----
ggplot(post_inla2) +
  geom_histogram(
    data = post_sdmTMB, aes(.value, after_stat(density)), bins = 30,
    fill = "blue", alpha = 0.2
  ) +
  geom_line(data = post_inla2, aes(x, y), inherit.aes = FALSE) +
  facet_wrap(~variable, scales = "free")


## ----sdm-predictions--------------------------------------------------------------------------------------------------------
nd <- replicate_df(qcs_grid, time_name = "year", time_values = unique(d$year))
pred <- predict(m, newdata = nd)


## ----maps, out.width="5.5in", fig.cap="Predicted occurence probabilities in time and space from the sdmTMB model.", fig.align='center'----
inverse_link <- m$family$linkinv
ggplot(pred) +
  geom_raster(aes(X, Y, fill = inverse_link(est))) +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  coord_fixed()


## ----inlabru-pred, out.width="5.5in", fig.cap="Predicted occurence probabilities in time and space from the inlabru model.", fig.align='center'----
nd_bru <- nd
nd_bru <- dplyr::left_join(nd_bru, yr_lu, by = "year")
sp::coordinates(nd_bru) <- c("X", "Y")
nd_bru <- as(nd_bru, "SpatialPixelsDataFrame")
pred_bru <- predict(m_bru, data = nd_bru, formula = ~ plogis(f + fac))

ggplot(as.data.frame(pred_bru)) +
  geom_raster(aes(X, Y, fill = mean)) +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  coord_fixed()


## ----inla-pred, out.width="5.5in"-------------------------------------------------------------------------------------------
# idx <- inla.stack.index(sdat, "stdata")$data
# d$inla_pred <- m_inla$summary.fitted.values[idx, "mean"]
#
# ggplot(d) +
#   geom_point(aes(X, Y, colour = inla_pred), size = 2) +
#   facet_wrap(~year) +
#   scale_colour_viridis_c() +
#   coord_fixed()




## ----pkg-versions, echo=TRUE------------------------------------------------------------------------------------------------
packageVersion("INLA")
packageVersion("inlabru")
packageVersion("mgcv")
packageVersion("spaMM")
packageVersion("sdmTMB")


## ----pkgs, warning=FALSE, message=FALSE, cache=FALSE, echo=FALSE------------------------------------------------------------
library("INLA")
library("inlabru")
library("ggplot2")
library("sdmTMB")
library("mgcv")
library("spaMM")


## ----mesh-timing, echo=TRUE-------------------------------------------------------------------------------------------------
max_edge <- 0.06
loc_bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
segm_bnd <- INLA::inla.mesh.segment(loc_bnd)
mesh <- INLA::inla.mesh.2d(
  boundary = segm_bnd,
  max.edge = c(max_edge, 0.2),
  offset = c(0.1, 0.05)
)


## ----mesh-vis-timing, fig.width=7, fig.asp=1.1, echo=FALSE, fig.cap="The meshes used in simulations from least to most vertices.", fig.pos='ht', fig.align='center'----
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


## ----sim-timing, echo=TRUE--------------------------------------------------------------------------------------------------
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


## ----sim-dat-plot-timing, out.width="5in", echo=FALSE, fig.cap="Example simulated dataset with a spatial random field.", fig.pos='ht', fig.align='center'----
ggplot(sim_dat, aes(X, Y, colour = observed, size = abs(observed))) +
  geom_point() +
  scale_color_gradient2() +
  coord_fixed() +
  theme_light()


## ----sdmTMBfit-timing, echo=TRUE--------------------------------------------------------------------------------------------
fit_sdmTMB <- sdmTMB(
  observed ~ a1,
  data = sim_dat,
  mesh = mesh_sdmTMB,
  family = gaussian(),
  priors = sdmTMBpriors(
    matern_s = pc_matern(range_gt = 0.05, sigma_lt = 2)
  )
)


## ----spaMMfit-timing, echo=TRUE, cache=TRUE---------------------------------------------------------------------------------
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


## ----inlabrufit-timing, echo=TRUE, cache=TRUE-------------------------------------------------------------------------------
# convert to sp SpatialPointsDataFrame first:
dat_sp <- sp::SpatialPointsDataFrame(
  cbind(sim_dat$X, sim_dat$Y),
  proj4string = sp::CRS(
    "+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000
+ +y_0=0 +datum=NAD83 +units=km +no_defs"
  ), data = sim_dat
)
components <- observed ~ -1 + Intercept(1) + a1 +
  spatrf(main = coordinates, model = spde)
like <- like(observed ~ Intercept + a1 + spatrf,
  family = "gaussian", data = dat_sp
)
fit_bru <- bru(
  like,
  components = components,
  options = bru_options(
    control.inla = list(int.strategy = "eb", strategy = "gaussian"),
    bru_max_iter = 1, num.threads = "1:1"
  )
)


## ----mgcv-spde-funcs-timing, echo=FALSE, eval=TRUE--------------------------------------------------------------------------
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
# End of code from
# Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic Partial
# Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020).
# https://doi.org/10.1007/s13253-019-00377-z
#
# Re-used here under a Creative Commons Attribution 4.0 International License
# http://creativecommons.org/licenses/by/4.0/
# -------------------------------------------------------------------------


## ----mgcvfit-timing, echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------
# define smooth.construct.spde.smooth.spec() and Predict.matrix.spde.smooth()
# from supplement of:
# Miller, D.L., Glennie, R. & Seaton, A.E. (2019). Understanding the Stochastic
# Partial Differential Equation approach to smoothing. Journal of Agricultural,
# Biological and Environmental Statistics.
# https://doi.org/10.1007/s13253-019-00377-z

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


## ----Pacific cod figures from main text, echo=TRUE-------------------------------------------------------------------------------------------------

library("sdmTMB")
library("dplyr")
library("ggplot2")
library("sf")
library("patchwork")
dir.create("figs", showWarnings = FALSE)

if (!require("rnaturalearthhires", quietly = TRUE)) {
  stop(
    "Please install 'rnaturalearthhires'.\n",
    "`remotes::install_github('ropensci/rnaturalearthhires')`"
  )
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
        colour = "grey70", size = 1
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
    # ggsave("figs/pcod-fig.pdf", width = 9, height = 7)

## ----Timing figure in main text, echo=TRUE-------------------------------------------------------------------------------------------------

library("sdmTMB")
library("dplyr")
library("ggplot2")
library("INLA")
library("inlabru")
library("mgcv")
library("spaMM")
dir.create("figs", showWarnings = FALSE)
dir.create("analysis", showWarnings = FALSE)

INLA::inla.setOption(num.threads = "1:1")

# source(here::here("analysis/mgcv_spde_smooth.R"))
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
# End of code from
# Miller, D.L., Glennie, R. & Seaton, A.E. Understanding the Stochastic Partial
# Differential Equation Approach to Smoothing. JABES 25, 1–16 (2020).
# https://doi.org/10.1007/s13253-019-00377-z
#
# Re-used here under a Creative Commons Attribution 4.0 International License
# http://creativecommons.org/licenses/by/4.0/
# -------------------------------------------------------------------------

simulate_dat <- function(n_obs = 100,
  cutoff = 0.1,
  range = 0.5,
  phi = 0.1,
  tweedie_p = 1.5,
  max.edge = 0.1,
  family = gaussian(),
  sigma_O = 0.2,
  iter = 1,
  plot = FALSE,
  seed = sample.int(.Machine$integer.max, 1L)) {
  set.seed(seed)
  predictor_dat <- data.frame(
    X = runif(n_obs), Y = runif(n_obs),
    a1 = rnorm(n_obs)
  )

  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- inla.mesh.segment(loc.bnd)

  me <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )

  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), mesh = me)
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    mesh = mesh,
    family = family,
    range = range,
    phi = phi,
    tweedie_p = tweedie_p,
    sigma_O = sigma_O,
    seed = seed,
    B = c(0.2, -0.4) # B0 = intercept, B1 = a1 slope
  )
  if (plot) {
    g <- ggplot(sim_dat, aes(X, Y, colour = observed)) +
      geom_point() +
      scale_color_gradient2()
    print(g)
  }
  list(mesh = mesh, dat = sim_dat)
}

sim_fit_time <- function(n_obs = 1000, cutoff = 0.1, iter = 1, phi = 0.3, tweedie_p = 1.5,
  seed = sample.int(.Machine$integer.max, 1), family = gaussian(), sigma_O = 0.2,
  max.edge = 0.1) {
  cat("max.edge:", max.edge, "\n")
  cat("n_obs:", n_obs, "\n")
  cat("iter:", iter, "\n")
  cat("simulating...\n")

  s <- simulate_dat(
    n_obs = n_obs, cutoff = cutoff, seed = seed,
    iter = iter, family = family, max.edge = max.edge, tweedie_p = tweedie_p, phi = phi,
    sigma_O = sigma_O
  )
  sim_dat <- s$dat
  mesh <- s$mesh
  times <- list()
  cat("fitting...\n")

  cat("sdmTMB\n")
  out <- system.time({
    fit <- sdmTMB(observed ~ a1,
      data = sim_dat, mesh = mesh, family = family,
      priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 0.05, sigma_lt = 2))
    )
  })
  times$sdmTMB <- out[["elapsed"]]

  dat_sp <- sp::SpatialPointsDataFrame(cbind(sim_dat$X, sim_dat$Y),
    proj4string = sp::CRS("+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000
+ +y_0=0 +datum=NAD83 +units=km +no_defs"), data = sim_dat
  )
  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd)
  mesh2 <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )
  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh2,
    prior.range = c(0.05, 0.05),
    prior.sigma = c(2, 0.05)
  )

  dat <- s$dat
  dat_sp <- sp::SpatialPointsDataFrame(cbind(dat$X, dat$Y),
    proj4string = sp::CRS("+proj=aea +lat_0=45 +lon_0=-126 +lat_1=50 +lat_2=58.5 +x_0=1000000
  + +y_0=0 +datum=NAD83 +units=km +no_defs"), data = dat
  )
  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd)
  mesh <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )

  spde <- INLA::inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(0.05, 0.05),
    prior.sigma = c(2, 0.05)
  )

  components <- observed ~ -1 + Intercept(1) + a1 +
    spatrf(main = coordinates, model = spde)

  dat_sp$observed <- dat_sp$observed * 10 # help convergence
  if (family$family == "gaussian") {
    like <- like(observed ~ Intercept + a1 + spatrf, family = "gaussian", data = dat_sp)
  } else if (family$family == "tweedie") {
    like <- like(observed ~ Intercept + a1 + spatrf, family = "tweedie", data = dat_sp)
  } else {
    stop("Family not found")
  }

  out <- system.time({
    tryCatch(
      {
        fit_bru <- bru(
          like,
          components = components,
          options = bru_options(
            control.inla = list(int.strategy = "eb", strategy = "gaussian"),
            bru_max_iter = 1, num.threads = "1:1"
          )
        )
      },
      error = function(e) NA
    )
  })
  if (!is.null(fit_bru$error)) {
    inla_eb_error <- TRUE
  } else {
    inla_eb_error <- FALSE
  }
  times$inlabru_eb <- out[["elapsed"]]

  cat("spaMM\n")

  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- INLA::inla.mesh.segment(loc.bnd)
  mesh2 <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )
  # note `<<-` below:
  # Avoids: `Hmmm... it looks like you should put some object(s) in control.HLfit$formula_env`
  spde_spaMM <<- INLA::inla.spde2.pcmatern(
    mesh = mesh2,
    prior.range = c(0.05, 0.05),
    prior.sigma = c(2, 0.05)
  )
  out <- system.time({
    fit_spaMM <- tryCatch(
      {
        fitme(observed ~ a1 + IMRF(1 | X + Y, model = spde_spaMM),
          family = family, data = sim_dat
        )
      },
      error = function(e) NA
    )
  })
  times$spaMM <- if (all(is.na(fit_spaMM))) NA else out[["elapsed"]]

  cat("mgcv disc.\n")
  loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
  segm.bnd <- inla.mesh.segment(loc.bnd)
  mesh2 <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )
  out <- system.time({
    fit <- tryCatch(
      {
        bam(
          observed ~ a1 + s(X, Y,
            bs = "spde", k = mesh2$n,
            xt = list(mesh = mesh2)
          ),
          data = sim_dat,
          family = family,
          method = "fREML",
          control = gam.control(scalePenalty = FALSE), discrete = TRUE
        )
      },
      error = function(e) NA
    )
  })
  times$mgcv_disc <- if (all(is.na(fit))) NA else out[["elapsed"]]

  out <- as_tibble(times)
  out$cutoff <- cutoff
  out$n_obs <- as.integer(n_obs)
  out$iter <- as.integer(iter)
  out$seed <- seed
  out$mesh_n <- as.integer(mesh$n)
  out$max.edge <- max.edge
  out$inla_eb_error <- inla_eb_error
  out
}

# fast run:
to_run <- expand.grid(
  n_obs = c(1000L),
  max.edge = c(0.06, 0.1, 0.2),
  iter = seq_len(2L)
)

# full (slow) run to match paper:
if (FALSE) {
  to_run <- expand.grid(
    n_obs = c(1000L),
    max.edge = c(0.06, 0.075, 0.1, 0.15, 0.2),
    iter = seq_len(50L)
  )
}

to_run$seed <- to_run$iter * 29212
nrow(to_run)

# don't run in parallel to ensure all timing is comparable:
system.time({
  out <- purrr::pmap_dfr(
    to_run,
    sim_fit_time
  )
})

# saveRDS(out, file = "analysis/timing-cache-parallel-2022-11-01.rds")
# out <- readRDS("analysis/timing-cache-parallel-2022-11-01.rds")

out_long <- tidyr::pivot_longer(
  out,
  sdmTMB:mgcv_disc,
  names_to = "model",
  values_to = "time"
)
clean_names <- tribble(
  ~model, ~model_clean,
  "sdmTMB", "sdmTMB",
  "inlabru_eb", "inlabru EB",
  "spaMM", "spaMM",
  "mgcv_disc", "mgcv::bam\n(discretize = TRUE) SPDE"
)
clean_names$model_clean <- factor(clean_names$model_clean, levels = clean_names$model_clean)
out_long <- left_join(out_long, clean_names)
out_long <- out_long %>%
  group_by(model_clean, cutoff) %>%
  mutate(mean_mesh_n = paste0("Mesh n = ", round(mean(mesh_n))))
out_long$mean_mesh_n <- as.factor(out_long$mean_mesh_n)
out_long$mean_mesh_n <- forcats::fct_reorder(out_long$mean_mesh_n, -out_long$cutoff)
out_long_sum <- group_by(out_long, model_clean, n_obs, cutoff, max.edge) %>%
  summarise(lwr = quantile(time, probs = 0.025), upr = quantile(time, probs = 0.975), time = mean(time), mean_mesh_n = mean(mesh_n))

cols <- RColorBrewer::brewer.pal(6, "Dark2")
colour_vect <- c(
  "sdmTMB" = cols[1],
  "inlabru EB" = cols[3],
  "spaMM" = cols[4],
  "mgcv::bam\n(discretize = TRUE) SPDE" = cols[2]
)

leg <- tribble(
  ~model_clean, ~mean_mesh_n, ~time,
  "sdmTMB", 173, 0.26,
  "inlabru EB", 173, 3.5,
  "spaMM", 500, 2.0,
  "mgcv::bam\n(discretize = TRUE) SPDE", 249, 20.5
)

g <- out_long_sum %>%
  filter(model_clean != "mgcv::bam\n(discretize = F) SPDE") %>%
  ggplot(aes(mean_mesh_n, time, colour = model_clean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model_clean),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(lwd = 0.7) +
  scale_y_log10() +
  scale_x_log10(breaks = c(250, 500, 1000)) +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(legend.position = c(0.4, 0.85)) +
  geom_text(aes(label = model_clean), data = leg, hjust = 0) +
  scale_fill_manual(values = colour_vect) +
  scale_colour_manual(values = colour_vect) +
  labs(y = "Time (s)", x = "Mesh nodes", colour = "Model", fill = "Model") +
  coord_cartesian(expand = FALSE, ylim = c(0.09, 45)) +
  theme(panel.spacing.x = unit(20, "pt"), legend.position = "none")
print(g)

.width <- 5
# ggsave("figs/timing-spatial.pdf", width = .width, height = .width / 1.3)

cat("inlabru to sdmTMB timing ratios\n")
x <- filter(
  out_long_sum,
  model_clean == "inlabru EB"
)$time /
  filter(out_long_sum, model_clean == "sdmTMB")$time
round(x, 1)

sessionInfo()
