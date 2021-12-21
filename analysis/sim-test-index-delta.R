library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
theme_set(ggsidekick::theme_sleek())
plan(multisession, workers = floor(availableCores() / 2) - 0)
# options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
N <- 625 # per year
time_steps <- 20
sigma_O <- 0
sigma_E <- 0.6
phi <- 9
tweedie_p <- 1.5 # Tweedie p
.range <- 70

# x <- seq(-1, 1, length.out = 25)
# y <- seq(-1, 1, length.out = 25)
x <- seq(340, 580, length.out = 25)
y <- seq(5500, 5800, length.out = 25)
loc <- expand.grid(x = x, y = y)
x <- rep(loc$x, time_steps)
y <- rep(loc$y, time_steps)
loc <- data.frame(x = x, y = y)
spde <- sdmTMB::make_mesh(loc, xy_cols = c("x", "y"), cutoff = 10)
plot(spde)
plot(spde$mesh, asp = 1)

sim_test_index <- function(i) {
  set.seed(i * 1029)
  betas <- rnorm(time_steps, 1, 1)
  dat <- data.frame(year = rep(seq_len(time_steps), each = N))
  s <- sdmTMB_simulate(
    formula = ~ 0 + as.factor(year),
    data = dat,
    mesh = spde,
    B = betas,
    phi = phi,
    range = .range,
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    seed = SEED * i,
    time = "year",
    family = tweedie(),
    tweedie_p = tweedie_p
  )
  s <- mutate(s, year = dat$year)
  s <- mutate(s, x = loc$x)
  s <- mutate(s, y = loc$y)
  s_sampled <- s %>% group_by(year) %>% slice_sample(n = 200L)

  if (FALSE) {
    ggplot(s, aes(x, y, fill = mu)) +
      geom_raster() +
      geom_point(aes(size = observed), data = s_sampled, pch = 21, fill = NA) +
      scale_fill_viridis_c() +
      scale_size_area() +
      facet_wrap(vars(year)) +
      coord_cartesian(expand = FALSE)
  }

  FieldConfig <- matrix(c(0, "IID", "IID", 0, "IID", "IID"),
    ncol = 2, nrow = 3,
    dimnames = list(
      c("Omega", "Epsilon", "Beta"),
      c("Component_1", "Component_2")
    )
  )
  RhoConfig <- c("Beta1" = 0, "Beta2" = 0, "Epsilon1" = 0, "Epsilon2" = 0)
  # remove replicate locations for each year and format for VAST
  grid_ll <- unique(as.data.frame(loc))
  grid_ll <- rename(grid_ll, X = x, Y = y)
  grid_ll$Y <- grid_ll$Y * 1000
  grid_ll$X <- grid_ll$X * 1000
  sp::coordinates(grid_ll) <- ~ X + Y
  sp::proj4string(grid_ll) <- sp::CRS("+proj=utm +zone=9")
  grid_ll <- as.data.frame(sp::spTransform(grid_ll, sp::CRS("+proj=longlat +datum=WGS84")))
  input_grid <- cbind(Lat = grid_ll$Y, Lon = grid_ll$X, Area_km2 = 1)

  s_ll <- s_sampled
  s_ll$Y <- s_ll$y * 1000
  s_ll$X <- s_ll$x * 1000
  sp::coordinates(s_ll) <- ~ X + Y
  sp::proj4string(s_ll) <- sp::CRS("+proj=utm +zone=9")
  s_ll <- as.data.frame(sp::spTransform(s_ll, sp::CRS("+proj=longlat +datum=WGS84")))
  s_ll <- rename(s_ll, Lon = "X", Lat = "Y")

  settings <- FishStatsUtils::make_settings(
    n_x = 125, # number of vertices in the SPDE mesh
    Region = "User",
    purpose = "index2", # index of abundance with Gamma for positive catches
    fine_scale = TRUE, # use bilinear interpolation from the INLA 'A' matrix
    zone = 9,
    FieldConfig = FieldConfig,
    RhoConfig = RhoConfig,
    ObsModel = c(2, 0), # conventional logit-linked delta-Gamma; c(10, 2) for Tweedie
    bias.correct = TRUE,
    use_anisotropy = FALSE,
    max_cells = Inf, # use all grid cells from the extrapolation grid
    knot_method = "grid" # or "samples"
  )
  dir.create(paste0(here::here("temp-"), i), showWarnings = FALSE)
  # effort is 1 when using CPUE instead of observed weight as the response:
  s_ll$effort <- 1
  s_ll <- as.data.frame(s_ll)
  library(VAST)
  index_vast <- tryCatch({
    m_vast <- FishStatsUtils::fit_model(
      settings = settings,
      Lat_i = s_ll[, "Lat"],
      Lon_i = s_ll[, "Lon"],
      t_i = s_ll[, "year"],
      b_i = s_ll[, "observed"],
      a_i = s_ll[, "effort"],
      input_grid = input_grid,
      working_dir = paste0(here::here("temp-"), i, "/"))
    suppressWarnings(out <- plot_biomass_index(m_vast, DirName = paste0(here::here("temp"), "/")))
    index_vast <- data.frame(
      year = unique(s_ll$year),
      est_vast = out$Table$Estimate,
      log_est_vast = log(out$Table$Estimate),
      se_vast = out$log_Index_ctl[1,,1,2],
      lwr_vast = log(out$Table$Estimate) + qnorm(0.025) * out$log_Index_ctl[1,,1,2],
      upr_vast = log(out$Table$Estimate) + qnorm(0.975) * out$log_Index_ctl[1,,1,2]
    )
    index_vast
  }, error = function(e) return(NULL))
  if (is.null(index_vast)) {
    index_vast <- structure(list(year = NA_real_, est_vast = NA_real_,
      lwr_vast = NA_real_,
      upr_vast = NA_real_, log_est_vast = NA_real_, se_vast = NA_real_),
      row.names = 1L, class = "data.frame")
  }

  if (!is.null(index_vast) && exists('m_vast')) {
    mesh <- sdmTMB::make_mesh(s_sampled, xy_cols = c("x", "y"), mesh = m_vast$spatial_list$MeshList$isotropic_mesh)
  } else {
    mesh <- sdmTMB::make_mesh(s_sampled, xy_cols = c("x", "y"), n_knots = 125)
  }
  # mesh <- make_mesh(s_sampled, xy_cols = c("x", "y"), cutoff = 0.2)
  # mesh <- sdmTMB::make_mesh(s_sampled, xy_cols = c("x", "y"), n_knots = 125)
  # plot(mesh)
  s_sampled$present <- ifelse(s_sampled$observed > 0, 1, 0)
  m1 <- tryCatch({sdmTMB(
    data = s_sampled,
    formula = present ~ 0 + as.factor(year),
    time = "year",
    mesh = mesh,
    spatial = "off",
    # silent = FALSE,
    family = binomial(),
  )}, error = function(e) return(NULL))

  if (max(m1$gradients) > 0.001) {
    m1 <- tryCatch({run_extra_optimization(m1,
      nlminb_loops = 0, newton_loops = 1)},
      error = function(e) return(NULL))
  }

  s_pos <- subset(s_sampled, present == 1)
  mesh_pos <- sdmTMB::make_mesh(s_pos, xy_cols = c("x", "y"), mesh = mesh$mesh)
  # plot(mesh_pos)

  m2 <- tryCatch({sdmTMB(
    data = s_pos,
    formula = observed ~ 0 + as.factor(year),
    time = "year",
    mesh = mesh_pos,
    spatial = "off",
    # silent = FALSE,
    family = Gamma(link = "log"),
  )}, error = function(e) return(NULL))

  if (max(m2$gradients) > 0.001) {
    m2 <- tryCatch({run_extra_optimization(m2,
      nlminb_loops = 0, newton_loops = 1)},
      error = function(e) return(NULL))
  }

  m3 <- tryCatch({sdmTMB(
    data = s_sampled,
    formula = observed ~ 0 + as.factor(year),
    time = "year",
    mesh = mesh,
    spatial = "off",
    # silent = FALSE,
    family = tweedie(),
  )}, error = function(e) return(NULL))

  if (max(m1$gradients) > 0.001) {
    m3 <- tryCatch({run_extra_optimization(m3,
      nlminb_loops = 0, newton_loops = 1)},
      error = function(e) return(NULL))
  }

  grid <- select(s, x, y, year)

  p_sims1 <- tryCatch({
    predict(m1, newdata = grid, sims = 800L)
  }, error = function(e) return(NULL))
  p_sims2 <- tryCatch({
    predict(m2, newdata = grid, sims = 800L)
  }, error = function(e) return(NULL))
  p_sims3 <- tryCatch({
    predict(m3, newdata = grid, sims = 800L)
  }, error = function(e) return(NULL))

  if (!is.null(p_sims1) && !is.null(p_sims2)) {
    p_bin_prob <- plogis(p_sims1)
    p_pos_exp <- exp(p_sims2)
    p_combined <- log(p_bin_prob * p_pos_exp)
    index_sims <- get_index_sims(p_combined)
    names(index_sims) <- paste0(names(index_sims), "_sim")
  } else {
    index_sims <- structure(list(year_sim = NA_real_, est_sim = NA_real_,
      lwr_sim = NA_real_,
      upr_sim = NA_real_, log_est_sim = NA_real_, se_sim = NA_real_),
      row.names = 1L, class = "data.frame")
  }

  if (!is.null(p_sims3)) {
    index_sims_tw <- get_index_sims(p_sims3)
    names(index_sims_tw) <- paste0(names(index_sims_tw), "_sim_tw")
  } else {
    index_sims_tw <- structure(list(est_sim_tw = NA_real_,
      lwr_sim_tw = NA_real_,
      upr_sim_tw = NA_real_, log_est_sim_tw = NA_real_, se_sim_tw = NA_real_),
      row.names = 1L, class = "data.frame")
  }

  true_index <- group_by(s, year) %>%
    summarise(true = sum(mu), .groups = "drop") %>%
    mutate(year = as.integer(year))

  if (FALSE) {
    ggplot(index_sims, aes(year_sim, est_sim, ymin = lwr_sim, ymax = upr_sim)) +
      geom_ribbon(fill = "grey70") +
      geom_line() +
      geom_line(aes(year, true), true_index, colour = "red",
        inherit.aes = FALSE) +
      geom_line(aes(year, est_vast), index_vast, colour = "blue",
        inherit.aes = FALSE)
  }

  out <- bind_cols(select(true_index, true), index_sims)
  out <- bind_cols(out, index_vast)
  out <- bind_cols(out, index_sims_tw)
  out$pdHess1 <- m1$sd_report$pdHess
  out$pdHess2 <- m2$sd_report$pdHess
  bind_cols(data.frame(iter = i), out)
}

# est <- purrr::map_dfr(seq_len(1L), sim_test_index)
est <- furrr::future_map_dfr(seq_len(100), sim_test_index,
  .options = furrr::furrr_options(seed = TRUE))

dir.create("data/generated", showWarnings = FALSE)
saveRDS(est, "data/generated/sim-test-delta-index-vast2.rds")
est <- readRDS("data/generated/sim-test-delta-index-vast2.rds")

est$lwr_vast <- exp(est$lwr_vast)
est$upr_vast <- exp(est$upr_vast)

# est %>% filter(max_gradient > 0.001)
# est %>% filter(bad_eig)
est %>% filter(is.na(se_sim)) %>% nrow()
est %>% filter(!is.na(se_sim)) %>% nrow()

keep <- est %>% filter(!is.na(se_sim), !is.na(se_vast)) %>% pull(iter) %>% unique()
keep
est %>% filter(is.na(se_sim)) %>% pull(iter) %>% unique()
est %>% filter(is.na(se_vast)) %>% pull(iter) %>% unique()

est$year <- est$year_sim

dplyr::filter(est, iter %in% sample(keep, 6)) %>%
# b$lwr_vast <- exp(b$lwr_vast)
# b$upr_vast <- exp(b$upr_vast)
# b %>%
  ggplot(aes(year, est_sim, ymin = lwr_sim, ymax = upr_sim)) +
  geom_ribbon(fill = "black", alpha = 0.4) +
  geom_line() +
  geom_ribbon(aes(y = est_vast, ymin = lwr_vast, ymax = upr_vast),
    fill = "blue", alpha = 0.4) +
  geom_ribbon(aes(y = est_sim_tw, ymin = lwr_sim_tw, ymax = upr_sim_tw),
    fill = "green", alpha = 0.2) +
  geom_line(aes(y = est_vast), colour = "blue") +
  geom_line(aes(y = est_sim_tw), colour = "green") +
  geom_line(aes(y = true), colour = "red") +
  facet_wrap(vars(iter), scales = "free_y") +
  scale_y_log10()

sum(is.na(est$est_sim))
sum(!is.na(est$est_sim))
sum(!is.na(est$est_vast))
est <- filter(est, !is.na(est_sim), !is.na(est_vast))
# est <- filter(est, !is.infinite(upr))

coverage <- dplyr::filter(est, iter %in% keep) %>%
  mutate(covered = lwr_sim < true & upr_sim > true) %>%
  summarise(coverage = mean(covered, na.rm = TRUE))
coverage

coverage <- dplyr::filter(est, iter %in% keep) %>%
  mutate(covered = lwr_vast < true & upr_vast > true) %>%
  summarise(coverage = mean(covered, na.rm = TRUE))
coverage

# # epsilon bias correction
# bias <- est %>%
#   mutate(ratio = est/true) %>%
#   summarise(ratio = mean(ratio, na.rm = TRUE))
# bias

# MVN simulation
bias <- dplyr::filter(est, iter %in% keep) %>%
  mutate(ratio = est_sim/true) %>%
  summarise(ratio = mean(ratio, na.rm = TRUE))
bias

bias <- dplyr::filter(est, iter %in% keep) %>%
  mutate(ratio = est_sim_tw/true) %>%
  summarise(ratio = mean(ratio, na.rm = TRUE))
bias

error <- dplyr::filter(est, iter %in% keep) %>%
  mutate(error = (true - est_sim)/true) %>%
  summarise(error = median(error))
error

bias <- dplyr::filter(est, iter %in% keep) %>%
  mutate(ratio = est_vast/true) %>%
  summarise(ratio = mean(ratio, na.rm = TRUE))
bias

error <- dplyr::filter(est, iter %in% keep) %>%
  mutate(error = (true - est_vast)/true) %>%
  summarise(error = median(error))
error

bias <- dplyr::filter(est, iter %in% keep) %>%
  mutate(ratio = est_vast/est_sim) %>%
  summarise(ratio = mean(ratio, na.rm = TRUE))
bias

se_diff <- dplyr::filter(est, iter %in% keep) %>%
  mutate(ratio = se_vast/se_sim) %>%
  summarise(ratio = mean(ratio, na.rm = TRUE))
se_diff

ci_diff <- dplyr::filter(est, iter %in% keep) %>%
  mutate(ratio = (upr_vast - lwr_vast)/(upr_sim - lwr_sim)) %>%
  summarise(ratio = mean(ratio, na.rm = TRUE))
ci_diff


# ggplot(ratios_long, aes(se)) + geom_boxplot() + coord_flip()
# ggplot(ratios_long, aes(lwr)) + geom_histogram()
# ggplot(ratios_long, aes(est)) + geom_boxplot()

plan(sequential)
