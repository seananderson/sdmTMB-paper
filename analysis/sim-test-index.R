library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
theme_set(ggsidekick::theme_sleek())
plan(multisession, workers = floor(availableCores() / 2))
# options(future.rng.onMisuse = "ignore")

SEED <- 1
set.seed(SEED)
N <- 625 # per year
time_steps <- 20
sigma_O <- 0.6
sigma_E <- 0.3
phi <- 6
tweedie_p <- 1.6 # Tweedie p
.range <- 1

x <- seq(-1, 1, length.out = 25)
y <- seq(-1, 1, length.out = 25)
loc <- expand.grid(x = x, y = y)
x <- rep(loc$x, time_steps)
y <- rep(loc$y, time_steps)
loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.2)
plot(spde)

sim_test_index <- function(i) {
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
  s_sampled <- s %>% group_by(year) %>% slice_sample(n = 25L)

  if (FALSE) {
    ggplot(s, aes(x, y, fill = mu)) +
      geom_raster() +
      geom_point(aes(size = observed), data = s_sampled, pch = 21, fill = NA) +
      scale_fill_viridis_c() +
      scale_size_area() +
      facet_wrap(vars(year)) +
      coord_cartesian(expand = FALSE)
  }

  mesh <- make_mesh(s_sampled, xy_cols = c("x", "y"), cutoff = 0.2)
  m <- tryCatch({sdmTMB(
    data = s_sampled,
    formula = observed ~ 0 + as.factor(year),
    time = "year",
    mesh = mesh,
    reml = FALSE,
    family = tweedie(),
  )}, error = function(e) return(NULL))

  if (max(m$gradients) > 0.001) {
    m <- tryCatch({run_extra_optimization(m,
      nlminb_loops = 0, newton_loops = 1)},
      error = function(e) return(NULL))
  }

  grid <- select(s, x, y, year)
  p <- predict(m, newdata = grid, return_tmb_object = TRUE)
  index_bc <- get_index(p, bias_correct = TRUE)

  index_nbc <- get_index(p, bias_correct = FALSE)
  names(index_nbc) <- paste0(names(index_nbc), "_nbc")

  p_sims <- tryCatch({
    predict(m, newdata = grid, sims = 1000L)
  }, error = function(e) return(NULL))
  if (!is.null(p_sims)) {
    index_sims <- get_index_sims(p_sims)
    names(index_sims) <- paste0(names(index_sims), "_sim")
  } else {
    index_sims <- structure(list(year_sim = NA_real_, est_sim = NA_real_,
      lwr_sim = NA_real_,
      upr_sim = NA_real_, log_est_sim = NA_real_, se_sim = NA_real_),
      row.names = 1L, class = "data.frame")
  }

  true_index <- group_by(s, year) %>%
    summarise(true = sum(mu), .groups = "drop") %>%
    mutate(year = as.integer(year))

  if (FALSE) {
    ggplot(index_bc, aes(year, est, ymin = lwr, ymax = upr)) +
      geom_ribbon(fill = "grey70") +
      geom_line() +
      geom_line(aes(year, true), true_index, colour = "red",
        inherit.aes = FALSE)
  }

  out <- bind_cols(select(true_index, true), index_bc)
  out <- bind_cols(out, index_sims)
  out <- bind_cols(out, index_nbc)
  bind_cols(data.frame(iter = i), out)
}

est <- furrr::future_map_dfr(seq_len(100L), sim_test_index,
  .options = furrr::furrr_options(seed = TRUE))

dir.create("data/generated", showWarnings = FALSE)
saveRDS(est, "data/generated/sim-test-index.rds")
est <- readRDS("data/generated/sim-test-index.rds")

est %>% filter(max_gradient > 0.001)

ggplot(est, aes(year, est, ymin = lwr, ymax = upr)) +
  geom_ribbon(fill = "grey70") +
  geom_line() +
  geom_line(aes(y = true), colour = "red") +
  facet_wrap(vars(iter), scales = "free_y") #+
# scale_y_log10()

ggplot(est, aes(year, est_sim, ymin = lwr_sim, ymax = upr_sim)) +
  geom_ribbon(fill = "grey70") +
  geom_line() +
  geom_line(aes(y = true), colour = "red") +
  facet_wrap(vars(iter), scales = "free_y") #+

coverage <- est %>%
  mutate(covered = lwr < true & upr > true) %>%
  summarise(coverage = mean(covered))
coverage

coverage <- est %>%
  mutate(covered = lwr_sim < true & upr_sim > true) %>%
  summarise(coverage = mean(covered))
coverage

error <- est %>%
  mutate(error = (true - est)/true) %>%
  summarise(error = median(error))
error

error <- est %>%
  mutate(error = (true - est_sim)/true) %>%
  summarise(error = median(error))
error

ratios <- est %>%
  mutate(
    se = se_sim / se,
    lwr = lwr_sim / lwr,
    upr = upr_sim/upr,
    est = est_sim/est
  ) %>%
  summarise(
    se = mean(se, na.rm = TRUE),
    lwr = mean(lwr, na.rm = TRUE),
    upr = mean(upr, na.rm = TRUE),
    est = mean(est, na.rm = TRUE)
  )
ratios

plan(sequential)
