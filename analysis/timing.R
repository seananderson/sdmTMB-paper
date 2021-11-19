library(sdmTMB)
library(dplyr)
library(ggplot2)
library(INLA)
library(inlabru)
# library(tictoc)
library(future)
plan(multisession, workers = 8L)

simulate_dat <- function(n_obs = 100,
                         cutoff = 0.1,
                         range = 0.5,
                         phi = 0.1,
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
  mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = cutoff)
  mesh$mesh$n
  sim_dat <- sdmTMB_simulate(
    formula = ~ 1 + a1,
    data = predictor_dat,
    mesh = mesh,
    family = family,
    range = range,
    phi = phi,
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
sim_fit_time <- function(n_obs = 1000, cutoff = 0.1, iter = 1, seed = sample.int(.Machine$integer.max, 1)) {
  cat("cutoff:", cutoff, "\n")
  cat("n_obs:", n_obs, "\n")
  cat("iter:", iter, "\n")
  cat("simulating...\n")
  s <- simulate_dat(n_obs = n_obs, cutoff = cutoff, seed = seed, iter = iter)
  sim_dat <- s$dat
  mesh <- s$mesh
  times <- list()

  cat("fitting...\n")
  cat("sdmTMB\n")
  out <- system.time({
    fit <- sdmTMB(observed ~ a1,
      data = sim_dat, mesh = mesh,
      priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 0.05, sigma_lt = 2))
    )
  })
  times$sdmTMB <- out[["elapsed"]]

  cat("sdmTMB normalize\n")
  out <- system.time({
    fit2 <- sdmTMB(observed ~ a1,
      data = sim_dat, mesh = mesh,
      priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 0.05, sigma_lt = 2)),
      control = sdmTMBcontrol(normalize = TRUE)
    )
  })
  times$sdmTMB_norm <- out[["elapsed"]]

  spde <- inla.spde2.pcmatern(
    mesh = mesh$mesh,
    prior.range = c(0.05, 0.05), # P(range < 10) = 0.05
    prior.sigma = c(2, 0.05)
  ) # P(sigma > 2) = 0.05
  components <- observed ~ a1 + spatrf(main = coordinates, model = spde)
  cat("inlabru\n")
  out <- system.time({
    fit_bru <- bru(
      components = components,
      family = "gaussian", data = sim_dat, options = bru_options(bru_verbose = FALSE)
    )
  })
  times$inlabru <- out[["elapsed"]]

  # bru_coefs <- fit_bru$summary.fixed[, c("mean", "sd")]
  # bru_coefs
  # fit_bru$summary.hyperpar[,c("mean", "sd")]
  # 1 / sqrt(fit_bru$summary.hyperpar[1,c("mean")]) # N(0, sigma^2)
  #

  cat("inlabru eb\n\n")
  out <- system.time({
    fit_bru2 <- bru(
      components = components,
      family = "gaussian", data = sim_dat,
      options = bru_options(bru_verbose = FALSE, control.inla = list(int.strategy = "eb"))
    )
  })
  times$inlabru_eb <- out[["elapsed"]]

  # bru_coefs <- fit_bru2$summary.fixed[, c("mean", "sd")]
  # bru_coefs
  # fit_bru2$summary.hyperpar[,c("mean", "sd")]
  # 1 / sqrt(fit_bru2$summary.hyperpar[1,c("mean")]) # N(0, sigma^2)

  out <- as_tibble(times)
  out$cutoff <- cutoff
  out$n_obs <- as.integer(n_obs)
  out$iter <- as.integer(iter)
  out$seed <- seed
  out$mesh_n <- as.integer(mesh$mesh$n)
  out
}

test <- sim_fit_time(n_obs = 100, cutoff = 0.5)
test

to_run <- expand.grid(
  n_obs = c(100, 200, 500, 1000, 2000, 5000),
  cutoff = c(0.2, 0.1, 0.05, 0.025),
  iter = seq_len(2L)
)
to_run$seed <- to_run$iter * 29172830

# out <- purrr::pmap_dfr(to_run, sim_fit_time)
out <- furrr::future_pmap_dfr(
  to_run,
  sim_fit_time,
  .progress = TRUE,
  .options = furrr::furrr_options(seed = TRUE)
)
saveRDS(out, file = "analysis/timing-cache.rds")
out <- readRDS("analysis/timing-cache.rds")

out_long <- tidyr::pivot_longer(
  out,
  sdmTMB:inlabru_eb,
  names_to = "model",
  values_to = "time"
)
clean_names <- tribble(
  ~model, ~model_clean,
  "sdmTMB", "sdmTMB",
  "sdmTMB_norm", "sdmTMB normalized",
  "inlabru", "inlabru",
  "inlabru_eb", "inlabru EB"
)
clean_names$model_clean <- factor(clean_names$model_clean, levels = clean_names$model_clean)
out_long <- left_join(out_long, clean_names)
out_long <- out_long %>%
  group_by(model_clean, cutoff) %>%
  mutate(mean_mesh_n = paste0("Mesh n = ", round(mean(mesh_n))))
out_long$mean_mesh_n <- as.factor(out_long$mean_mesh_n)
out_long$mean_mesh_n <- forcats::fct_reorder(out_long$mean_mesh_n, -out_long$cutoff)
out_long <- filter(out_long, n_obs <= 5000)
g <- ggplot(out_long, aes(n_obs, time, colour = model_clean)) +
  geom_line(lwd = 0.7) +
  scale_y_log10() +
  ggsidekick::theme_sleek() +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  scale_x_log10(breaks = unique(out_long$n_obs)) +
  # scale_x_continuous(breaks = unique(out_long$n_obs)) +
  facet_wrap(vars(mean_mesh_n), nrow = 1L) +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  labs(y = "Time (s)", x = "Number of observations", colour = "Model") +
  coord_cartesian(expand = FALSE) +
  theme(panel.spacing.x = unit(20, "pt"))
ggsave("figs/timing.pdf", width = 10, height = 3.25)
