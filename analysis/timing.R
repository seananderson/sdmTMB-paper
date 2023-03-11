# Timing figure -----------------------------------------------------------

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
    fit <- sdmTMB(
      observed ~ a1,
      data = sim_dat,
      mesh = mesh,
      family = family,
      control = sdmTMBcontrol(newton_loops = 0L),
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

saveRDS(out, file = "analysis/timing-cache-parallel-2023-03-10.rds")
out <- readRDS("analysis/timing-cache-parallel-2023-03-10.rds")

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

leg <- tribble(
  ~model_clean, ~mean_mesh_n, ~time,
  "sdmTMB", 290, 0.40,
  "inlabru EB", 155, 3.8,
  "spaMM", 442, 1.95,
  "mgcv::bam\n(discretize = TRUE) SPDE", 155, 5
)

g <- out_long_sum %>%
  ggplot(aes(mean_mesh_n, time, colour = model_clean)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = model_clean),
    alpha = 0.2,
    colour = NA
  ) +
  geom_line(lwd = 0.7) +
  scale_x_log10(breaks = c(250, 500, 1000)) +
  theme_sleek() +
  theme(panel.grid.major = element_line(colour = "grey90")) +
  theme(legend.position = c(0.4, 0.85)) +
  geom_text(aes(label = model_clean), data = leg, hjust = 0) +
  scale_fill_manual(values = colour_vect) +
  scale_colour_manual(values = colour_vect) +
  labs(y = "Time (s)", x = "Mesh nodes", colour = "Model", fill = "Model") +
  coord_cartesian(expand = FALSE, ylim = c(0, 5.9)) +
  theme(panel.spacing.x = unit(20, "pt"), legend.position = "none")
print(g)

.width <- 5
ggsave("figs/timing-spatial-2023-01-24-xkcd.pdf", width = .width, height = .width / 1.3)

cat("inlabru to sdmTMB timing ratios\n")
x <- filter(
  out_long_sum,
  model_clean == "inlabru EB"
)$time /
  filter(out_long_sum, model_clean == "sdmTMB")$time
round(x, 1)

sessionInfo()
