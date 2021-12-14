library(dplyr)
library(testthat)
library(ggplot2)
library(sdmTMB)
library(future)
theme_set(ggsidekick::theme_sleek())
plan(multisession, workers = floor(availableCores() / 2))

SEED <- 1
set.seed(SEED)
x <- runif(800, -1, 1)
y <- runif(800, -1, 1)
N <- length(x)
time_steps <- 1
betas <- c(0.5, 0.7)
sigma_O <- 0.4
phi <- 0.5
.range <- 1
X <- model.matrix(~ x1, data.frame(x1 = rnorm(N * time_steps)))

true <- tribble(
  ~variable, ~true_value,
  "b0.est", betas[1],
  "b1.est", betas[2],
  "phi.est", phi,
  "range.est", .range,
  "sigma_O.est", sigma_O
)

loc <- data.frame(x = x, y = y)
spde <- make_mesh(loc, xy_cols = c("x", "y"), cutoff = 0.1)
plot(spde)

families <- list(
  Beta(),
  Gamma(link = "log"),
  # Gamma(link = "inverse"),
  binomial(),
  gaussian(),
  lognormal(),
  nbinom2(),
  nbinom1(),
  # truncated_nbinom1(),
  # truncated_nbinom2()
  poisson(),
  student(df = 3),
  tweedie()
)

fake_NAs <- function(iter, family, link) {
  structure(list(iter = iter, family = family, link = link,
    max_gradient = NA_real_, b0.est = NA_real_, b0.lwr = NA_real_,
    b0.upr = NA_real_, b1.est = NA_real_, b1.lwr = NA_real_,
    b1.upr = NA_real_, phi.est = NA_real_, phi.lwr = NA_real_,
    phi.upr = NA_real_, range.est = NA_real_, range.lwr = NA_real_,
    range.upr = NA_real_, sigma_O.est = NA_real_, sigma_O.lwr = NA_real_,
    sigma_O.upr = NA_real_), row.names = 1L, class = "data.frame")
}

est <- purrr::map_dfr(families, function(.fam) {
  cat(.fam$family, "\n")
  furrr::future_map_dfr(seq_len(8 * 16), function(i) {
  # furrr::future_map_dfr(seq_len(1), function(i) {
  # purrr::map_dfr(seq_len(1), function(i) {
    if (.fam$family == "Beta") {
      phi <- phi * 5 # small phi's can cause ~0.0 and ~1.0
    }
    if (.fam$link == "inverse") {
      betas[1] <- 5 # needs to keep all mu(i) > 0
    }
    s <- sdmTMB_sim(
      x = x, y = y, mesh = spde, X = X, sigma_V = c(0, 0),
      betas = betas, time_steps = time_steps,
      phi = phi, range = .range, sigma_O = sigma_O, sigma_E = 0, tweedie_p = 1.3,
      seed = SEED * i, family = .fam, df = 3
    )
    # ggplot(s, aes(x, y, colour = mu)) +
    #   geom_point() +
    #   scale_colour_gradient2() +
    #   facet_wrap(vars(time))
    if (.fam$family == "Beta") { # small phi's can cause ~0.0 and ~1.0
      s$observed[s$observed > 0.9999] <- 0.9999
      s$observed[s$observed < 0.0001] <- 0.0001
    }
    m <- tryCatch({sdmTMB(
      data = s, formula = observed ~ x1,
      time = "time", mesh = spde, reml = FALSE,
      spatiotemporal = "off", family = .fam
    )}, error = function(e) return(NULL))

    if (max(m$gradients) > 0.001) {
      m <- tryCatch({run_extra_optimization(m,
        nlminb_loops = 0, newton_loops = 1)},
        error = function(e) return(NULL))
    }
    if (is.null(m)) return(fake_NAs(i, .fam$family, .fam$link))

    # if (.fam$family == "truncated_nbinom1") browser()
    est <- tidy(m, conf.int = TRUE)
    est_ran <- tidy(m, "ran_pars", conf.int = TRUE)
    est <- bind_rows(est, est_ran)
    if (.fam$family == "binomial" || .fam$family == "poisson") {
      est <- bind_rows(est, data.frame(term = "phi", estimate = NA))
    }
    data.frame(
      iter = i,
      family = .fam$family,
      link = .fam$link,
      max_gradient = max(m$gradients),
      b0.est = est[est$term == "(Intercept)", "estimate"],
      b0.lwr = est[est$term == "(Intercept)", "conf.low"],
      b0.upr = est[est$term == "(Intercept)", "conf.high"],
      b1.est = est[est$term == "x1", "estimate"],
      b1.lwr = est[est$term == "x1", "conf.low"],
      b1.upr = est[est$term == "x1", "conf.high"],
      phi.est = est[est$term == "phi", "estimate"],
      phi.lwr = est[est$term == "phi", "conf.low"],
      phi.upr = est[est$term == "phi", "conf.high"],
      range.est = est[est$term == "range", "estimate"],
      range.lwr = est[est$term == "range", "conf.low"],
      range.upr = est[est$term == "range", "conf.high"],
      sigma_O.est = est[est$term == "sigma_O", "estimate"],
      sigma_O.lwr = est[est$term == "sigma_O", "conf.low"],
      sigma_O.upr = est[est$term == "sigma_O", "conf.high"],
      stringsAsFactors = FALSE
    )
  }, .options = furrr::furrr_options(seed = TRUE))
  # })
})

dir.create("data/generated", showWarnings = FALSE)
saveRDS(est, "data/generated/sim-test-family.rds")
est <- readRDS("data/generated/sim-test-family.rds")
est <- est %>% mutate(
  phi.lwr = if_else(family == "Beta", phi.lwr / 5, phi.lwr),
  phi.est = if_else(family == "Beta", phi.est / 5, phi.est),
  phi.upr = if_else(family == "Beta", phi.upr / 5, phi.upr)
)

est %>% filter(max_gradient > 0.002)
est %>% filter(sigma_O.est > 2)

est <- est %>% filter(max_gradient < 0.002) %>%
  filter(sigma_O.est < 2)

meds <- est %>% reshape2::melt() %>%
  right_join(true, by = "variable") %>%
  filter(!(family == "binomial" & variable == "phi.est")) %>%
  filter(!(family == "poisson" & variable == "phi.est")) %>%
  mutate(family_link = paste(family, link)) %>%
  group_by(family_link, variable) %>%
  mutate(median_value = median(value))

est %>%
  reshape2::melt() %>%
  right_join(true, by = "variable") %>%
  filter(!(family == "binomial" & variable == "phi.est")) %>%
  filter(!(family == "poisson" & variable == "phi.est")) %>%
  ggplot(aes(value)) +
  facet_grid(vars(paste(family, link)), vars(variable), scales = "free") +
  geom_histogram(bins = 15) +
  geom_vline(aes(xintercept = true_value), colour = "red") +
  geom_vline(aes(xintercept = median_value), colour = "blue", data = meds)
ggsave("figs/sim-test-families-pars.pdf", width = 11, height = 11)

coverage <- est %>%
  filter(family != "poisson") %>%
  filter(family != "binomial") %>%
  group_by(family, link) %>%
  mutate(covered = phi.lwr < phi & phi.upr > phi) %>%
  summarise(coverage = mean(covered))
coverage

worm_plot <- function(dat, lwr, est, upr, true, title) {
  median_est <- pull(dat, {{est}}) %>% median()
  dat %>%
    mutate(sim = seq_len(n())) %>%
    mutate(covered = {{lwr}} < true & {{upr}} > true) %>%
    ggplot(aes({{est}}, sim, xmin = {{lwr}}, xmax = {{upr}}, shape = covered)) +
    geom_vline(xintercept = true, col = "red", lty = 2) +
    geom_point(size = 1.5) +
    geom_linerange(size = 0.3) +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 21)) +
    guides(shape = "none") +
    ggtitle(title) +
    labs(x = "Parameter value") +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5))
}

est %>%
  filter(family != "poisson") %>%
  filter(family != "binomial") %>%
  worm_plot(phi.lwr, phi.est, phi.upr, phi, expression(phi)) +
  coord_cartesian(xlim = c(phi-0.2, phi+0.2)) +
  facet_wrap(~paste(family, link), scales = "free_y", ncol = 5)
ggsave("figs/sim-test-families-phi.pdf", width = 12, height = 7)

est %>%
  worm_plot(b1.lwr, b1.est, b1.upr, betas[2], "b1") +
  facet_wrap(~paste(family, link), scales = "free_y", ncol = 5) +
  coord_cartesian(xlim = c(0.5, 0.9))
ggsave("figs/sim-test-families-b1.pdf", width = 12, height = 7)

est %>%
  filter(range.upr < 1000) %>%
  worm_plot(range.lwr, range.est, range.upr, .range, "range") +
  facet_wrap(~paste(family, link), scales = "free_y", ncol = 5) +
  scale_x_log10() +
  coord_cartesian(xlim = c(0.1, 10))
ggsave("figs/sim-test-families-sigmaO.pdf", width = 12, height = 7)

est %>%
  worm_plot(sigma_O.lwr, sigma_O.est, sigma_O.upr, sigma_O, "sigma_O") +
  facet_wrap(~paste(family, link), scales = "free", ncol = 5) +
  scale_x_log10()
ggsave("figs/sim-test-families-range.pdf", width = 12, height = 7)

coverage <- est %>%
  group_by(family, link) %>%
  mutate(covered = b1.lwr < betas[2] & b1.upr > betas[2]) %>%
  summarise(coverage = mean(covered))
coverage

coverage <- est %>%
  group_by(family, link) %>%
  mutate(covered = range.lwr < .range & range.upr > .range) %>%
  summarise(coverage = mean(covered, na.rm = TRUE))
coverage

coverage <- est %>%
  group_by(family, link) %>%
  mutate(covered = sigma_O.lwr < sigma_O & sigma_O.upr > sigma_O) %>%
  summarise(coverage = mean(covered, na.rm = TRUE))
coverage

plan(sequential)
