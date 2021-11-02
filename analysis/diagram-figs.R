library(sdmTMB)
library(dplyr)
library(ggplot2)

predictor_dat <- expand.grid(X = seq(0, 1, length.out = 100), Y = seq(0, 1, length.out = 100))
predictor_dat <- bind_rows(
  mutate(predictor_dat, year = 1L),
  mutate(predictor_dat, year = 2L),
  mutate(predictor_dat, year = 3L),
  mutate(predictor_dat, year = 4L)
  )
mesh <- make_mesh(predictor_dat, xy_cols = c("X", "Y"), cutoff = 0.03)
plot(mesh$mesh)

head(predictor_dat)
nrow(predictor_dat)

d <- sdmTMB_simulate(
  formula = ~ 1,
  data = predictor_dat,
  time = "year",
  mesh = mesh,
  family = gaussian(link = "identity"),
  range = 0.3,
  phi = 0.1,
  sigma_O = 0.2,
  sigma_E = 0.2,
  seed = 123,
  B = c(0)
)

dir.create("figs", showWarnings = FALSE)
ggplot(filter(d, year == 1L), aes(X, Y, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/omega_s.png", width = 3, height = 3)

ggplot(filter(d, year == 1L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() + scale_fill_gradient2() + theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/eps_st1.png", width = 3, height = 3)

ggplot(filter(d, year == 2L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() + scale_fill_gradient2() + theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/eps_st2.png", width = 3, height = 3)

ggplot(filter(d, year == 3L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() + scale_fill_gradient2() + theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/eps_st3.png", width = 3, height = 3)

ggplot(filter(d, year == 4L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() + scale_fill_gradient2() + theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/zeta_s.png", width = 3, height = 3)


spde <- make_mesh(pcod, xy_cols = c("X", "Y"), cutoff = 10)
plot(spde)
m <- sdmTMB(density ~ 1, time_varying = ~ 0 + depth_scaled + depth_scaled2, data = pcod, mesh = spde, time = "year", spatial = "on", spatiotemporal = "off")
#
# m

x <- seq(-1, 1, length.out = 100)
x2 <- x^2


dat <- data.frame(x = x, x2 = x2, y = -x2*4)


set.seed(299)
devs <- rnorm(10, sd = 0.4)
d <- data.frame(dat, group = "a")
d <- bind_rows(d, mutate(dat, x = x - devs[1], group = "b"))
d <- bind_rows(d, mutate(dat, x = x - devs[2], group = "c"))
d <- bind_rows(d, mutate(dat, x = x - devs[3], group = "d"))
d <- bind_rows(d, mutate(dat, x = x - devs[4], group = "e"))
d <- bind_rows(d, mutate(dat, x = x - devs[5], group = "f"))

ggplot(d, aes(x, exp(y), colour = group)) + geom_line(lwd = 2.5) +
  theme_void() +
  coord_cartesian(ylim = c(min(exp(d$y)) + 0.001, max(exp(d$y)))) +
  scale_colour_viridis_d() +
  theme(legend.position = "none")

ggsave("figs/time-varying.pdf", width = 3.5, height = 3)
