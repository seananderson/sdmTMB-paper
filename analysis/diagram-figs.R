library(sdmTMB)
library(dplyr)
library(ggplot2)

## spatial fields
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


## Main effects plots

### linear
sens <- readRDS("all-sensor-data-processed.rds") %>% select(year, X, Y,
                                                            temp = temperature_c,
                                                            do = do_mlpl) %>% unique()

dat <- left_join(pcod, sens) %>% filter(year > 2008)
dat <- na.omit(dat)
hist(dat$do)
plot(dat$do~dat$depth)
plot(dat$do~dat$temp)

dat <- mutate(dat,
              do_mean = mean(log(do), na.rm = TRUE),
              do_sd = sd(log(do), na.rm = TRUE),
              do_scaled = (log(do) - do_mean[1]) / do_sd[1]
)

spde <- make_mesh(dat, xy_cols = c("X", "Y"), cutoff = 10)
plot(spde)

nd <- data.frame(
  do_scaled = seq(min(dat$do_scaled) + 0.2,
                  max(dat$do_scaled) - 0.2, length.out = 100),
  year = 2015 # a chosen year
)

m1 <- sdmTMB(present ~ 0 + as.factor(year) + (do_scaled), data = dat, mesh = spde, time = "year",
            family = binomial(link = "logit"),
            spatial = "on", spatiotemporal = "IID")
# m1
p1 <- predict(m1, newdata = nd, se_fit = TRUE, re_form = NA)

(pp1 <- p1 %>%
    ggplot(., aes(do_scaled, est, ymin = (est - 1.96 * est_se), ymax = (est + 1.96 * est_se))) +
    geom_line() +
    geom_ribbon(alpha = 0.1) +
    theme_void())

### gam

m2 <- sdmTMB(present ~ 0 + as.factor(year) + s(do_scaled), data = dat, mesh = spde, time = "year",
             family = binomial(link = "logit"),
             spatial = "on", spatiotemporal = "IID")
# m2
p2 <- predict(m2, newdata = nd, se_fit = TRUE, re_form = NA)

(pp2 <- p2 %>%
  ggplot(., aes(do_scaled, est, ymin = (est - 1.96 * est_se), ymax = (est + 1.96 * est_se))) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  theme_void())

### breakpoint

m3 <- sdmTMB(present ~ 0 + as.factor(year) + breakpt(do_scaled),
             data = dat, mesh = spde2, time = "year",
             family = binomial(link = "logit"),
             spatial = "on", spatiotemporal = "IID")
# m3
p3 <- predict(m3, newdata = nd, se_fit = TRUE, re_form = NA)

(pp3 <- p3 %>%
  ggplot(., aes(do_scaled, est, ymin = (est - 1.96 * est_se), ymax = (est + 1.96 * est_se))) +
  geom_line() +
  geom_ribbon(alpha = 0.1) +
  theme_void())

pp1 + pp2 + pp3 + patchwork::plot_layout(nrow = 1)

ggsave("figs/fixed-effects-do.pdf", width = 1.5, height = 0.8)


## time-varying

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



