# Various small illustration plots for Figure 1

library(sdmTMB)
library(dplyr)
library(ggplot2)

x <- seq(-3.5, 3.5, length.out = 300)
d <- data.frame(x = x, y = dnorm(x))
ggplot(d, aes(x, y)) +
  geom_line(lwd = 2) +
  theme_void()
ggsave("figs/dnorm.pdf", width = 1.5, height = 1.5)

r <- MASS::rnegbin(1e5, mu = 2, theta = 10)
y <- table(r)
d <- data.frame(x = as.numeric(names(y)), y = as.numeric(y))
ggplot(d, aes(x, y)) +
  geom_col(width = 0.5) +
  theme_void()
ggsave("figs/nb.pdf", width = 1.5, height = 1.5)

# spatial fields

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
  formula = ~1,
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

DPI <- 25
dir.create("figs", showWarnings = FALSE)
ggplot(filter(d, year == 1L), aes(X, Y, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/omega_s.png", width = 3, height = 3, dpi = DPI)

ggplot(filter(d, year == 1L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  scale_fill_gradient2() +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/eps_st1.png", width = 3, height = 3, dpi = DPI)

ggplot(filter(d, year == 2L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  scale_fill_gradient2() +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/eps_st2.png", width = 3, height = 3, dpi = DPI)

ggplot(filter(d, year == 3L), aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  scale_fill_gradient2() +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/eps_st3.png", width = 3, height = 3, dpi = DPI)

.d <- filter(d, year == 4L)
.d$epsilon_st <- .d$epsilon_st - mean(.d$epsilon_st)
ggplot(.d, aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  scale_fill_gradient2(low = scales::muted("purple"), high = scales::muted("green")) +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/zeta_s.png", width = 3, height = 3, dpi = DPI)


ggplot(.d, aes(X, Y, fill = epsilon_st)) +
  geom_raster() +
  scale_fill_viridis_c() +
  theme_void() +
  theme(legend.position = "none") +
  coord_fixed()
ggsave("figs/zeta_s_viridis.png", width = 3, height = 3, dpi = DPI)

# Main effects plots

x <- seq(0, 1, length.out = 100)
y <- x * 2
plot(x, y)
d1 <- data.frame(x = x, y = y, type = "a - linear")

x <- seq(-2, 2, length.out = 100)
y <- -x^2
plot(x, exp(y))
d2 <- data.frame(x = x, y = exp(y), type = "b - smooth")

# $x < b_{1}$, $s(x) = x \cdot b_{0}$, and for $x > b_{1}$, $s(x) = b_{1} \cdot b_{0}$.

bkpt <- function(x, b0, b1) {
  if (x < b1) y <- x * b0
  if (x >= b1) y <- b1 * b0
  y
}

x <- seq(0, 2, length.out = 100)
y <- sapply(x, bkpt, b0 = 2, b1 = 1)
plot(x, y)
d3 <- data.frame(x = x, y = y, type = "c - bkpt")

d <- bind_rows(d1, d2, d3)

ggplot(filter(d, type == "a - linear"), aes(x, y)) +
  geom_line(lwd = 2.7) +
  facet_wrap(~type, scales = "free") +
  theme_void() +
  theme(strip.text = element_blank())
ggsave("figs/main-effects-a.pdf", width = 2.7, height = 2.7)
ggplot(filter(d, type == "b - smooth"), aes(x, y)) +
  geom_line(lwd = 2.7) +
  facet_wrap(~type, scales = "free") +
  theme_void() +
  theme(strip.text = element_blank())
ggsave("figs/main-effects-b.pdf", width = 2.7, height = 2.7)
ggplot(filter(d, type == "c - bkpt"), aes(x, y)) +
  geom_line(lwd = 2.7) +
  facet_wrap(~type, scales = "free") +
  theme_void() +
  theme(strip.text = element_blank())
ggsave("figs/main-effects-c.pdf", width = 2.7, height = 2.7)

# time-varying effects

x <- seq(-1, 1, length.out = 100)
x2 <- x^2

dat <- data.frame(x = x, x2 = x2, y = -x2 * 4)

set.seed(299)
devs <- rnorm(10, sd = 0.4)
d <- data.frame(dat, group = "a")
d <- bind_rows(d, mutate(dat, x = x - devs[1], group = "b"))
d <- bind_rows(d, mutate(dat, x = x - devs[2], group = "c"))
d <- bind_rows(d, mutate(dat, x = x - devs[3], group = "d"))
d <- bind_rows(d, mutate(dat, x = x - devs[4], group = "e"))
d <- bind_rows(d, mutate(dat, x = x - devs[5], group = "f"))

ggplot(d, aes(x, exp(y), colour = group)) +
  geom_line(lwd = 2.5) +
  theme_void() +
  coord_cartesian(ylim = c(min(exp(d$y)) + 0.001, max(exp(d$y)))) +
  scale_colour_viridis_d() +
  theme(legend.position = "none")

ggsave("figs/time-varying.pdf", width = 3.5, height = 3)

# time-varying depth

# of for pcod
d3 <- pcod

spde3 <- make_mesh(d3, xy_cols = c("X", "Y"), cutoff = 10)
plot(spde3)

m5 <- sdmTMB(density ~ 0 + as.factor(year),
  time_varying = ~ 0 + depth_scaled + depth_scaled2,
  time = "year",
  data = d3, mesh = spde3,
  family = tweedie(link = "log"),
  spatial = "on",
  spatiotemporal = "IID"
)
m5

nd5 <- expand.grid(
  depth_scaled = seq(min(d3$depth_scaled) + 0.2,
    max(d3$depth_scaled) - 0.2,
    length.out = 100
  ),
  year = unique(d3$year)
)
nd5$depth_scaled2 <- nd5$depth_scaled^2

p6 <- predict(m5, newdata = nd5, se_fit = TRUE, re_form = NA)

(pp6 <- p6 %>%
  ggplot(., aes(depth_scaled, exp(est), colour = year, group = as.factor(year))) +
  geom_line(size = 0.8, alpha = 0.6) +
  scale_color_viridis_c() +
  guides(colour = "none") +
  theme_void())

ggsave("figs/time-varying-pcod.pdf", width = 1.5, height = 1.5)
