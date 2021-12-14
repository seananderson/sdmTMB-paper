library(INLA)
library(sdmTMB)

d <- pcod_2011
yr_lu <- data.frame(year = unique(d$year), i_year = seq_along(unique(d$year)))
nyear <- length(unique(d$year))
d <- dplyr::left_join(d, yr_lu)
mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 7)
plot(mesh)
tictoc::tic("INLA")
spde <- inla.spde2.pcmatern(mesh = mesh$mesh,
  prior.range = c(5, 0.05), # P(range < 5) = 0.05
  prior.sigma = c(5, 0.05)) # P(sigma > 5) = 0.05
iset <- inla.spde.make.index('i', n.spde = spde$n.spde,
  n.group = nyear)
A <- inla.spde.make.A(mesh = mesh$mesh,
  loc = cbind(d$X, d$Y), group = d$i_year)
sdat <- inla.stack(
  data = list(y = d$present),
  A = list(A, 1),
  effects = list(iset, year_factor = as.factor(d$year)),
  tag = 'stdata')
# Pr(cor > 0) = 0.9:
h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))
# h.spec <- list(theta = list(prior = 'normal', param = c(0, 1)))
formulae <- y ~ 0 + year_factor + f(i, model = spde, group = i.group,
  control.group = list(model = 'ar1', hyper = h.spec))
# PC prior on the autoreg. param.
# prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
m_inla <- inla(formulae,
  data = inla.stack.data(sdat),
  control.predictor = list(compute = TRUE, A = inla.stack.A(sdat)),
  family = "binomial",
  # control.family = list(hyper = list(theta = prec.prior)),
  control.fixed = list(expand.factor.strategy = 'inla', mean = 0, prec = 1/(100*100)),
  verbose = FALSE,
  control.inla = list(int.strategy = "eb", strategy = "gaussian"),
  control.compute = list(config = TRUE)
  )
post_inla <- inla.posterior.sample(result = m_inla, n = 400L)
# head(post_inla)
tictoc::toc()
summary(m_inla)

# Priors:
# names(inla.models()$prior)
# inla.models()$prior
# inla.set.control.fixed.default()[c("mean.intercept", "prec.intercept", "mean", "prec")]
# https://www.paulamoraga.com/book-geospatial/sec-inla.html

tictoc::tic("sdmTMB")
m <- sdmTMB(present ~ 0 + as.factor(year), data = d, mesh = mesh, time = "year",
  spatiotemporal = "AR1", spatial = "off", family = binomial(link = "logit"),
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 5, range_prob = 0.05, sigma_lt = 5, sigma_prob = 0.05),
    ar1_rho = normal(0, 1),
    b = normal(0, 100)
  ))
post_sdmTMB <- gather_sims(m, n_sims = 2000L)
tictoc::toc()
m

head(post_sdmTMB)

library(dplyr)
library(ggplot2)
year_factor2011_inla <- unlist(lapply(post_inla, function(s) s$latent["year_factor2017:1", 1]))
year_factor2011_tmb <- filter(post_sdmTMB, .variable == "as.factor.year.2017") %>% pull(.value)

ggplot(bind_rows(data.frame(type = "INLA", post = year_factor2011_inla),
  data.frame(type = "sdmTMB", post = year_factor2011_tmb)), aes(post, colour = type)) +
  geom_density()

plot(inla.smarginal(m_inla$marginals.fixed$year_factor2011), type = "l")
hist(post_sdmTMB)

# inla.smarginal(m_inla$marginals.hyperpar$`Range for i`)

inla.qmarginal(0.025, m_inla$marginals.fixed$year_factor2011)
inla.qmarginal(0.5, m_inla$marginals.fixed$year_factor2011)
inla.qmarginal(0.975, m_inla$marginals.fixed$year_factor2011)

inla.qmarginal(0.025, m_inla$marginals.fixed$year_factor2013)
inla.qmarginal(0.5, m_inla$marginals.fixed$year_factor2013)
inla.qmarginal(0.975, m_inla$marginals.fixed$year_factor2013)

inla.qmarginal(0.025, m_inla$marginals.hyperpar$`Range for i`)
inla.qmarginal(0.5, m_inla$marginals.hyperpar$`Range for i`)
inla.qmarginal(0.975, m_inla$marginals.hyperpar$`Range for i`)

inla.qmarginal(0.025, m_inla$marginals.hyperpar$`Stdev for i`)
inla.qmarginal(0.5, m_inla$marginals.hyperpar$`Stdev for i`)
inla.qmarginal(0.975, m_inla$marginals.hyperpar$`Stdev for i`)

inla.qmarginal(0.025, m_inla$marginals.hyperpar$`GroupRho for i`)
inla.qmarginal(0.5, m_inla$marginals.hyperpar$`GroupRho for i`)
inla.qmarginal(0.975, m_inla$marginals.hyperpar$`GroupRho for i`)

tidy(m, conf.int = TRUE)
tidy(m, "ran_pars", conf.int = TRUE)

x <- inla.smarginal(m_inla$marginals.fixed$year_factor2013)
b <- filter(post_sdmTMB, .variable == "as.factor.year.2013") %>% pull(.value)
g2 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = b), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 0.1, fill = "blue", alpha = 0.2) +
  ggtitle("Year 2013")

x <- inla.smarginal(m_inla$marginals.fixed$year_factor2011)
b <- filter(post_sdmTMB, .variable == "as.factor.year.2011") %>% pull(.value)
g1 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = b), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 0.1, fill = "blue", alpha = 0.2) +
  ggtitle("Year 2011")

x <- inla.smarginal(m_inla$marginals.fixed$year_factor2015)
b <- filter(post_sdmTMB, .variable == "as.factor.year.2015") %>% pull(.value)
g2015 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = b), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 0.1, fill = "blue", alpha = 0.2) +
  ggtitle("Year 2015")

x <- inla.smarginal(m_inla$marginals.fixed$year_factor2017)
b <- filter(post_sdmTMB, .variable == "as.factor.year.2017") %>% pull(.value)
g2017 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = b), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 0.1, fill = "blue", alpha = 0.2) +
  ggtitle("Year 2017")

x <- inla.smarginal(m_inla$marginals.hyperpar$`Range for i`)
range_post <- filter(post_sdmTMB, .variable == "range") %>% pull(.value)
g3 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = range_post), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 1.5, fill = "blue", alpha = 0.2) +
  ggtitle("Range")

x <- inla.smarginal(m_inla$marginals.hyperpar$`Stdev for i`)
sigma_E_post <- filter(post_sdmTMB, .variable == "sigma_E") %>% pull(.value)
g4 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = sigma_E_post), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 0.1, fill = "blue", alpha = 0.2) +
  ggtitle("Random Field SD")

x <- inla.smarginal(m_inla$marginals.hyperpar$`GroupRho for i`)
rho_post <- filter(post_sdmTMB, .variable == "ar1_rho") %>% pull(.value)
g5 <- ggplot(as.data.frame(x), aes(x, y)) + geom_line() +
  geom_histogram(data = data.frame(x = rho_post), mapping = aes(x, after_stat(density)),
    inherit.aes = FALSE, binwidth = 0.02, fill = "blue", alpha = 0.2) +
  ggtitle("AR1 rho")

cowplot::plot_grid(g1, g2, g2015, g2017, g3, g4, g5)

# ------------------


# library(inlabru)
# library(tidyverse)
# data <- data.frame(name_idx = 1:100) %>%
#   mutate(name_val = rnorm(100),
#     mu = 100 * exp(name_val),
#     y = rpois(100, mu))
# comp <- ~ Intercept(1) + name(name_idx, model = "iid", mapper = bru_mapper_index(100))
# fit <- bru(comp, formula = y ~ Intercept + name, family = "poisson", data = data)
# newdata <- data.frame(name_idx = 1:200)
# gen <- generate(fit,
#   data = newdata,
#   formula = ~ exp(Intercept + name_eval(name_idx)),
#   n.samples = 1)
# ggplot(data, aes(name_idx, mu)) +
#   geom_point(aes(col = "Observed")) +
#   geom_point(aes(col = "Generated"), data = cbind(newdata, mu = gen[,1])) +
#   scale_y_log10()
