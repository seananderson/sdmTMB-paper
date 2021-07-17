library(INLA)
library(sdmTMB)

d <- pcod_2011
yr_lu <- data.frame(year = unique(d$year), i_year = seq_along(unique(d$year)))
nyear <- length(unique(d$year))
d <- dplyr::left_join(d, yr_lu)
mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 5)
plot(mesh)
tictoc::tic("INLA")
spde <- inla.spde2.pcmatern(mesh = mesh$mesh,
  prior.range = c(5, 0.05), # P(range < 0.05) = 0.01
  prior.sigma = c(5, 0.05)) # P(sigma > 1) = 0.01
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
# h.spec <- list(theta = list(prior = 'pccor1', param = c(0, 0.9)))
h.spec <- list(theta = list(prior = 'normal', param = c(0, 1)))
formulae <- y ~ 0 + year_factor + f(i, model = spde, group = i.group,
  control.group = list(model = 'ar1', hyper = h.spec))
# PC prior on the autoreg. param.
# prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
m_inla <- inla(formulae,
  data = inla.stack.data(sdat),
  control.predictor = list(compute = TRUE, A = inla.stack.A(sdat)),
  family = "binomial",
  # control.family = list(hyper = list(theta = prec.prior)),
  control.fixed = list(expand.factor.strategy = 'inla', mean = 0, prec = 0.0001),
  verbose = FALSE,
  control.compute = list(config = TRUE)
  )
# post_inla <- inla.posterior.sample(result = m_inla, n = 30L)
# a <- unlist(lapply(post_inla, function(s) s$latent["(Intercept)", 1]))
# head(post_inla)
tictoc::toc()
summary(m_inla)

# Priors:
# names(inla.models()$prior)
# inla.models()$prior
# inla.set.control.fixed.default()[c("mean.intercept", "prec.intercept", "mean", "prec")]
# https://www.paulamoraga.com/book-geospatial/sec-inla.html


tictoc::tic("sdmTMB")
m <- sdmTMB(present ~ 0 + as.factor(year), data = d, spde = mesh, time = "year",
  fields = "AR1", include_spatial = FALSE, family = binomial(link = "logit"),
  reml = FALSE,
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 5, range_prob = 0.05, sigma_lt = 5, sigma_prob = 0.05),
    ar1_rho = normal(0, 1)
  ))
post_sdmTMB <- gather_sims(m, n_sims = 400L)
tictoc::toc()
m

head(post_sdmTMB)

plot(inla.smarginal(m_inla$marginals.fixed$year_factor2011), type = "l")

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

# ------------------

