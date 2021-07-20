library(glmmTMB)
library(sdmTMB)

# all pcod2011 data with Matern far too slow with glmmTMB (45 minutes, not really converged)
# try a version on just 2017 as an example

(pos <- numFactor(d$X, d$Y))
parseNumLevels(levels(pos))
# d$times <- numFactor(d$year)
# levels(d$times)
d$pos <- pos
d$group <- factor(rep(1, nrow(d)))

tictoc::tic("glmmTMB")
d <- pcod
d2017 <- subset(d, year == 2017)
d2017$pos <- numFactor(d2017$X, d2017$Y)
d2017$group <- factor(rep(1, nrow(d2017)))
m_glmmTMB <- glmmTMB(present ~ 1 + mat(pos + 0 | group), data = d2017,
  family = binomial(link = "logit"), verbose = TRUE,
  map = list(theta = factor(c(1, 2, NA))),
  start = list(theta = c(log(1.2), log(30), log(1))) # SD, range, smoothness v (fixed at 1 usually)
)
tictoc::toc()

mesh <- make_mesh(d2017, c("X", "Y"), cutoff = 1.5)
mesh$mesh$n
m <- sdmTMB(present ~ 1, data = d2017, spde = mesh,
  family = binomial(link = "logit"), silent = FALSE
  # priors = sdmTMBpriors(
  #   matern_st = pc_matern(range_gt = 5, range_prob = 0.05, sigma_lt = 5, sigma_prob = 0.05),
  #   ar1_rho = normal(0, 1),
  #   b = normal(0, 100)
)

summary(m_glmmTMB)
# confint(m_glmmTMB, "theta")
m_glmmTMB$sdr
# exp(theta(0)) = SD?
# exp(theta(1)) /* range */,
# exp(theta(2)) /* smoothness */) );
m

tidy(m, "ran_pars", conf.int = TRUE)
s <- as.list(m_glmmTMB$sdr, "Estimate")
exp(s$theta)
m
# https://kaskr.github.io/adcomp/group__special__functions.html#ga4bdabb9b3e59feb7c552b30c46c8959c


library(INLA)

tictoc::tic("INLA")
d <- d2017
d$i_year <- 1L
mesh_inla <- inla.mesh.2d(cbind(d$X, d$Y), max.edge = c(10, 20))
plot(mesh_inla)
mesh <- mesh_inla
spde <- inla.spde2.pcmatern(mesh = mesh, # mesh$mesh,
  prior.range = c(5, 0.05), # P(range < 5) = 0.05
  prior.sigma = c(5, 0.05)) # P(sigma > 5) = 0.05
nyear <- 1
iset <- inla.spde.make.index('i', n.spde = spde$n.spde,
  n.group = nyear)
A <- inla.spde.make.A(mesh = mesh, #mesh$mesh,
  loc = cbind(d$X, d$Y), group = d$i_year)
sdat <- inla.stack(
  data = list(y = d$present),
  A = list(A, 1),
  effects = list(iset, Intercept = rep(1L, nrow(d))),
  tag = 'stdata')
formulae <- y ~ Intercept + f(i, model = spde, group = i.group,
  control.group = list(model = 'iid'))
# PC prior on the autoreg. param.
# prec.prior <- list(prior = 'pc.prec', param = c(1, 0.01))
m_inla2 <- inla(formulae,
  data = inla.stack.data(sdat),
  control.predictor = list(compute = TRUE, A = inla.stack.A(sdat)),
  family = "binomial",
  control.fixed = list(expand.factor.strategy = 'inla', mean = 0, prec = 1/(100*100)),
  verbose = F,
  control.compute = list(config = TRUE)
)
tictoc::toc()
summary(m_inla)
summary(m_inla2)
