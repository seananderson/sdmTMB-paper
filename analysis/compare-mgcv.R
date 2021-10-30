library(sdmTMB)
library(mgcv)

d <- pcod_2011
mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 5)
mesh$mesh$n
tictoc::tic("mgcv")
d$year_factor <- as.factor(d$year)
m_mgcv <- gam(present ~ 0 + year_factor + s(X, Y, by = year_factor, k = 60), # FIXME: k has huge influence on timing...
  data = d, family = binomial(link = "logit"))
# m_mgcv
# mgcv::plot.gam(m_mgcv)
tictoc::toc()
p1 <- predict(m_mgcv)

tictoc::tic("sdmTMB")
m <- sdmTMB(present ~ 0 + as.factor(year), data = d, spde = mesh, time = "year",
  fields = "AR1", include_spatial = FALSE, family = binomial(link = "logit"),
  priors = sdmTMBpriors(
    matern_st = pc_matern(range_gt = 5, range_prob = 0.05, sigma_lt = 5, sigma_prob = 0.05),
    ar1_rho = normal(0, 1),
    b = normal(0, 100)
  ))
tictoc::toc()
p2 <- predict(m)

plot(p1, p2$est);abline(a = 0, b = 1)


mesh$mesh$n

mesh2 <- inla.mesh.2d(
  loc = d[,c("X", "Y")],
  # boundary = boundary,
  # offset = c(-0.05, -0.05),
  # offset = c(-0.05, -0.05),
  max.n = 100,
  # max.n.strict = 5000,
  cutoff = 5
  # max.edge = c(100, 500)
)
plot(mesh2)

library(INLA)
fit <- gam(present ~ 0 + as.factor(year) + s(X, Y, bs = "spde", k = mesh2$n,
  xt = list(mesh = mesh2)),
  data = d,
  family= binomial(),
  method = "ML",
  control =  gam.control(scalePenalty = FALSE))
