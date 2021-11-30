library(sdmTMB)
library(dplyr)
library(ggplot2)
library(INLA)
library(spaMM)

# with sample sim data
sim_dat <- readRDS("simdat.rds")

max.edge = 0.1
loc.bnd <- matrix(c(0, 0, 1, 0, 1, 1, 0, 1), 4, 2, byrow = TRUE)
segm.bnd <- inla.mesh.segment(loc.bnd)


me <- INLA::inla.mesh.2d(
    boundary = segm.bnd,
    max.edge = c(max.edge, 0.2),
    offset = c(0.1, 0.05)
  )

mesh <- make_mesh(sim_dat, xy_cols = c("X", "Y"), mesh = me)

fit_sdmTMB  <- sdmTMB(observed ~ a1,
                data = sim_dat, mesh = mesh, family = gaussian(),
                priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 0.05, sigma_lt = 2)))


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

fit_spaMM <- fitme(observed~
                     a1 + IMRF(1|X+Y, model=spde),
                   family=gaussian, data=sim_dat)
fit_sdmTMB
fit_spaMM
