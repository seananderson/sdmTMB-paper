library(INLA)
library(geoR)
library(sdmTMB)
library(matlab) # just for timing
library(inlabru)

data(parana)
# from inla book, https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldatatheory.html#spatial-modeling-of-rainfall-in-paran%C3%A1-brazil
parana = as.data.frame(parana)
coo <- as.matrix(parana[,c("east","north")])
bnd <- inla.nonconvex.hull(coo)
meshb <- inla.mesh.2d(
  boundary = bnd, offset = c(50, 100),
  cutoff = 1, max.edge = c(30, 60)
)

# use sdmTMB to fit the model
fit = sdmTMB(data ~ 1,
             data=parana,
             spde = make_mesh(data = parana, c("east","north"),cutoff=0.1))

# use sdmTMB to fit the model, using custom mesh above
fit = sdmTMB(data ~ 1,
             data=parana,
             spde = make_mesh(data = parana, c("east","north"),
                                        mesh=meshb))

# Add PC priors
fit = sdmTMB(data ~ 1,
             data=parana,
             spde = make_mesh(data = parana, c("east","north"),
                              mesh=meshb),
             priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 5, range_prob = 0.05,
                                                        sigma_lt = 10, sigma_prob = 0.05)))


# use inlabru to fit the same model/mesh
matern <- inla.spde2.pcmatern(meshb, prior.sigma=c(5,0.05),
                              prior.range=c(10,0.05))
components <- data ~ field(main=coordinates,model=matern)

# add coordinates for dataframe
parana <- as.data.frame(parana)
coordinates(parana) = c("east","north")

fit2 <- bru(components, data = parana)
