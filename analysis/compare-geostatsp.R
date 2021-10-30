library(geostatsp)
data("swissRain")
names(swissRain)
names(swissAltitude)

library(sdmTMB)
library(ggplot2)
parana <- geoR::parana
d <- data.frame(rain = parana$data, x = parana$coords[, 1], y = parana$coords[, 2])
head(d)
rain_mesh <- make_mesh(d, xy_cols = c("x", "y"), cutoff = 10)
rain_mesh$mesh$n
ggplot(d, aes(x, y, size = rain)) +
  geom_point() +
  scale_size_area()
plot(rain_mesh)

parana <- as.data.frame(parana)
coo <- as.matrix(parana[, c("east", "north")])
bnd <- inla.nonconvex.hull(coo)
meshb <- inla.mesh.2d(
  boundary = bnd, offset = c(50, 100),
  cutoff = 5, max.edge = c(75, 150)
)
plot(meshb)

rain_mesh <- make_mesh(d, xy_cols = c("x", "y"), mesh = meshb)
tictoc::tic()
fit <- sdmTMB(rain ~ 1,
  data = d, family = Gamma(link = "log"),
  spde = rain_mesh, silent = TRUE
)
tictoc::toc()
fit

bb <- bbox(parana$border)
x <- seq(bb[1, "min"] - 1, bb[1, "max"] + 1, length.out = 50)
y <- seq(bb[2, "min"] - 1, bb[2, "max"] + 1, length.out = 50)
coop <- as.matrix(expand.grid(x, y))
ind <- point.in.polygon(
  coop[, 1], coop[, 2],
  parana$border[, 1], parana$border[, 2]
)
coop <- coop[which(ind == 1), ]
coop <- as.data.frame(coop)
names(coop) <- c("x", "y")

p <- predict(fit, newdata = coop)

ggplot(p, aes(x, y, fill = exp(est))) +
  geom_raster() +
  scale_fill_viridis_c(option = "C") +
  geom_point(aes(size = rain), data = d)

psims <- predict(fit, newdata = coop, sims = 1000)
p$cv <- apply(psims, 1, function(x) sd(exp(x)) / mean(exp(x)))

ggplot(p, aes(x, y, fill = cv)) +
  geom_raster() +
  scale_fill_viridis_c(option = "B")


library(inlabru)
library(INLA)

# use inlabru to fit the same model/mesh
matern <- inla.spde2.pcmatern(meshb,
  prior.sigma = c(5, 0.05),
  prior.range = c(10, 0.05)
)
components <- data ~ field(main = coordinates, model = matern)

# add coordinates for dataframe
parana <- as.data.frame(parana)
coordinates(parana) <- c("east", "north")

tictoc::tic()
fit2 <- bru(components, data = parana)
tictoc::toc()
