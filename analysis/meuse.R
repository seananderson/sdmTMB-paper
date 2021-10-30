library("spatstat")
library("sp")
library("maptools")
library(INLA)

library("gstat")
data(meuse)

coordinates(meuse) <- ~ x + y
proj4string(meuse) <- CRS("+init=epsg:28992")

# Code from gstat
data(meuse.grid)
coordinates(meuse.grid) <- ~ x + y
proj4string(meuse.grid) <- CRS("+init=epsg:28992")
gridded(meuse.grid) <- TRUE

meuse.bdy <- unionSpatialPolygons(
  as(meuse.grid, "SpatialPolygons"), rep(1, length(meuse.grid))
)

pts <- meuse.bdy@polygons[[1]]@Polygons[[1]]@coords
pts <- pts
mesh <- inla.mesh.2d(
  loc.domain = pts, max.edge = c(300, 600),
  offset = c(100, 250)
)
plot(mesh, asp = 1, main = "")
lines(pts, col = 3, lwd = 2)

pts2 <- pts / 1000
mesh2 <- inla.mesh.2d(
  loc.domain = pts2, max.edge = c(.300, .600),
  offset = c(.100, .250)
)

# Create SPDE
meuse.spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
A.meuse <- inla.spde.make.A(mesh = mesh, loc = coordinates(meuse))
s.index <- inla.spde.make.index(
  name = "spatial.field",
  n.spde = meuse.spde$n.spde
)

meuse.stack <- inla.stack(
  data = list(zinc = meuse$zinc),
  A = list(A.meuse, 1),
  effects = list(
    c(s.index, list(Intercept = 1)),
    list(dist = meuse$dist)
  ),
  tag = "meuse.data"
)

A.pred <- inla.spde.make.A(mesh = mesh, loc = coordinates(meuse.grid))
meuse.stack.pred <- inla.stack(
  data = list(zinc = NA),
  A = list(A.pred, 1),
  effects = list(
    c(s.index, list(Intercept = 1)),
    list(dist = meuse.grid$dist)
  ),
  tag = "meuse.pred"
)

# Join stack
join.stack <- inla.stack(meuse.stack, meuse.stack.pred)

# Fit model
form <- log(zinc) ~ -1 + Intercept + dist + f(spatial.field, model = spde)

tictoc::tic()
m1 <- inla(form,
  data = inla.stack.data(join.stack, spde = meuse.spde),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(join.stack), compute = TRUE),
  control.compute = list(cpo = FALSE, dic = FALSE)
)
tictoc::toc()

# Get predicted data on grid
index.pred <- inla.stack.index(join.stack, "meuse.pred")$data
meuse.grid$zinc.spde <- m1$summary.fitted.values[index.pred, "mean"]
meuse.grid$zinc.spde.sd <- m1$summary.fitted.values[index.pred, "sd"]


# Summary of results
summary(m1)

g_inla <- ggplot(as.data.frame(meuse.grid), aes(x, y, fill = zinc.spde)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(fill = "Estimate") +
  ggtitle("INLA") +
  coord_fixed()

g_inla_sd <- ggplot(as.data.frame(meuse.grid), aes(x, y, fill = zinc.spde.sd)) +
  geom_raster() +
  scale_fill_viridis_c(option = "A") +
  labs(fill = "SD") +
  ggtitle("INLA") +
  coord_fixed()

#############

library(sdmTMB)
d <- cbind(
  as.data.frame(coordinates(meuse)),
  data.frame(zinc = meuse$zinc, dist = meuse$dist)
)
d$x <- d$x / 1000
d$y <- d$y / 1000
# data(meuse.grid)
# coordinates(meuse.grid)= ~x+y
nd <- as.data.frame(coordinates(meuse.grid))
nd$x <- nd$x / 1000
nd$y <- nd$y / 1000
nd$dist <- meuse.grid$dist

mesh3 <- make_mesh(d, xy_cols = c("x", "y"), mesh = mesh2)
# plot(mesh3)
tictoc::tic()
m <- sdmTMB(log(zinc) ~ dist, data = d, mesh = mesh3)
tictoc::toc()

# p <- predict(m, newdata = nd)
psims <- predict(m, newdata = nd, sims = 500)
p$sd <- apply(psims, 1, sd)
p$mean <- apply(psims, 1, mean)


g_sdmTMB <- ggplot(p, aes(x, y, fill = mean)) +
  geom_raster() +
  scale_fill_viridis_c() +
  labs(fill = "Estimate") +
  ggtitle("sdmTMB") +
  coord_fixed()
g_sdmTMB_sd <- ggplot(p, aes(x, y, fill = sd)) +
  geom_raster() +
  scale_fill_viridis_c(option = "A") +
  labs(fill = "SD") +
  ggtitle("sdmTMB") +
  coord_fixed()

cowplot::plot_grid(g_inla, g_sdmTMB, g_inla_sd, g_sdmTMB_sd)
