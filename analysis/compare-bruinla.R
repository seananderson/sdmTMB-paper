# from MEE inlabru paper

library(INLA)
library(inlabru)
library(RColorBrewer)
data(gorillas)

gnests <- gorillas$nests
mesh <- gorillas$mesh
boundary <- gorillas$boundary
gcounts <- gorillas$plotsample$counts
plots <- gorillas$plotsample$plots
plotnests <- gorillas$plotsample$nests

fig.plots <- ggplot() +
  gg(boundary) +
  gg(plots) +
  gg(gnests, size=0.07,color = "blue") +
  gg(plotnests,size=0.07,color = "red") +
  geom_label(aes(label=gcounts$count, x=gcounts$x+0.35, y=gcounts$y-0.3),
    fill = "white", size = 1.5, label.padding = unit(0.1, "lines")) +
  coord_fixed() + xlab("x") + ylab("y")
fig.plots

# Inference using Poisson count likelihood
# Define random variables

constr <- FALSE
mdl <- count ~ spat(main = coordinates,
  model = inla.spde2.matern(mesh, constr = constr)) +
  Intercept

# The `constr = TRUE` parameter imposes an integrate-to-zero constraint on the random field to make the model more identifiable.
# Fit the model

fit <- bru(components = mdl, family = "poisson", data = gcounts,
  options = list(E = gcounts$exposure), formula = ~ spat + Intercept)
print(fit)

# Results
# Predict the intensity

pxl <- pixels(mesh, mask = boundary, nx = 500, ny = 400)
dens <- predict(fit, pxl, ~ exp(spat+Intercept))

# Plot it

fig.brudens <- ggplot() +
  gg(dens) +
  gg(boundary, alpha = 0) +
  gg(plots) +
  scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"),
    limits = c(0,max(dens$mean))) +
  gg(gcounts, mapping = aes(size = count)) +
  scale_size_area(max_size = 2)
fig.brudens

# sdmTMB

library(sdmTMB)
d <- as.data.frame(gcounts)
# mesh <- make_mesh(d, xy_cols = c("x", "y"), )
mesh2 <- make_mesh(d, xy_cols = c("x", "y"), mesh = gorillas$mesh)
d$offset <- log(d$exposure)
fit2 <- sdmTMB(count ~ 1 + offset, family = poisson(link = "log"),
  data = d, silent = FALSE, mesh = mesh2)

nd <- as.data.frame(pxl)
nd$offset <- mean(d$offset)

p <- predict(fit2, newdata = nd)
head(p)

ggplot(p, aes(x, y, fill = exp(est))) + geom_raster() +
  scale_fill_viridis_c()
