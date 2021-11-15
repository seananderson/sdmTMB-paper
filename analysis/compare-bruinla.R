# from MEE inlabru paper

library(INLA)
library(inlabru)
library(RColorBrewer)
library(sdmTMB)

d <- pcod_2011

# mesh in inla
# mesh <- inla.mesh.2d(loc = d[,c("X","Y")], cutoff = 5)
mesh <- make_mesh(d, xy_cols = c("X", "Y"), cutoff = 5)
plot(mesh)
# spde with same priors as in the INLA example
spde <- inla.spde2.pcmatern(mesh = mesh$mesh,
                            prior.range = c(10, 0.05), # P(range < 10) = 0.05
                            prior.sigma = c(2, 0.05)) # P(sigma > 2) = 0.05

# fit sdmtmb model
fit_sdmtmb = sdmTMB(present ~ 0 + as.factor(year), data = d, spde = mesh, time = "year",
              fields = "IID", include_spatial = FALSE, family = binomial(link = "logit"),
              priors = sdmTMBpriors(
                matern_st = pc_matern(range_gt = 10, range_prob = 0.05,
                                      sigma_lt = 2, sigma_prob = 0.05)
              ))

constr <- FALSE
# set this up so year effects are independent fixed effects. This first approach is similar
# to this thread: https://githubmemory.com/repo/inlabru-org/inlabru/issues/87
components_1 <- present ~ 0 +
  fac(main = year, model = "factor_full") +
      #hyper = list(prec = list(fixed = TRUE))) +
  spatrf(main = coordinates,model = spde, group=year, ngroup=4,
         control.group = list(model = "iid"))

# note that 'd' not passed in -- b/c inlabru doesn't recognize tbl_df etc
# sdmTMB estimates (code below) are
# as.factor(year)2011    -0.69    0.45
# as.factor(year)2013     0.48    0.45
# as.factor(year)2015     0.18    0.45
# as.factor(year)2017    -0.65    0.45
fit_1 <- bru(components = components_1,
           family = "binomial", data = pcod_2011)
#print(fit)

# look at effects
fit_1$summary.random$fac


# Second option is to use model.matrix() and have fixed effects
m = lm(present ~ -1 + as.factor(year), data = d)
model_mat = model.matrix(m)
colnames(model_mat) = paste0("x",seq(2011,2017,2))

pcod_2011 = cbind(pcod_2011, model_mat)

components_2 <- present ~ -1 + x2011 + x2013 + x2015 + x2017 +
  #hyper = list(prec = list(fixed = TRUE))) +
  spatrf(main = coordinates,model = spde, group=year, ngroup=4,
         control.group = list(model = "iid"))

fit_2 <- bru(components = components_2,
             family = "binomial", data = pcod_2011)
fit_2$summary.fixed

# check that n random effects is the same as the dimension = knots
dim(fit_2$summary.random[[1]])/4

# tictoc::tic("mgcv")
#
#
# m_mgcv <- gam(present ~ 0 + year_factor + s(X, Y, by = year_factor, k = 60), # FIXME: k has huge influence on timing...
#               data = d, family = binomial(link = "logit"))

# gnests <- gorillas$nests
# mesh <- gorillas$mesh
# boundary <- gorillas$boundary
# gcounts <- gorillas$plotsample$counts
# plots <- gorillas$plotsample$plots
# plotnests <- gorillas$plotsample$nests
#
# fig.plots <- ggplot() +
#   gg(boundary) +
#   gg(plots) +
#   gg(gnests, size=0.07,color = "blue") +
#   gg(plotnests,size=0.07,color = "red") +
#   geom_label(aes(label=gcounts$count, x=gcounts$x+0.35, y=gcounts$y-0.3),
#     fill = "white", size = 1.5, label.padding = unit(0.1, "lines")) +
#   coord_fixed() + xlab("x") + ylab("y")
# fig.plots
#
# # Inference using Poisson count likelihood
# # Define random variables
#
# constr <- FALSE
# mdl <- count ~ spat(main = coordinates,
#   model = inla.spde2.matern(mesh, constr = constr)) +
#   Intercept
#
# # The `constr = TRUE` parameter imposes an integrate-to-zero constraint on the random field to make the model more identifiable.
# # Fit the model
#
# fit <- bru(components = mdl, family = "poisson", data = gcounts,
#   options = list(E = gcounts$exposure), formula = ~ spat + Intercept)
# print(fit)
#
# # Results
# # Predict the intensity
#
# pxl <- pixels(mesh, mask = boundary, nx = 500, ny = 400)
# dens <- predict(fit, pxl, ~ exp(spat+Intercept))
#
# # Plot it
#
# fig.brudens <- ggplot() +
#   gg(dens) +
#   gg(boundary, alpha = 0) +
#   gg(plots) +
#   scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"),
#     limits = c(0,max(dens$mean))) +
#   gg(gcounts, mapping = aes(size = count)) +
#   scale_size_area(max_size = 2)
# fig.brudens
#
# # sdmTMB
#
# library(sdmTMB)
# d <- as.data.frame(gcounts)
# # mesh <- make_mesh(d, xy_cols = c("x", "y"), )
# mesh2 <- make_mesh(d, xy_cols = c("x", "y"), mesh = gorillas$mesh)
# d$offset <- log(d$exposure)
# fit2 <- sdmTMB(count ~ 1 + offset, family = poisson(link = "log"),
#   data = d, silent = FALSE, mesh = mesh2)
#
# nd <- as.data.frame(pxl)
# nd$offset <- mean(d$offset)
#
# p <- predict(fit2, newdata = nd)
# head(p)
#
# ggplot(p, aes(x, y, fill = exp(est))) + geom_raster() +
#   scale_fill_viridis_c()
