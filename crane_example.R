devtools::install_github("cmjt/lgcpSPDE")
library(lgcpSPDE)
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/ecog.03771
# https://github.com/cmjt/examples/blob/master/species_distribution.md
data(cranes)
locs <- as.matrix(cranes[,2:3]) ## matrix of wetland epicentre locations
mesh <- inla.mesh.2d(loc = locs, cutoff = 0.15, max.edge = c(0.15,2))

covariates <- data.frame(Area_sc = scale(cranes$Area), PA_ratio_sc = scale(cranes$PA_ratio),
                         Wet_density_nosea_sc = scale(cranes$Wet_density_buf_NoSea),
                         Urb_density_nosea_sc = scale(cranes$Urb_density_buf_NoSea))
cranes <- cbind(cranes, covariates)

# use sdmTMB to fit the model -- spatial only
cranes$Year = as.factor(cranes$Year)
fit1 = sdmTMB(mark ~ -1 + Year + Area_sc + PA_ratio_sc,
             data=cranes,
             time="Year",
             spde = make_mesh(data = cranes, c("Lon","Lat"),
                              mesh=mesh), family = binomial(),
             spatial_only = TRUE)

# fit second model with different fields by year
fit2 = sdmTMB(mark ~ -1 + Year + Area_sc + PA_ratio_sc,
             data=cranes,
             time="Year",
             spde = make_mesh(data = cranes, c("Lon","Lat"),
                              mesh=mesh), family = binomial(),
             include_spatial = FALSE)

round(AIC(fit1) - AIC(fit2),2)
