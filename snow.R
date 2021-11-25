# owls

# library(rgdal)
library(sp)
library(spdep)
library(sf)
library(tidyverse)

#### get clean snowy owl data ####

sdat <- readRDS("owls/SNOW_data.rds")

# length(unique(sdat$CBCID))
ggplot(sdat, aes(year, count)) + geom_point() + geom_smooth()


#### spatially varying model ####

library(sdmTMB)

mesh <- make_mesh(sdat, xy_cols = c("X", "Y"), cutoff = 1.5)
plot(mesh)
mesh$mesh$n

m0 <- sdmTMB(count ~ 1 + nao + (1|year_f),
             time = "year",
             spatial_varying = ~ 0 + nao,
             family = poisson(link = "log"),
             # control = sdmTMBcontrol(normalize = TRUE),
             spatial = "on", spatiotemporal = "IID",
             data = sdat, mesh = mesh)
m0
saveRDS(m0, "owls/snow_w_main_effect_0_150km.rds")

m1 <- sdmTMB(count ~ 1 + nao + (1|year_f),
            time = "year",
            spatial_varying = ~ 0 + nao,
            family = nbinom1(link = "log"),
            # control = sdmTMBcontrol(normalize = TRUE),
            spatial = "on", spatiotemporal = "IID",
            data = sdat, mesh = mesh)
m1
saveRDS(m1, "owls/snow_w_main_effect_1_150km.rds")


m <- sdmTMB(count ~ 1 + nao + (1|year_f),
              time = "year",
              # time_varying = ~ 1,
              spatial_varying = ~ 0 + nao,
              family = nbinom2(link = "log"),
              # control = sdmTMBcontrol(normalize = TRUE),
              spatial = "on", spatiotemporal = "IID",
              # silent = F,
              reml = T,
              data = sdat, mesh = mesh)
m

saveRDS(m, "owls/snow_w_main_effect_150km_reml.rds")

# mf <- sdmTMB(count ~ 1 + nao + (1|year_f),
#               time = "year",
#               # time_varying = ~ 1,
#               spatial_varying = ~ 0 + nao,
#               family = nbinom2(link = "log"),
#               control = sdmTMBcontrol(normalize = TRUE),
#               share_range = F, # tried without priors first and results used to inform priors
#               spatial = "on", spatiotemporal = "IID",
#               # silent = F,
#               data = sdat, mesh = mesh)
# mf
#
# saveRDS(mf, "owls/snow_w_diff_ranges_150km.rds")


# mf3 <- sdmTMB(count ~ 1 + nao + (1|year_f),
#             time = "year",
#             # time_varying = ~ 1,
#             spatial_varying = ~ 0 + nao,
#             family = nbinom2(link = "log"),
#             # control = sdmTMBcontrol(normalize = TRUE),
#             share_range = F, # tried without priors first and results used to inform priors
#             # this helps estimate ranges by constraining complexity/preventing ranges getting too small
#             priors = sdmTMBpriors(matern_s = pc_matern(range_gt = 2, sigma_lt = 2),
#                                   matern_st = pc_matern(range_gt = 10, sigma_lt = 1)),
#             spatial = "on", spatiotemporal = "IID",
#             # silent = F,
#             reml = T,
#             data = sdat, mesh = mesh)
# mf3
#
# saveRDS(mf3, "owls/snow_w_reml_priors_150km.rds")

## without filtering for TotalSpecies == T when after 1997
# mr <- readRDS("owls/snow_w_main_effect.rds")
# mr

m2 <- readRDS( "owls/snow_w_main_effect_150km.rds")

AIC(m0, m1, m2)

m <- mf3

m <- readRDS( "owls/snow_w_main_effect_150km_reml.rds")

tidy(mf2, "fixed", conf.int = T)
tidy(m, "ran_pars", conf.int = T)


# check residuals
sdat$residuals <- residuals(m)
qqnorm(sdat$residuals);abline(a = 0, b = 1)

sdat$residuals0 <- residuals(m0)
qqnorm(sdat$residuals0)

sdat$residuals1 <- residuals(m1)
qqnorm(sdat$residuals1);abline(a = 0, b = 1)



library(rstan)
library(tmbstan)
s <- tmbstan(m$tmb_obj, iter = 200, chains = 1, warmup = 199)
pstan <- predict(m, tmbstan_model = s)
resid <- as.numeric(pstan) - m$data$observed
# or do any PIT residuals randomization here...
qqnorm(resid);qqline(resid) # should be right

# # try model with SOI for comparison -- a little worse
# m2 <- sdmTMB(count ~ 1 + soi, time = "year",
#              spatial_varying = ~ 0 + nao,
#              family = nbinom2(link = "log"),
#              data = sdat, mesh = mesh,
#              spatial = "on", spatiotemporal = "IID")
# m2
# tidy(m2, "ran_pars", conf.int = T)


## predict for sampling points
#
p <- predict(m)
p <- p %>% mutate(
  X10 = X, Y10 = Y,
  X = X * 100000,
  Y = Y * 100000
)


# mean(p$count)
#
# ggplot(p) + geom_point(aes(X,Y, colour = zeta_s))
# ggplot(p) + geom_point(aes(X,Y, colour = omega_s))
# ggplot(p) + geom_point(aes(X,Y, colour = epsilon_st)) + facet_wrap(~year)
# ggplot(p) + geom_point(aes(X,Y, colour = est)) + facet_wrap(~year)

# ggplot(p) + geom_point(aes(nao,count), alpha = 0.3)

# # check pattern in raw data
# plow <- filter(p, X > 0 & nao < 0)
# phi <- filter(p, X > 0 & nao > 0)
# mean(plow$count)
# mean(phi$count)
#
# plow <- filter(p, X < 0 & nao < 0)
# phi <- filter(p, X < 0 & nao > 0)
# mean(plow$count)
# mean(phi$count)




## --- get coastlines and lakes
library("rnaturalearth")
library("rnaturalearthdata")
coast <- ne_coastline(scale = "medium", returnclass = "sf")
coast2 <- coast %>% sf::st_transform(crs = proj)

land <- ne_countries(scale = "medium", returnclass = "sf")
land2  <- land  %>% sf::st_transform(crs = proj)

## --- from http://www.naturalearthdata.com/downloads/10m-physical-vectors/
lakes <- sf::st_read("owls/ne_10m_lakes")
## --- Only the largest lakes
lakes <- lakes[lakes$scalerank==0,]
# lakes <- lakes %>% sf::st_transform(crs = sf::st_crs(cbc_na_grid))

lakes <- lakes %>% sf::st_transform(crs = proj)

# make maps

# Make average predictions
# horizonal scale in bottom left
# exp(zeta)
# log10 scale
# name

# p <- droplevels(p)
p_proj <- p %>% mutate(x = X, y = Y) %>% st_as_sf(., coords = c("x", "y"), crs = proj)

# p_mean <- p %>% group_by(CBCID) %>% summarise(X = mean(X), Y= mean(Y),est = mean(est), zeta_s = mean(zeta_s))

b <- tidy(m)
p_mean <- p %>% group_by(CBCID) %>% summarise(
  X = mean(X), Y= mean(Y),
  mean_est_count = exp(mean(est)),
  nao_effect = exp(b[b$term == "nao", 2] + mean(zeta_s))
)

# mean(p_mean$zeta_s)
ggplot(data = p_proj) +
  geom_sf(data = land2, fill = "white", colour = "white",lwd = 0.35) +
  geom_sf(data = lakes, colour= "gray23", fill = "grey90", lwd = 0.35) +
  geom_point(
    # data = filter(p, year == 2020),
    data = p_mean,
             aes(X,Y, colour = nao_effect, size = mean_est_count), alpha = 0.5) +
  geom_sf(data = coast2, colour = "gray23", fill = NA, lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray23", fill = NA, lwd = 0.35) +
  coord_sf(xlim = c(min(p_proj$X)-50000, max(p_proj$X)-50000), ylim = c(min(p_proj$Y), max(p_proj$Y)))+
  # scale_fill_viridis_c() +
  scale_colour_viridis_c(
    # trans = "log", #limits =c(0.85, 1.45),
                         # breaks = c(1.0, 1.4, 1.8),
                         #breaks = scales::trans_breaks("log10", function(x) 10^x),
          guide = guide_colourbar(direction = "horizontal", title.vjust = 1,
                                  title.position = "top", label.position = "bottom")) +
  guides(size = "none") +
  labs(x= "Longitude", y = "Latitude", colour = "NAO effect\non Snowy Owl count") +
  theme_bw() + theme(
    legend.title = element_text(size= 9, hjust = 0),
    legend.key.height = unit(0.2, "cm"),
    panel.background = element_rect(fill = "grey90", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA) # get rid of legend bg
    , legend.box.background = element_rect(fill = "transparent", colour = NA), # get rid of legend panel bg
    legend.position = c(0.25,0.13),
    # legend.position = "bottom",
    axis.title = element_blank())

ggsave("owls/nao_effect_w_main_effect_150_reml.png", width = 5.5, height = 3.1)
# ggsave("snow_nao_zeta_s_grid.png", width = 6, height = 3)

ggplot(data = p_proj) +
  geom_sf(data = land2, fill = "white", colour = "white",lwd = 0.35) +
  geom_sf(data = lakes, colour= "gray23", fill = "grey90", lwd = 0.35) +
  geom_point(data = p, aes(X,Y, colour = omega_s)) +
  geom_sf(data = coast2, colour = "gray23", fill = NA, lwd = 0.35) +
  geom_sf(data = lakes, colour = "gray23", fill = NA, lwd = 0.35) +
  coord_sf(xlim = c(min(p_proj$X)-50000, max(p_proj$X)-50000), ylim = c(min(p_proj$Y), max(p_proj$Y)))+
  scale_colour_viridis_c(
    # option = "turbo", begin = 0.2,
    guide = guide_colourbar(direction = "horizontal", title.vjust = 1,
                            title.position = "top", label.position = "bottom")) +
  labs(x= "Longitude", y = "Latitude") +
  theme_bw() + theme(
    legend.title = element_text(size= 9, hjust = 0),
    legend.key.height = unit(0.2, "cm"),
    panel.background = element_rect(fill = "grey90", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    legend.position = c(0.25,0.13),
    # legend.position = "bottom",
    axis.title = element_blank())

ggsave("snow_nao_omega_150_reml.png", width = 5.5, height = 3.1)

ggplot(data = filter(p_proj, year > 1995)) +
  # geom_sf(data = land2, fill = "white", colour = "gray23", lwd = 0.35) +
  # geom_sf(data = lakes, fill = "white", colour = "gray23", lwd = 0.35) +
  geom_point(data = filter(p, year > 1995), aes(X, Y, colour = exp(est)), size = 0.4, alpha = 0.5) +
  # geom_sf(data = coast2, colour = "gray23", fill = NA, lwd = 0.35) +
  # geom_sf(data = lakes, colour = "gray23", fill = NA, lwd = 0.35) +
  coord_sf(xlim = c(min(p_proj$X)-50000, max(p_proj$X)-50000), ylim = c(min(p_proj$Y), max(p_proj$Y)))+
  scale_colour_viridis_c(trans="sqrt", na.value = "yellow",
                         limits = c(0, quantile(exp(p_proj$est), 0.995))) +
  labs(x = "Longitude", y = "Latitude", colour = "Count") +
  ggsidekick::theme_sleek() + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                    axis.title = element_blank()) +
  # theme(strip.text = element_blank(), #remove strip text
  #      strip.background = element_blank()) + #remove strip rectangles
  # geom_text(data = pred, aes(label = year, x = -15*100000, y = 0)) + # add titles using geom_text()
  facet_wrap(~year_f, drop=TRUE)

ggsave("snow_nao_annual_p_est_post1996_150.png", width = 9, height = 6.5)

glimpse(p)
