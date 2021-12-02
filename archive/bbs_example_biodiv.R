library(bbsBayes)
library(sdmTMB)
library(ggplot2)
library(raster)
library(viridis)

fetch_bbs_data() # takes a few minutes, run once
# break BBS data into lon-lat 1 deg squares
stratified_data <- stratify(by = "latlong")
# Pull out species of interest, e.g. Smith et al. PLoS
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130768
route_strat = stratified_data$route_strat

# Biodiversity example
dat = route_strat

dat = dplyr::filter(dat,Year >= 1970, St_Abrev%in%c("AK","YT","NT")==FALSE)

# mesh$n ~ 858
mesh = make_mesh(data = dat,
                 xy_cols = c("Longitude","Latitude"),
                 cutoff = 1)
# get jday info for phenology
dat$date = paste0(dat$Month,"-",dat$Day,"-",dat$Year)
dat$date = lubridate::parse_date_time(dat$date,orders="mdy")
dat$jday = lubridate::yday(dat$date)
# fit model
fit <- sdmTMB(TotalSpp ~ -1 + s(Year,jday,k=10),
              spde = mesh,
              time="Year",
              spatial_only = TRUE,
              data=dat)

# build prediction grid for map
states <- getData(country="USA", level=1)
provinces <- getData(country="Canada", level=1)
# build a raster that's on 1/4 deg scale
r <- raster(ncol = 180*10, nrow = 360*10)
r_usa <- rasterize(states, r)
r_canada <- rasterize(provinces, r)
# convert to points
r2p_usa = rasterToPoints(r_usa)
r2p_canada = rasterToPoints(r_canada)
# filter out AK and phillipines / hawaii
r2p_usa = dplyr::filter(as.data.frame(r2p_usa), y < 50, x < 0)
# bind both together
r2p_na = rbind(as.data.frame(r2p_canada), r2p_usa)

pred_grid = data.frame(Latitude=r2p_na$y, Longitude=r2p_na$x)
pred_grid$Year = as.integer(2000)
pred_grid$jday=170
# make predictions
pred <- predict(fit,pred_grid)

# remove points that don't fall near surveys (about 50%)
pred$Latitude_1deg = floor(pred$Latitude)
pred$Longitude_1deg = floor(pred$Longitude)
dat$Latitude_1deg = floor(dat$Latitude)
dat$Longitude_1deg = floor(dat$Longitude)
dat$cell = paste0(dat$Latitude_1deg, " ", dat$Longitude_1deg)
pred$cell = paste0(pred$Latitude_1deg, " ", pred$Longitude_1deg)
pred = dplyr::filter(pred, cell %in% unique(dat$cell))

g1 = ggplot(pred,aes(Longitude,Latitude)) +
  geom_tile(aes(fill=est)) +
  scale_fill_viridis(end=0.8) +
  #scale_fill_gradient2(low="red",high="blue",midpoint=max(pred$est,na.rm=TRUE)/2,na.value="grey90") +
  ggtitle("Total BBS Species") +
  theme_bw() + ylim(c(25,56))

jpeg("sdmTMB_bbs_diversity_map_hires.jpeg")
g1
dev.off()

pred <- predict(fit,newdata)
g2 = ggplot(pred,aes(Year,jday)) +
  geom_tile(aes(fill=est)) +
  ggtitle("Total BBS Species")

g3 = dplyr::filter(pred, jday %in% c(150, 170, 190)) %>%
  ggplot(aes(Year, est, group=jday,col=jday)) +
  geom_line() + ylab("Species")

pdf("bbs_example_biodiv.pdf")
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()
