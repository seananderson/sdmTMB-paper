library(bbsBayes)

fetch_bbs_data() # takes a few minutes, run once
# break BBS data into lon-lat 1 deg squares
stratified_data <- stratify(by = "latlong")
# Pull out species of interest, e.g. Smith et al. PLoS
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0130768
route_strat = stratified_data$route_strat
#Olive-sided Flycatcher in Quebec/BCR 12
indx <- grep("Olive-sided Flycatcher", stratified_data$species_strat$english)

# Pull out info from all routes
all_routes <- stratified_data$route_strat
# join data
stratified_data$bird_strat <- dplyr::filter(stratified_data$bird_strat,
                                            AOU == as.numeric(stratified_data$species_strat$aou[indx]))
# join in route info
dat <- dplyr::left_join(stratified_data$route_strat, stratified_data$bird_strat)

# grab routes for prediction later
on_routes = dplyr::filter(stratified_data$route_strat, St_Abrev=="ON") %>%
  dplyr::group_by(Route) %>%
  dplyr::summarize(Latitude=Latitude[1],Longitude=Longitude[1])

# pick out canadian data
dat <- dplyr::filter(dat, St_Abrev %in% c("AB", "BC","MB","NB","NS","NT",
                                          "NU","SK",
                                          "ON","QC","YT"))

# StopTotal is the total for all stops
dat=dplyr::select(dat, Year,StopTotal,
                  Latitude,Longitude,Month,Day,ObsN)


# filter based on longitude
dat = dplyr::filter(dat,Year >= 1970)

library(sdmTMB)
mesh = make_mesh(data = dat,
                 xy_cols = c("Longitude","Latitude"),
                 cutoff = 1)

dat$date = paste0(dat$Month,"-",dat$Day,"-",dat$Year)
dat$date = lubridate::parse_date_time(dat$date,orders="mdy")
dat$jday = lubridate::yday(dat$date)
dat$StopTotal[which(is.na(dat$StopTotal))] <- 0
dat$present = ifelse(dat$StopTotal>0,1,0)

ggplot(dat, aes(Longitude, Latitude,col=present)) +
  geom_point() +
  facet_wrap(~Year)
# full model with smooths -- does better than fixed
# effects or modeling them separate

fit <- sdmTMB(present ~ -1 + s(Year,jday,k=10),
                   spde = mesh,
                   time="Year",
                   family=binomial(),
                   spatial_only = TRUE,
                   data=dat)

newdata = expand.grid(Route=on_routes$Route,
                      Year=as.integer("2000"),
                      jday=170)
newdata <- dplyr::left_join(newdata,on_routes)
pred <- predict(fit,newdata)

# first panel could be something like a map
pred$p = plogis(pred$est)
g1 = ggplot(pred,aes(Longitude,Latitude,col=p)) +
  geom_point()

# make second prediction to show
newdata = expand.grid(Latitude = mean(on_routes$Latitude),
                      Longitude=mean(on_routes$Longitude),
                      Year=as.integer(unique(dat$Year)),
                      jday=150:190)
pred <- predict(fit,newdata)
pred$p = plogis(pred$est)
g2 = ggplot(pred,aes(Year,jday)) +
  geom_raster(aes(fill=p)) +
  ylab("Day of year")
pdf("bbs_example.pdf")
gridExtra::grid.arrange(g1,g2,nrow=1)
dev.off()


# Biodiversity example
dat = route_strat

dat = dplyr::filter(dat,Year >= 1970, St_Abrev%in%c("AK","YT","NT")==FALSE)

library(sdmTMB)
mesh = make_mesh(data = dat,
                 xy_cols = c("Longitude","Latitude"),
                 cutoff = 1)

dat$date = paste0(dat$Month,"-",dat$Day,"-",dat$Year)
dat$date = lubridate::parse_date_time(dat$date,orders="mdy")
dat$jday = lubridate::yday(dat$date)


fit <- sdmTMB(TotalSpp ~ -1 + s(Year,jday,k=10),
              spde = mesh,
              time="Year",
              spatial_only = TRUE,
              data=dat)

# come up with grids for prediction
step=1
pred_grid = expand.grid(Longitude = seq(round(min(dat$Longitude)),
                                        round(max(dat$Longitude)), by = step),
                        Latitude = seq(round(min(dat$Latitude)),
                                        round(max(dat$Latitude))), by = step)
floor_dat = dplyr::mutate(dat,
                          Latitude = floor(Latitude/step)*step,
                          Longitude = floor(Longitude/step)*step,
                          cell = paste0(Longitude,Latitude)) %>%
  dplyr::group_by(cell) %>%
  dplyr::summarise(Latitude = Latitude[1], Longitude = Longitude[1])
floor_dat$obs = 1

pred_grid = dplyr::left_join(pred_grid, floor_dat)
pred_grid$Year = as.integer(2000)
pred_grid$jday=170

pred <- predict(fit,pred_grid)
pred_grid$est = pred$est
pred_grid$est = pred_grid$est * pred_grid$obs
g1 = ggplot(pred_grid,aes(Longitude,Latitude)) +
  geom_raster(aes(fill=est)) +
  scale_fill_gradient2(low="red",high="blue",midpoint=max(pred_grid$est,na.rm=TRUE)/2,na.value="grey90") +
  ggtitle("Total BBS Species")

newdata = expand.grid("Year"=as.integer(unique(dat$Year)),
                      "jday"=150:190,
                      "Longitude"=mean(dat$Longitude),
                      "Latitude"=mean(dat$Latitude))
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
