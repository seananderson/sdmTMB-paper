# owls

library(readxl)
library(rgdal)
library(sp)
library(spdep)
library(sf)
library(tidyverse)

#### get owl data ####

d <- read_excel("owls/1-121-CBC_Count_History_Report.xlsx") %>%
  select(CBCID = Abbrev, Latitude, Longitude, surveyID = Count_yr, TotalSpecies) %>%
  filter(surveyID > 79) %>%
  # filter(TotalSpecies == TRUE) %>% # remove survey years not conducted?
  mutate(year = ifelse(surveyID <= 100, paste0(19, surveyID-1),
                       ifelse(surveyID <= 110, paste0(200, surveyID-101), paste0(20, surveyID-101))))

a <- read_excel("owls/PE-CBC_data/PE-CBC_Circle_Species_Report_SQL_updated.xlsx") %>%
  select(CBCID = Abbrev, Species = COM_NAME, Latitude, Longitude, surveyID = Count_yr,
         count = how_many, SppCount = TotalSpecies) %>%
  mutate(year = ifelse(surveyID <= 100, paste0(19, surveyID-1),
                       ifelse(surveyID <= 110, paste0(200, surveyID-101), paste0(20, surveyID-101)))) %>%
  filter(Longitude > -125 & Longitude < -60 & Latitude > 32 & Latitude < 55)


e <- read_excel("owls/PE-CBC_data/PE-CBC_Effort_Report_SQL_updated-1.xlsx") %>%
  select(CBCID = Abbrev, surveyID = Count_yr, Field_counters) %>%
  filter(surveyID > 82) %>%
  mutate(year = ifelse(surveyID <= 100, paste0(19, surveyID-1),
                       ifelse(surveyID <= 110, paste0(200, surveyID-101), paste0(20, surveyID-101))),
         year = as.integer(year))

s <- filter(a, Species == "Snowy Owl")
# select locations with at least one record (includes 0 that represent count week presence)
s_loc <- s %>% select(CBCID) %>%
  # keep only circles with records in more than 3 years
  group_by(CBCID) %>% mutate(n = n()) %>% filter(n > 0)

# get all years with records for these circles
ds <- filter(d, CBCID %in% c(unique(s_loc$CBCID)))

# merge on owl data
s <- left_join(ds, s) %>%
  mutate(lat = Latitude, lon = Longitude)
s$count[is.na(s$count)] <- 0 # fill in 0s
s <- st_as_sf(s, coords = c("lon", "lat"))
st_crs(s) <- 4326

# use Albers projection for X Y coordinates
## same as grid used in https://doi.org/10.1002/ecs2.2707
## https://github.com/tmeeha/inlaSVCBC
# cbc_na_grid <- rgdal::readOGR(dsn="owls/grid", layer="cbc_na_grid")
# proj <- raster::crs(cbc_na_grid)
# CRS arguments:
proj <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"


sxy <- s %>% sf::st_transform(crs = proj) %>%
  sf::st_coordinates() %>%
  as.data.frame()
s <- bind_cols(s, sxy)

# select relevant variables and scale X and Y for mesh
# use unit of 100 kms
sdat <- s %>% select(X, Y, count, year, CBCID, Latitude, Longitude,TotalSpecies, SppCount) %>%
  mutate(year = as.integer(year),
         year_f = as.factor(year),
         X = X / 100000,
         Y = Y / 100000
  )

# remove geometry
st_geometry(sdat) <- NULL

# check data
sdat %>% ggplot(.) + geom_point(aes(X,Y), alpha = 0.2)


check <- left_join(sdat,e)

## facet map of all snowy owl data (slow to build)
# ggplot(s) +
#   geom_sf(aes(size = count, colour = count), alpha = 0.5) +
#   coord_sf(xlim = c(-125,-60), ylim = c(30,55)) +
#   scale_colour_viridis_c("sqrt", option = "C", direction = -1) +
#   labs(fill = "Count") +
#   facet_wrap(~year) +
#   theme_light()

#### get climate data ####
# Climatic Research Unit, University of East Anglia
# https://crudata.uea.ac.uk/cru/data/nao/
nao <- read.table("owls/nao.dat", header=F)
colnames(nao) <- c("year", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, "nao")
nao[nao == -99.99] <- NA

# https://www.cpc.ncep.noaa.gov/data/indices/soi
# https://crudata.uea.ac.uk/cru/data/soi/soi.dat
soi <- read.table("owls/soi.dat", header=F)
colnames(soi) <- c("year", 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, "soi")
soi[soi == -99.99] <- NA

nao <- nao %>% select(year, nao)
soi <- soi %>% select(year, soi)


# add to dataframe
sdat <- left_join(sdat, nao) %>% left_join(., soi)

ggplot(sdat) + geom_point(aes(year, nao))+
  geom_point(aes(year, soi), colour = "red") + ylab("index")

hist(sdat$nao)


sdat %>% ggplot(.) + geom_point(aes(X,Y, colour = TotalSpecies), alpha = 0.2) + facet_wrap(~year)

# could filter for True after 1997?
# sdat <- sdat %>% filter(year < 1998 | TotalSpecies == TRUE)

# remove suspect rows in the post-1998 era that don't have a total species number
# prior to 1996 all the Canadian ones lack species totals, but still have SNOW counts
sdat <- sdat %>% filter(year < 1998 | TotalSpecies == TRUE) %>% select(-SppCount, -TotalSpecies)

saveRDS(sdat, "owls/SNOW_data.rds")



#### make a prediction grid ####

cbc_na_grid <- rgdal::readOGR(dsn="owls/grid", layer="cbc_na_grid")
cbc_na_grid$na_id <- cbc_na_grid$id; head(cbc_na_grid@data)
# plot(cbc_na_grid)

# select only grid cells relevant to this species
s_sf <- s %>% sf::st_transform(crs = sf::st_crs(cbc_na_grid))
cbc_na_sf <- st_as_sf(cbc_na_grid)
keep <- st_intersection(cbc_na_sf, s_sf)
cbc_na_grid_s <- cbc_na_sf[unlist(keep),]
# plot(cbc_na_grid_s)

# select and label varibles to match data
s_grid <- cbc_na_grid_s %>% select(id, X = centroid_x, Y = centroid_y)
st_geometry(s_grid) <- NULL
ggplot(s_grid) + geom_point(aes(X, Y), alpha = 0.1)
hist(s_grid$Y)

# removes mystery rows in far north
s_grid <- s_grid %>% unique() %>% filter(Y < 2200000)
hist(s_grid$Y)
ggplot(s_grid) + geom_point(aes(X, Y), alpha = 0.1)

s_grid_yrs <- expand.grid(
  id = unique(s_grid$id),
  year = as.integer(unique(s$year))
)

# convert to unit of 100 kms
s_grid <- left_join(s_grid_yrs, s_grid) %>% mutate(
  X = X / 100000,
  Y = Y / 100000
)

ggplot(s_grid) + geom_point(aes(X, Y), alpha = 0.1)
# # add climate to grid
# s_grid <- left_join(s_grid, nao) %>% left_join(., soi) %>%
#   mutate(nao_sd = sd(nao),
#          nao_scaled = nao/nao_sd,
#          soi_sd = sd(soi),
#          soi_scaled = soi/soi_sd)

saveRDS("owls/snow_grid.rds")

