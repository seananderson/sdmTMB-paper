#' Spatio-temporal modelling
#'=========================================================
#+results="hide",warning=FALSE,message=FALSE
# rm(list=ls())

library(inlabru)
init.tutorial()
library(rgdal)
library(ggmap) # for fancy plotting
library(INLA)
library(RColorBrewer)

#' Make a shortcut to a nicer colour scale:
#+results="hide",warning=FALSE,message=FALSE
colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}

#' Read in the data
#+results="hide"

load("Nigeria/terror_data.RData")
terror.data<-data[data$country=='NGA',]
#' look at it
#+warning=FALSE,message=FALSE
head(terror.data)
class(terror.data)

studyarea<-studyarea[studyarea$sov_a3=='NGA',]
data<-data[data$country=='NGA',]
bdry <- inla.sp2segment(studyarea)
bdry$loc <- inla.mesh.map(bdry$loc)

#' this is quite a rough mesh; used here to avoid long running times
#' during the practical
#+warning=FALSE,message=FALSE
mesh3<-inla.mesh.2d(boundary=bdry, max.edge=c(80,6000)/180,cutoff=80/180,
                    crs=CRS(proj4string(studyarea))) #nv=330
plot(mesh3)

terror.data <- terror.data[terror.data$iyear>=2010 & terror.data$iyear<=2015,]
round.to <- function(x, b) {round(x/b)*b}
terror.data$originalnkill<-round.to(terror.data$nkill, 1)#some observations have non-integer values
terror.data$originallethal<-ifelse(terror.data$originalnkill==0,0,1)
resp.bern<- terror.data$originallethal

#' The data are stored as a `data.frame` and the columns holding the coorinates are
#' called `longitude` and `latitude`. Using this information we can convert the data into a
#' `SpatialPointsDataFrame` object
#+warning=FALSE,message=FALSE
coordinates(terror.data) = c("longitude","latitude")

#' The object `terror.data` now is in the proper format:
#+warning=FALSE,message=FALSE

class(terror.data)

#' and the coordinates of each point are equivalent to the previous `Lat` and `Long`:

head(terror.data)

#' Moreover, the names of the coordinates are still available using

coordnames(terror.data)

#' There is one more thing to take care of. So far we have only told the object which 
#' columns to take the coordinates from. However, it is still unaware of how to interpret these
#' coordinates since it does not know which coordinate reference system was used 
#' when the data were stored:

proj4string(terror.data)

#' The data are stored as coordinates in longitude and latitude so the 
#' corresponding proj4string is `+proj=longlat`
proj4string(terror.data) = "+proj=longlat"

#' Finally, we make the mesh live in the same CRS:
#+results="hide",warning=FALSE,message=FALSE

mesh3$crs <- inla.CRS("+proj=longlat")

#' Plot the sampling locations and the mesh
#' 
#+results="hide",warning=FALSE,message=FALSE

ggplot() + gg(mesh3) + gg(terror.data)

#' ... on a map:
#' 
#+results="hide",warning=FALSE,message=FALSE

gmap(terror.data) + gg(mesh3) + gg(terror.data)

#' We will now set up a spatio-temporal model. In order to do so, we
#' must define a temporal index (must be an integer, starting at 1)
#' and the number of discrete time points we want to model (the first
#' line is currently needed because of a bug; in the future it will
#' either be unnecessary or replaced by something else):
##+results="hide",warning=FALSE,message=FALSE

terror.data$Year <- terror.data$iyear - min(terror.data$iyear) + 1
Year <- terror.data$Year
nyear <- length(unique(terror.data$Year))

#' Note that `Year` holds values from 1 to 7 instead of the actual year, 
#' is it really is just an index!
#'
#' The spatio-temporal model is set up as follows. The `alpha=1.5`
#' parameter specifies an exponential covariance model.
##+results="hide",warning=FALSE,message=FALSE

cmp.terror = originallethal ~ Intercept+ myspde(map = coordinates, 
                                    group = Year,
                                    ngroup = nyear,
                                    model = inla.spde2.pcmatern(mesh3, alpha=1.5,
                                                                prior.range=c(0.1,0.01),
                                                                prior.sigma=c(1,0.01)),
                                    mesh = mesh3,
                                    control.group=list(model="ar1"))


#' Using this model we run a Binomial regression with `bru`
#' 
#+results="hide",warning=FALSE,message=FALSE

bru.terror = bru(cmp.terror,             
            family = "binomial",
             data= terror.data)

#' Check the resulting parameter estimates
#' 
#+results='hide',warning=FALSE,message=FALSE
summary(bru.terror)

#' ##########################
#' Predict log intensity using a predefined set of points 
#' (mesh locations and years 1 to 6)
#' 
#+results="hide",warning=FALSE,message=FALSE

df <- pixels(mesh3)
df1 <- SpatialPixelsDataFrame(df, data.frame(Year=rep(1, nrow(coordinates(df)))),
                              proj4string=CRS(proj4string(mesh3)))
df2 <- SpatialPixelsDataFrame(df, data.frame(Year=rep(2, nrow(coordinates(df)))),
                              proj4string=CRS(proj4string(mesh3)))
df3 <- SpatialPixelsDataFrame(df, data.frame(Year=rep(3, nrow(coordinates(df)))),
                              proj4string=CRS(proj4string(mesh3)))
df4 <- SpatialPixelsDataFrame(df, data.frame(Year=rep(4, nrow(coordinates(df)))),
                              proj4string=CRS(proj4string(mesh3)))
df5 <- SpatialPixelsDataFrame(df, data.frame(Year=rep(5, nrow(coordinates(df)))),
                              proj4string=CRS(proj4string(mesh3)))
df6 <- SpatialPixelsDataFrame(df, data.frame(Year=rep(6, nrow(coordinates(df)))),
                              proj4string=CRS(proj4string(mesh3)))

logint1 = predict(bru.terror, df1, ~ (myspde + Intercept))
logint2 = predict(bru.terror, df2, ~ (myspde + Intercept))
logint3 = predict(bru.terror, df3, ~ (myspde + Intercept))
logint4 = predict(bru.terror, df4, ~ (myspde + Intercept))
logint5 = predict(bru.terror, df5, ~ (myspde + Intercept))
logint6 = predict(bru.terror, df6, ~ (myspde + Intercept))

csc = scale_fill_gradientn(colours = brewer.pal(9,"YlOrRd"), limits=c(-4, 4))

# Note: gg by default picks the first variable in the output, but that's not
#       necessarily the mean; in this case we need to specify it:
multiplot(ggplot() + gg(logint1, mapping = aes_string(fill = "mean")) + csc,
          ggplot() + gg(logint2, mapping = aes_string(fill = "mean")) + csc,
          ggplot() + gg(logint3, mapping = aes_string(fill = "mean")) + csc,
          ggplot() + gg(logint4, mapping = aes_string(fill = "mean")) + csc,
          ggplot() + gg(logint5, mapping = aes_string(fill = "mean")) + csc,
          ggplot() + gg(logint6, mapping = aes_string(fill = "mean")) + csc,
          cols = 3, layout= matrix(c(1,2,3,4,5,6), nrow=3, byrow=TRUE))


