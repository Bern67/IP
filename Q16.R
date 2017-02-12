## IP field flow 2016

## ------- Data -----------
q_dta <- read.csv('Q16.csv')
gps <- read.csv("ip_gps.csv")

attach(q_dta)
q_dta$q <- ((ave(Depth_1+Depth_2+Depth_3))*wcw*Distance_m)/Time_s #Calculate CMS and add to flow data
detach(q_dta)

library(dplyr)
q_reach <- q_dta%>%
  group_by(ReachID)%>%
  summarise(ave_q = mean(q))
q_dt <-merge(q_reach,gps)
rm(q_dta)
##
names(q_reach)[1] <- "ip"
names(gps)[2] <-"ip"
q_dt <- merge(gps,q_reach)
rm(gps,q_reach)

## Add redd density, run DFA script to line 27
redd <- read.csv("redd.csv")
redd <- redd[,-1]

#ip_fsh <- read.csv("ip_fsh.csv")
#ip_fsh <-ip_fsh[,-1]
#q_dt <- merge(q_dt,ip_fsh)
#rm(ip_fsh)

## Maps using shapefile
library(maptools) # for shapefiles
## WGS 1984 projection
crk<-readShapeLines("D:/R/KetaIP16/Shapefile/str.shp",
  proj4string=CRS("+proj=longlat"), verbose = T
)
summary(crk)
basin<-readShapePoly("D:/R/KetaIP16/Shapefile/ip_basin.shp",
                   proj4string=CRS("+proj=longlat"), verbose = T
)
summary(basin)

plot(basin[basin$Shape_Area=="secondary",],add=T,lwd=2,col="grey",axes=F)#check projection

## ------ Plots -------
## Flow(Q)
par(mfrow=c(1,1))
plot(crk,lwd=1,col="lightblue",axes=T)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=q_dt$ave_q, inches=0.2, add=T, lwd=2)
#par(new=T)
#plot(q_dt$POINT_Y~q_dt$POINT_X, type='n',ylab="Lat",xlab="Lon",axes=F)
title(main="2016 Flow (cms)")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)


## Plot density by lat long location

par(mfrow=c(2,2))
plot(crk,lwd=1,col="lightblue",axes=TRUE)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$pr15, inches=0.2, add=T, lwd=2)
par(new=T)
plot(q_dt$POINT_Y~q_dt$POINT_X, type='n',ylab="Lat",xlab="Lon",axes=F)
text(-135.51,58.09,paste("Downstream"))
text(-135.3, 57.955, paste("Upstream"))
text(-135.13, 57.95, paste("2015"))
title(main="Pink Redd Density")

plot(crk,lwd=1,col="lightblue",axes=T)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$pr16, inches=0.2, add=T, lwd=2)
par(new=T)
plot(q_dt$POINT_Y~q_dt$POINT_X, type='n',ylab="Lat",xlab="Lon", axes=F)
text(-135.13, 57.95, paste("2016"))

plot(crk,lwd=1,col="lightblue",axes=T)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$cr15, inches=0.2, add=T, lwd=2)
par(new=T)
plot(q_dt$POINT_Y~q_dt$POINT_X, type='n',ylab="Lat",xlab="Lon",axes=F)
text(-135.13, 57.95, paste("2015"))
title(main="Chum Redd Density")

plot(crk,lwd=1,col="lightblue",axes=T)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$cr16, inches=0.2, add=T, lwd=2)
par(new=T)
plot(q_dt$POINT_Y~q_dt$POINT_X,type='n',ylab="Lat",xlab="Lon",axes=F)
text(-135.13, 57.95, paste("2016"))



#************************

## Examples

nc_SP <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1],
                       proj4string=CRS("+proj=longlat +ellps=clrk66"))
## Not run:
pls <- slot(nc_SP, "polygons")
pls_new <- lapply(pls, checkPolygonsHoles)
nc_SP <- SpatialPolygonsDataFrame(SpatialPolygons(pls_new,
                                                  proj4string=CRS(proj4string(nc_SP))), data=as(nc_SP, "data.frame"))
## End(Not run)
try1 <- dotsInPolys(nc_SP, as.integer(nc_SP$SID74))
plot(nc_SP, axes=TRUE)
plot(try1, add=TRUE, pch=18, col="red")
try2 <- dotsInPolys(nc_SP, as.integer(nc_SP$SID74), f="regular")
plot(nc_SP, axes=TRUE)
plot(try2, add=TRUE, pch=18, col="red")


