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
library(maptools)
## In ArcMap set dataframe projection to WGS 1984 before exporting shapefile
## so that all shapefiles are in same projection.
crk<-readShapeLines("D:/R/KetaIP/Shapefile/str.shp", 
                    proj4string=CRS("+proj=longlat"), verbose = T)
summary(crk) # Shows shapfile attributes and lat/lon extent

basin<-readShapePoly("D:/R/KetaIP/Shapefile/ip_study_basin.shp", 
                     proj4string=CRS("+proj=longlat"), verbose = T)
summary(basin)

## ------ Plots -------
## Flow(Q)
par(bg="grey87",fg="grey")
plot(basin,lwd=2,bg="white",fg="grey67",axes=T)
plot(crk,lwd=1,col="lightblue",add=T,axes=F)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=q_dt$ave_q, fg = "green",inches=0.2, add=T, lwd=2)
title(main="2016 Flow (cms)",xlab='Longitude', ylab = "Latitude")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)




## Plot density by lat long location
par(mfrow=c(2,2),mai=c(.75,.75,.75,.75),bg="white", fg="grey80")
plot(basin,lwd=2,bg="white",fg="grey67",axes=T)
plot(crk,lwd=1, add=T, col="lightblue",axes=F)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$pr15, fg = "grey17", inches=0.15, add=T, lwd=2)
text(-135.6,58.1,paste("Downstream"), col = "grey17")
text(-135.3, 57.92, paste("Upstream"), col="grey17")
text(-135.1, 57.92, paste("2015"), font=2, col="grey17")
#mtext("Pink Redd Density",side=3,line=-3,cex=1,outer=TRUE, col="grey17")
title(main="Pink Redd Density", ylab="Latitude",xlab="Longitude")

plot(basin,lwd=2,bg="white",fg="grey67",axes=T)
plot(crk,lwd=1,add=T, col="lightblue",axes=F)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$pr16, fg = "grey17", inches=0.15, add=T, lwd=2)
text(-135.1, 57.92, paste("2016"), font=2, col="grey17")
title(ylab="Latitude",xlab="Longitude")

plot(basin,lwd=2,bg="white",fg="grey67",axes=T)
plot(crk,lwd=1,add=T,col="lightblue",axes=F)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$cr15, fg = "grey17", inches=0.15, add=T, lwd=2)
text(-135.1, 57.92, col="grey17", font=2, paste("2015"))
#mtext("Chum Redd Density",side=1,line=-17,cex=1,outer=TRUE, col="grey17")
title(main="Chum redd density", ylab="Latitude",xlab="Longitude")

plot(basin,lwd=2,bg="white",fg="grey67",axes=T)
plot(crk,lwd=1,add=T,col="lightblue",axes=T)
symbols(q_dt$POINT_Y~q_dt$POINT_X,circles=redd$cr16, fg = "grey17", inches=0.15, add=T, lwd=2)
text(-135.1, 57.92, paste("2016"), font=2, col="grey17")
title(ylab="Latitude",xlab="Longitude")



## Proportion of redd for each species per reach
redd$pp15 <- 1-(redd$prn15/(redd$prn15+redd$crn15))
redd$pc15 <- 1-(redd$crn15/(redd$prn15+redd$crn15))
redd$pp16 <- 1-(redd$prn16/(redd$prn16+redd$crn16))
redd$pc16 <- 1-(redd$crn16/(redd$prn16+redd$crn16))


library(plotrix)
floating.pie(q_dt$POINT_X[13:49],q_dt$POINT_Y[13:49],x=c(redd$pp15[13:49],4),radius=4,col=c("green","red"))



