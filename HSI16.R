## HSI curves

## Predictors that are best correlated (or linear combination) with field data, and with redd density.
## Chum IP persistent variable - Random Forest & PCA

#**************
##------ Data Setup ----------
redd <- read.csv("redd.csv")
r <-redd[,c(2,21:24)] #redd density
## Density plot of pink 15 redds
hist(r$pr15, xlab = "Pin redds 15", ylab = "Frequency",
        probability = TRUE, main = "Gaussian kernel",
        border = "gray")
lines(density(r$pr15, width = 12), lwd = 2)
rug(r$pr15)


hab <- read.csv("D:/R/IP16/hnl2.csv")
h <- hab[,c("ip","GRADIENT","MEANANNCMS","VWI_Floor","AREA_SQKM")]
h$Slope <- h$GRADIENT*100# %
rm(hab,redd)

pk <- stack(r[,2:3])#all chum density
cm <- stack(r[4:5])#all pink density
q <- stack(h[,c(3,3)])#flow
g <- stack(h[,c(2,2)])#gradient
vc <- stack(h[,c(4,4)])#valley confinement
ws <- stack(h[,c(5,5)])
## Rename header
names(q)[names(q) == 'values'] <- 'Mean_Q'
names(g)[names(g) == 'values'] <- 'Gradient'
names(vc)[names(vc) == 'values'] <- 'VWI'
names(ws)[names(ws) == 'values'] <- 'ws'
h2 <- cbind(q,g,vc,ws)#Stacked habitat for both years
rm(q,g,vc,ws)
#**********************
## Outliers removed-Pink
h_p <- h2[-c(1,10,46),]
pk <- pk[-c(1,10,46),]
p_ip <- cbind(pk,h_p)
p_ip <- p_ip[,-c(2,4,6,8,10)]
# Scale redd densities (0-1); Standardize denstiy
p_ip$ip <- ((p_ip$values) - min(p_ip$values))/diff(range(p_ip$values))
rm(pk,h_p)
## Outliers removed-Chum
h_c <- h2[-c(3),]
cm <- cm[-c(3),]
c_ip <- cbind(cm,h_c)
c_ip <- c_ip[,-c(2,4,6,8,10)]
# Scale redd densities (0-1); Standardize denstiy
c_ip$ip <-  ((c_ip$values) - min(c_ip$values))/diff(range(c_ip$values))
rm(cm,h_c)
#r$c_ip <- ((r$cr16) - min(r$cr16))/diff(range(r$cr16))

## HSI Curves
#--------- HSI Curves --------------

## Chum HSI for both 2015 & 2016
par(mfrow=c(2,3))
plot(c_ip$Gradient, c_ip$ip, xaxt = "n", xlab='',ylab='',xlim=c(0,.06))
par(new=T)
x <- c(0,.0030,.030,.05,.05)#variables for NetMap, 5% for PNW (Sheer 2009, WDFW 2009)
y <- c(1,1,0,0,0)#variables for NetMap
plot(x,y, type='l',lwd=2, col='darkgreen',xlab=expression(paste('Channel gradient ', (m/m))),ylab='IP Index')
#grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
#     lwd = par("lwd"), equilogs = TRUE)
## Gradient was observed to be more a barrier to Cm than Pk, as indicated by adults present.

plot(c_ip$Mean_Q, c_ip$ip, main='Chum HSI', xaxt = "n", xlab='',ylab='',xlim=c(0,15))
par(new=T)
x <- c(0,.25,5.5,11.5,15)
y <- c(0,1,1,.45,.45)
plot(x,y, type='l',lwd=2,col='darkgreen',xlab=expression(paste('Mean annual stream flow ', (m^3/s))), ylab='IP Index')
## Calibrated VWI?
plot(c_ip$VWI, c_ip$ip, xaxt = "n", xlab='',ylab='',xlim=c(1,12))#Based on random forest & RSR (Moore 2002)
par(new=T)
## 2.18 & 8.82 from VWI analysis, used all anadromous channels <5% in hnl2
x <- c(1,3.18,8.5,10,12)# Based on resource selection ratio analysis
y <- c(.3,.3,1,1,1)
plot(x,y, type='l',lwd=2,col='darkgreen',xlim=c(1,12),ylim=c(0,1),xlab='Callibrated valley width index',ylab='IP Index')


## Pink HSI for both 2015 & 2016
plot(p_ip$Gradient, p_ip$ip, xaxt = "n", xlab='',ylab='',xlim=c(0,.06))
par(new=T)
x <- c(0,.011,.011,.046,.06)#variables for NetMap, & adult Pk 2015 P/A survey (Grad.< 4.6%)
y <- c(1,1,1,0,0)
plot(x,y, type='l',lwd=2, col='darkgreen',xlab=expression(paste('Channel gradient ', (m/m))),ylab='IP Index')
## slope for pinks is slightly better because we have both IP survey & distribution

plot(p_ip$Mean_Q, p_ip$ip, main='Pink HSI', xaxt = "n", xlab='',ylab='',xlim=c(0,15))
par(new=T)
x <- c(0,.15,.5,11.5,15)
y <- c(0,1,1,.25,.25)
plot(x,y, type='l',lwd=2,col='darkgreen',xlab=expression(paste('Mean annual stream flow ', (m^3/s))),ylab='IP Index')
## 2.18 & 8.82 from VWI analysis, used all anadromous channels <5% in hnl2

plot(p_ip$VWI, p_ip$ip, xaxt = "n", xlab='',ylab='',xlim=c(1,12))
par(new=T)
x <- c(1,2.18,7.5,12)# Based on resource selection ratio analysis
y <- c(.1,.1,1,1)
plot(x,y, type='l', xlim=c(1,12), ylim=c(0,1), lwd=2, col='darkgreen',xlab='Calibrated valley-width index',ylab='IP Index')

## Use binary recursive partisioning, or moving average for thresholds.  The focus of the HSC should be
## on the main bulk of the data, not one outlier.
## Or, use 75th percentile to remove extreme outliers to seperate from lower values

## HSI values 1.0 (excellent), 0.75 (good), 0.50 (average), &  0.25 (below average) are from
## (USFWS 1981)

## Spawning flow (Bjornn 1991,p.91;)
## Flow velocity comparison below


#***********************************
## Sample size required for HSC; accurate representation of fish use of habitat
sd <- sd(redd$pr15)
z <- 1.96
E <- 1.5# target sd from the mean
(n <- (z^2)*(sd^2/E^2))

#------ Flow Velocity (m/s)  -------------
hab <- read.csv("D:/R/IP16/hnl2.csv")
redd <- read.csv("redd.csv")

qv <- hab[,17]
rd <- redd[,c(21:24)]
rm(hab,redd)

par(mfrow=c(2,2))
names<-colnames(rd) # 2015
for (i in 1:ncol(rd)){
  plot(rd[,i]~qv,main=paste("Boxplot of",names[i]),xlab="Groups", ylab=paste(names[i]))
}

#--------------
## Relative Population for Game Creek (For low and high flow years)
GC_ip <- read.csv('GC_ip.csv')
GC_ip$ip <- as.numeric(GC_ip$ip)
(GCpop <- aggregate(cbind(IP_Pink,IP_Chum)~ip,GC_ip,mean))#Mean for each ip reach by species
rm(GC_ip)

