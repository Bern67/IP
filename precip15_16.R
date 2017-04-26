## Precipitation analysis comparing 2015 and 2016 to historical precipitation and snow
## Hoonah precip data link: 


#prcip <- read.csv("precip15_16.csv")
prcp <- read.csv("prcpHoonah.csv")# in mm and celcius
d_p <- read.csv("day_precip.csv")

prcp$p_cm <- prcp$PRCP/10 # Convert from mm to cm
prcp$s_cm <-prcp$SNOW/10

## change date to Rdate
library(zoo)# needed when only year & month
## Monthly precip
prcp$Rdate <- as.Date(as.yearmon(prcp$DATE))
prcp$yr <- strftime(prcp$Rdate,"%y")
prcp$mo <- strftime(prcp$Rdate,"%m")
## Daily precip
d_p$Rdate<- strptime(as.character(d_p$DATE),"%d-%b-%y")
d_p$yr <- strftime(d_p$Rdate,"%y")
d_p$mo <- strftime(d_p$Rdate,"%m")
d_p$d <- strftime(d_p$Rdate,"%d")
d_p$j <- strftime(d_p$Rdate,"%j")
pday <- na.omit(d_p)
rm(d_p)

boxplot(pday$PRCP~pday$mo)


aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,max)
aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,min)
aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,median)
aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,length)# number of observations
mo <- aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,mean)
sum(mo$PRCP)# average yearly precip in mm

p15 <- prcp[prcp$mo==15,"PRCP"]# for use in below plots
p16 <- prcp[prcp$mo==16,"PRCP"]# only precip for 2016
s15 <- prcp[prcp$mo==15,"SNOW"]
s16 <- prcp[prcp$mo==16,"SNOW"]
m <- c(1:12)
prs <- as.data.frame(cbind(m,p15,s15,p16,s16))
rm(p15,p16,s15,s16)

with(prs,#compare snow with precipitation for fall vs spring high water events
     {plot(p15~m,type="l",col="blue", ylim=range(s15), xlab="", ylab="")# or xaxt="n", yaxt="n", will supress the entire axes
       par(new=T)
       plot(s15~m,type="l", col="blue", lty=2, ylim=range(s15), xlab="", ylab="")
       par(new=T)
       plot(p16~m,type="l", col="red", ylim=range(s15), xlab="", ylab="")
       par(new=T)
       plot(s16~m, type="l", col="red", lty=2, ylim=range(s15), xlab="", ylab="")
      title(main = "Precip & Snow", xlab="Depth (mm)", ylab="Month")
     })


## Precip plot
with(prcp,
     {season2 <-terrain.colors(12, alpha = 1)
    temp <- c(12,12,12,12,12,12,1,1,1,12,12,12)
    par(mar=c(4,4.5,2,1)+0.1)
    plot(as.factor(mo),PRCP,col=season2[temp],xlim=c(1,12),ylim=c(0,460), xaxt = "n", yaxt = 'n')
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    par(new=T)
    plot(prcp[prcp$yr==15,"mo"],prcp[prcp$yr==15,"PRCP"],xlim=c(1,12),ylim=c(0,460), cex=1.25,pch=16,col="grey17",xlab = '', ylab='')#Plot x,y
    #lines(prcp[prcp$yr==15,"mo"],prcfggfp[prcp$yr==15,"PRCP"], col='red')
    par(new=T)
    plot(prcp[prcp$yr==16,"mo"],prcp[prcp$yr==16,"PRCP"],xlim=c(1,12),ylim=c(0,460),cex=1.25,pch=17,col="grey17", ylab=expression(paste('Precipitation ', (mm/mo))),
     xlab='Month', main="")# Historical Precipitation (1941 - 2016)
    #lines(prcp[prcp$yr==16,"mo"],prcp[prcp$yr==16,"PRCP"], col='blue')
    legend("topleft",c("2015","2016"),col=c("grey17","grey17"),pch=c(16,17),cex=1.25,inset=.01)
})

## 1941 to 2016 (~35yrs) of precip 


library(ggplot2)
pr <- ggplot(prcp,aes(mo,PRCP))+
  geom_boxplot() + geom_jitter(width = 0.2,colour= "darkgreen")
pr

sn <- ggplot(prcp,aes(mo,SNOW))+
  geom_boxplot() + geom_jitter(width = 0.2,colour= "darkgreen")
sn

tmp <- ggplot(prcp,aes(mo,TAVG))+
  geom_boxplot() + geom_jitter(width = 0.2,colour= "darkgreen")
tmp



