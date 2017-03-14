## Precipitation analysis comparing 2015 and 2016 to historical precipitation and snow
## Hoonah precip data link: 


#prcip <- read.csv("precip15_16.csv")
prcp <- read.csv("prcpHoonah.csv")# in mm and celcius
str(prcp)
prcp$p_cm <- prcp$PRCP/10 # Convert from mm to cm
prcp$s_cm <-prcp$SNOW/10

## change date to Rdate
library(zoo)# needed when only year & month
prcp$Rdate <- as.Date(as.yearmon(prcp$DATE))
prcp$yr <- strftime(prcp$Rdate,"%y")
prcp$mo <- strftime(prcp$Rdate,"%m")

aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,max)
aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,min)
aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,median)
aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,length)# number of observations
mo <- aggregate(cbind(PRCP,SNOW,TAVG,TMAX,TMIN)~mo,prcp,mean)
sum(mo$PRCP)# average yearly precip in mm

#prcp[prcp$yr==15,"PRCP"]# for use in below plots
#prcp[prcp$yr==16,"PRCP"]# only precip for 2016

## Precip plot
with(prcp,
     {season2 <-terrain.colors(12, alpha = 1)
    temp <- c(12,12,12,12,12,12,1,1,4,12,12,12)
    plot(as.factor(mo),PRCP,col=season2[temp],xlim=c(1,12),ylim=c(0,460), xaxt = "n", yaxt = 'n')
    par(new=T)
    plot(prcp[prcp$yr==15,"mo"],prcp[prcp$yr==15,"PRCP"],xlim=c(1,12),ylim=c(0,460), pch=".",xlab = '', ylab='')#Plot x,y
    lines(prcp[prcp$yr==15,"mo"],prcp[prcp$yr==15,"PRCP"], col='red')
    par(new=T)
    plot(prcp[prcp$yr==16,"mo"],prcp[prcp$yr==16,"PRCP"],xlim=c(1,12),ylim=c(0,460),pch=".", ylab=expression(paste('Precipitation ', (mm/mo))),
     xlab='Month', main="Historical Precipitation (1941 - 2016)")#Plot x,y
    lines(prcp[prcp$yr==16,"mo"],prcp[prcp$yr==16,"PRCP"], col='blue')
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE)
    legend("topleft",c("2015","2016"),col=c("red","blue"),pch=c(".","."), lty=c(1,1),cex=1,inset=.01)
})

## 1941 to 2016 (~35yrs) of precip 


library(ggplot2)
b <- ggplot(prcp,aes(mo,PRCP))+
  geom_boxplot() + geom_jitter(width = 0.2,colour= "darkgreen")
b


