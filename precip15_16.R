## Precipitation analysis comparing 2015 and 2016 to historical precipitation and snow
## Hoonah precip data link: 


prcip <- read.csv("precip15_16.csv")
str(prcip)
prcip$p_cm <- prcip$PRCP/10 # Convert from mm to cm
prcip$s_cm <-prcip$SNOW/10

library(dplyr)
pcp <- prcip%>% # Use the prcipa object,
  filter(year > 2014)%>% # to filter by years < 2014,
  group_by(year, month) %>% # and group by year and month,
#  summarise(avg = mean(PRCP)%>% # then summarise by mean,
  summarise(prcp = first(PRCP),
            snow = first(SNOW),
            p_cm = first(p_cm),
            s_cm = first(s_cm))%>%
    arrange(year, month, prcp) # and arrange by means from small to large

## 1941 to 2016 (~35yrs) of precip prcipa
#bp <- boxplot(PRCP~month, prcipa=prcip,col="lightgray", horizontal=F,xlab="Month", ylab="Precipitation (mm)")



## Only 2015 & 2016 precip
sn15 <- pcp%>%
  filter(year == "2015")
sn16 <- pcp%>%
  filter(year == "2016")

#season <- heat.colors(12)
prcip$month <-factor(prcip$month) #Change month to factor for plot catagory
season2 <-terrain.colors(12, alpha = 1)
temp <- c(12,12,12,12,12,12,1,1,4,12,12,12)
plot(prcip$month, prcip$p_cm,col=season2[temp],xlim=c(1,12),ylim=c(0,46), xaxt = "n", yaxt = 'n')
prcip$month <-integer(prcip$month) #Change month from factor back to integer
par(new=T)
plot(sn15$month,sn15$p_cm,xlim=c(1,12),ylim=c(0,46), pch=17,xlab = '', ylab='')#Plot x,y
lines(sn15$month,sn15$p_cm, col='red')
par(new=T)
plot(sn16$month,sn16$p_cm,xlim=c(1,12),ylim=c(0,46),pch=16, ylab=expression(paste('Precipitation ', (cm^3/mo))),
     xlab='Month', main="HNFP Precipitation 1941 - 2016")#Plot x,y
lines(sn16$month,sn16$p_cm, col='blue')


#points(7,15.4,type='p', pch=8,col='blue',cex=2,lwd=2)
#points(8,17.1,type='p', pch=8,col='blue',cex=2,lwd=2)
#points(7,5.3,type='p', pch=8,col='red',cex=2,lwd=2)
#points(8,5.9,type='p', pch=8,col='red',cex=2,lwd=2)




library(ggplot2)
b <- ggplot(prcip,aes(month,PRCP))+
  geom_boxplot() + geom_jitter(width = 0.2,colour= "darkgreen")
b


