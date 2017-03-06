## BootTrib - Selection ratio for mainstem vs tributary with bootstrap CI
## for mainstem vs trib redd counts

###------ Data ---------
habcat <- read.csv("ip_habcat.csv")
t <- habcat[,c(1,4)] ## Trib vs mainstem catagories
rm(habcat)
print(table(t$trib))#table of group counts

#Redd counts
redd <-read.csv("redd.csv")
redd <- redd[,-1]
fc <- redd[,c(1,13:19)]
str(fc)
fc <- merge(fc,t) #add tributary class
rm(t,redd)
head(fc)

##----- Bootstrap --------
library(boot)
##Pink 2015 Tributary (t) - Mainstem (m)
set.seed(167)
t_p15 <- function(data,indices){
  d = data[indices, ]
  o_t = (d$prn15[d$trib=="t"]/sum(d$prn15))# for each CC study reach
  p_t = (d$ha15[d$trib=="t"]/sum(d$ha15))
  o_m = (d$prn15[d$trib=="m"]/sum(d$prn15))# for each UV study reach
  p_m = (d$ha15[d$trib=="m"]/sum(d$ha15))
  
  wi_t = mean(o_t/p_t)
  wi_m = mean(o_m/p_m)
  wi_df = mean((o_t/p_t)-(o_m/p_m)) #difference between wi's
  wi = c(wi_t,wi_m,wi_df)
  return(wi)
}
p15 <- boot(data = fc, statistic =t_p15, strata=fc$trib,R = 3000)
p15
plot(p15,index=1)
plot(p15,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci1<-boot.ci(p15, index=1,type='bca',conf = 0.975) 
ci2<-boot.ci(p15, index=2,type='bca',conf = 0.975)
ci1$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci2$bca[1,c(4:5)]# for UV
(m1 <- mean(p15$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m2 <- mean(p15$t[,2]))#for UV

t.test(p15$t[,1],p15$t[,2])

##Chum 2015 t/m
set.seed(1267)
t_c15 <- function(data,indices){
  d = data[indices, ]
  o_t = (d$crn15[d$trib=="t"]/sum(d$crn15))# for each CC study reach
  p_t = (d$ha15[d$trib=="t"]/sum(d$ha15))
  o_m = (d$crn15[d$trib=="m"]/sum(d$crn15))# for each UV study reach
  p_m = (d$ha15[d$trib=="m"]/sum(d$ha15))
  
  wi_t = mean(o_t/p_t)
  wi_m = mean(o_m/p_m)
  wi_df = mean((o_m/p_m)-(o_t/p_t)) #difference between wi's
  wi = c(wi_t,wi_m,wi_df)
  return(wi)
}
c15 <- boot(data = fc, statistic =t_c15, strata=fc$trib,R = 3000)
c15
plot(c15,index=1)
plot(c15,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci3<-boot.ci(c15, index=1,type="bca",conf = 0.975) 
ci4<-boot.ci(c15, index=2,type="bca",conf = 0.975)
ci3$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci4$bca[1,c(4:5)]# for UV
(m3 <- mean(c15$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m4 <- mean(c15$t[,2]))#for UV

t.test(c15$t[,1],c15$t[,2])

##Pink 2016 t/m
set.seed(167)
t_p16 <- function(data,indices){
  d = data[indices, ]
  o_t = (d$prn16[d$trib=="t"]/sum(d$prn16))# for each CC study reach
  p_t = (d$ha16[d$trib=="t"]/sum(d$ha16))
  o_m = (d$prn16[d$trib=="m"]/sum(d$prn16))# for each UV study reach
  p_m = (d$ha16[d$trib=="m"]/sum(d$ha16))
  
  wi_t = mean(o_t/p_t)
  wi_m = mean(o_m/p_m)
  wi_df = mean((o_m/p_m)-(o_t/p_t)) #difference between wi's
  wi = c(wi_t,wi_m,wi_df)
  return(wi)
}
p16 <- boot(data = fc, statistic =t_p16, strata=fc$trib,R = 3000)
p16
plot(p16,index=1)
plot(p16,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci5<-boot.ci(p16, index=1,type="bca",conf = 0.975) 
ci6<-boot.ci(p16, index=2,type="bca",conf = 0.975)
ci5$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci6$bca[1,c(4:5)]# for UV
(m5 <- mean(p16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m6 <- mean(p16$t[,2]))#for UV

t.test(p16$t[,1],p16$t[,2])

##Chum 2016 t/m
set.seed(167)
t_c16 <- function(data,indices){
  d = data[indices, ]
  o_t = (d$crn16[d$trib=="t"]/sum(d$crn16))# for each CC study reach
  p_t = (d$ha16[d$trib=="t"]/sum(d$ha16))
  o_m = (d$crn16[d$trib=="m"]/sum(d$crn16))# for each UV study reach
  p_m = (d$ha16[d$trib=="m"]/sum(d$ha16))
  
  wi_t = mean(o_t/p_t)
  wi_m = mean(o_m/p_m)
  wi_df = mean((o_m/p_m)-(o_t/p_t)) #difference between wi's
  wi = c(wi_t,wi_m,wi_df)
  return(wi)
}
c16 <- boot(data = fc, statistic =t_c16, strata=fc$trib,R = 3000)
c16
plot(c16,index=1)
plot(c16,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci7<-boot.ci(c16, index=1,type="bca",conf = 0.975) 
ci8<-boot.ci(c16, index=2,type="bca",conf = 0.975)
ci7$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci8$bca[1,c(4:5)]# for UV
(m7 <- mean(c16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m8 <- mean(c16$t[,2]))#for UV

t.test(c16$t[,1],c16$t[,2])

##---- Create DF for plot of mean proportion with CI ----

bt1 <- data.frame(sp=c("p","p","c","c","p","p","c","c"), yr=c(15,15,15,15,16,16,16,16),
                 trib=c("t","m","t","m","t","m","t","m"),
                 s_y=c("p15","p15","c15","c15","p16","p16","c16","c16"),
                 id=c("p15t","p15m","p16t","p16m","c15t","c15m","c16t","c16m"),
                 w_ave=c(m1,m2,m3,m4,m5,m6,m7,m8),lci=c(ci1$bca[1,c(4)],ci2$bca[1,c(4)],
                 ci3$bca[1,c(4)],ci4$bca[1,c(4)],ci5$bca[1,c(4)],ci6$bca[1,c(4)],ci7$bca[1,
                 c(4)],ci8$bca[1,c(4)]), uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],
                 ci4$bca[1,c(5)],ci5$bca[1,c(5)],ci6$bca[1,c(5)],ci7$bca[1,c(5)],ci8$bca[1,c(5)]))

bt1 <- with(bt1,bt1[order(sp),])

##Calcuate standard deviation from CI
bt1$lsd <- abs(bt1$w_ave-bt1$lci)
bt1$usd <- abs(bt1$w_ave-bt1$uci)
bt1$yr <-as.factor(bt1$yr)

## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt1$w_stnd <- ((bt1$w_ave) - min(bt1$w_ave))/diff(range(bt1$w_ave))

## ------ BW Plot of selection ratio w Bonferroni CI ------

x <- c(1:8)
with(bt1,plot(x, w_ave,xaxt = "n",
     ylim=range(c(w_ave-lsd, w_ave+usd)),
     pch=c(16:17),bg=par(trib), cex=1.2,xlab="Species/Year", ylab=expression("Selection ratio  "(w[i])),
     main="Tributary vs Mainstem Channels"# with standard error bars
))
#Draw arrows with length of sd, & horizontal bar with tip length of .1
with(bt1,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
abline(1,0,lty=2,lwd=1,col="black")
abline(v = c(4.5), col = "grey", lty = 1)
xl <- c("C15","C16","P15","P16")
axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
with(bt1,legend("topright",c("Tributary","Mainstem"),pch=c(16:17),bg=par(trib), cex=1,inset=.01))

#dev.copy(png,"Fig/rsr_trib.png",width=800,height=600);#Copies plot to file as seen in plot window, then turns off
#dev.off()

## Pinks selected tributary habitat proportionaly higher than mainstem habitat in 2015 & 2016.  
## Chum avoided tributaries in both 2015 & 2016.
## Due to the samll sample size for tributaries the variance was high resulting in a large CI.
## Even though not significant, the pattern in the data suggests that pink select tributary habitat
## proportionly more than chum.

## During years when flows are low, pink salmon that use tributaries will be impacted at a higher rate than chum.
## Thus, during high flows pink are more likely to use tributary habitat than are chum.
## The data suggest that due to global warming there will be less snow spring melt in coastal Alaska.  This 
## This may reduce tribitary flows impacting pinks more than chum.


