## Selection Ratio by confinement (Manly 2002) w Bootstrap CI
## Reference Btutorial.pdf presentation for example; Zotero

## Need to run one-way ANOVA to on ranked data to test sig dif between species and channel unit features among years

#***************************
##--- 1) VWI (Wi) Data ------
#***************************
hnl2 <- read.csv("D:/R/IP16/hnl2.csv", header=T)
vc <- hnl2[,c('ip', 'ValCnstrnt','VWI_Floor')]

## Channel confinement groups (Grant 1995,Moore 2002), based on VWI calibrated
vc$c <- with(vc,ifelse(vc$VWI_Floor <= 3.09, "CC",#Constrained Canyon; from VWI analysis
                              ifelse(vc$VWI_Floor > 3.09, "UV", "")))#Unconstrained Valley
## Need to take the median of CC and UV for threshold HSI values (Burnett 2007, pg 69)
rm(hnl2)
print(table(vc$c))#table of group counts
## Sample size too small, combine T & UV
c <- vc$c

## Netmap confinement default for hnl2 is CC < 4, Tran 4-6, & UV > 6 (Hall 2007)

#Redd counts
redd <-read.csv("redd.csv")
redd <- redd[,-1]
fc <- redd[,c(1,13:19)]
str(fc)
fc <- cbind(fc,c) #add contrained class
head(fc)

#hist(fc$hadif) # Some reaches had side channel in 2015, and 2016, & some did not.
#plot(fc$ha15~fc$ha16)


#-------Bootstrap----------------
library(boot)
##Pink 2015 CC/UV
set.seed(167)
c_p15 <- function(data,indices){
  d = data[indices, ]
  o_cc = (d$prn15[d$c=="CC"]/sum(d$prn15))# for each CC study reach
  p_cc = (d$ha15[d$c=="CC"]/sum(d$ha15))
  o_uv = (d$prn15[d$c=="UV"]/sum(d$prn15))# for each UV study reach
  p_uv = (d$ha15[d$c=="UV"]/sum(d$ha15))
  wi_cc = mean(o_cc/p_cc)
  wi_uv = mean(o_uv/p_uv)
  wi_df = mean((o_uv/p_uv)-(o_cc/p_cc)) #difference between wi's
  wi = c(wi_cc,wi_uv,wi_df)
  return(wi)
}
p15 <- boot(data = fc, statistic =c_p15, strata=fc$c,R = 3000)
p15
plot(p15,index=1)
plot(p15,index=2)

# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci1<-boot.ci(p15, index=1,type='bca',conf = 0.975)
ci1$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci2<-boot.ci(p15, index=2,type='bca',conf = 0.975)
ci2$bca[1,c(4:5)]# for UV
(m1 <- mean(p15$t[,1]))#Bootstrap mean corrected for bias, report this mean (CC)
(m2 <- mean(p15$t[,2]))#for UV

t.test(p15$t[,1],p15$t[,2])#independant 2-group t.test

##Chum 2015 CC/UV
set.seed(167)
c_c15 <- function(data,indices){
  d = data[indices, ]
  o_cc = (d$crn15[d$c=="CC"]/sum(d$crn15))# for each CC study reach
  p_cc = (d$ha15[d$c=="CC"]/sum(d$ha15))
  o_uv = (d$crn15[d$c=="UV"]/sum(d$crn15))# for each UV study reach
  p_uv = (d$ha15[d$c=="UV"]/sum(d$ha15))
  wi_cc = mean(o_cc/p_cc)
  wi_uv = mean(o_uv/p_uv)
  wi_df = mean((o_uv/p_uv)-(o_cc/p_cc)) #difference between wi's
  wi = c(wi_cc,wi_uv,wi_df)
  return(wi)
}
c15 <- boot(data = fc, statistic =c_c15, strata=fc$c,R = 3000)
c15
plot(c15,index=1)
plot(c15,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci3<-boot.ci(c15, index=1,type="bca",conf = 0.975)
ci3$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci4<-boot.ci(c15, index=2,type="bca",conf = 0.975)
ci4$bca[1,c(4:5)]# for UV
m3 <- mean(c15$t[,1])#Bootstrap mean corrected for bias, report this mean
m4 <- mean(c15$t[,2])#for UV

t.test(c15$t[,1],c15$t[,2])

##Pink 2016 CC/UV
set.seed(167)
c_p16 <- function(data,indices){
  d = data[indices, ]
  o_cc = (d$prn16[d$c=="CC"]/sum(d$prn16))# for each CC study reach
  p_cc = (d$ha16[d$c=="CC"]/sum(d$ha16))
  o_uv = (d$prn16[d$c=="UV"]/sum(d$prn16))# for each UV study reach
  p_uv = (d$ha16[d$c=="UV"]/sum(d$ha16))
  
  wi_cc = mean(o_cc/p_cc)
  wi_uv = mean(o_uv/p_uv)
  wi_df = mean((o_uv/p_uv)-(o_cc/p_cc)) #difference between wi's
  wi = c(wi_cc,wi_uv,wi_df)
  return(wi)
}
p16 <- boot(data = fc, statistic =c_p16, strata=fc$c,R = 3000)
p16
plot(p16,index=1)
plot(p16,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci5<-boot.ci(p16, index=1,type="bca",conf = 0.975)
ci5$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci6<-boot.ci(p16, index=2,type="bca",conf = 0.975)
ci6$bca[1,c(4:5)]# for UV
(m5 <- mean(p16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m6 <- mean(p16$t[,2]))#for UV

t.test(p16$t[,1],p16$t[,2])

##Chum 2016 CC/UV
set.seed(167)
c_c16 <- function(data,indices){
  d = data[indices, ]
  o_cc = (d$crn16[d$c=="CC"]/sum(d$crn16))# for each CC study reach
  p_cc = (d$ha16[d$c=="CC"]/sum(d$ha16))
  o_uv = (d$crn16[d$c=="UV"]/sum(d$crn16))# for each UV study reach
  p_uv = (d$ha16[d$c=="UV"]/sum(d$ha16))
  wi_cc = mean(o_cc/p_cc)
  wi_uv = mean(o_uv/p_uv)
  wi_df = mean((o_uv/p_uv)-(o_cc/p_cc)) #difference between wi's
  wi = c(wi_cc,wi_uv,wi_df)
  return(wi)
}
c16 <- boot(data = fc, statistic =c_c16, strata=fc$c,R = 3000)
c16
plot(c16,index=1)
plot(c16,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci7<-boot.ci(c16, index=1,type="bca",conf = 0.975)
ci7$bca[1,c(4:5)]# BCa 97.5%CI, for CC
ci8<-boot.ci(c16, index=2,type="bca",conf = 0.975)
ci8$bca[1,c(4:5)]# for UV
(m7 <- mean(c16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m8 <- mean(c16$t[,2]))#for UV

t.test(c16$t[,1],c16$t[,2])

##------- Create DF for plot of mean proportion with CI ------

bt <- data.frame(sp=c("p","p","c","c","p","p","c","c"), yr=c(15,15,15,15,16,16,16,16),
                 con=c("cc","uv","cc","uv","cc","uv","cc","uv"),
                 s_y=c("p15","p15","c15","c15","p16","p16","c16","c16"),
                 id=c("p15cc","p15uv","c15cc","c15uv","p16cc","p16uv","c16cc","c16uv"),
                 w_ave=c(m1,m2,m3,m4,m5,m6,m7,m8),lci=c(ci1$bca[1,c(4)],ci2$bca[1,c(4)],
                 ci3$bca[1,c(4)],ci4$bca[1,c(4)],ci5$bca[1,c(4)],ci6$bca[1,c(4)],ci7$bca[1,
                 c(4)],ci8$bca[1,c(4)]), uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],
                 ci4$bca[1,c(5)],ci5$bca[1,c(5)],ci6$bca[1,c(5)],ci7$bca[1,c(5)],ci8$bca[1,c(5)]))
bt <- with(bt,bt[order(sp),])

##Calcuate standard deviation from CI
bt$lsd <- abs(bt$w_ave-bt$lci)
bt$usd <- abs(bt$w_ave-bt$uci)
bt$yr <-as.factor(bt$yr)
## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt$w_stnd <- ((bt$w_ave) - min(bt$w_ave))/diff(range(bt$w_ave))

#------- BW Plot of selection ratio w Bonferroni CI -------

x <- c(1:8)
rsr_c <- with(bt,
              {plot(x, w_ave,xaxt = "n", ylim=range(c(w_ave-lsd, w_ave+usd)),
     pch=c(16:17),bg=par(con), cex=1.2, col=c("black","black"),xlab="Species/Year", 
     ylab=expression("Selection ratio  "(w[i])),
     main="Calibrated constrained vs unconstrained channels")# std.dev error bars
#Draw arrows with length of sd, & horizontal bar with tip length of .1
    with(bt,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
    abline(1,0,lty=2,lwd=1,col="black")
    abline( v = c(4.5), col = "grey", lty = 1)
    xl <- c("C15","C16","P15","P16")
    axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
    legend("topright",c("CC","UV"),pch=c(16:17),bg=par(con), cex=1,inset=.01) 
              })
## Chum & pink salmon significantly avoided constrained channels in 2015 & 2016, while both selected 
## for unconstrained valley channels in 2015 at a higher proportion than available; and both species 
## selected UV in proportion to availability during 2016.
## Important to use constraint to isolate species difference in spawning habitat.
## Important to note that constraint is a better predictor than geology due to scale

## Low flow limit resource selection for most resource types


#****************************
##--- 2) Geology Data --------
#****************************
habcat <- read.csv("ip_habcat.csv")

print(table(habcat$geol))#geology groups
g <-habcat[,c(1,5)]

#Redd counts
redd <-read.csv("redd.csv")
redd <- redd[,-1]
fc <- redd[,c(1,13:19)]
str(fc)
fc <- merge(fc,g) #add contrained class
head(fc)

fc <- subset(fc,fc$geol!="Sc" & fc$geol!="Ss")#remove Sc & Ss geology types, too few
rm(g,redd)
print(table(fc$geol))

##----- Bootstrap --------
library(boot)
##Pink 2015 geol (Dv) - (Qs)
set.seed(167)
g_p15 <- function(data,i){
  d = data[i, ]
  o_Dv = (d$prn15[d$geol=="Dv"]/sum(d$prn15))# for each CC study reach
  p_Dv = (d$ha15[d$geol=="Dv"]/sum(d$ha15))
  o_Qs = (d$prn15[d$geol=="Qs"]/sum(d$prn15))# for each UV study reach
  p_Qs = (d$ha15[d$geol=="Qs"]/sum(d$ha15))
  wi_Dv = mean(o_Dv/p_Dv)
  wi_Qs = mean(o_Qs/p_Qs)
  wi_df = mean((o_Dv/p_Dv)-(o_Qs/p_Qs)) #difference between wi's
  wi = c(wi_Dv,wi_Qs,wi_df)
  return(wi)
}
p15 <- boot(data = fc, statistic =g_p15, strata=fc$geol,R = 3000)
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

##Chum 2015 
set.seed(167)
g_c15 <- function(data,i){
  d = data[i, ]
  o_Dv = (d$crn15[d$geol=="Dv"]/sum(d$crn15))# for each CC study reach
  p_Dv = (d$ha15[d$geol=="Dv"]/sum(d$ha15))
  o_Qs = (d$crn15[d$geol=="Qs"]/sum(d$crn15))# for each UV study reach
  p_Qs = (d$ha15[d$geol=="Qs"]/sum(d$ha15))
  wi_Dv = mean(o_Dv/p_Dv)
  wi_Qs = mean(o_Qs/p_Qs)
  wi_df = mean((o_Qs/p_Qs)-(o_Dv/p_Dv)) #difference between wi's
  wi = c(wi_Dv,wi_Qs,wi_df)
  return(wi)
}
c15 <- boot(data = fc, statistic =g_c15, strata=fc$geol,R = 3000)
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
g_p16 <- function(data,i){
  d = data[i, ]
  o_Dv = (d$prn16[d$geol=="Dv"]/sum(d$prn16))# for each CC study reach
  p_Dv = (d$ha16[d$geol=="Dv"]/sum(d$ha16))
  o_Qs = (d$prn16[d$geol=="Qs"]/sum(d$prn16))# for each UV study reach
  p_Qs = (d$ha16[d$geol=="Qs"]/sum(d$ha16))
  
  wi_Dv = mean(o_Dv/p_Dv)
  wi_Qs = mean(o_Qs/p_Qs)
  wi_df = mean((o_Qs/p_Qs)-(o_Dv/p_Dv)) #difference between wi's
  wi = c(wi_Dv,wi_Qs,wi_df)
  return(wi)
}
p16 <- boot(data = fc, statistic =g_p16, strata=fc$geol,R = 3000)
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
g_c16 <- function(data,i){
  d = data[i, ]
  o_Dv = (d$crn16[d$geol=="Dv"]/sum(d$crn16))# for each CC study reach
  p_Dv = (d$ha16[d$geol=="Dv"]/sum(d$ha16))
  o_Qs = (d$crn16[d$geol=="Qs"]/sum(d$crn16))# for each UV study reach
  p_Qs = (d$ha16[d$geol=="Qs"]/sum(d$ha16))
  wi_Dv = mean(o_Dv/p_Dv)
  wi_Qs = mean(o_Qs/p_Qs)
  wi_df = mean((o_Qs/p_Qs)-(o_Dv/p_Dv)) #difference between wi's
  wi = c(wi_Dv,wi_Qs,wi_df)
  return(wi)
}
(c16 <- boot(data = fc, statistic =g_c16, strata=fc$geol,R = 3000))
plot(c16,index=1)
plot(c16,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
(ci7<-boot.ci(c16, index=1,type="bca",conf = 0.975))
ci7$bca[1,c(4:5)]# BCa 97.5%CI, for CC
(ci8<-boot.ci(c16, index=2,type="bca",conf = 0.975))
ci8$bca[1,c(4:5)]# for UV
(m7 <- mean(c16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m8 <- mean(c16$t[,2]))#for UV

t.test(c16$t[,1],c16$t[,2])

##---- Create DF for plot of mean proportion with CI ----

bt2 <- data.frame(sp=c("p","p","c","c","p","p","c","c"), yr=c(15,15,15,15,16,16,16,16),
                  geol=c("Dv","Qs","Dv","Qs","Dv","Qs","Dv","Qs"),
                  s_y=c("p15","p15","c15","c15","p16","p16","c16","c16"),
                  id=c("p15Dv","p15Qs","c15Dv","c15Qs","p16Dv","p16Qs","c16Dv","c16Qs"),
                  w_ave=c(m1,m2,m3,m4,m5,m6,m7,m8),lci=c(ci1$bca[1,c(4)],ci2$bca[1,c(4)],
                 ci3$bca[1,c(4)],ci4$bca[1,c(4)],ci5$bca[1,c(4)],ci6$bca[1,c(4)],ci7$bca[1,
                 c(4)],ci8$bca[1,c(4)]), uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],                                                                                                                                                 ci4$bca[1,c(5)],ci5$bca[1,c(5)],ci6$bca[1,c(5)],ci7$bca[1,c(5)],ci8$bca[1,c(5)]))

bt2 <- with(bt2,bt2[order(sp),])

##Calcuate standard deviation from CI
bt2$lsd <- abs(bt2$w_ave-bt2$lci)
bt2$usd <- abs(bt2$w_ave-bt2$uci)
bt2$yr <-as.factor(bt2$yr)
## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt2$w_stnd <- ((bt2$w_ave) - min(bt2$w_ave))/diff(range(bt2$w_ave))

#------- Geology Plot of selection ratio w Bonferroni CI-------

rsr_g <- with(bt2,
     {x <- c(1:8)
      plot(x, w_ave,xaxt = "n",
              ylim=range(c(w_ave-lsd, w_ave+usd)),
              pch=c(16:17),bg=par(geol), cex=1.2,xlab="Species/Year", ylab=expression("Selection ratio "(w[i])), 
              main='Glacial alluvium (Qs) vs volcanic diorite (Dv)')
#Draw arrows with length of sd, & horizontal bar with tip length of .1
with(bt2,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
abline(1,0,lty=2,lwd=1,col="black")
abline(v = c(4.5), col = "grey", lty = 1)
xl <- c("C15","C16","P15","P16")
axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
with(bt2,legend("topright",c("Dv","Qs"),pch=c(16:17),bg=par(geol), cex=1,inset=.01))
})

## Dv - Karheen Formation volcanic and hypabyssal rocks; dolerite, similar to basalt, but cooled slower (Game Creek Canyon)
## Qs - Surficial deposits, undifferentiated (Alluvial deposits)

## Pinks significantly selected for Qs geology in 2015 & 2016, and avoided Dv geology during 2015.
## Chum significantly avoided Dv geology in 2016.
## Dv is the harder basalt, and the Qs is the glacial alluviam deposits

#***************************
##--- 3) Tributary Data ---------
## BootTrib - Selection ratio for mainstem vs tributary with bootstrap CI
## for mainstem vs trib redd counts
#habcat <- read.csv("ip_habcat.csv")
t <- habcat[,c(1,4)] ## Trib vs mainstem catagories
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

bt3 <- data.frame(sp=c("p","p","c","c","p","p","c","c"), yr=c(15,15,15,15,16,16,16,16),
                  trib=c("t","m","t","m","t","m","t","m"),
                  s_y=c("p15","p15","c15","c15","p16","p16","c16","c16"),
                  id=c("p15t","p15m","p16t","p16m","c15t","c15m","c16t","c16m"),
                  w_ave=c(m1,m2,m3,m4,m5,m6,m7,m8),lci=c(ci1$bca[1,c(4)],ci2$bca[1,c(4)],
                 ci3$bca[1,c(4)],ci4$bca[1,c(4)],ci5$bca[1,c(4)],ci6$bca[1,c(4)],ci7$bca[1,
                 c(4)],ci8$bca[1,c(4)]), uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],
                 ci4$bca[1,c(5)],ci5$bca[1,c(5)],ci6$bca[1,c(5)],ci7$bca[1,c(5)],ci8$bca[1,c(5)]))

bt3 <- with(bt3,bt3[order(sp),])

##Calcuate standard deviation from CI
bt3$lsd <- abs(bt3$w_ave-bt3$lci)
bt3$usd <- abs(bt3$w_ave-bt3$uci)
bt3$yr <-as.factor(bt3$yr)

## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt3$w_stnd <- ((bt3$w_ave) - min(bt3$w_ave))/diff(range(bt3$w_ave))

## ------ BW Plot of selection ratio w Bonferroni CI ------


rsr_t <- with(bt3,
      {x <- c(1:8)
      plot(x, w_ave,xaxt = "n",
      ylim=range(c(w_ave-lsd, w_ave+usd)),
      pch=c(16:17),bg=par(trib), cex=1.2,xlab="Species/Year", ylab=expression("Selection ratio  "(w[i])),
      main="Tributary vs Mainstem Channels")# with standard error bars
      #Draw arrows with length of sd, & horizontal bar with tip length of .1
      with(bt3,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
      abline(1,0,lty=2,lwd=1,col="black")
      abline(v = c(4.5), col = "grey", lty = 1)
      xl <- c("C15","C16","P15","P16")
      axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
      with(bt3,legend("topright",c("Tributary","Mainstem"),pch=c(16:17),bg=par(trib), cex=1,inset=.01))
})

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


#*********************************
##--- 4) Habitat Size Data ---------
#*********************************
## Resourse Selection Ratio - bootstrapped channel size & gradient (small/large)
## **** Threshold based on PCA exploratory analysis
#habcat <- read.csv("ip_habcat.csv")

print(table(habcat$size))#habitat size groups
s <-habcat[,c(1,7)]
hnl2 <- read.csv("D:/R/IP16/hnl2.csv", header=T)
h <- hnl2[,c(7,11,15)]
rm(hnl2)
hs <- cbind(h,s)
## Habitat small high gradient/large low gradient threshold: based on PCA PC1 
## Use min flow and max gradient from large, or max flow and min gradient for small as threshold
aggregate(cbind(GRADIENT,MEANANNCMS)~size,mean, data=hs)
aggregate(cbind(GRADIENT,MEANANNCMS)~size,median, data=hs)
aggregate(cbind(GRADIENT,MEANANNCMS)~size,max, data=hs)# 1% gradient, 1.9 cms
aggregate(cbind(GRADIENT,MEANANNCMS)~size,min, data=hs)
median(with(hs,GRADIENT))
median(with(hs,MEANANNCMS))

#Redd counts
redd <-read.csv("redd.csv")
redd <- redd[,-1]
rc <- redd[,c(1,13:19)]#redd counts & habitat
str(rc)
rc <- merge(rc,s) #add contrained class
head(rc)


##--- Bootstrap --------
library(boot)
##Pink 2015 size (s) - (l)
set.seed(167)
s_p15 <- function(data,i){
  d = data[i, ]
  o_s = (d$prn15[d$size=="s"]/sum(d$prn15))# 
  p_s = (d$ha15[d$size=="s"]/sum(d$ha15))
  o_m = (d$prn15[d$size=="m"]/sum(d$prn15))# 
  p_m = (d$ha15[d$size=="m"]/sum(d$ha15))  
  o_l = (d$prn15[d$size=="l"]/sum(d$prn15))# 
  p_l = (d$ha15[d$size=="l"]/sum(d$ha15))
  wi_s = mean(o_s/p_s)
  wi_m = mean(o_m/p_m)
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_s/p_s)-(o_l/p_l)) #difference between wi's
  wi = c(wi_s,wi_m,wi_l,wi_df)
  return(wi)
}
p15 <- boot(data = rc, statistic =s_p15, strata=rc$size,R = 3000)
p15
plot(p15,index=1)# small 
plot(p15,index=2)# medium
plot(p15,index=3)# Large
#plot(p15,index=4)# difference between small-large
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci1<-boot.ci(p15, index=1,type='bca',conf = 0.983) 
ci2<-boot.ci(p15, index=2,type='bca',conf = 0.983)
ci3 <- boot.ci(p15, index=3,type='bca',conf = 0.983)
ci1$bca[1,c(4:5)]# BCa 97.5%CI, small
ci2$bca[1,c(4:5)]# medium
ci3$bca[1,c(4:5)]# large
(m1 <- mean(p15$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m2 <- mean(p15$t[,2]))#for med
(m2 <- mean(p15$t[,3]))#for large


##Chum 2015 
set.seed(167)
s_c15 <- function(data,i){
  d = data[i, ]
  o_s = (d$crn15[d$size=="s"]/sum(d$crn15))# 
  p_s = (d$ha15[d$size=="s"]/sum(d$ha15))
  o_m = (d$crn15[d$size=="m"]/sum(d$crn15))# 
  p_m = (d$ha15[d$size=="m"]/sum(d$ha15))
  o_l = (d$crn15[d$size=="l"]/sum(d$crn15))# 
  p_l = (d$ha15[d$size=="l"]/sum(d$ha15))
  wi_s = mean(o_s/p_s)
  wi_m = mean(o_m/p_m)  
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_l/p_l)-(o_s/p_s)) #difference between wi's
  wi = c(wi_s,wi_m,wi_l,wi_df)
  return(wi)
}
c15 <- boot(data = rc, statistic =s_c15, strata=rc$size,R = 3000)
c15
plot(c15,index=1)# Small
plot(c15,index=2)# med
plot(c15,index=3)# large
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci4<-boot.ci(c15, index=1,type="bca",conf = 0.983) 
ci5<-boot.ci(c15, index=2,type="bca",conf = 0.983)
ci6<-boot.ci(c15, index=3,type="bca",conf = 0.983)
ci4$bca[1,c(4:5)]# BCa 97.5%CI, for small
ci5$bca[1,c(4:5)]# for med
ci5$bca[1,c(4:5)]# large
(m4 <- mean(c15$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m5 <- mean(c15$t[,2]))#for med
(m6 <- mean(c15$t[,3]))#for large


##Pink 2016 t/m
set.seed(167)
s_p16 <- function(data,i){
  d = data[i, ]
  o_s = (d$prn16[d$size=="s"]/sum(d$prn16))# 
  p_s = (d$ha16[d$size=="s"]/sum(d$ha16))
  o_m = (d$prn16[d$size=="m"]/sum(d$prn16))# 
  p_m = (d$ha16[d$size=="m"]/sum(d$ha16))  
  o_l = (d$prn16[d$size=="l"]/sum(d$prn16))# 
  p_l = (d$ha16[d$size=="l"]/sum(d$ha16))
  
  wi_s = mean(o_s/p_s)
  wi_m = mean(o_m/p_m)
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_l/p_l)-(o_s/p_s)) #difference between wi's
  wi = c(wi_s,wi_m,wi_l,wi_df)
  return(wi)
}
p16 <- boot(data = rc, statistic =s_p16, strata=rc$size,R = 3000)
p16
plot(p16,index=1)
plot(p16,index=2)
plot(p16,index=3)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
ci7<-boot.ci(p16, index=1,type="bca",conf = 0.983) 
ci8<-boot.ci(p16, index=2,type="bca",conf = 0.983)
ci9<-boot.ci(p16, index=3,type="bca",conf = 0.983)
ci7$bca[1,c(4:5)]# BCa 97.5%CI
ci8$bca[1,c(4:5)]# 
ci9$bca[1,c(4:5)]# 
(m7 <- mean(p16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m8 <- mean(p16$t[,2]))# for med
(m9 <- mean(p16$t[,2]))# large


##Chum 2016 t/m
set.seed(167)
s_c16 <- function(data,i){
  d = data[i, ]
  o_s = (d$crn16[d$size=="s"]/sum(d$crn16))#
  p_s = (d$ha16[d$size=="s"]/sum(d$ha16))
  o_m = (d$crn16[d$size=="m"]/sum(d$crn16))#
  p_m = (d$ha16[d$size=="m"]/sum(d$ha16))
  o_l = (d$crn16[d$size=="l"]/sum(d$crn16))# for each UV study reach
  p_l = (d$ha16[d$size=="l"]/sum(d$ha16))
  wi_s = mean(o_s/p_s)
  wi_m = mean(o_m/p_m)  
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_l/p_l)-(o_s/p_s)) #difference between wi's
  wi = c(wi_s,wi_m,wi_l,wi_df)
  return(wi)
}
(c16 <- boot(data = rc, statistic =s_c16, strata=rc$size,R = 3000))
plot(c16,index=1)
plot(c16,index=2)
# Bonferroni adj. Alpha=0.1/(2I), I=#Groups,& BCa is corrected CI preffered by statisticians (Manley pg.58)
(ci10<-boot.ci(c16, index=1,type="bca",conf = 0.983))
ci10$bca[1,c(4:5)]# BCa 97.5%CI
(ci11<-boot.ci(c16, index=2,type="bca",conf = 0.983))
ci11$bca[1,c(4:5)]# med
(ci12<-boot.ci(c16, index=3,type="bca",conf = 0.983))
ci12$bca[1,c(4:5)]# large

(m10 <- mean(c16$t[,1]))#Bootstrap mean corrected for bias, report this mean
(m11 <- mean(c16$t[,2]))#for med
(m12 <- mean(c16$t[,3]))#for med

##---- Create DF for plot of mean proportion with CI ----

bt4 <- data.frame(sp=c("p","p","p","c","c","c","p","p","p","c","c","c"), yr=c(15,15,15,15,15,15,16,16,16,16,16,16),
       size=c("s","m","l","s","m","l","s","m","l","s","m","l"),
       s_y=c("p15","p15","p15","c15","c15","c15","p16","p16","p16","c16","c16","c16"),
       id=c("p15s","p16m","p15l","c15s","c15m","c15l","p16s","p16m","p16l","c16s","c16m","c16l"),
       w_ave=c(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12),lci=c(ci1$bca[1,c(4)],ci2$bca[1,c(4)],ci3$bca[1,c(4)],
    ci4$bca[1,c(4)],ci5$bca[1,c(4)],ci6$bca[1,c(4)],ci7$bca[1,c(4)],ci8$bca[1,c(4)],ci9$bca[1,c(4)],
    ci10$bca[1,c(4)], ci11$bca[1,c(4)], ci12$bca[1,c(4)]),
    uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],ci4$bca[1,c(5)],ci5$bca[1,c(5)], ci6$bca[1,c(5)],ci7$bca[1,c(5)],
    ci8$bca[1,c(5)],ci9$bca[1,c(5)],ci10$bca[1,c(5)],ci11$bca[1,c(5)],ci12$bca[1,c(5)]))

(bt4 <- with(bt4,bt4[order(sp),]))

##Calcuate standard deviation from CI
bt4$lsd <- abs(bt4$w_ave-bt4$lci)
bt4$usd <- abs(bt4$w_ave-bt4$uci)
bt4$yr <-as.factor(bt4$yr)
## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt4$w_stnd <- ((bt4$w_ave) - min(bt4$w_ave))/diff(range(bt4$w_ave))

#------- size Plot of selection ratio w Bonferroni CI-------
## Habitat size based on original basin area GRTS selected based on even 1/3 frequency of all reaches

rsr_hs <- with(bt4,
        {x <- c(1:12)
          plot(x, w_ave,xaxt = "n",ylim=range(c(w_ave-lsd, w_ave+usd)), pch=c(16:17,15),bg=par(size), 
               cex=1.2,xlab="Species/Year", ylab=expression("Selection ratio "(w[i])), 
               main='Small, medium, & large habitat')
          #Draw arrows with length of sd, & horizontal bar with tip length of .1
          with(bt4,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
          abline(1,0,lty=2,lwd=1,col="black")
          abline(v = c(6.5), col = "grey", lty = 1)
          xl <- c("C15","C16","P15","P16")
          axis(1, at=c(2,5,8,11), lwd=2,col.axis="Blue",labels=xl)
         # axis(1, at=c(3.5,6.5,9.5), lwd=2,tcl=.5,col.axis="red",xlab = "")
          with(bt4,legend("topright",c("Small","Medium","Large"),pch=c(16:17,15),bg=par(size), cex=1,inset=.01))
})


## During high average flows in 2015 both chum and pink select all size habitat in proportion to availability, with pink selecting
## small habitat at a slightly higher propotion than chum.
## During low average spawning flows in 2016 chum and pink significantly avoided small habitat and selected
## medium and large in proportion to availability.
## The large variance in small habitat for 2015 is due to 

##--- 5) Multi-Plots --------

par(mfrow=c(2,2))
matrix(2,2,2,2)
plot(rsr_c,rsr_g,rsr_t,rsr_hs)

