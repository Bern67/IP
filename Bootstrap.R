## Selection Ratio by confinement (Manly 2002) w Bootstrap CI
## Reference Btutorial.pdf presentation for example; located in Zotera
## Create channel morphology group analysis

## Need to run one-way ANOVA to on ranked data to test sig dif between species and channel unit features among years

##------ VWI Wi Data ------
hnl2 <- read.csv("D:/R/IP16/hnl2.csv", header=T)
vc <- hnl2[,c('ip', 'ValCnstrnt','VWI_Floor')]

## Channel confinement groups (Grant 1995,Moore 2002), based on VWI
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
rsr_c <- with(bt,plot(x, w_ave,xaxt = "n", ylim=range(c(w_ave-lsd, w_ave+usd)),
     pch=c(16:17),bg=par(con), cex=1.2, col=c("black","black"),xlab="Species/Year", 
     ylab=expression("Selection ratio  "(w[i])),
     main="Calibrated constrained vs unconstrained channels")# std.dev error bars
)
#Draw arrows with length of sd, & horizontal bar with tip length of .1
with(bt,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
abline(1,0,lty=2,lwd=1,col="black")
abline( v = c(4.5), col = "grey", lty = 1)
xl <- c("C15","C16","P15","P16")
axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
with(bt,legend("topright",c("CC","UV"),pch=c(16:17),bg=par(con), cex=1,inset=.01)) 

## Chum & pink salmon significantly avoided constrained channels in 2015 & 2016, while both selected 
## for unconstrained valley channels in 2015 at a higher proportion than available; and both species 
## selected UV in proportion to availability during 2016.
## Important to use constraint to isolate species difference in spawning habitat.

## Low flow limit resource selection for most resource types