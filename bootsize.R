## Resourse Selection Ratio - bootstrapped channel size & gradient (small high grad/large low grad)

## **** do medium along with small and large
###------ Data ---------
habcat <- read.csv("ip_habcat.csv")

print(table(habcat$size))#habitat size groups
s <-habcat[,c(1,7)]
rm(habcat)
hnl2 <- read.csv("D:/R/IP16/hnl2.csv", header=T)
h <- hnl2[,c(7,11,15)]
rm(hnl2)
hs <- cbind(h,s)
## Habitat small high gradient/large low gradient threshold: based on PCA PC1 
aggregate(cbind(GRADIENT,MEANANNCMS)~size,mean, data=hs)
aggregate(cbind(GRADIENT,MEANANNCMS)~size,max, data=hs)# 1% gradient, 1.9 cms
aggregate(cbind(GRADIENT,MEANANNCMS)~size,min, data=hs)
mean(with(hs,GRADIENT))
mean(with(hs,MEANANNCMS))

#Redd counts
redd <-read.csv("redd.csv")
redd <- redd[,-1]
rc <- redd[,c(1,13:19)]#redd counts & habitat
str(rc)
rc <- merge(rc,s) #add contrained class
head(rc)


##----- Bootstrap --------
library(boot)
##Pink 2015 size (s) - (l)
set.seed(167)
g_p15 <- function(data,i){
  d = data[i, ]
  o_s = (d$prn15[d$size=="s"]/sum(d$prn15))# for each CC study reach
  p_s = (d$ha15[d$size=="s"]/sum(d$ha15))
  o_l = (d$prn15[d$size=="l"]/sum(d$prn15))# for each UV study reach
  p_l = (d$ha15[d$size=="l"]/sum(d$ha15))
  wi_s = mean(o_s/p_s)
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_s/p_s)-(o_l/p_l)) #difference between wi's
  wi = c(wi_s,wi_l,wi_df)
  return(wi)
}
p15 <- boot(data = rc, statistic =g_p15, strata=rc$size,R = 3000)
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
  o_s = (d$crn15[d$size=="s"]/sum(d$crn15))# for each CC study reach
  p_s = (d$ha15[d$size=="s"]/sum(d$ha15))
  o_l = (d$crn15[d$size=="l"]/sum(d$crn15))# for each UV study reach
  p_l = (d$ha15[d$size=="l"]/sum(d$ha15))
  wi_s = mean(o_s/p_s)
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_l/p_l)-(o_s/p_s)) #difference between wi's
  wi = c(wi_s,wi_l,wi_df)
  return(wi)
}
c15 <- boot(data = rc, statistic =g_c15, strata=rc$size,R = 3000)
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
  o_s = (d$prn16[d$size=="s"]/sum(d$prn16))# for each CC study reach
  p_s = (d$ha16[d$size=="s"]/sum(d$ha16))
  o_l = (d$prn16[d$size=="l"]/sum(d$prn16))# for each UV study reach
  p_l = (d$ha16[d$size=="l"]/sum(d$ha16))
  
  wi_s = mean(o_s/p_s)
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_l/p_l)-(o_s/p_s)) #difference between wi's
  wi = c(wi_s,wi_l,wi_df)
  return(wi)
}
p16 <- boot(data = rc, statistic =g_p16, strata=rc$size,R = 3000)
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
  o_s = (d$crn16[d$size=="s"]/sum(d$crn16))# for each CC study reach
  p_s = (d$ha16[d$size=="s"]/sum(d$ha16))
  o_l = (d$crn16[d$size=="l"]/sum(d$crn16))# for each UV study reach
  p_l = (d$ha16[d$size=="l"]/sum(d$ha16))
  wi_s = mean(o_s/p_s)
  wi_l = mean(o_l/p_l)
  wi_df = mean((o_l/p_l)-(o_s/p_s)) #difference between wi's
  wi = c(wi_s,wi_l,wi_df)
  return(wi)
}
(c16 <- boot(data = rc, statistic =g_c16, strata=rc$size,R = 3000))
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

bt4 <- data.frame(sp=c("p","p","c","c","p","p","c","c"), yr=c(15,15,15,15,16,16,16,16),
                  size=c("s","l","s","l","s","l","s","l"),
                  s_y=c("p15","p15","c15","c15","p16","p16","c16","c16"),
                  id=c("p15s","p15l","c15s","c15l","p16s","p16l","c16s","c16l"),
                  w_ave=c(m1,m2,m3,m4,m5,m6,m7,m8),lci=c(ci1$bca[1,c(4)],ci2$bca[1,c(4)],
                  ci3$bca[1,c(4)],ci4$bca[1,c(4)],ci5$bca[1,c(4)],ci6$bca[1,c(4)],ci7$bca[1,
                  c(4)],ci8$bca[1,c(4)]), uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],
                  ci4$bca[1,c(5)],ci5$bca[1,c(5)],ci6$bca[1,c(5)],ci7$bca[1,c(5)],ci8$bca[1,c(5)]))

bt4 <- with(bt4,bt4[order(sp),])

##Calcuate standard deviation from CI
bt4$lsd <- abs(bt4$w_ave-bt4$lci)
bt4$usd <- abs(bt4$w_ave-bt4$uci)
bt4$yr <-as.factor(bt4$yr)
## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt4$w_stnd <- ((bt4$w_ave) - min(bt4$w_ave))/diff(range(bt4$w_ave))

#------- size Plot of selection ratio w Bonferroni CI-------

x <- c(1:8)
with(bt4,plot(x, w_ave,xaxt = "n",
              ylim=range(c(w_ave-lsd, w_ave+usd)),
              pch=c(16:17),bg=par(size), cex=1.2,xlab="Species/Year", ylab=expression("Selection ratio "(w[i])), 
              main='Large low grad. (llg) vs Small high grad. (shg)'
))
#Draw arrows with length of sd, & horizontal bar with tip length of .1
with(bt4,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
abline(1,0,lty=2,lwd=1,col="black")
abline(v = c(4.5), col = "grey", lty = 1)
xl <- c("C15","C16","P15","P16")
axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
with(bt4,legend("topright",c("SHG","LLG"),pch=c(16:17),bg=par(size), cex=1,inset=.01))

## Chum salmon significantly selected for large low gradient habitat in 2015, and significantly avoided 
## small high gradient habitat in 2016.  Pink salmon selected both habitat types in proportion to availability
## during both years.  Although, the data sugests that pinks also avoided small high gradient habitat
## more during 2016 low water year than during 2015 high water year.
