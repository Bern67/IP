## Selection ratio for Geology with bootsrap CI


###------ Data ---------
habcat <- read.csv("ip_habcat.csv")

print(table(habcat$geol))#geology groups
g <-habcat[,c(1,5)]
rm(habcat)

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
      c(4)],ci8$bca[1,c(4)]), uci=c(ci1$bca[1,c(5)],ci2$bca[1,c(5)],ci3$bca[1,c(5)],
      ci4$bca[1,c(5)],ci5$bca[1,c(5)],ci6$bca[1,c(5)],ci7$bca[1,c(5)],ci8$bca[1,c(5)]))

bt2 <- with(bt2,bt2[order(sp),])

##Calcuate standard deviation from CI
bt2$lsd <- abs(bt2$w_ave-bt2$lci)
bt2$usd <- abs(bt2$w_ave-bt2$uci)
bt2$yr <-as.factor(bt2$yr)
## Calcualte standardized index for comparison, pg 53,table 4.1 Manly
bt2$w_stnd <- ((bt2$w_ave) - min(bt2$w_ave))/diff(range(bt2$w_ave))

#------- Geology Plot of selection ratio w Bonferroni CI-------

x <- c(1:8)
with(bt2,plot(x, w_ave,xaxt = "n",
     ylim=range(c(w_ave-lsd, w_ave+usd)),
     pch=c(16:17),bg=par(geol), cex=1.2,xlab="Species/Year", ylab=expression("Selection ratio "(w[i])), 
     main='Glacial alluvium (Qs) vs volcanic diorite (Dv)'
))
#Draw arrows with length of sd, & horizontal bar with tip length of .1
with(bt2,arrows(x, w_ave-lsd, x, w_ave+usd, length=0.1, angle=90, code=3))
abline(1,0,lty=2,lwd=1,col="black")
abline(v = c(4.5), col = "grey", lty = 1)
xl <- c("C15","C16","P15","P16")
axis(1, at=c(1.5,3.5,5.5,7.5), labels=xl)
with(bt2,legend("topright",c("Dv","Qs"),pch=c(16:17),bg=par(geol), cex=1,inset=.01))

## Dv - Karheen Formation volcanic and hypabyssal rocks; dolerite, similar to basalt, but cooled slower (Game Creek Canyon)
## Qs - Surficial deposits, undifferentiated (Alluvial deposits)

## Pinks significantly selected for Qs geology in 2015 & 2016, and avoided Dv geology during 2015.
## Chum significantly avoided Dv geology in 2016.
## Dv is the harder basalt, and the Qs is the glacial alluviam deposits

