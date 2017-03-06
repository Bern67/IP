## Field vs hnl2 analysis:
## Including VWI, BFW, MAQ, D50, BFD

##------ VWI Analysis ------
## Look at the correlation between field and hnl2 VC
hb <- read.csv("D:/R/IP16/field_demV1.CSV")
hab <- hb[,c(1,5,8:12,15,17,24)]#BFW, VW, and VWI
hab <- hab[c(-25,-43),]#site 63,GarT2 outliers for VWI
plot(hnl2_Xs_bfw~hnl2_bfw,data=hab)#hnl2 Xs vs reach BFW
plot(hnl2_bfw~f2_bfw,hab)
plot(hnl2_Xs_vw~hnl2_vw,data=hab) # hnl2 Xs vs reach VW
plot(hnl2_Xs_vw~opt_vw,data=hab)#Optical measure vs hnl2_ Xs VW
boxplot(hab)

## Plot of Field 2 BFW >12m & hnl2 BFW 
#****************************
f2bfw <- na.omit(hab[,c(8,10)])
bfw <- na.omit(hab[,c("hnl2_Xs_bfw",'f_bfw')])

plot(hnl2_bfw~f2_bfw, f2bfw, pch=17,col='grey17', xlim=range(bfw$f_bfw), 
     ylim=c(0,42), xlab='',ylab = '')
par(new=T)# Pay attention to the bfw measure method (what scale?)

## ---- BFW Nonlinear model for field Xs & hnl2 ----
#*****************************************
with(bfw,
     {(bfw_mod <- nls(hnl2_Xs_bfw~a*f_bfw^b, start=list(a=.0362, b=.74)))
       (summary(bfw_mod))
       plot(hnl2_Xs_bfw~f_bfw,pch=21,col='black', cex=1, bg='white', 
            xlab = (""), ylab = (""), xlim=range(bfw$f_bfw),ylim=c(0,42))
       abline(0,1, lty=2, col='grey60')#1:1 line
       par(new=T)
       curve(4.33198*x^0.49750, xlim=range(f_bfw), ylim=c(0,42),
             col="darkgreen",lwd= 2,xlab="Field BFW (m)", ylab="DEM model BFW (m)")
       grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
            lwd = par("lwd"), equilogs = TRUE)
       
       ## R-sqr
       (RSS.p <- sum(residuals(bfw_mod)^2))  # Residual sum of squares
       (TSS <- sum((hnl2_Xs_bfw - mean(hnl2_Xs_bfw))^2))  # Total sum of squares
       (cd <- round(1 - (RSS.p/TSS),2)) # Coefficient of determination;R-squared
       (cf <- round(coef(bfw_mod), 2)) # Model coefficient (a,b)
       ## sign check to avoid having plus followed by minus for negative coefficients
       eq <- paste0("DEM BFW (m) = ", cf[1],'* Field BFW (m)^', cf[2])
       ## printing of the equation
       mtext(eq, side=3, line=-2,cex=.85)#side: (1=bottom, 2=left, 3=top, 4=right)
       mtext(as.expression(substitute(R^2==cd,list(cd=cd))),3,-3.5,cex=.85)
       legend("topright",c("F15","F16"),pch=c(21,17), col=c("black","grey17"),inset=.01) })


## VW Data; opt_vw was optical estimate from DEM
hab$f_vwi <- hab$opt_vw/hab$f_bfw
hab$d_vwi <- hab$hnl2_Xs_vw/hab$hnl2_Xs_bfw
vw <- na.omit(hab[,c('xs','hnl2_Xs_vw','opt_vw','d_vwi','f_vwi')]) #opt_vw is filter for anadromous
# Only reaches with anadromous fish access used (n=35 -> 33 with outliers removed)

## ---- Valley width plot (hnl2 v optical opt) ----
#************************************
with(vw,
     {plot(hnl2_Xs_vw~opt_vw, xlim=range(opt_vw), ylim=c(0,424),
     xlab='Optical valley width (m)',ylab='Model valley width (m)')
abline(0,1, lty=2,col="grey60")#1:1 line
(vw_mod <- lm(hnl2_Xs_vw~opt_vw,data=vw))
summary(vw_mod)
clip(min(x=min(opt_vw)), max(x=max(opt_vw)), min(y=min(hnl2_Xs_vw)), max(y=max(hnl2_Xs_vw)))
abline(vw_mod, lwd=2.5,col='darkolivegreen')
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)

## rounded coefficients for better output
names(summary(vw_mod))
cf <- round(coef(vw_mod), 2) 
(cd <- round(summary(vw_mod)$r.squared,2))#coefficient of determination
## sign check to avoid having plus followed by minus for negative coefficients
eq <- paste0("Model valley width (m) = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " Optical valley width (m) ")
## printing of the equation
mtext(eq, side=3, line=-2,cex=.85)#side: (1=bottom, 2=left, 3=top, 4=right)
mtext(as.expression(substitute(R^2==cd,list(cd=cd))),3,-3.5,cex=.85)

opar <- par(mfrow=c(2,2), oma = c(0, 0, 1.1, 0))
(plot(vw_mod))
par(opar)
anova(vw_mod) })


## ---- Valley Width Index **(VWI)**; Constraint Plot; hnl2 Xs vs optical opt ----
#************************
with(vw,
     {plot(d_vwi~f_vwi, col='black', xlim=range(f_vwi), ylim=c(1,18.7),
     xlab="Optical field VWI",ylab="Model VWI")
abline(0,1, lty=2, col='grey60')#1:1 line
(vwi_mod <- lm(d_vwi~f_vwi,vw))
summary(vwi_mod)
summary(vwi_mod)$r.squared
clip(min(x=min(f_vwi)), max(x=max(f_vwi)), min(y=min(d_vwi) ), max(y=max(d_vwi)))
abline(vwi_mod ,lty=1,lwd=2,col="darkgreen")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
(cf1 <- round(coef(vwi_mod), 2))
(cd1 <- round(summary(vwi_mod)$r.squared,2))#coefficient of determination (r^2)
eq1 <- paste0("Model VWI = ", cf1[1],
              ifelse(sign(cf1[2])==1, " + ", " - "), abs(cf1[2]), " (Optical field VWI)")
## printing of the equation
mtext(eq1, side=3, line=-2,cex=.85)#side: (1=bottom, 2=left, 3=top, 4=right)
mtext(as.expression(substitute(R^2==cd1,list(cd1=cd1))),3,-3.25,cex=.85)

opar <- par(mfrow=c(2,2), oma = c(0, 0, 1.1, 0))
(plot(vwi_mod))
par(opar)
anova(vwi_mod) 
(vc <- round(coef(vwi_mod)[2]*2.9,2)+round(coef(vwi_mod)[1],2))# Used in bootstrap Wi for VC
})
## For HSI16 HSI VWI calibration: 2.9 = 3.09 (Grant 1995, Moore 2002, Hall 2007, pg. 789)
## Hall 2007 used 3.8 for 20-m DEM, we calibrated with 3.09 for LiDAR

## ---- Analysis & Summary Tables ----
## HNFP Ch & Pk anadromous stream data <= 5% gradient for entire project area
allhab <- read.csv('hnl2_anad.csv')
vwi <- allhab[,c('VWI_Floor','GRADIENT','MEANANNCMS','WIDTH_M')]
rm(allhab)
boxplot(vwi)
hist(vwi)
summary(vwi)

vwi$con <- as.factor(with(vwi,ifelse(vwi$VWI_Floor <= 3.09, "CC",#Constrained Canyon; from VWI analysis
                       ifelse(vwi$VWI_Floor > 3.09, "UV", ""))))#Unconstrained Valley
quantile(vwi$VWI_Floor)
tapply(vwi$VWI_Floor,vwi$con,median)# What is the median of VWI for CC & UV?
tapply(vwi$VWI_Floor,vwi$con,var)# Variance

aggregate(VWI_Floor~con,vwi,mean)
aggregate(VWI_Floor~con,vwi,median)

xtabs(~con,vwi)#how many CC & UV sites for all stream < 6% gradient?
round((prop.table(xtabs(~con,vwi))*100),2)#what percent of sites were CC & UV?
chisq.test(prop.table(xtabs(~con,vwi)))#Is the proporiton of CC sig. dif. from UV?
#warning message for small sample size



#**************************
## ---- Channel Depth (BFD) ----
bfd <- na.omit(hb[,c('f_bfd','hnl2_bfd')])
with(bfd,
     {plot(hnl2_bfd~f_bfd, col='black', xlim=range(bfd$f_bfd), ylim=c(0.05,.81),
     xlab="",ylab="")
abline(0,1, lty=2, col='grey60')#1:1 line
(bfd_mod <-nls(hnl2_bfd~a*f_bfd^b, bfd, start=list(a=.0362, b=.74)))
(summary(bfd_mod))
par(new=T)
curve(0.75909*x^0.27979, xlim=range(f_bfd), ylim=c(0.05,.81),
      col="darkgreen",lwd= 2,xlab="Field BFD (m)", ylab="DEM hnl2 BFD (m)")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
## R-sqr
(RSS.p <- sum(residuals(bfd_mod)^2))  # Residual sum of squares
(TSS <- sum((bfd$hnl2_bfd - mean(bfd$hnl2_bfd))^2))  # Total sum of squares
(cd <- round(1 - (RSS.p/TSS),2)) # Coefficient of determination;R-squared
(cf <- round(coef(bfd_mod), 2)) # Model coefficient (a,b)
eq <- paste0("DEM BFD (m) = ", cf[1],'* Field BFD (m)^', cf[2])
## printing of the equation
mtext(eq, side=3, line=-2,cex=.85)#side: (1=bottom, 2=left, 3=top, 4=right)
mtext(as.expression(substitute(R^2==cd,list(cd=cd))),3,-3.5,cex=.85)
})

#***************************************
## ---- D50 - Not a good correlation ----
d50 <- na.omit(hb[,c('f_d50','hnl2_d50')])
with(d50,
     {plot(hnl2_d50~f_d50, col='black', xlim=range(5,165), ylim=c(5,165),
     xlab="Field D50 (mm)",ylab="DEM D50 (mm)")
abline(0,1, lty=2, col='grey60')#1:1 line
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)})


#*********************
## ---- MAQ & BFW ----
hnl2 <- read.csv('D:/R/IP16/hnl2.csv')
q_w1 <- hnl2[,c("WIDTH_M",'MEANANNCMS')]
anad <- read.csv('hnl2_anad.csv')
q_w <- anad[,c('WIDTH_M','MEANANNCMS')]
rm(anad)
with(q_w,
     {plot(WIDTH_M~MEANANNCMS, col='grey', xlim=c(0,12), ylim=c(0,30),
     xlab="DEM Q (cms)",ylab="DEM BFW (cms)")
#abline(0,1, lty=2, col='grey60')#1:1 line
par(new=T)# redd dots are IP sites
plot(q_w1$WIDTH_M~q_w1$MEANANNCMS, q_w1,col='red', pch=16, cex=.75,xlim=c(0,12), ylim=c(0,30),
     xlab="DEM Q (cms)",ylab="DEM BFW (cms)")
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)
## Different curves most likely a reflection of the different hnl2 layers in NetMap
## Can run a non-linear regression on data to obtain model used, but there are multiple models})
})

