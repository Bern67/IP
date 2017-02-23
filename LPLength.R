### Bern Romey, DEM Channel length vs slope, 5/17/2016

####  (~:|:~) Global Settings (~:|:~)  #### 
op <- par(mfrow=c(1,1)) # set the default plot grid to 1X1
source("cor.matrix.r")
library(ggplot2)
dta <- read.csv("D:/R/IP16/field_demV1.CSV")

##### Data Analysis #####
grad <- dta[,c("f_slope","hnl2_slope")]
grad <- na.omit(grad)
#slope$f_slope <- slope$f_slope*100
#slope$hnl2_slope <- slope$hnl2_slope*100

cor.matrix(grad)
cor.matrix(log(grad+1))

# names(data)[names(data)=="s30m"] <- "s30ml" #rename header
dta$dif <- dta$hnli_lpl-dta$f_lpl # Length difference between LiDAR stream and field length
dta$s30dif <- dta$hnli_s30m-dta$f_lpl  # Length difference between LiDAR 30m smoothed length and field length
dta$prop <- (1-(dta$f_lpl/dta$hnli_lpl)) # LiDAR stream length proportion longer than field length
dta$pct <- dta$prop*100
dta$s30p <- 1-(dta$f_lpl/dta$hnli_s30m) # LiDAR 30m smoothed proportion longer than field length
dta$s30pct <- dta$s30p*100
str(dta)
head(dta)
summary(dta)


## Matrix of hnl2 & field 2015 variables, NA removed
#dt1 <-na.omit(dta)
dt1 <- na.omit(dta[,c(2:11)])
cor.matrix(dt1)
cor.matrix(log(dt1+1))

hist(with(dta,f_bfw))

## Create size catagories based on bfw histogram above
dta$size <- with(dta,ifelse(dta$f_bfw < 6, "s",
            ifelse(dta$f_bfw > 6 & dta$f_bfw <= 20, "m",
            ifelse(dta$f_bfw > 20, "l",""))))

library(dplyr)
sdt <- filter(dta, size == 's')
mdt <- filter(dta, size == 'm')
ldt <- filter(dta, size == 'l')

median(sdt$prop)
median(mdt$prop)
median(ldt$prop)
median(sdt$s30p) # negative because reduced to shorter length than field measured length
median(mdt$s30p)
median(ldt$s30p)

####  :~) Plots (~: ####
library(ggplot2)
## slope 
ps3 <- ggplot(log(grad+1), aes(f_slope, hnl2_slope)) +
  geom_point() +
#  lims(x = c(0, .2), y = c(0, .2)) +
  geom_smooth(method = "lm", se = TRUE) +
  geom_abline(intercept = 0, slope = 1, linetype = 4) +
  labs(list(title = "", x = "Ln(Field gradient)", y = "Ln(Model gradient)"))
ps3

##R-plot of gradient; same as ps3 plot
#*************************************
plot(log(grad$hnl2_slope+1)~log(grad$f_slope+1),ylab="Ln(Model gradient)",xlab='Ln(Field gradient)',
     ylim=c(0,.2),xlim=c(0,.2))
abline(0,1,lty=2, col="grey60")
(smod <- lm(hnl2_slope~f_slope,data=log(grad+1)))
summary(log(grad+1))
clip(min(x=0.001568), max(x=0.156303), min(y=0.001024), max(y=0.175825))
abline(smod, lwd=2.5,col='darkolivegreen')
## rounded coefficients for better output
cf <- round(coef(smod), 2) 
cd <- round(summary(smod)$r.squared,2)#coefficient of determination
## sign check to avoid having plus followed by minus for negative coefficients
eq <- paste0("Ln(model gradient) = ", cf[1],
             ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), " Ln(Field gradient)")
## printing of the equation
mtext(eq, side=3, outer =F,line=-2,cex=.85)#side: (1=bottom, 2=left, 3=top, 4=right)
mtext(as.expression(substitute(R^2==cd,list(cd=cd))),3,-3.5,cex=.85)


## Boxplot of proportional difference
p1 <- ggplot(dta, aes(size, prop)) +
    geom_boxplot() +  geom_jitter(width = 0.2) +
    geom_hline(yintercept = 0, colour = 'grey') +
    labs(list(title = "", x = "Field Measured Channel Width (m)", y = "LiDAR Length Increase (Proportion)"))
p1

## Boxplot of Smooth 30m proportional difference
s30p <- ggplot(dta, aes(size, s30p)) +
    geom_boxplot() +  geom_jitter(width = 0.2) +
    geom_hline(yintercept = 0, colour = 'grey') +
    labs(list(title = "", x = "Field Measured Channel Width (m)", y = "LiDAR 30m Smooth Length Increase (Proportion)"))
s30p

## Field channel length increase vs field bfw (by size)
p2 <- ggplot(dta, aes(f_bfw, dif, colour=size)) +
    geom_point() +
    geom_smooth(se = TRUE, method = "lm") +
    labs(list(title = "", x = "Field Measured Channel Width (m)", y = "LiDAR Stream Length Increase (m)"))
p2
## Synthetic channel lenght increase vs field bfw
p3 <- ggplot(dta, aes(f_bfw, dif)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    labs(list(title = "", x = "Field Measured Channel Width (m)", y = "LiDAR Stream Length Increase (m)"))
p3

# hnl2: DEM length vs field length
p4 <- ggplot(dta, aes(f_lpl, hnl2_lpl, label = dta$size)) +
  geom_point() +
  lims(x = c(0, 1200), y = c(0, 1200)) +
  geom_smooth(method = "lm", se = T, colour ='darkgreen') +
  geom_abline(intercept = 0, slope = 1, linetype = 4) +
#  geom_point() + geom_text(hjust = 0, nudge_x = 0.0015) +
  annotate("text", x = 750, y = 500, label = "y = 1.19x - 20.74") +
  annotate("text", x = 750, y = 450, label = "R-sqr. = 0.988") +
  labs(list(title = "", x = "Field measured channel length (m)", y = "Hnl2 Model channel length (m)"))
p4

## hnli: DEM length vs field length
p5 <- ggplot(dta, aes(f_lpl, hnli_lpl)) +
    geom_point() +
    lims(x = c(0, 1200), y = c(0, 1200)) +
    geom_smooth(method = "lm", se = T, colour ='darkgreen') +
    geom_abline(intercept = 0, slope = 1, linetype = 4) +
    annotate("text", x = 750, y = 500, label = "y = 1.32x - 24.36") +
    annotate("text", x = 750, y = 450, label = "R-sqr. = 0.987") +
    labs(list(title = "", x = "Field Measured Channel length (m)", y = "Hnli DEM Channel Length (m)"))
p5

s30 <- ggplot(dta, aes(f_lpl, hnli_s30m)) +
    geom_point() +
    lims(x = c(0, 950), y = c(0, 950)) +
    geom_smooth(method = "lm", se = T) +
    geom_abline(intercept = 0, slope = 1, linetype = 4) +
    labs(list(title = "", x = "Field Measured Channel length (m)", y = "LiDAR 30m Smooth Stream Length (m)"))
s30


#### Field measured mid channel LP with range finder on July 11-15, 2016 ####
require(ggplot2)
lp <- na.omit(data[ ,-c(12,13)])
## Matrix of hnl2 & 2016 variables (saved as m2.png)
cor.matrix(lp[,c(3,4,8,9,12:14)])

#write.csv(lp,"LP.csv") #Save data used for analysis as LP.csv
# lp$dif <-abs(lp$f2_lpl-lp$f_lpl)
lp$f2dif <-lp$f2_lpl-lp$f2_lpl_map # difference between f2 field measured and GPS ditance between points
lp$f2_lpldif <- lp$hnl2_lpl-lp$f2_lpl 
lp$f_lpldif <-lp$hnl2_lpl-lp$f_lpl
lp$f2_mdd <- lp$hnl2_lpl-lp$f2_lpl_map

lpl <- lp[ , c(1,3,12,13)]
colnames(lpl)
colnames(lpl) <-c("xs", "HC_lpl","RF_lpl","GPS_lpl")

#### Percent predicted hnl2 stream length is longer than field length
x1 <- sum(lpl$HC_lpl) #total length of streams mesured with hipchain
x2 <- sum(lpl$RF_lpl) #total length with rangefinder
x3 <- sum(lpl$GPS_lpl) #total length with GPS
x4 <- sum(lp$hnl2_lpl) # total synthetic length

x5 <- mean(lpl$HC_lpl) #meanlength with hipchain
x6 <- mean(lpl$RF_lpl) #mean length with rangefinder
x7 <- mean(lpl$GPS_lpl) #mean length with GPS
x8 <- mean(lp$hnl2_lpl) #mean synthetic length

x9 <- x4-x1 # hnl2 difference from hipchain
x10 <- x4-x2
x11 <- x4-x3
x12 <- x4-x4

x13 <- (1-(x1/x4))*100 # Percent longer hnl2 channel length
x14 <- (1-(x2/x4))*100
x15 <- (1-(x3/x4))*100
x16 <- (1-(x4/x4))*100

lp_tbl <- matrix(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16),
                 nrow = 4, ncol = 4, byrow = T)

dimnames(lp_tbl) = list(c("Channel length: total(m)", "Channel length: mean(m)","Hnl2-field: difference (m)", "Hnl2 % longer"),         # row names 
                        c("Hipchain", "Rangefinder", "GPS", "Hnl2")) # column names 
lp_tbl <- as.data.frame(lp_tbl)
# colnames(lp_tbl)[1] <- "HC" #Change header for a specific column
write.csv(lp_tbl,"lp_tbl.csv")


f2dif <-tidyr::gather(lpl,"group", "dif",2:4)
bp <- ggplot(f2dif, aes(group,dif))+
geom_boxplot(fill = "darkgrey", colour = "#006400") + geom_jitter(width = 0.1)+
#  geom_boxplot() + geom_jitter(aes(colour = f2dif$xs))+
labs(list(title = "", x = "Measurement Type", y = "Channel Length (m)"))
png("D:/R/IP16/Figs/LPbox.png", width=5, height=5, units="in", res=300) #Save plot to folder
plot(bp) #plot
dev.off() #turn off the dev
bp # plot name


summary(lpl)
ln_lpl <- log(lpl[,c(2:4)]) #normalize data
boxplot(ln_lpl)

 ## Is there a sig diff between  hipchain, rangefinder, and GPS map distance from hnl2 distance?
# One Way Anova (Completely Randomized Design)
fit <- aov(dif ~ group, data=f2dif) 
summary(fit) # P-value is > .05, not significant different
par(mfrow=c(2,2))
plot(fit)

f2 <-ggplot(lp, aes(f2_lpl, f_lpl,label = lp$xs)) +
  geom_point() +
#   lims(x = c(0, 1200), y = c(0, 1200)) +
  geom_point() + geom_text(hjust = 0, nudge_x = 2.5, size=3.5) +
  geom_smooth(method = "lm", se = T, colour ='darkgreen') +
  geom_abline(intercept = 0, slope = 1, linetype = 4) +
  annotate("text", x = 750, y = 550, label = "y = 0.96x + 21.2") +
  annotate("text", x = 750, y = 520, label = "R-sqr. = 0.980") +
  labs(list(title = "", x = "Rangefinder Channel length (m)", y = "Hipchain Channel Length (m)"))
png("D:/R/IP16/Figs/LPcal.png", width=5, height=5, units="in", res=300) #Save plot to folder
plot(f2)
dev.off()
f2

# Longitudinal profile with rangefinder, July 11-15
mod2 <- lm(lp$f_lpl~lp$f2_lpl)
summary(mod2)
plot(mod2)

lp1 <-ggplot(lp, aes(f2_lpl, hnl2_lpl,label = lp$xs)) +
  geom_point() +
  #   lims(x = c(0, 1200), y = c(0, 1200)) +
#  geom_point() + geom_text(hjust = 0, nudge_x = 2.5, size=3.5) +
  geom_smooth(method = "lm", se = T, colour ='darkgreen') +
  geom_abline(intercept = 0, slope = 1, linetype = 4) +
  annotate("text", x = 750, y = 550, label = "y = 0.80x + 29.5") +
  annotate("text", x = 750, y = 520, label = "R-sqr. = 0.991") +
  labs(list(title = "", x = "Rangefinder Channel length (m)", y = "Predicted (hnl2) Channel Length (m)"))
png("D:/R/IP16/Figs/lplot.png", width=5, height=5, units="in", res=300) #Save plot to folder
plot(lp1)
dev.off()
lp1

mod_lp <- lm(lp$f2_lpl~lp$hnl2_lpl)
summary(mod_lp)
par(mfro=c(2,2))
plot(mod_lp)
par(op)


# Plot of length difference between rangefinder and hnl2 compared with field BFW
f2a <-ggplot(lp, aes(f2_bfw,f2_lpldif,label = lp$xs)) +
  geom_point() +
  #   lims(x = c(0, 1200), y = c(0, 1200)) +
#  geom_point() + geom_text(hjust = 0, nudge_x = 0.0015) +
  geom_smooth(method = "lm", se = T, colour ='darkgreen') +
#  geom_abline(intercept = 0, slope = 1, linetype = 4) +
  annotate("text", x = 30, y = 65, label = "y = 0.11x + 12.9", col="black") +
  annotate("text", x = 30, y = 55, label = "R-sqr. = 0.672", col="black") +
  labs(list(title = "", y = "Length Difference (hnl2-field,m)", x = "Field Bank Full Width (m)"))
png("D:/R/IP16/Figs/LPcal2.png", width=5, height=5, units="in", res=300) #Save plot to folder
plot(f2a)
dev.off()
f2a

## Is the length error correlated with BFW?
mod2a <- lm(lp$f2_bfw~lp$f2_lpldif)
summary(mod2a)
par(mfrow=c(2,2))
plot(mod2a)


#### BFW comparison ####
#***********************

dt <- read.csv("D:/R/IP16/field_demV1.CSV")#Xs data
bw <- dt[,c(5,10,23)]
rm(dt)
bw1 <- na.omit(bw[,c(2,3)]) # field 2 and hnl2
bw2 <- na.omit(bw[,c(1,2)]) # field & hnl2
## Field2 & hnl2 BFW: multiple width measurements per channel, mean
plot(with(bw1,hnl2_Xs_bfw~f2_bfw), xlim=c(10,40),ylim=c(10,40),
     xlab='Field 2 BFW (m)', ylab='hnl2 BFW (m)')
abline(0,1,lty=2, col="grey60")
(w_mod <- lm(hnl2_Xs_bfw~f2_bfw,data=bw1))
clip(min(x=min(bw1$f2_bfw)), max(x=max(bw1$f2_bfw)), min(y=min(bw1$hnl2_Xs_bfw)), max(y=max(bw1$hnl2_Xs_bfw)))
abline(w_mod, lwd=2.5,col='darkolivegreen')
summary(w_mod)
##Field & hnl2 BFW: Only one width per channel
plot(with(bw2,hnl2_Xs_bfw~f_bfw), xlim=c(0,45),ylim=c(0,45),
     xlab='Field BFW (m)', ylab='hnl2 BFW (m)')
abline(0,1,lty=2, col="grey60")
(w_mod2 <- lm(hnl2_Xs_bfw~f_bfw,data=bw2))
clip(min(x=min(bw2$f_bfw)), max(x=max(bw2$f_bfw)), min(y=min(bw2$hnl2_Xs_bfw)), max(y=max(bw2$hnl2_Xs_bfw)))
abline(w_mod2, lwd=2.5,col='darkolivegreen')
summary(w_mod2)
summary(w_mod2)$r.squared
par(mfrow=c(2,2))
plot(w_mod2)#Equal variance, but not linear

## BFW Nonlinear model for field Xs & hnl2
#*****************************************
bfw_mod <-nls(hnl2_Xs_bfw~a*f_bfw^b, bw2, start=list(a=.0362, b=.74))
(summary(bfw_mod))
plot(with(bw2,hnl2_Xs_bfw~f_bfw), xlab = (""), ylab = (""))
abline(0,1, lty=2, col='grey60')#1:1 line
par(new=T)
curve(4.33198*x^0.49750, xlim=range(bw2$f_bfw), ylim=range(bw2$hnl2_Xs_bfw),
      col="darkgreen",lwd= 2,xlab="Field BFW (m)", ylab="Hnl2 BFW (m)")
## R-sqr
(RSS.p <- sum(residuals(bfw_mod)^2))  # Residual sum of squares
(TSS <- sum((bw2$hnl2_Xs_bfw - mean(bw2$hnl2_Xs_bfw))^2))  # Total sum of squares
(cd <- round(1 - (RSS.p/TSS),2))  # Coefficient of determination;R-squared
(cf <- round(coef(bfw_mod), 2)) #Model coefficient (a,b)
## sign check to avoid having plus followed by minus for negative coefficients
eq <- paste0("DEM channel width (m) = ", cf[1],'* Field channel width (m)^', cf[2])
## printing of the equation
mtext(eq, side=3, line=-2,cex=.85)#side: (1=bottom, 2=left, 3=top, 4=right)
mtext(as.expression(substitute(r^2==cd,list(cd=cd))),3,-3.5,cex=.85)


####  (~:|:~) LinModels (~:|:~)  ####
cor.matrix(dta[-c(39,48) ,-c(1,17)])       

library(MVN) # load MVN package
uniNorm(dta[-c(39,48) ,c(2:16)], type = "SW", desc = T) # summary statistics and Shipiro-Wilks normality test


mod <- lm(dta$hnl2_lpl~dta$f_lpl) # hnl2
summary(mod)
par(mfrow=c(2,2))
plot(mod)
par(op)

mod1 <- lm(dta$hnli_lpl~dta$f_lpl) # hnli
summary(mod1)
par(mfrow=c(2,2))
plot(mod1)
par(op)


### Assumption (~a~).  Random errors have a mean of zero, linearity.  Residuals should be normaly distributed about 
### the zero line for the residuals vs fitted (y-hat) plot.  
### Assumption (~b~).  Random error terms uncorrelated. If correlated, then multi-collinearity exists (two varaibles have same informatoin).  Use VIF test.
### Assumption (~c~). Random error terms have a constant variance (homoscedasticity). Scale-locatin plot should be uniform from 
### left to right.  Slight horn shape indicates heteroscedasticity.  
### Assumption (~d~).  Random error tersm follow a normal distribution.  Residuals Q-Q plot; all observatins fall on or near the line 
### Cooks distance shows what values are outliers and may bias the model.  0.5 is a concern, > 1 should be removed because it inflates the RSS

# paired t-test
t.test(dta$dif,dta$s30dif,paired=TRUE) # 

