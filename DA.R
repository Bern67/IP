## DFA for 2015 & 2016

## Assumptions of DFA:
## 1) Homoscedasticity - equal variance,
## 2) Multivariate normality,
## 3) absence of multicollinearity,
## 4) groups are predefined and no observation can be in two groups


#************************
## hnl2 habitat data ####
#************************
hnl2 <- read.csv("D:/R/IP16/hnl2.csv", header=T)
ip <- hnl2$ip

source("cor.matrix.r")#Linear combinations
cor.matrix(scale(sqrt(hnl2[,c("SRC_DIST","DEPTH_M","WIDTH_M","AREA_SQKM","MEANANNCMS")]))) #PC1
cor.matrix(hnl2[,c("StrmPow","FlowVel","Shear","d50", "GRADIENT")]) #PC1
cor.matrix(scale(sqrt(hnl2[,c("ValCnstrnt","VWI_Floor","VAL_WIDTH","FP_WIDTH")]))) #PC2
cor.matrix(hnl2[,c("GEP_DEL","GEP","BFQ","GEP_Cum")])#

h <-hnl2[,-c(2,5,8,11:12,14,17,19,21,25,26:28)] #Only variables with no multicollinearity; r<.07
#h <- hnl2[,-1]
rm(hnl2)

library(ResourceSelection)
kdepairs(scale(log(h[,-1]+1)))

#### ************************
## Redd density vs hnl2  ####  
#### ************************
redd <-read.csv("redd.csv")
redd <- redd[,-1]
r <- redd[,c(1,9,12)]

##Need to merge redd data with habitat for pink2015
library(Hmisc) ## Just needs to be signif correlated, reference Flitcroft 2014
####  Correlation matrix & also includes p-values matrix for significance; p<.05 signif.
#rcorr(as.matrix(p15_hnl2), type="spearman") #can use pearson if data is normal, or spearman if not
# If majority are not normal, use spearman.- subjective


#### **********************
## ALL DFA MODEL VARIABLES ####
#### **********************
head(redd) # all 49 study reaches
x <- as.factor(r[,"da15"]) #independent variable (factor); 2015
x1 <- as.factor(r[,"da16"])#independent variable (factor); 2016
y <- as.matrix(h[,-1]) #multivariate response variable matrix; hnl2 synthetic habitat


#### ***************************
## HNL2 MANOVA & DFA Assumption analysis ####
#### ***************************
##Run one-way MANOVA - Is there a significant difference between density groups
mod<-manova(y~x) # 2015
summary.manova(mod) #MANOVA table, compare it with univariate ANOVA table
summary.manova(mod, test = "Wilks")
mod1<-manova(y~x1) # 2016
summary.manova(mod1) #MANOVA table, compare it with univariate ANOVA table
summary.manova(mod1, test = "Wilks")

#### ***************************************
### Check for the two most important assumptions for MANOVA
##  ~~ 1. Multi-Normality assumption ~~
#### ***************************************
#Graphic assessment of multi-normality: chi-square plot written by Everitt
#it is similar to a q-q plot in univariate ANOVA
# residuals need to be normally distributed, the original data does not
chisplot <- function(y) {
  if (!is.matrix(y)) stop("x is not a matrix")
  ### determine dimensions
  n <- nrow(y)
  p <- ncol(y)
  #
  xbar <- apply(y, 2, mean)
  S <- var(y)
  S <- solve(S)
  index <- (1:n)/(n+1)
  #
  xcent <- t(t(y) - xbar)
  di <- apply(xcent, 1, function(y,S) y %*% S %*% y,S)
  #
  quant <- qchisq(index,p)
  plot(quant, sort(di), ylab = "Ordered distances",
       xlab = "Chi-square quantile", lwd=2,pch=1)
  return (di)
}
chisplot(y)
chisplot(resid(mod))# the dataset must be in matrix format
abline(a=0,b=1)

## Run a normality test 
library(mvnormtest)
mshapiro.test(t(y)) #data has to be a matrix and has to be transposed
## If p < .05, then sig dif than normal distribution

##Individual normality
par(mfrow=c(4,4))
names<-colnames(y)
for (i in 1:ncol(y)){
  qqnorm(y[,i],main=paste("Q-Q plot of",names[i]))
  qqline(y[,i])
}

##refit manova with sqrt-transformed response variables
(mod.sq<-manova(sqrt(y)~x)) #2015
summary.manova(mod.sq)
mod1.sq<-manova(sqrt(y)~x1) #2016
summary.manova(mod1.sq)
## They are still significantly different

##test the normality assumption again for multinormality
par(mfrow=c(2,1))
chisplot(resid(mod.sq))
abline(a=0,b=1)
mshapiro.test(t(resid(mod.sq))) # 2015
chisplot(resid(mod1.sq))
abline(a=0,b=1)
mshapiro.test(t(resid(mod1.sq))) # 2016 One site (12, row 8) need to be removed for normality
par(mfrow=c(1,1))

## Refit manova excluding 1 datapoint (8), for 2016
redd.48<-r[-c(8),]
x.48<-as.factor(redd.48[,"da16"]) # 
y.48<-as.matrix(h[-c(8),]) #
y.48 <- y.48[,-1]
mod.48<-manova(sqrt(y.48)~x.48)
summary.manova(mod.48)
##test the normality assumption again
chisplot(resid(mod.48))
abline(a=0,b=1)
mshapiro.test(t(resid(mod.48)))  #small size issue!

par(mfrow=c(2,1))
chisplot(resid(mod.sq)) #2015
abline(a=0,b=1)
chisplot(resid(mod.48)) #2016 with removed site 12, only site with pa
abline(a=0,b=1)

#### ***************************************
##  ~~ 2.Test homogeneity of covariance matrices ~~
#### ***************************************
## First look at boxplots for each variable
par(mfrow=c(3,4))
names<-colnames(y) # 2015
for (i in 1:ncol(y)){
  boxplot(y[,i]~x,main=paste("Boxplot of",names[i]),xlab="Groups", ylab=paste(names[i]))
}

#sqrt-transform
par(mfrow=c(3,4))
names<-colnames(y) #2015
for (i in 1:ncol(y)){
  boxplot(log(y[,i]+1)~x,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]), col="wheat3")
}

#sqrt-transform y.48 2016
par(mfrow=c(3,4))
names<-colnames(y.48)
for (i in 1:ncol(y.48)){
  boxplot(log(y.48[,i]+1)~x.48,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]), col="darkgreen")
}

#### ***************************************
## Each species group: c15 & 16, p15 & 16

## 2015 chum P/A,log-transform
par(mfrow=c(3,4))
names<-colnames(y)
for (i in 1:ncol(y)){
  boxplot(log(y[,i]+1)~redd$cpa15,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]), col="wheat3")
}## HSI: SRC_Dist,Grad, W,D, MAQ,Flow_vel, StrmPow,VW,GEP,d50

##2016 chum P/A,log-transform
par(mfrow=c(3,4))
names<-colnames(y)
for (i in 1:ncol(y)){
  boxplot(log(y[,i]+1)~redd$cpa16,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]), col="lightgrey")
}## HSI: SRC_Dist,Grad, W,D,MAQ,flow_vel,StrmPow, FPW,d50

##2015 pink H/L, log-transform
par(mfrow=c(3,4))
names<-colnames(y)
for (i in 1:ncol(y)){
  boxplot(log(y[,i]+1)~redd$ppa15,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]), col="wheat3")
}## MNANPRC,... Not a binomial distribution ->  Use linear regression, then compare predictors to 2016

##2016 pink p/a, log-transform
par(mfrow=c(3,4))
names<-colnames(y)
for (i in 1:ncol(y)){
  boxplot(log(y[,i]+1)~redd$ppa16,main=paste("Boxplot of log",names[i]),xlab="Groups", ylab=paste("log",names[i]), col="lightgrey")
}## Elev, src_dist, grad, fitElev, W,D, MAQ,flowVel,d50, FPW-marginal

#### ***************************************
## 2. Test homogeneity of covariance matrices; WITHIN GROUP EQUAL VARIANCE
#### ***************************************

## The test is developed by Anderson and is similar to Levene's univariate test for equal variance
## ***Intervals that do not overlap the vertical dashed line are significantly different***

#2015 & 2016 redd density for hnl2 habitat.

library(vegan) #2016 redd density groups - only one pa, remove
par(mfrow=c(1,2))
m<-betadisper(d=dist(y),group=redd$da16, type="centroid") 
m
m.HSD<-TukeyHSD(m)  
plot(m.HSD) 
#repeat the same test with square root tranformed 2016 data
m.sq<-betadisper(dist(sqrt(y)),redd$da16, type="centroid") 
m.sq
m.sq.HSD<-TukeyHSD(m.sq) 
plot(m.sq.HSD, sub="sqrt 2016")

library(vegan) #2015 redd density groups
par(mfrow=c(2,2))
m1<-betadisper(dist(y),redd$da15, type="centroid") #calculate average distance to its group centroid for each group
m1.HSD<-TukeyHSD(m1) #run TukeyHSD pair-wise test difference in dispersion between groups 
plot(m1.HSD) #graphic display of pair-wise comparisons of difference in dispersion between groups
#repeat the same test with sqrt data from 2015
m1.sq<-betadisper(dist(sqrt(y)),redd$da15, type="centroid") #calculate average distance to its group centroid for each group
m1.sq
m1.sq.HSD<-TukeyHSD(m1.sq) #run TukeyHSD pair-wise test difference in dispersion between groups 
plot(m1.sq.HSD, sub="sqrt (hnl2) 2015") #graphic display of pair-wise comparisons of difference in dispersion between groups
## Groups have same variance after log transformation
## y.48 - Removed one problem reach - only one pa for 2016
library(vegan) #Required functions for testing are in the 'vegan' Package
m2<-betadisper(dist(y.48),redd.48$da16, type="centroid") #calculate average distance to its group centroid for each group
m2.HSD<-TukeyHSD(m2) #run TukeyHSD pair-wise test difference in dispersion between groups 
plot(m2.HSD) #graphic display of pair-wise comparisons of difference in dispersion between groups
#repeat the same test with log-transformed data
m2.sq<-betadisper(dist(sqrt(y.48)),redd.48$da16, type="centroid") #calculate average distance to its group centroid for each group
m2.sq
m2.sq.HSD<-TukeyHSD(m2.sq) #run TukeyHSD pair-wise test difference in dispersion between groups 
plot(m2.sq.HSD,sub="y.48 log (hnl2) 2016") 
par(mfrow=c(1,1))

#### ***************************************
### Linear DFA for redd density & hnl2 habitat ####
#### ***************************************

## ~~**~~ 2015 ~~**~~ ##

library(MASS) #"lda" function is in the MASS package
## x is grouping variable, y is data matrix
(dfa<-lda(x~scale(log(y+1)))) #run linear Discriminant Function Analysis
## Proportion of trace is eigenvalue. Or, is the percentage seperation achieved by each discriminant function.
## The first discriminant function has a good separation between the groups (64%), & 2nd DF improves seperation by 19% seperation

##Each variable in the eigenvector is the loading on the DF
plot(dfa)
names(dfa)
## Copy dataframe to clipboard in Excel format
ld <- as.data.frame(round(dfa$scaling,2)) 
write.table(ld,"clipboard",sep="\t",col.names=NA) # In Excel, press Ctrl+V 

#show a confusion table (classification results using the DFA model)
# Assess the accuracy of the prediction - percent correct for each category of redd density
dfa.hat <- predict(dfa)
dfa.hat
tb <-table(x,dfa.hat$class)
tb ## The effectiveness of DFA in classifying the groups, by comparing the predicted
## to actual group assignments.
diag(prop.table(tb,1))# if near 1, overfit
sum(tb[row(tb)==col(tb)])/sum(tb)
# total percent correct
#Confusion table: resubstitution error (predictive accuracy),
# is how well the samples are classified into groups.

library(klaR)
g<-greedy.wilks(scale(sqrt(y)),x) #Stepwise variable selection for group classification 
## based on first variable that seperates groups most, then on Wilks lambda criterion.
g

gg <-stepclass(scale(sqrt(y)),x,"lda") #cross-validation selection; includes all, & excludes all
#?stepclass
gg


## ~~** 2016 **~~ ##

## had to remove reach 12, only reach with a pa redd catagory, was outlier
library(MASS) #"lda" function is in the MASS package
dfa1<-lda(x.48~scale(log(y.48+1))) #run linear Discriminant Function Analysis
dfa1 #outputs

## Save loadings, copy to Excel clipboard
ld1 <- as.data.frame(round(dfa1$scaling,2)) #Scaling = Loading = DA coefficients (need to standardize for DA coefficient)
write.table(ld1,"clipboard",sep="\t",col.names=NA) # In Excel, press Ctrl+V 

#show a confusion table (classification results using the DFA model)
(dfa1.hat <- predict(dfa1))
dfa1.hat
tb1 <-table(x.48,dfa1.hat$class)
tb1 ## The effectiveness of DFA in classifying the groups, by comparing the predicted
## to actual group assignments.
sum(tb1[row(tb1)==col(tb1)])/sum(tb1) #Observed on top, preditcted on left
#Confusion table: resubstitution error (predictive accuracy),
# is how well the samples are classified into groups.

## Plot both 2015 and 2016 together
##2015
plot(dfa.hat$x[ ,1], dfa.hat$x[ ,2], pch=as.numeric(x), col=c(1,2,3,4)[x], lwd=2,
     xlab="LD1 (74/84%)", ylab="LD2 (20/16%)", main = "2015 HNL2", xlim = c(-6.5,5), ylim=c(-3.5,4.5)) #eigenvalues for each axis
text(dfa.hat$x[ ,1], dfa.hat$x[ ,2], labels=ip, cex=0.75,adj=c(-.5,-.25), offset=0.5,col="grey")
#legend("topright", c("HA","HP","LA","LP"), pch=c(1,2,3,4), col=c(1,2,3,4), inset=.005, cex = .75) 
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
library(vegan)
ordihull(dfa, x, scaling=3, col='grey',label=F)
ordispider(dfa, x, col='darkolivegreen', lty=3, label = T)
ordiellipse(dfa, x, kind="se", conf=0.95, lwd=2,
            draw = "polygon", col=1:4, border=1:4, alpha=63)
par(new=T)
##2016
plot(dfa1.hat$x[ ,1], dfa1.hat$x[ ,2],pch=as.numeric(x.48), col=c(1,3,4)[x.48], lwd=2,
     xlab="", ylab="", main = "2016 HNL2", xlim = c(-6.5,5), ylim=c(-3.5,4.5)) #eigenvalues for each axis
text(dfa1.hat$x[ ,1], dfa1.hat$x[ ,2], labels=redd.48$ip, cex=0.75,adj=c(-.5,-.25), offset=0.5,col="grey")
ordispider(dfa1,groups=x.48, col="red")
ordihull(dfa1, x.48, scaling=3, col='grey',label=F)
ordispider(dfa1, x.48, col='darkolivegreen', lty=3, label = T)
ordiellipse(dfa1, x.48, kind="se", conf=0.95, lwd=2,
            draw = "polygon", col=1:4, border=1:4, alpha=63)

## HP and pp are large, low gradient unconfined channels (floodplain).  LP & ap are large, confined channels
## LA & aa are high gradient, small, confined channels.  


library(klaR)
g1<-greedy.wilks(scale(sqrt(y)),grouping=x1) # Performs a stepwise forward variable selection for classification
## using the Wilk's Lambda criterion.  The model starts with the variable which seperates the groups most.
## Further variables are then selected based on the Wilk's lambda criterion.
g1

gr1 <-stepclass(scale(sqrt(y.48)),x.48,"lda") #cross-validation selection; includes all, & excludes all
#?stepclass
gr1
