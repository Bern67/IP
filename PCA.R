## PCA

## Master redd.csv file created in DFA script, line 132

## Assumptions for equal variance and normality not met for all habitat data
## Thus, exploratory analysis
## Weather to use variance or co-variance matrix:
## If all variable are the same unit of measure, use covaraince matrix. If all variable are different
## units of measure, need to standardize variables and use correlation matrix if the data is somewhat multinormal
## Default is correlation matrix, so only need to scale (standardize to Z-score)

## The habitat data is the response, and the overlaid redd data is the predictor variable

## PCA:
## Eigenvector is the sum of the loadings on each axis; PC1 is an eigenvector


##-----Know your data ----
hnl2 <- read.csv("D:/R/IP16/hnl2.csv", header=T)
r <- read.csv("redd.csv")
r <- r[,-1]

r$plh15 <- as.factor(with(r,ifelse(r$pr15 <= 3.5, 'l',#low; based on GRADIENT of 0.016,RPRT
           ifelse(r$pr15 >3.5 , 'h',""))))#high

par(mfrow=c(2,2))
boxplot(hnl2[,-1], main="Raw data")
boxplot(sqrt(hnl2[,-1]), main="sqrt transformed")
boxplot(scale(hnl2[,-1],center=T), main="Standardized & Centered")
boxplot(scale(sqrt(hnl2[,-1])),center=F, main="sqrt & scaled to Z-score")

h_c15 <- hnl2
h_c15$c15 <- r$crd15

round(cov(h_c15[,-1]),2)# Covariance matrix: variance & covariance (covariance is diagonal)
##Correlation matrix - all diagonal elements are equal to 1
round(var(scale(h_c15[,-1])),2)#standardized variance matrix = correlation matrix
round(cor(h_c15[,-1]),2)#correlation matrix; default method="pearson" -> parametric

## Assumptions for PCA:
## Multivariate normality:
cp <- as.matrix(hnl2[,-1])
chisplot <- function(cp) {
  if (!is.matrix(cp)) stop("x is not a matrix")
  ### determine dimensions
  n <- nrow(cp)
  p <- ncol(cp)
  #
  xbar <- apply(cp, 2, mean)
  S <- var(cp)
  S <- solve(S)
  index <- (1:n)/(n+1)
  #
  xcent <- t(t(cp) - xbar)
  di <- apply(xcent, 1, function(cp,S) cp %*% S %*% cp,S)
  #
  quant <- qchisq(index,p)
  plot(quant, sort(di), ylab = "Ordered distances",
       xlab = "Chi-square quantile", lwd=2,pch=1)
  return (di)
}
chisplot(cp)#Multinormal plot
abline(a=0,b=1)
## Multinormal Shapiro test
library(mvnormtest)
mshapiro.test(t(hnl2[,-1]+1)) #data has to be a matrix and has to be transposed
## If p < .05, then sig dif than normal distribution

library(MVN)# reference MVN.pdf
(mnt <- mardiaTest(sqrt(hnl2[,-1]),qqplot=T))

## Boxplot Look at each variable distribution
par(mfrow=c(3,4))
names<-colnames(cp) # 2015
for (i in 1:ncol(cp)){
  boxplot(sqrt(cp)[,i],main=paste("Boxplot of",names[i]),xlab="Groups", ylab=paste(names[i]))
}

## Histogram Look at each variable distribution
par(mfrow=c(3,4))
names<-colnames(cp) # 2015
for (i in 1:ncol(cp)){
  hist((cp)[,i],main=paste("histogram of",names[i]),xlab="Groups", ylab=paste(names[i]))
}

##Different multi-normality tests from the MVN package
ip <- hnl2$ip
## ******** REmoved colinear and redundant variabl *********
h <- hnl2[,-c(1,6,8,10)]#No ip number column, & FROM_DIST is a redundant variable; see line 111
names(h)

##Multinorm test
mardiaTest(h,qqplot=T)# Mardia (1970) using squared Mahalanobis distance
hzTest(h, qqplot = F)# Distance bwtween distribution functions
roystonTest(h,qqplot=F)# Shapiro-Wilk/Shapiro-Francia statistic

uniPlot(h,type="qqplot")#univariate
uniPlot(h,type="histogram")#univariate
uniNorm(h,type="SW",desc=F)#Shapiro-Wilk normality test on each variable
uniNorm(log(h+1),type="SW",desc=F)
uniNorm(sqrt(h),type="SW",desc=F)# Looks like sqrt normalizes the most habitat variables
## Normal w/o t.: WIDTH_M, BFQ
## log t.: OUT_DIST, SRC_DIST, FROM_DIST, FlowVel, FP_WIDTH, VAL_WIDTH, VWI_Floor, d50, Shear
## Sqrt t.: AREA_SQKM, ELEV_M, FitElev, MEANANNCMS, p_trib, GEP

##Multivariate outlier detection
# Mahalanobis distance
par(mfrow=c(1,2))
res2 <- mvOutlier(h, qqplot = TRUE, method = "quan")
# Adjusted Mahalanobis distance
res2 <- mvOutlier(h, qqplot = TRUE, method = "adj.quan")#New dataset with declared outliers removed

source("cor.matrix.r")#Linear combinations
cor.matrix(scale(sqrt(h[,c("OUT_DIST","SRC_DIST","DEPTH_M","WIDTH_M","AREA_SQKM","MEANANNCMS")]))) #PC1+; Meananncms
cor.matrix(h[,c("StrmPow","FlowVel","Shear","d50", "GRADIENT")]) #PC1-; Gradient
cor.matrix(scale(sqrt(h[,c("ValCnstrnt","VWI_Floor","VAL_WIDTH","FP_WIDTH")]))) #PC2-; ValCnstrnt
cor.matrix(h[,c("GEP_DEL","GEP","BFQ","GEP_Cum")])#PC2+; BFQ or GEP
## OUT_DIST = Distance downstream to estuary
## SRC_DIST = Distance upstream of reach
## FROM_DIST = Distance donwnstream to confluance, or estuary for mainstem, thus removed -> same as OUT_DIST for all but 10 sites

## From Distance & Out Distance are same;is distance to confluance node, this explains the multiple correlation lines with SRC_DIST

##*******************
## Transformations
##*******************

### Since all reach habitat variables are continuous, we will use log transformation

library(mvnormtest)
mshapiro.test(t(h)) #data has to be a matrix and has to be transposed
## If p < .05, then sig dif than normal distribution

#Continuous predictor (log = ln;natural log, log10 = log base 10)
hb <- h[,c('WIDTH_M',"BFQ")]
h_ln <- log(h[,-c(7,14)]+1)# all others log transformed

th <-  cbind(h_ln,hb)#transformed habitat data
rm(h_ln,hb)
library(MVN)
roystonTest(th,qqplot=T)#Close enought to multivariate normal for exploratory PCA
hzTest(th,qqplot=T)
mardiaTest(th,qqplot=T)
(unrm <- uniNorm(th,type="SW",desc=T))# summary statistics and Shipiro-Wilks normality test
write.table(unrm$`Descriptive Statistics`,"clipboard",sep="\t",col.names=NA) # In Excel, press Ctrl+V 

mvOutlier(th, qqplot = TRUE, method = "adj.quan")


#**********************************************
## What variables have the highest correlation?
#**********************************************
mhc <- function(th,numtoreport)#Most Highly Correlated (MHC) function
{
  # find the correlations
  cormatrix <- cor(h,method='pearson')# pearson(default), kendall, spearman (nonparametric)
  # set the correlations on the diagonal or lower triangle to zero,
  # so they will not be reported as the highest ones:
  diag(cormatrix) <- 0
  cormatrix[lower.tri(cormatrix)] <- 0
  # flatten the matrix into a dataframe for easy sorting
  fm <- as.data.frame(as.table(cormatrix))
  # assign human-friendly names
  names(fm) <- c("First.Variable", "Second.Variable","Correlation")
  # sort and print the top n correlations
  head(fm[order(abs(fm$Correlation),decreasing=T),],n=numtoreport)
}
mhc(h[2:27], 27)#Show top 20 correlations


## ||||||||||||||||||
##-----*~*   PCA  *~*---------------
## ||||||||||||||||||
require(MASS) #loads the PCA package
rownames(th)<-ip#changes row number to site number for biplot site number (must be unique)
(summary(pca <- princomp(scale(th), cor=T))) #creates a PC matrix using the correlation matrix; scales & centers different variables
signif(pca$sdev^2,4)#eigenvalues = sd^2

## Biplot and assumption plots
with(th,
     {biplot(pca, main = "", arrow.len=0.05,expand=0.9, cex=c(.75,.75), col=c("grey","darkgreen"),
            xlab="PC1 (42.1% variance)", ylab="PC2 (20.9% variance)")
     abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
     abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
#     mtext('Unconstrained  <----->  Constrained', side=4, line=-1)
#     mtext('High gradient,small  <----->  Low gradient, large', side=3, line=2.2)
     ## PC1 is habitat size (flow), and a channel gradient (slope)
     ## PC2 is valley constraint (VWI), and sediment transport (GEP)
     ## Flow and gradient are opposite sign indicating the best predictor seperation between A/P
     ## Refer to Baltz & Moyle 1993 PCA in california for methods and results 
     st <- r[,'st']
     rownames(th)<-st
     (summary(pca <- princomp(scale(th)))) #creates a PC matrix using the c
     biplot(pca, main = "Biplot", arrow.len=0.05,expand=.95, cex=c(.75,.75), col=c("darkgreen","grey"),
            xlab="PC1 (42.1%)", ylab="PC2 (20.9%)")
     abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
     abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
     
     # biplot: Scale for sites(PC matrix-pca$scores) on top, scale for variables (vectors-loadings) along bottom
     round(pca$scores[ ,c(1:4)],2) # PC Matrix (eigenvalues); use to plot PC1 vs PC2
     (pc <- round(pca$loadings[,c(1:3)],3)) # The matrix of variable loadings (i.e, colums of Eigenvectors)
     write.table(pc,"clipboard",sep="\t",col.names=NA) # In Excel, press Ctrl+V 
     
     library(vegan)
     screeplot(pca, bstick = T, npcs = 8, bst.col = "darkgreen", ylab = "Variance") #inertia = variance in PCA
     ## The most obvious change in slope occurs between PC1& PC2 ("elbow"), thus first two retained
     ## Or, a minimum of 50% of the variance
     
     source("broken.stick.r")
     broken.stick(3)#Keep increasing until last Comp. < proportion of variance
     summary(pca)# compare to proportion of variane
     #Expected first eigenvalue (PC variance) by random is 0.61, if the 1st eigenvalue from your data > 0.61, keep it
     #Repeat the above until the eigenvalues from your data is < the eigenvalues generated by random
     ## The broken stick model indicated that the PCA was reduced to 3 principal components that make up
     (tv <- 0.3992822+0.1930861+0.1345968)
     ## a total variance of 73%.
     
     ## Shepard diagram; compare euclidean distances
     ## Are the distances among samples well preserved in the reduced space? 
     euc <-dist(log(hnl2[,-1])+1)
     euc.1 <-dist(pca$scores[,c(1,2,3)])
     plot(euc,euc.1,main="PC=3",col="darkgreen", xlab="Distance in Multidimentional space", 
          ylab="Distance in Reduced space") # Shepard diagram
})

## Bubbleplot of Pink redd density for 2015
plot(pca$scores[ ,1], pca$scores[ ,2], type="p",xlab="PCA1 (42.1%)", ylab="PCA2 (20.9%)")
symbols(pca$scores[ ,1], pca$scores[ ,2],circles=r$pr15,inches=.5, fg = "white", 
        bg = "grey", add=T, lwd=2)  
abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
title(main="Pink 2015 Redd Density")
mtext('Unconstrained  <----->  Constrained', side=4, line=1)

## Cm & pk 2015-16 low and high density vs principal component
with(r,
     {plh15 <- as.factor(plh15) #Need to change character to factor for pch to work as numeric
     plot(pca$scores[,1], pca$scores[,2], pch=c(17,16)[plh15],col=c("darkgreen","darkgrey")[r$plh15], 
          lwd=2,cex=1.2,xlab="PCA 1 (42.1%)", ylab="PCA 2 (20.9%)")
     legend("topright",c("L","H"),pch=c(16,17),cex=1,col=c("grey","darkgreen"),inset=.01) 
     title(main="2015 Pink L/H Density")
     abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
     abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
     
     ## Redd density
     par(mfrow=c(2,2))
     plot(pr15,pca$scores[,1],pch=16,col="blue",ylab="PC 1",xlab = "Pink redd density 2015")
     plot(pr15,pca$scores[,2],pch=16,col="blue",ylab = "PC 2", xlab="Pink redd density 2015")
     plot(cr15,pca$scores[,1],pch=17,col="darkgreen",ylab = "PC 1", xlab="Chum redd density 2015")
     plot(cr15,pca$scores[,2],pch=17,col="darkgreen",ylab = "PC 2", xlab="Chum redd density 2015")
     ## Both chum and pink redd density incresed toward medium size, low gradient habitat (PC 1), &
     ## unconfined channels (PC2) - High ave flows
     ## Difference is starting low densities for both species for PC 1; chum primarily low gradient large hab
     ## pinks evenly for both large and small habitat.
     
     par(mfrow=c(2,2))
     plot(pr16,pca$scores[,1],pch=16,col="blue",ylab = "PC 1", xlab="Pink redd density 2016")
     plot(pr16,pca$scores[,2],pch=16,col="blue",ylab = "PC 2", xlab="Pink redd density 2016")
     plot(cr16,pca$scores[,1],pch=17,col="darkgreen",ylab = "PC 1", xlab="Chum redd density 2016")
     plot(cr16,pca$scores[,2],pch=17,col="darkgreen",ylab = "PC 2", xlab="Chum redd density 2016")
     ## Both chum and pink redd density increased toward larger low gradienet habitat (PC1) &
     ## unconfined channels (PC 2) - Low ave flows
     
     ## Look at densities for PC3 also
     par(mfrow=c(2,2))
     plot(pr16,pca$scores[,1],pch=16,col="blue",ylab = "PC 1", xlab="Pink redd density 2016")
     plot(pr16,pca$scores[,3],pch=16,col="blue",ylab = "PC 3", xlab="Pink redd density 2016")
     plot(cr16,pca$scores[,1],pch=17,col="darkgreen",ylab = "PC 1", xlab="Chum redd density 2016")
     plot(cr16,pca$scores[,3],pch=17,col="darkgreen",ylab = "PC 3", xlab="Chum redd density 2016")
     
})


## Pink analysis:
## 2015-high flow: The data sugests that the higher densities associated with all channels that converge to unconfined
## 2016-low flow: The data sugests that the higher densities converge to large unconfined channel

## Chum analysis:
## 2015-high flow: The pattern in the data suggests that densities are primarily larger channels converging to unconfined
## 2016-low flow:  The data suggests that densites transition from moderate to large habitat & confined to unconfined



## ----- Pink 2015 P/A -----
with(r,
     {par(mfrow=c(2,2))
       ppa15 <- as.factor(ppa15) #Need to change character to factor for pch to work as numeric
       plot(pca$scores[,1], pca$scores[,2], pch=c(16,17)[ppa15], col=c("darkgrey","darkgreen")[ppa15], 
            lwd=2, cex=1.2, xlab="PCA 1", ylab="PCA 2")
       legend("topright",c("A","P"), pch=c(16,17), cex=1.1, col=c("grey","darkgreen"),inset=.01) 
       title(main="Pink 2015")
       abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
       abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
       ## Pink 2016 P/A
       ppa16 <- as.factor(ppa16) #Need to change character to factor for pch to work as numeric
       plot(pca$scores[,1], pca$scores[,2], pch=c(16,17)[ppa16],col=c("grey","darkgreen")[ppa16], 
            lwd=2,cex=1.2,xlab="PCA 1", ylab="PCA 2")
       legend("topright",c("A","P"),pch=c(16,17),cex=1.1, col=c("grey","darkgreen"),inset=.01) 
       title(main="Pink 2016")
       #mtext('Unconstrained  <----->  Constrained', side=4, line=.75,cex=.75)
       abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
       abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
       ## Chum 2015 P/A
       cpa15 <- as.factor(cpa15) #Need to change character to factor for pch to work as numeric
       plot(pca$scores[,1], pca$scores[,2], pch=c(16,17)[cpa15], col=c("grey","darkgreen")[cpa15], 
            lwd=2,cex=1.2, xlab="PCA 1", ylab="PCA 2")
       legend("topright",c("A","P"),pch=c(16,17),cex=1.1,col=c("grey","darkgreen"),inset=.01) 
       title(main="Chum 2015")
       abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
       abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
       ## Chum 2016 P/A
       cpa16 <- as.factor(cpa16) #Need to change character to factor for pch to work as numeric
       plot(pca$scores[,1], pca$scores[,2], pch=c(16,17)[cpa16],col=c("grey","darkgreen")[cpa16], 
            lwd=2,cex=1.2,xlab="PCA 1", ylab="PCA 2")
       legend("topright",c("A","P"),pch=c(16,17),cex=1.1,col=c("grey","darkgreen"),inset=.01) 
       title(main="Chum 2016")
       abline(v=(seq(0,100,25)), col="lightgray", lty="dotted")
       abline(h=(seq(0,100,25)), col="lightgray", lty="dotted")
})


## The Pink plots show how flows influanced fish distribution during both years.  In 2015, pinks were able
## to access more habitat in the small category.  In 2016, they moved into the larger habitat due to low flow.
## The Chum data is telling; It stays relatively the same with the exception of a few small tribs,
## regardless of the large difference in flow between years.



## ------- Same plots as above using ggplot ------
## Chum 15
library(ggplot2)
ggplot(fsh, aes(pca$scores[ ,1], pca$scores[ ,2])) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(colour = factor(cpa15), shape = factor(cpa15)), size = 4) +
  labs(list(title = "", x = "PCA1 (40%)", y = "PCA2 (19%)"))
## Chum 16
ggplot(fsh, aes(pca$scores[ ,1], pca$scores[ ,2])) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(colour = factor(cpa16), shape = factor(cpa16)), size = 4) +
  labs(list(title = "", x = "PCA1 (40%)", y = "PCA2 (19%)"))
## Pink 15
ggplot(fsh, aes(pca$scores[ ,1], pca$scores[ ,2])) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(colour = factor(plh15), shape = factor(plh15)), size = 4) +
  labs(list(title = "", x = "PCA1 (40%)", y = "PCA2 (19%)"))
## Pink square root density 2015
ggplot(fsh, aes(pca$scores[ ,1], pca$scores[ ,2])) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(colour = sqrt(pr15)), size=5)+ scale_colour_gradient(low = "grey", high = "darkgreen") +
  labs(list(title = "", x = "PCA1 (40%)", y = "PCA2 (19%)"))

## Pink 16
ggplot(fsh, aes(pca$scores[ ,1], pca$scores[ ,2])) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_point(aes(colour = factor(ppa16), shape = factor(ppa16)), size = 4) +
  labs(list(title = "", x = "PCA1 (39%)", y = "PCA2 (18%)"))





