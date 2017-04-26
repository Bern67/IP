## REGRESSION/CLASSIFICATION (RPART) TREE - Non-parametric
## & RANDOM FOREST


## The objective is to select the most appropriate predictors for IP geometric mean model (nonlinear), not a
## lm, logit, or DFA model. 

#-----RPART description------
## Reference page 770 of The R Book for tree model details
## PDF - An introductio to recursive partitioning using the RPART routines

## Regression & classificaiton trees are a form of **THRESHOLD ANALYSIS**, or cluster analysis
## Classification tree (factor/catagorical) for Presence/Absence explanatory variable
## Thresholds are based on variance partisioning of the predictors for regression, thus the
## value of a split is defined as the reduction in residual sum of squares.
## Have a tendancy to over-interpret the data (overfit), thus need to prune back
## Uses forward selection of predictors
## Scaling predictors does not matter! Thus, imune to outliers in predictors! no need for transformation

## ** What to REPORT for tree model: 
## 1) final tree model: decision tree
## 2) Cross-validation: how to finalize tree model
## 3) R-sqr
## 4) Misclassification rate (for classification tree): Gini index = total classification error variance for all classes
## 5) Alternative models

## ASSUMPTIONS:  Also refer to Lecture 14 from ESM 567, slide 12
## Non-parametric, no model specification; No model equation produced, just an importance of predictors
## Hierarchical relationship among predictors
## Uses both continuous numerical (regression), and catagorical (classification) predictors
## Branch with largest average class goes to right
#------RandomForest Description -------
## ***** RandomForest looks for thresholds in the data that give low deviance.*****
## In other words, for either classificaiton or regression RF tree, the goal is to maximize homogeneiety
## of response variable for each predictor, with the most important predictor having greatest homogeneity
## Classificaiton uses proporiton, regression uses RSS, & MSE=RSS/n-1 for OOB error

## Random Forest; what to report:
## 1) Model used for analysis
## 2) 
## 3) 
## ***** Variable selection criteria for IP *****:
## Is there field data that coroborates/correlates?
## Does the variable have compounded error?  PCA-Use lowest order of error
## Is the variable biologily important?  Literature


##Data
#----------Explanatory variable------------
## Explanatory variable is redd desity, redd presence/absence
redd <- read.csv("redd.csv")
r <-redd[,c(21:24)] # redd density
hnl2 <- read.csv("D:/R/IP16/hnl2.csv")
## Only variables for within reach, independant, not linear combination,
## and ecologicly important.  Plotted and compared Person correlations in PCA  script.
h <- hnl2[,c("OUT_DIST","MNANPRC_M","MEANANNCMS","BFQ","p_trib","GRADIENT","GRAD_D","FP_WIDTH","VWI_Floor","GEP",
             "GEP_Cum","SINUOSITY","ELEV_M")]

hcat <-read.csv("ip_habcat.csv")
rm(hnl2)

#pr15 <-r[,"pr15"]
#p15d <-cbind(pr15,h1)# data for normal dist of pink 2015 redd density
#rm(pr15)

## Habitat summary table
library(MVN)# reference MVN.pdf
uniNorm(h1, type='SW',desc=T)
uniPlot(h1, type='qqplot')# creates univariate Q-Q plots
uniPlot(h1, type='histogram')# creates univariate histograms
uniPlot(h1, type='box')# Single boxplot with center & normalized
#uniPlot(h1, type='scatter')# same as cor.matrix


## Logistic Presence/absence for cm & pk 2015 & 2016
## Need to first create presence/absence (0/1) for each species

r$ppa15 <- as.factor(with(redd,ifelse(ppa15 =="p", 1,
                                      ifelse(ppa15=="a",0,""))))
r$pr16 <- as.factor(with(r,ifelse(r$pr16 == 0, 0,#absent
                                   ifelse(r$pr16 >0 , 1,""))))#present
r$cr15 <- as.factor(with(r,ifelse(r$cr15 == 0, 0,#absent
                                   ifelse(r$cr15 >0 , 1,""))))#present
r$cr16 <- as.factor(with(r,ifelse(r$cr16 == 0, 0,#absent
                                   ifelse(r$cr16 >0 , 1,""))))#present
ppa15 <- r$ppa15
p15 <- cbind(ppa15,h)
pr16 <-r[,"pr16"]
p16 <- cbind(pr16,h)
cr15 <- r[,"cr15"]
c15 <- cbind(cr15,h)
cr16 <- r[,"cr16"]
c16 <- cbind(cr16,h)
rm(ppa15,cr15,cr16,pr16)


#---------- RF -----------
## RANDOM FOREST VARIABLE SELECTION 
#**********************************

## Pink 2015 Regression Random Forest predictor importance
library(randomForest)
(13/3)# number of 'mtry' regression predictors to use for each random forest
## log transformation response variable to reduce heteroscedasticity (Breiman 2008,p586)
#tuneRF(p15d[,-1],p15d$pr15,ntree=10000,stepFactor = 2,improve=0.05,trace = F,plot=T, doBest=T)
#round(mean(p15.rf$rsq),3)#Pseudo R-sqr.

## Pink 2015 Classificaiton RF
(sqrt(13))
set.seed(167)
(p15.rf <- randomForest(ppa15~.,data = p15, mtry=4, importance=T, corr.bias=F,proximity=F, 
                        do.trace=1000, ntree=5000))# view results

## The %IncMSE is the variance of the predictor, IncNode Purity is the RSS for the predictor
## Variable is important if > absolute value of lowest predictor importance value (%IncMSE for regression, MDA for classification)
varImpPlot(p15.rf,sort=T,n.var=min(28,nrow(p15.rf$importance)),
           main = "Variable Importance (Pink 2015)") # 
#abline(v=3,lty=2,col="grey")#Predictor cutoff-Elbow where index is no longer reduced by a factor of 2X previous diff
## Before model interpretaiton, make sure error stabalizes
plot(p15.rf)# Cross validaiton error; Numbers of trees set to ensure that every input row gets predicted at least a few times.

##RF plot of variable importance
p <- as.data.frame(p15.rf$importance)
p$n <- rownames(p)
p <- p[rev(order(p$IncNodePurity)),]#In reverse order for plot
plot(p$IncNodePurity,type = "h", 
     ylab="Mean Decrease Gini",main = paste("Predictor Importance (Pink 2015)"))
text(x=c(1:28),p$IncNodePurity,p$n,cex=.5)

## Partial dependence plots
par(mfrow=c(2,2),oma=c(1,1,1,0))
partialPlot(p15.rf,ppa15,GRADIENT)# presence of pink redd is associated with out dist > 12
abline(v=12.5,lty=2,col="grey")
partialPlot(p15.rf,p15d,VWI_Floor)#The presence of pink redds is associated with VWI > 4
abline(v=40,lty=2,col="grey")
partialPlot(p15.rf,p15d,GRADIENT)# redd density quickly decreases when slope > ~ 2.4%
abline(v=0.022,lty=2,col="grey")
partialPlot(p15.rf,p15d,MEANANNCMS)# redd density highest when cms is < 2.2 cms
mtext("Logit probability presence?", side=2, outer=TRUE, line=-1) 
## Valley width, VWI, and Val constraint all have positive relationships with redd density. 
## Such that an increase in valley constraint result in an increase in redd density


##-------- Pink 16 Random Forest classification predictor importance
#*******************************
(sqrt(13))# Number of 'mtry' predictors to use for each random forest in classificaiton
set.seed(67)
(p16.rf <- randomForest(pr16~.,data = p16, mtry=4, do.trace=500,
                        importance=T, ntree=5000))# view results
varImpPlot(p16.rf,sort=T,n.var=min(28,nrow(p16.rf$importance)), cex=1,#use grey scale
           main = "Pink 2016")
abline(v=2.85,lty=2,col="grey")
plot(p16.rf,type="l")#Cross validation error; How many random trees to use for good performance
## Error rate or MSE

## Partial dependence plots; predicted outcome using X, after averaging out all other predictors
par(mfrow=c(2,2),oma=c(1,1,1,1))#oma creates margin size for lines of text
partialPlot(p16.rf,p16,FitElev)# redd presence increases at 25, max out at 75(randomForest.pdf,p.12)
partialPlot(p16.rf,p16,GRADIENT)
partialPlot(p16.rf,p16,VWI_Floor)
partialPlot(p16.rf,p16,MEANANNCMS)
mtext("Logit probability presence", side=2, outer=TRUE, line=-1) 

#round(importance(p16.rf),2) # importance of each predictor

#-------- Chum 15 Random Forest classification predictor importance
##*****************************
library(randomForest)
set.seed(67)
(c15.rf <- randomForest(cr15~.,data = c15, mtry=4, do.trace=500, 
                        importance=T, ntree=5000))# view results
varImpPlot(c15.rf,sort=T,n.var=min(28,nrow(c15.rf$importance)),
           main = "Variable Importance (Chum 2015")
abline(v=2.7,lty=2,col="grey")
plot(c15.rf)#Cross validation error

## Partial dependence plots
par(mfrow=c(2,1),oma=c(1,1,1,0))
#partialPlot(c15.rf,c15,Shear)
#abline(v=50,lty=2,col="grey")
partialPlot(c15.rf,c15,GRADIENT)
abline(v=0.011,lty=2,col="grey")
partialPlot(c15.rf,c15,MEANANNCMS)
abline(v=3,lty=2,col="grey")
partialPlot(c15.rf,c15,VWI_Floor)# Use to verify constrained channel difference (randomForest.pdf,p.12)
mtext("Logit probability presence", side=2, outer=TRUE, line=-1) 
partialPlot(c15.rf,c15,WIDTH_M)


#round(importance(c15.rf),2) # importance of each predictor
##RF plot of variable importance
p2 <- as.data.frame(c15.rf$importance)
p2$n <- rownames(p2)
p2 <- p2[rev(order(p2$MeanDecreaseGini)),]
plot(p2$MeanDecreaseGini,type = "h", 
     ylab="Mean Decrease Gini",main = paste("Predictor Importance (Chum 2015)"))
text(x=c(1:28),p2$MeanDecreaseGini,p2$n,cex=.5)


## Chum 2016 Random Forest classification predctor importance
#********************************************
library(randomForest)# Needed for confusion matrix and error rate
set.seed(167)
(c16.rf <- randomForest(cr16~.,data = c16, mtry=4, do.trace=500, proximity=T, 
                        importance=T, ntree=5000))# view results
varImpPlot(c16.rf,sort=T,n.var=min(28,nrow(c16.rf$importance)),#significant predictors of redd density
           main = "Variable Importance (Chum 2016)")# Biased toward correlated predictors, compare with party below
plot(c16.rf)

(1-c16.rf$confusion[2,3])*100#sensitivity = percent present correctly classified (Cutler 2007)
(1-c16.rf$confusion[1,3])*100#specifity = percent absent correctly classified (Cutler 2007)
#Percent Correclty Classified (PCC) = (1-OOB error rate)*100 (Cutler 2007)

## Look at variable importance:
#round(importance(c16.rf),2) # importance of each predictor
##RF plot of variable importance
p3 <- as.data.frame(c16.rf$importance)
p3$n <- rownames(p3)
p3 <- p3[rev(order(p3$MeanDecreaseGini)),]#In reverse order for plot
plot(p3$MeanDecreaseGini,type = "h", 
     ylab="Mean Decrease Gini", main = paste("Chum 2016 variable importance"))
text(x=c(1:28),p3$MeanDecreaseGini,p3$n,cex=.5)

## Partial dependence plots
par(mfrow=c(2,1),oma=c(1,1,1,0))
#partialPlot(c16.rf,c16,SRC_DIST)
#abline(v=c(6,12),lty=2,col="grey")#
partialPlot(c16.rf,c16,GRADIENT)
abline(v=c(.006,.014),lty=2,col="grey")
partialPlot(c16.rf,c16,MEANANNCMS)
abline(v=c(1.5,5.6),lty=2,col="grey")

partialPlot(c16.rf,c16,VWI_Floor)# Use to verify constrained channel difference (randomForest.pdf,p.12)
abline(v=c(3.1,9),lty=2,col="grey")
mtext("Logit probability of presence", side=2, outer=TRUE, line=-1) 


## ******** Variable imporance plots for all P/A 
pk15 <- p15.rf$importance[,2]
pk16 <- p16.rf$importance[,4]
cm15 <- c15.rf$importance[,4]
cm16 <- c16.rf$importance[,4]
v_imp <- as.data.frame(cbind(pk15,pk16,cm15,cm16))
rm(pk15,pk16,cm15,cm16)
v_imp$n <- rownames(v_imp)

par(mfrow=c(2,2))
p1 <- v_imp[rev(order(v_imp$pk15)),]# set the order
plot(p1$pk15,type = "h", xlim=c(0,14), 
     ylab="Mean decrese Gini index", xlab="Variable", main = paste("Pink 2015"))
text(x=c(1:13),p1$pk15,p1$n,cex=.7)

p2 <- v_imp[rev(order(v_imp$pk16)),]# set the order
plot(p2$pk16,type = "h", xlim=c(0,14), 
     ylab="Mean decrease Gini index",xlab="Variable", main = paste("Pink 2016"))
text(x=c(1:13),p2$pk16,p2$n,cex=.7)

c1 <- v_imp[rev(order(v_imp$cm15)),]
plot(c1$cm15,type = "h", xlim=c(0,14), 
     ylab="Mean decrease Gini index",xlab="Variable", main = paste("Chum 2015"))
text(x=c(1:13),c1$cm15,c1$n,cex=.7)

c2 <- v_imp[rev(order(v_imp$cm16)),]
plot(c2$cm16,type = "h", xlim=c(0,14), 
     ylab="Mean decrease Gini index",xlab="Variable", main = paste("Chum 2016"))
text(x=c(1:13),c2$cm16,c2$n,cex=.7)

#----------------------------------------------------
## Used for unsupervised learning when no class labels (ex. Present/Absent)
## Do MDS on 1 - proximity:
rf4.mds <- cmdscale(1 - c16.rf$proximity, eig=TRUE)
op <- par(pty="s")
pairs(cbind(c16[,c('SRC_DIST','MEANANNCMS','GRADIENT','VWI_Floor')], rf4.mds$points), cex=0.9, gap=0,#
      col=c("orange", "darkgreen")[as.numeric(hcat$c)],
      main="Chum16 Data: Predictors and MDS of Proximity Based on RandomForest")
par(op)
print(rf4.mds$GOF)


