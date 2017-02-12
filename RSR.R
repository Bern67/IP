## Manly  Resource Selection Ratio (RSR) for catagorical proportions
## Refer to adehabitatHS.pdf, section 3.4
## Uses ChiSq analysis and assumes a normal distribution.
# Not used for this analysis, distribution not normal
## Compare to Boodtrap results


redd <- read.csv("redd.csv")
r <-redd[,c(2,14:19)]
rm(redd)
hc <- read.csv("ip_habcat.csv")
cl <- hc[,c(1,3:5)]#class
cn <-merge(r,cl)
rm(cl,hc,r)
c <-cn[,"c"]
(tc <-table(c))

##Neither log nor sqrt normalize
boxplot(log(cn[,c(2:5)]+1))
boxplot(sqrt(cn[,c(2:5)]))

library(MVN) # load MVN package
uniNorm(sqrt(cn[,c(2:5)]), type = "SW", desc = T) # summary statistics and Shipiro-Wilks normality test
##Pink15 normal dist. after sqrt transformation

##---- Constrained: Create species and habitat class matrix for resourse selection ratio model ----
## Number of redds in each class
c15 <- as.matrix(xtabs(crn15~c, data=cn))#creates a contingency table: number of chum in each confinement class
c16 <- as.matrix(xtabs(crn16~c,cn)) #needs to be matrix, not table
p15 <- as.matrix(xtabs(prn15~c,cn))
p16 <- as.matrix(xtabs(prn16~c,cn))
##Ratio of Selected habitat in each class/ total habitat in class
h15 <- as.matrix(xtabs(ha15~c, cn))
av15 <- as.matrix(h15/sum(h15))
h16 <- as.matrix(xtabs(ha16~c,cn))
av16 <- as.matrix(h16/sum(h16))
  

## Habitat Selection for Constrained vs Unconstrained
## Chum 15
library(adehabitatHS)
?widesI 
(Wi <- widesI(c15,av15, avknown = F, alpha = 0.1)) #Needs to be vectors
names(Wi)
par(mfrow=c(2,2))
plot(Wi)#significant plot with CI
## Chum 16
(Wi1 <- widesI(c16,av16, avknown = F,alpha = 0.1))#avknown = True if population available proportion is known
par(mfrow=c(2,2))
plot(Wi1)#significant

## Pink 15
(Wi2 <- widesI(p15,av15,alpha = 0.1))
par(mfrow=c(2,2))
plot(Wi2)
## Pink 16
(Wi3 <- widesI(p16,av16,alpha = 0.1))
par(mfrow=c(2,2))
plot(Wi3)


##---- Trib: Create species and habitat class matrix for resourse selection ratio model ----
## Number of redds in each class
c15_t <- as.matrix(xtabs(crn15~trib, data=cn))#creates a contingency table: number of chum in each confinement class
c16_t <- as.matrix(xtabs(crn16~trib,cn)) #needs to be matrix, not table
p15_t <- as.matrix(xtabs(prn15~trib,cn))
p16_t <- as.matrix(xtabs(prn16~trib,cn))
##Ratio of Selected habitat in each class/ total habitat in class
h15_t <- as.matrix(xtabs(ha15~trib, cn))
av15_t <- as.matrix(h15/sum(h15))
h16_t <- as.matrix(xtabs(ha16~trib,cn))
av16_t <- as.matrix(h16/sum(h16))

(Wi_t <- widesI(c15_t,av15_t,alpha = 0.1))
(Wi_t1 <- widesI(c16_t,av16_t, avknown = F, alpha = 0.1))
(Wi_t2 <- widesI(p15_t,av15_t,alpha = 0.1))
(Wi_t3 <- widesI(p16_t,av16_t,alpha = 0.1))


##---- Geol: Create species and habitat class matrix for resourse selection ratio model ----
## Number of redds in each class
cn_g <- cn[-c(13,17,30,40),]
c15_g <- as.matrix(xtabs(crn15~geol, data=cn_g))#creates a contingency table: number of chum in each confinement class
c16_g <- as.matrix(xtabs(crn16~geol,cn_g)) #needs to be matrix, not table
p15_g <- as.matrix(xtabs(prn15~geol,cn_g))
p16_g <- as.matrix(xtabs(prn16~geol,cn_g))
##Ratio of Selected habitat in each class/ total habitat in class
h15_g <- as.matrix(xtabs(ha15~geol, cn_g))
av15_g <- as.matrix(h15/sum(h15))
h16_g <- as.matrix(xtabs(ha16~geol,cn_g))
av16_g <- as.matrix(h16/sum(h16))

(Wi_g <- widesI(c15_g,av15_g,alpha = 0.1))
(Wi_g1 <- widesI(c16_g,av16_g, avknown = F, alpha = 0.1))
(Wi_g2 <- widesI(p15_g,av15_g,alpha = 0.1))
(Wi_g3 <- widesI(p16_g,av16_g,alpha = 0.1))

