## Redd simulation
## https://sites.google.com/site/rforfishandwildlifegrads/home/week-9
## For sample size needed to detect sig dif

## Data- pilot study
redd <- read.csv("redd.csv")
r <-redd[,c(21:24)] #redd density
hnl2 <- read.csv("D:/R/IP16/hnl2.csv")
h <-hnl2[,-c(1,2,4:5,8,11:12,14,17,19:22,25,27:28)] #Only variables with no multicollinearity
rm(hnl2)
## 2015 pink linear model for redd density variables vs hnl2
pr15 <-r[,"pr15"]
p15d <-cbind(pr15,h)
rm(pr15)

summary(p15fm <- glm(log(pr15+1)~BFQ+GRAD_D+GEP,data=p15d))#Compare with simulated

## First simulate covariates
n <- 50# HOW LARGE A SAMPLE SIZE IS NEEDED?
bfq <- runif(n, 5, 33)#range of bfq
grad <- runif(n, 0.0045, 0.046)
B0 <- 4.09773  # THESE VALUES FROM PILOT DATA: Intercept
B1 <- -0.06708  # THESE VALUES FROM PILOT DATA: Slope 1
B2 <- -27.92894  # THESE VALUES FROM PILOT DATA: Slope 2
B3 <- -4.02222
sim_density <- B0 + B1 * bfq + B2 * grad + rnorm(n, mean = 0, sd = 1)
## Another way of adding error for different types of distributions
#pred <- B0 + B1 * bfq + B2 * grad + rnorm(n, mean = 0, sd = 1)
#sim_density <- rnorm(n, mean = pred, sd = 1)
## Plot the simulation
par(mfrow = c(3, 1), mar = c(4, 3, 0, 0), oma = c(1, 2, 1, 1))
plot(bfq, sim_density, ylab = "", las = 1, xlab = "BFQ")
plot(grad, sim_density, ylab = "", las = 1, xlab = "Gradient")
mylab <- expression(paste("Density (Pink Redds 100 meter "^-2, ")", sep = ""))
mtext(side = 2, mylab, outer = TRUE)

#Model the simulation
sim_dat <- data.frame(bfq = bfq, grad = grad, density = sim_density)
summary(fit <- lm(density ~ bfq + grad, sim_dat))
par(mfrow=c(2,2))
plot(fit)

## Extract standard errors and calculate coefficient of variation (CV)
stderrs <- sqrt(diag(vcov(fit)))
ests <- coef(fit)
abs(stderrs/ests) * 100#CV-relative standard deviation; dispersion from the mean predicted line
#Run multiple times to see how stable the coefficients are; including coefficient of determination (R^2)


