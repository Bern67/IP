## Chum length for Game Creek during 2015 - From Foot survey for otolith sub-contracted to ADFG
## E. Eric Knudsen, Ph.D., Manager, PWSSC Hatchery-Wild Program - eericknudsen@gmail.com

## Pink lengths from unknown stream in study area (Dickerson et al. 2002)

## Calculations from Riebe et al. 2014.

## ---- Chum salmon ----
cm_l <- read.csv("ChumLength.csv")# from Eric Knudsen
(cl <- aggregate(MEHLength~Sex,cm_l,mean))# table of male and female

(L_f <- cl[1,2])# Length of female chum (mm)
## A_redd = 3.3[L/600]^2.3; where L=fish length (in mm)
(A <- round(3.3*(L_f/600)^2.3,2))# chum redd m^2
## D_t = 115[L/600]^-0.62
(D <- round(115*(L_f/600)^0.62,2))# Threshold particle size (D_t) for chum

## ---- Pink salmon (Dickerson et al. 2002  ----
(pl <- mean(c(426,414,439,390)))# pink fish length; 1997-2000
(A_p <- round(3.3*(pl/600)^2.3,2))# pink redd m^2
(D_p <- round(115*(pl/600)^0.62,2))# pink threshold particle size (D_t)

## For each reach we surveyd, redd size and movable particle size threshold (D_t) was used to identify fish species
## for redds with no fish present at the time of survey.  


## For the redd calibration for estimated area of 0.36, refer to Wright 2011, pg 31, second paragraph

r1 <- 18
r2 <- 14
1-(r2/r1)
