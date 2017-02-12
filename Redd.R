##  BRomey.  11Sept2016

## Summer 2016 Redd desity

####  (~:|:~) Global Settings (~:|:~)  #### 
op <- par(mfrow=c(1,1)) # set the default plot grid to 1X1
source("cor.matrix.r")
library(ggplot2)
library(dplyr)
dt <- read.csv("ip16.CSV") #Added rotation in excel using filter based on Jack & Doug starting survey on 19th
summary(dt)
str(dt)


####  Average width of all reaches & density per 50 m
## SCPA is Side Channel Present/Absent
dt$wcw[is.na(dt$wcw)] <- 0 # Change NA's to zero
dta <- dt[,c("ReachID","rotation","Side_Channel","R_length_m","PinkR","ChumR","PpRedd","ChmRedd","wcw")]
rm(dt)
head(dta)
## Side channel y/n
dta$sc <- with(dta,ifelse(dta$Side_Channel == ("SChR"), "y",
                          ifelse(dta$Side_Channel == "SChL", "y",
                            ifelse(dta$Side_Channel == "", "n",""))))

width <- dplyr::filter(dta, wcw > 0) %>%
  filter(sc == "n") ## No side channels, main channel only

w <- width %>%
  group_by(ReachID) %>%
  summarise(aveW = mean(wcw)) %>%
  arrange(ReachID)
dta1 <- merge(dta, w)
rm(width,w)

scw <- dplyr::filter(dta, wcw > 0) %>% # Side Channel Width only
  filter(sc == "y")

scw$sc_area <- scw$wcw*scw$R_length_m # Side channel area
w2 <- scw %>%
  group_by(ReachID) %>%
  summarise(sca = sum(sc_area)) #Total side channel area for reach

dta1 <- merge(dta1,w2,all=T)
dta1$sca[is.na(dta1$sca)] <- 0 # Change NA's to zero
rm(scw,w2)

## Total counts & density for Pink & Chum per 250m reach, including side channel area

dta1$PinkR[is.na(dta1$PinkR)] <- 0 # Change NA's to zero
dta1$ChumR[is.na(dta1$ChumR)] <- 0
dta1$PpRedd[is.na(dta1$PpRedd)] <- 0
dta1$ChmRedd[is.na(dta1$ChmRedd)] <- 0

dta1$pr16 <- dta1$PinkR+dta1$PpRedd
dta1$cr16 <- dta1$ChumR+dta1$ChmRedd

require(dplyr)
pp<- dta1 %>% # Pink redds
  group_by(ReachID, rotation) %>%
  summarise(avW = mean(aveW), #ave channel width
            sca = mean(sca), # side channel area for reach
            pr16 = sum(pr16)) #total redds per 250 reach
pp$hab_area <-(pp$avW*250)+pp$sca # redd habitat area for reach and side channel
pp$pr16_den <- (pp$pr16)/pp$hab_area # redd density for 250 m reach (all rotations)

cc<- dta1%>%  # Chum redds
  group_by(ReachID, rotation)%>%
  summarise(avW = mean(aveW), #ave channel width
            sca = mean(sca), # side channel area for reach
            cr16 = sum(cr16)) #total redds per 250 reach
cc$hab_area <-(cc$avW*250)+cc$sca # redd habitat area for reach and side channel
cc$cr16_den <- (cc$cr16)/cc$hab_area # redd density for 250 m reach (all rotations)

#### ~(1)~ ISOLATE PEAK REDD DENSITY N REDD COUNTS FOR EACH SPECIES PER IP REACH #####

prd <- pp %>% # Pink peak redd density
  group_by(ReachID) %>%
  summarise(prd16 = max(pr16_den))  
crd <- cc %>% # Chum peak redd density
  group_by(ReachID) %>%
  summarise(crd16 = max(cr16_den)) 
prn <- pp %>% # Pink peak redd count
  group_by(ReachID) %>%
  summarise(prn16 = max(pr16))  
crn <- cc %>% # Chum peak redd count
  group_by(ReachID) %>%
  summarise(chn16 = max(cr16))
ha <- cc%>%
  group_by(ReachID)%>%
             summarise(ha16 = mean(hab_area))

## Reach habitat area 2016
write.csv(ha,"ra16.csv")  #Write reach area 2016 file to directory

rden <-  dplyr::full_join(prd,crd,by = "ReachID") # Join Chum and Pink peak redd density by reach ID
redd_n <- dplyr::full_join(prn,crn,by = "ReachID") #join redd count by reach ID
ip16 <- dplyr::full_join(rden,redd_n,by = "ReachID")

rm(prd,crd,prn,crn,rden,redd_n)
names(ip16)[1] <- "ip"  ## peak redd density and redd count for each IP reach

## 2016 Redd density
write.csv(ip16,"redd16.csv")  #Write redd density file to directory

par(mfrow=c(2,1))
hist(ip16$prd16)
hist(ip16$crd16)
par(op)





