## D50 binary catagory for Resource Selection Function (RSF)
## Range of D50 by species taken from Montgomery (2004)

h <- read.csv("D:/R/IP16/hnl2.csv")
h <- h[,c("ip","d50")]

## Pink D50 suitable (suit) vs unsuitable (unui); 25th to 75th %tile (6-10mm)
h$p_D50 <- with(h,ifelse(h$d50 >= 5 & h$d50 <=10, "su",
                        ifelse(h$d50 < 5 | h$d50 >10, "un",
                        ifelse(h$d50 == 0, "error",""))))

## Chum D50 suitable (suit) vs unsuitable (unui); 25th to 75th %tile (13-40mm)
h$c_D50 <- with(h,ifelse(h$d50 >= 13 & h$d50 <= 40, "su", #if inside a range, use &
                         ifelse(h$d50 < 13 | h$d50 >40, "un",#if outside a range, use OR (|)
                                ifelse(h$d50 == 0, "error",""))))


## These D50 ranges for chum and pink salmon do not represent what was observed in the project area per 250m reach.
## The proplem may be associated with the reach mean D50 not representing what is found at each spawning area.