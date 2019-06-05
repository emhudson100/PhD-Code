#setwd("/Users/emilyhudson/Dropbox/PhD/Rabies_Model/Distance_Kernel_Analysis_3/")
library(geosphere)
library(proj4)
library(sp)
library(adehabitatHR)
library(mc2d)
library(parallel)
library(foreach)
library(doParallel)
options(digits = 7)
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

distance <- seq(10,1000, by = 10)

#study_dur is a csv of all dog id's, community names, coordinates of the dogs houses,
#start time fo GPS and end time fo GPS
study_dur <- read.csv("example.csv")


projected_locations <- study_dur[, c(4,3)]

loc_Zone54 <- project(projected_locations, "+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
study_dur$Zone54_Lon <- loc_Zone54$x                                              # add projected lon to the table
study_dur$Zone54_Lat <- loc_Zone54$y

#x <- "Umagico_21"




sapply (study_dur[,1], function (x) {
  
  fileName <- paste (x, "clean", sep="_")
  dat <- read.csv(paste("YOUR_DOG_DATASET",fileName,".csv",sep=""),h=T)
  
  # add a new column of distance between a location to the previous one
  
  house <-as.matrix(cbind(study_dur[which(study_dur$collars == x),4],study_dur[which(study_dur$collars == x),3]))
  points <-as.matrix(cbind(dat$Longitude,dat$Latitude))
  
  dat$Distance.from.home <- distHaversine(house,points,r=6378137) 
  
  dat$bearings <- bearing(cbind(study_dur[which(study_dur$collars == x),4],study_dur[which(study_dur$collars == x),3]), 
                          cbind(dat$Longitude,dat$Latitude))
  
  dat$bearings_radians <- deg2rad(dat$bearings)
  
  
  dat$ang.dist <- dat$Distance.from.home/6378137
  
  
  assign (paste(x), dat, envir = .GlobalEnv) 
  
})


study_dur$category <- as.character(study_dur$collars)

#separate stay-at-home dogs
sah <- study_dur[which(study_dur$category == "sah"),]
  
sah.pairs <- c()

sapply(sah[,1], function(x){
  
  dogs <- data.frame(rep(paste(x),8))
  dogs$second <- as.character(sample(sah[-which(sah$collars == x),1], 8,replace = F))
  
  sah.pairs <<- rbind(sah.pairs,dogs)
  
})

no_cores <- 8
cl<-makeForkCluster(no_cores)
registerDoParallel(cl)


sah.probability <- foreach( z = distance, .combine = "cbind") %:% 
  foreach(x = 1:nrow(sah.pairs), .combine = "rbind") %dopar% {

    start.lat <- deg2rad(-10) # these are arbitrary - can choose any lat/long
    start.long <- deg2rad(142)# these are arbitrary - can choose any lat/long
    reloc_ang.dist <- z/6378137
    reloc_bearings_radians <- deg2rad(sample(0:360,1))
    
    reloc_lat <- asin((sin(start.lat)*cos(reloc_ang.dist))+(cos(start.lat)*sin(reloc_ang.dist)*cos(reloc_bearings_radians)))
    
    reloc_long <- start.long + atan2(sin(reloc_bearings_radians)*sin(reloc_ang.dist)*cos(start.lat),
                                     cos(reloc_ang.dist)-(sin(start.lat)*sin(reloc_lat)))
    reloc_lat <-rad2deg(reloc_lat) 
    reloc_long <- rad2deg(reloc_long)
    
    dogA <- get(as.character(sah.pairs[x,1]))
    dogB <- get(as.character(sah.pairs[x,2]))
    # 
    ###get new geolocations for house A
    
    start.latA <- deg2rad(-10)# these are arbitrary - can choose any lat/long
    start.lonA <- deg2rad(142)# these are arbitrary - can choose any lat/long
    
    latA <- asin((sin(start.latA)*cos(dogA$ang.dist))+(cos(start.latA)*sin(dogA$ang.dist)*cos(dogA$bearings_radians)))
    
    longA <- start.lonA + atan2(sin(dogA$bearings_radians)*sin(dogA$ang.dist)*cos(start.latA),
                                cos(dogA$ang.dist)-(sin(start.latA)*sin(latA)))
    
    newlatA <-rad2deg(latA) 
    newlongA <- rad2deg(longA)
    ID <- rep("dogA", length(newlatA))
    newdogA <- data.frame(ID,newlatA,newlongA)
    
    colnames(newdogA) <- c("ID","Lat", "Long")
    
    ###get new geolocations for house B
    
    start.latB <-deg2rad(reloc_lat) 
    start.lonB <- deg2rad(reloc_long) 
    
    latB <- asin((sin(start.latB)*cos(dogB$ang.dist))+(cos(start.latB)*sin(dogB$ang.dist)*cos(dogB$bearings_radians)))
    
    longB <- start.lonB + atan2(sin(dogB$bearings_radians)*sin(dogB$ang.dist)*cos(start.latB),
                                cos(dogB$ang.dist)-(sin(start.latB)*sin(latB)))
    
    
    newlatB <-rad2deg(latB) 
    newlongB <- rad2deg(longB)
    ID <- rep("dogB", length(newlatB))
    newdogB <- data.frame(ID,newlatB,newlongB)
    
    colnames(newdogB) <- c("ID","Lat", "Long")
    
    locs <- rbind(newdogA,newdogB)
    
    projected_locs <- locs[, c("Long","Lat")]
    
    projected_overlap <- project(projected_locs, "+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
    locs$Zone54_Lon <- projected_overlap$x                                              # add projected lon to the table
    locs$ Zone54_Lat<- projected_overlap$y
    
    coordinates(locs) <- locs[,c("Zone54_Lon","Zone54_Lat")]
    
    
    overlap <- kerneloverlap(locs[,1], grid=500, method="PHR", percent=95, conditional=TRUE)
    overlap[2,1]*overlap[1,2]
    
  }

stopCluster(cl)

write.csv(sah.probability, "sah_sah_all_probs.csv", row.names = F)

sah <- as.data.frame(sah.probability)
lower <- c()
upper <- c()
medians <- c()
sapply(1:length(sah), function (x){
  
  low <- quantile(sah[,x],c(0.025), na.rm = T)
  lower <<- c(lower,low)
  up <- quantile(sah[,x],c(0.975),na.rm = T)
  upper <<- c(upper,up)
  med <- median(sah[,x],na.rm = T )
  medians <<- c(medians,med)
  
})

sah.prob <- data.frame(medians,lower,upper)
colnames(sah.prob) <- c("Median", "Lower", "Upper")
sah.prob$Distance <- distance

write.csv(sah.prob, "sah_kernel.csv", row.names = F) 