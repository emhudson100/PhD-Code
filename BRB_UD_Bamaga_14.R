
library(adehabitatLT)
library (adehabitatHR)
library(matrixStats)
#setwd ("G:/VPHI/Epi/Projects/66_RabiesModelling_USyd (D?rr)/Emily_PhD/GPS_work/GPS_data")


HR_summary <- c()

sapply (c(1:14), function(period) {

# period = time duration in days for which we want to calcuate the HR

dat <- read.csv("Bamaga_14/Bamaga_14_clean.csv")
dat$TIME <- strptime(as.character(dat$TIME), "%Y-%m-%d %H:%M:%S")
dat$TIME <- as.POSIXct(dat$TIME)

#create dataframe (like study_dur) with start and end dates based on random start day and +24 hours for end
#When changing collar remember to change end date from 24 hours before the end so every point chosen will have 24 hours
#worth of data - change iteration number based on how many sets of days you can get/power of computer

# define the earliest possible end point
earliest_end <- min(dat$TIME [which(dat$TIME >= dat$TIME[1] + period*24*60*60)])
# set the nb_iteration to 500 or to the max possible if less than 500 are available
nb_iteration <- min(500, nrow(dat)-which(dat$TIME == min(dat$TIME [which(dat$TIME >= dat$TIME[1] + period*24*60*60)])))

id <- paste("subset", (seq(1,nb_iteration,by=1)), sep = "")    
end <- sample(dat$TIME [which(dat$TIME >= dat$TIME[1] + period*24*60*60)],nb_iteration, replace = F) #last time point that will give at least 24 hours
start <- end - period*24*60*60
sub.df <- data.frame(id,start,end)

#create subsets of data based on start and end dates
sapply (sub.df[,1], function(x){
  minTime <- paste(sub.df$start[which(sub.df$id == x)])
  maxTime <- paste(sub.df$end[which(sub.df$id ==x)])
  dat <- subset(dat, dat$TIME >= minTime & dat$TIME <= maxTime)
  dat$id <- x
#  fileName <- paste (x)
  assign (paste(x), dat, envir = .GlobalEnv) 
})

#create column in sub.df for number of locations in each subset sub.df = subset dataframe
sub.df$nb_obs <- sapply (sub.df[,1], function (x) {
  nrow(eval(parse(text=paste(x))))
})

#create column in sub.df for average timediff between each point in each subset
sub.df$ave_timediff_min <- sapply (sub.df[,1], function (x) {
  mean(eval(parse(text=paste(x)))$TIMEDIFF/60,na.rm=T)
})

#create a column for min timediff
sub.df$min_timediff_min <-sapply (sub.df[,1], function (x) {
  min(eval(parse(text=paste(x)))$TIMEDIFF/60,na.rm=T)
})

#create a column for max timediff
sub.df$max_timediff_min <-sapply (sub.df[,1], function (x) {
  max(eval(parse(text=paste(x)))$TIMEDIFF/60,na.rm=T)
})

#save sub.df
write.csv(sub.df, file=paste("Bamaga_14/B14_subsetSummary_period",period,".csv",sep=""), row.names=F)

# bind all subsets into one table
all_sub <- data.frame()
sapply (sub.df[,1], function (x) {
  current_sub <- eval(parse(text=paste(x)))
  all_sub <<- rbind(all_sub, current_sub)
})

# preparing data for BRB analysis
da <- strptime(as.character(all_sub[,"TIME"]), "%Y-%m-%d %H:%M:%S")
da <- as.POSIXct(da)
tr_sub <- as.ltraj(xy = all_sub[,c("Zone54_Lon","Zone54_Lat")], date = da, id=all_sub[,"id"])

#Tmax = 90 mins as GPS units take point every 15 mins, Lmin = 2* 22 (average dist between points) based on 15min static test
#hmin = 20 (average distance from centroid) based on 15min static test
BRB_D <- BRB.D(tr_sub, Tmax=90*60, Lmin=44)
BRB_UD <- BRB(tr_sub, D = BRB_D, Tmax = 90*60, Lmin = 44, hmin=20, type = "UD", 
              filtershort= F, extent=1.5)


sapply(1:nb_iteration, function (x) {
  current_sub <- paste(sub.df[x,1])
  for (i in c(50,95)) {
    Name <- paste("BRB_UD_level",i,"_",current_sub,sep="")
    kernUD_i <- getverticeshr(BRB_UD[[x]], percent=i)                                               
    assign (Name, kernUD_i, envir=.GlobalEnv) }
  #kernUD_99 <- getverticeshr(BRB_UD[[x]], percent=99)  # percent 100 is not possible
})

BRB_HR_levels <- c()
sapply(c(50,95), function(j) {
  
  BRB_level_j_HR <- c()
  sapply (sub.df[,1], function (x) {
    cur_sub <- paste (x)
    cur_file_name <- paste("BRB_UD_level",j,"_", cur_sub, sep="")
    cur_BRB_level_j_HR <- eval(parse(text=cur_file_name))[[2]]
    BRB_level_j_HR <<- rbind (BRB_level_j_HR, cur_BRB_level_j_HR)
  })
  
  BRB_HR_levels <<- cbind (BRB_HR_levels, BRB_level_j_HR)
})

rownames(BRB_HR_levels) <- paste(sub.df[,1])
colnames(BRB_HR_levels) <- c("50%", "95%") #only care about 50% and 95%
write.csv(BRB_HR_levels, file=paste("Bamaga_14/HR_levels_",period,".csv",sep=""))

# summary measures for all the iterations 
means <- colMeans(BRB_HR_levels) #means for both 50% and 95% from all subset UD's calculated
median_50 <- median(BRB_HR_levels[,1])
median_95 <- median(BRB_HR_levels[,2])
st.dev <- colSds(BRB_HR_levels)
sem <- st.dev/sqrt(nrow(BRB_HR_levels))
lowerCI <- means-1.96*sem
upperCI <- means+1.96*sem
quant_2.5_50 <- quantile(BRB_HR_levels[,1], probs=0.025)
quant_97.5_50 <- quantile(BRB_HR_levels[,1], probs=0.975)
quant_2.5_95 <- quantile(BRB_HR_levels[,2], probs=0.025)
quant_97.5_95 <- quantile(BRB_HR_levels[,2], probs=0.975)
nb_subsets <- nb_iteration

# summary table for period 
HR_summary_period <- rbind(means, cbind(median_50,median_95), st.dev, sem, lowerCI, upperCI, cbind(quant_2.5_50, quant_2.5_95), 
                    cbind(quant_97.5_50, quant_97.5_95), nb_subsets)
rownames(HR_summary_period) <- c("means", "median", "st.dev", "sem", "lowerCI", "upperCI","quant2.5","quant97.5","nb_subsets")
colnames(HR_summary_period) <- c(paste("50%_",period,sep=""), paste("95%_",period,sep="")) #only care about 50% and 95%

HR_summary <<- cbind(HR_summary, HR_summary_period)

})


write.csv(HR_summary, file= "Bamaga_14/HR_summary.csv")

#save.image(file="Bamaga_14.Rdata")

