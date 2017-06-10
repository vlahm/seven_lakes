#consumers

dirs<-c("1_mirror_fish","2_mirror_fish","3_C_caddis","4_C_caddis","6_milk_caddis","7_mirror_caddis",
        "8_morg_caddis","9_morg_caddis","10_noname_caddis","11_y015_caddis","12_y025_caddis",
        "13_Z_caddis","90_C_cal","92_O_cal","108_clear_caddis", "94_Z_cal", "85_mirror_cal", 
        "86_y025_cal", "87_y015_cal", "80_clear_chaob", "96_L_cal", "81_clear_clad", "82_clear_cal", 
        "83_morg_cal", "84_morg_clad", "93_O_clad", "105_milk_cal", "106_milk_chaob", "97_L_clad", "91_C_clad", "95_Z_clad")


#producers

dirs<-c("17_O_filAlgae", "32_O_peri", "29_morg_peri", "34_y015_peri", "41_C_pom", "44_clear_pom", "48_milk_pom", "51_mirror_pom",
        "54_morg_pom", "64_Z_pom", "71_O_soil", "73_plant", "74_plant", "75_plant", "76_plant", "77_plant", "109_Z_moss")

FAmerger<-function(dirs.list)
{
  #list that will contain each data frame
  reports<-list()
  
  #vector X will be used to merge and sort each data frame
  X<-as.data.frame(1:98)
  colnames(X)<-"X"
  
  #loop through directories and perform data filtering/extraction
  for(i in 1:length(dirs.list))
  {
    data<-read.csv(paste("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/GC readouts/producers/",dirs.list[i],"/REPORT01.CSV",sep=""))
    data<-data[!is.na(data$X),]
    data<-merge(X,data, all.y=T)
    data<-data[,c(1,3,6)]
    colnames(data)<-c("Peak", paste("Reten_",dirs.list[i], sep=""), paste("Area_",dirs.list[i],sep=""))
    reports[[i]]<-data
  }
  
  #bind all reports together
  merged<-Reduce(function (...) {merge(..., all=T)}, reports)
  
  #make vector of total appearances of each FA
  FA.frequency<-vector(length=length(merged[,1]))
  for(i in 1:length(merged[,1]))
  {
    FA.frequency[i]<-sum(sapply(na.omit(as.numeric(as.vector(merged[i,(2:length(merged[1,]))]))),is.numeric))/2
  }
  
  #bind vector to "merged"
  newCol<-(length(merged[1,])+1)
  merged[[newCol]]<-FA.frequency
  colnames(merged)[newCol]<-"total_appearances"
  
  #add informative prependix
  prependix<-read.csv("C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/Misc/prependix.csv")
  merged<-merge(merged,prependix,all=T)
  
  return(merged)
}

x<-FAmerger(dirs)

#for producers
write.csv(x, file="C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/producer_FA_data.csv")

#for consumers
write.csv(x, file="C:/Users/Mike/Desktop/Grad/Projects/Thesis/Seven Lakes Project 2014/Data/FA/consumer_FA_data.csv")
