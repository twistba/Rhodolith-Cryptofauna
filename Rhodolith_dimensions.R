# --- Crypotofuana Rhodolith paper - NIWA ---
## Author:       Brenton Twist ##
## Date Created: 2020-06-04 


# This script contains code for Running analyses on data collected on cryptofauna from two rhodolith species in northen NZ
# This script examines differences between sites for the same rhodolith species  -
# looking at differences in Rhodolith dimensions

# --- Packages ---
library(tidyverse)

# --- Load data ---
#   The following code loads the dataset

wd<-getwd()# gets working directory 
dim<-read.csv(paste(wd,"/data/rhodolith_dimensions.csv",sep=""))


## --- Summarising data ----

# SE function
se <- function(x, na.rm=FALSE){
  if(na.rm){
    se <- sd(x, na.rm=TRUE)/sqrt(length(na.omit(x)))
  } else {
    se <- sd(x)/sqrt(length(x))
  }
  return(se)
}


grouping1 <- group_by(dim, species_site)
grouping2 <- group_by(dim, species)


volume1 <- summarise(grouping1, value=as.factor("Volume"),mean=mean(volume), 
                     se=se(volume),n=sum(volume/volume))
names(volume1)[names(volume1) == "species_site"] <- "species"

volume2 <- summarise(grouping2, value=as.factor("Volume"),mean=mean(volume), 
                     se=se(volume),n=sum(volume/volume))


free.space1 <- summarise(grouping1,value=as.factor("Free_space"), mean=mean(volume_free_space),
                        se=se(volume_free_space),n=sum(volume_free_space/volume_free_space))
names(free.space1)[names(free.space1) == "species_site"] <- "species"

free.space2 <- summarise(grouping2,value=as.factor("Free_space"), mean=mean(volume_free_space),
                         se=se(volume_free_space),n=sum(volume_free_space/volume_free_space))


tip1 <- summarise(grouping1,value=as.factor("Number_of_tips"), mean=mean(tip_number, na.rm=TRUE),
                 se=se(tip_number,na.rm=TRUE),n=sum(tip_number/tip_number,na.rm=TRUE))
names(tip1)[names(tip1) == "species_site"] <- "species"

tip2 <- summarise(grouping2,value=as.factor("Number_of_tips"), mean=mean(tip_number, na.rm=TRUE),
                  se=se(tip_number,na.rm=TRUE),n=sum(tip_number/tip_number,na.rm=TRUE))


dimensions_summary<-rbind(volume1,volume2,free.space1,free.space2,tip1,tip2)
dimensions_summary

write.csv(dimensions_summary,"./output/rhodolith_dimensions_summary.csv")
