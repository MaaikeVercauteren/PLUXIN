########################
#libraries
########################

library(dplyr)
library(ggplot2)
library(tidyverse)
library(openxlsx)

########################
#datasets
########################
data<-read.csv("Final analysis/Final dataset/data_full_quality.csv")



################################################
## Check number of samples in each category
################################################
dim(data) #should be 62236 rows and 45 columns

#how many samples in each category: 
#check how many samples in each category contains plastics: 
#macro_spot_wat
unique_macro_spot_wat<-data%>%
  filter(Matrix == "Water",`Sample.Type` == "Macroplastics (> 5 mm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_macro_spot_wat$`Unique.Sample.Identifier`)
# 72 samples 

#macro_spot_sed
unique_macro_spot_sed<-data%>%
  filter(Matrix == "Sediment",`Sample.Type` == "Macroplastics (> 5 mm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_macro_spot_sed$`Unique.Sample.Identifier`)
# 61 samples

#micro_spot_wat
unique_micro_spot_wat<-data%>%
  filter(Matrix == "Water",`Sample.Type` == "Microplastics (> 100 µm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_micro_spot_wat$`Unique.Sample.Identifier`)
#64 samples

#micro_spot_sed
unique_micro_spot_sed<-data%>%
  filter(Matrix == "Sediment",`Sample.Type` == "Microplastics (> 100 µm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_micro_spot_sed$`Unique.Sample.Identifier`)
# 27 samples 


########################
#data split
########################

##Sediment
data_full_SED<-data%>%
  filter(`Matrix` == "Sediment")
##Check
unique(data_full_SED$Matrix)


#water
data_full_WAT<-data%>%
  filter(`Matrix` == "Water")
##Check
unique(data_full_WAT$Matrix)

##subsets based on matrix and sample type (do contain both spot and tidal data!)
data_full_SED_micro<-data_full_SED%>%
  filter(`Sample.Type` == "Microplastics (> 100 µm)")
length(unique(data_full_SED_micro$`Unique.Sample.Identifier`))

data_full_WAT_micro<- data_full_WAT%>%
  filter(`Sample.Type` == "Microplastics (> 100 µm)") 
length(unique(data_full_WAT_micro$`Unique.Sample.Identifier`))

data_full_SED_macro<-data_full_SED%>%
  filter(`Sample.Type` == "Macroplastics (> 5 mm)")
length(unique(data_full_SED_macro$`Unique.Sample.Identifier`))

data_full_WAT_macro<- data_full_WAT%>%
  filter(`Sample.Type` == "Macroplastics (> 5 mm)") 
length(unique(data_full_WAT_macro$`Unique.Sample.Identifier`))



############################################################
######Dataset creation spot sampling
############################################################

unique(data_full_WAT$`Type.of.Campaign`) #tidal and spot sampling
unique(data_full_SED$`Type.of.Campaign`) #only spot sampling

#data_full_SED_micro ==> no need to change as only spot samples are collected
#data_full_SED_macro ==> no need to change as only spot samples are collected



#data_full_WAT_micro_spot
data_full_WAT_micro_spot<-data_full_WAT_micro%>%
  filter(data_full_WAT_micro$`Type.of.Campaign` == "Spotsampling")
length(unique(data_full_WAT_micro_spot$ `Unique.Sample.Identifier`))
#64 samples

#data_full_WAT_macro_spot
data_full_WAT_macro_spot<-data_full_WAT_macro%>%
  filter(data_full_WAT_macro$`Type.of.Campaign` == "Spotsampling")
length(unique(data_full_WAT_macro_spot$ `Unique.Sample.Identifier`))
#72 samples 


length(unique(data_full_SED_macro$ `Unique.Sample.Identifier`))#61 samples
length(unique(data_full_SED_micro$ `Unique.Sample.Identifier`))#27 samples




########################
#save datasets
########################
write.csv(data_full_SED_micro, "Final analysis/Final dataset/Datasplits/data_full_SED_micro.csv")

write.csv(data_full_WAT_micro_spot, "Final analysis/Final dataset/Datasplits/data_full_WAT_micro_spot.csv")

write.csv(data_full_SED_macro, "Final analysis/Final dataset/Datasplits/data_full_SED_macro.csv")

write.csv(data_full_WAT_macro_spot, "Final analysis/Final dataset/Datasplits/data_full_WAT_macro_spot.csv")
