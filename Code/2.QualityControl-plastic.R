########################
#libraries
########################
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(tidyr)


########################
#dataset
########################
data_full_qual<-read.csv("Final analysis/Final dataset/data_full_cleaned.csv")


################################################
## Check number of samples in each category
################################################
dim(data_full_qual)

#how many samples in each category: 
#check how many samples in each category contains plastics: 
#macro_spot_wat
unique_macro_spot_wat<-data_full_qual%>%
  filter(Matrix == "Water",`Sample.Type` == "Macroplastics (> 5 mm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_macro_spot_wat$`Unique.Sample.Identifier`)
# 72 samples 

#macro_spot_sed
unique_macro_spot_sed<-data_full_qual%>%
  filter(Matrix == "Sediment",`Sample.Type` == "Macroplastics (> 5 mm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_macro_spot_sed$`Unique.Sample.Identifier`)
# 61 samples

#micro_spot_wat
unique_micro_spot_wat<-data_full_qual%>%
  filter(Matrix == "Water",`Sample.Type` == "Microplastics (> 100 µm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_micro_spot_wat$`Unique.Sample.Identifier`)
#64 samples

#micro_spot_sed
unique_micro_spot_sed<-data_full_qual%>%
  filter(Matrix == "Sediment",`Sample.Type` == "Microplastics (> 100 µm)", `Type.of.Campaign` == "Spotsampling")%>%
  distinct(`Unique.Sample.Identifier`)
length(unique_micro_spot_sed$`Unique.Sample.Identifier`)
# 27 samples 


########################################
# Extended quality control of the data
########################################


#1. levels of Replicate (sediment)
unique(data_full_qual$`Replicate..sediment.`)
sum(is.na(data_full_qual$`Replicate..sediment.`))
data_full_qual$`Replicate..sediment.`[data_full_qual$`Replicate..sediment.` =="NA"]<-NA #as "NA" is considered a text value, it might be better to change it to NA (missing value)



##2.polymers identified ==> for all analysis Polymer_NERC should be used.
unique(data_full_qual$Polymer_NERC)
sum(is.na(data_full_qual$Polymer_NERC))#120 samples without plastic have NA as polymer type




#3.terminology color
unique(data_full_qual$Color) 
unique(data_full_qual$Color_NERC)
sum(is.na(data_full_qual$Color_NERC)) # a lot of micropalstics do not have color

unique(subset(data_full_qual, `Sample.Type` == "Macroplastics (> 5 mm)")$Color_NERC)
unique(subset(data_full_qual, `Sample.Type` == "Microplastics (> 100 µm)")$Color_NERC)
#For microplastics: this should be NA as FTIR doesn't measure color; some of the Microplastic samples were measured with the ATR ==> there a color is registered
#For macroplastics: this follows NERC vocabulary #but pellets are not all registered so color might be missing
sum(is.na(subset(data_full_qual, `Sample.Type` == "Microplastics (> 100 µm)")$Color_NERC))
sum(is.na(subset(data_full_qual, `Sample.Type` == "Macroplastics (> 5 mm)")$Color_NERC))




#4. matrix
unique(data_full_qual$Matrix)
#water and sediment ==> correct



#5.project ==> not in dataset anymore



#6.depth 
summary(data_full_qual$`Depth.Sample..m.`)
unique(data_full_qual$`Depth.Sample..m.`)

#As this deals with depth of the sample ==> this should be 1 for water samples and > 1 for sediment samples
unique(subset(data_full_qual, Matrix == "Water")$`Depth.Sample..m.`)
unique(subset(data_full_qual, Matrix == "Sediment")$`Depth.Sample..m.`)
#Correct



#7. Check volume for errors (extremely large or small values)
summary(subset(data_full_qual, Matrix=="Water")$`Volume..L.`)
#no NA's
#between 2092 and 46710 L ==> very high variabilty  

ggplot(subset(data_full_qual, Matrix=="Water"), aes(x=`Volume..L.`))+
  geom_histogram()+
  theme_minimal()
#Correlation Volume and Time of trawl?
ggplot(data_full_qual, aes(x=`Volume..L.`, y=Diff_Time))+
  geom_point()+
  geom_smooth(method = "lm", se=T)+
  theme_minimal()
#small (visual) correlation with the time of trawl 
summary(subset(data_full_qual, Matrix=="Sediment")$`Volume..L.`)
#all NA which is good



#8. Check dry weight for errors (extremely large or small values)
summary(subset(data_full_qual, Matrix=="Sediment")$`DW..Sediment..kg.`)
#from 0.250 to 17 kg==> very large vaiability (but mean is 1.404 so few outliers)
#outliers are known but no reason to suspect that it was wrong  (although really lot of sediment...)
#12775 NA's ==> a lot ==> should be checked 
#==> for microplastics pooled DW should be used which is not added in the sample overview ==> so NA's could be ok  ==> look into how to add this in the datafile!! 

ggplot(subset(data_full_qual, Matrix=="Sediment"), aes(x=`DW..Sediment..kg.`))+
  geom_histogram()+
  theme_minimal()

summary(subset(data_full_qual, Matrix=="Water")$`DW..Sediment..kg.`)
#only NA's which is good




#9. Check length and width for extreme values

summary(data_full_qual$`Lenght..mm.`)

ggplot(data_full_qual, aes(x= Matrix, y=`Lenght..mm.`)) +
  geom_point()+
  theme_minimal()
#length between 0.1 and 327 mm; 327 seems like an outlier however, not clear where the threshold would be 
#326 +120 NA's which is not that much based on >60 000 particles

summary(data_full_qual$`Width..mm.`)

ggplot(data_full_qual, aes(x= Matrix, y=`Width..mm.`)) +
  geom_point()+
  theme_minimal()
#width between 0.01 and 149 which seems logic
#948 NA's which is ok since macro only have length. 



#10. No T2 for sediment samples 
summary(subset(data_full_qual, Matrix=="Sediment")$`T2..STOP.`)
##all sediment samples are NA which is good 


#11. Consistency was checked before so no problems



#12. NA values; should be ok (maybe some descripancy due to missing data in sample_metadata upon joinin, but should have no effect)


#14. Microplastic can never have NA/blank as polymer type except for the 2 samples without plastic
unique(subset(data_full_qual, `Sample.Type` == "Microplastics (> 100 µm)")$Polymer)
sum(is.na(subset(data_full_qual, `Sample.Type` == "Microplastics (> 100 µm)")$Polymer)) #1 undefined/unidentified plastic type
unique(subset(data_full_qual, `Sample.Type` == "Microplastics (> 100 µm)")$Polymer_NERC)
sum(is.na(subset(data_full_qual, `Sample.Type` == "Microplastics (> 100 µm)")$Polymer_NERC)) 
##correct


#15.  Length of manta trawls + sediment trawls should have T2 as NA
unique(subset(data_full_qual, Matrix=="Sediment")$`T2..STOP.`)
#all NA's which is good



#16. length of manta trawls should be logic

summary(data_full_qual$Diff_Time)

ggplot(subset(data_full_qual, Matrix=="Water"), aes(x= Matrix, y=Diff_Time)) +
  geom_point()+
  theme_minimal()
unique(subset(data_full_qual, Diff_Time > 15)$`Unique.Sample.Identifier`)
#one sample (MAC_W_033) with 20 minutes ==> seems like an outlier; other all below 15 minutes 




#17. Check max dimension for extreme values
summary(data_full_qual$`Maximal.dimension..mm.`)

ggplot(data_full_qual, aes(x= Matrix, y=`Maximal.dimension..mm.`)) +
  geom_point()+
  theme_minimal()
#same range as length



#18. Check weight for extreme values
summary(data_full_qual$`Weight..mg.`)

ggplot(data_full_qual, aes(x= Matrix, y=`Weight..mg.`)) +
  geom_point()+ geom_boxplot()+
  theme_minimal() 
#minimum weight 0 ==> not possible==> some pellets have weigth 0 ==> should be NA
data_full_qual$`Weight..mg.`[data_full_qual$`Weight..mg.` == 0] <- NA

#max weight 21453 mg ==> seems like outlier ==> remove??? 
unique(subset(data_full_qual, `Weight..mg.` > 20000)$`Unique.Particle.Identifier`)




#19. Check sampling location & sample area
unique(data_full_qual$`Sampling.Location`)
unique(data_full_qual$`Sampling.Area`)
#ok



#20.check sample type
unique(data_full_qual$`Sample.Type`)
#ok



#21. Check shorter names of sampling locations
unique(data_full_qual$`Sampling.Location.specific`)
#looks ok, no NA's



#22. check sampling area general

unique(data_full_qual$`Sampling.Area.general`)
#looks ok, no NA's



#23. number of campaigns #expect 51 
length(unique(data_full_qual$`Code.Campaign`)) 
#==> correct


########################
#actions
########################
#remove T2 and Diff_time of MAC_W_033

data_full_qual$`T2 (STOP)`[data_full_qual$`Unique.Sample.Identifier` == "MAC_W_033"]<-NA
data_full_qual$`Diff_time`[data_full_qual$`Unique.Sample.Identifier` == "MAC_W_033"]<-NA

#remove weight of particle with weight 21453 mg (MAC_W_029_1)
data_full_qual$`Weight..mg.`[data_full_qual$`Unique.Particle.Identifier` == "MAC_W_029_1"]<-NA


########################
#save datasets
########################
dim(data_full_qual)
write.csv(data_full_qual, "Final analysis/Final dataset/data_full_quality.csv")
