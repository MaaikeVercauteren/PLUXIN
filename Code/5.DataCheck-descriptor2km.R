########################
#libraries
########################

library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(stringr)
library(ggsci) #voor kleurpalet grafieken
library(lubridate)
library(rlang)
library(scales)


#setting encoding to UTF-8 to avoid some problems with 'µ'
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

##lay-out grafieken
gglayer_theme<-list(
  theme_classic(), 
  scale_fill_npg(),
  scale_color_npg(),
  theme(legend.position = "top", text= element_text(size=15)))


########################
#datasets
########################

##1. Location specific descriptors + Anthropogenic descriptors 
LocSpecific_metadata<-as.data.frame(read_xlsx("Final analysis/Raw data/FINAL_spatial-2km.xlsx"))

##2. population densities (of each sample and in a 2km radius)
Pop_metadata<-as.data.frame(read_xlsx("Final analysis/Raw data/FINAL_PopulationIDs-2km.xlsx"))
unique(Pop_metadata$`Unique Sample Identifier`)


##3. Temporal descriptors
TempSpecific_metadata<-as.data.frame(read_xlsx("Final analysis/Raw data/Weather1.xlsx"))


##4. Sample metadata
sample_metadata<-as.data.frame(read_xlsx("Final analysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Sampling overview"))




########################
#dataset clean up
########################
##LocSpecific_metadata can be merged immediately with the sample_metadata dataset
#Pop_metadata can be merged immediately with the sample_metadata dataset

#TempSpecific_metadata needs a bit of clean-up as it cannot be merged like this
#removal of unncessary columns
colnames(TempSpecific_metadata)

TempSpecific_metadata<-TempSpecific_metadata%>%
  select(-c("...1","set","code.x","datestamp.x","y.x","x.x","Country.x", "code.y","y.y", "x.y","Country.y","datestamp.y"))
TempSpecific_metadata$tot_precip<-as.numeric(TempSpecific_metadata$tot_precip)
TempSpecific_metadata$mean_temp<-as.numeric(TempSpecific_metadata$mean_temp)
TempSpecific_metadata$mean_windspeed<-as.numeric(TempSpecific_metadata$mean_windspeed)
TempSpecific_metadata$mean_winddirection<-as.numeric(TempSpecific_metadata$mean_winddirection)
TempSpecific_metadata$mean_pressure<-as.numeric(TempSpecific_metadata$mean_pressure)
TempSpecific_metadata$mean_cloudiness<-as.numeric(TempSpecific_metadata$mean_cloudiness)

nrow(TempSpecific_metadata)
#For each sample (unique sample identifier) data on day of sampling 
#for each sample a combination of the information of 3 days before
#percipitation: total of first three timepoints & mean added 
#rest of the variables: mean 

TempSpecific_metadata<-TempSpecific_metadata%>%
  group_by(`USI`) %>%
  mutate_at(vars(1:6), .funs= list(mean = ~mean(., na.rm=T)))%>%
  mutate(tot_precip_TOT=sum(tot_precip, na.rm=T))%>%
  filter(TIME == 0)%>%
  rename(`Unique Sample Identifier` = `USI`)%>%
  rename(`Sampling Location` = `NAME`)
nrow(TempSpecific_metadata)






########################
#merge datasets
########################

Descriptor<-sample_metadata%>%
  left_join(LocSpecific_metadata, by = "Sampling Location")%>%
  left_join(Pop_metadata, by = "Unique Sample Identifier")%>%
  left_join(TempSpecific_metadata, by = "Unique Sample Identifier")




########################
#check and clean datasets
########################


##for Location specific descriptors
unique(Descriptor$`Sampling Location`[is.na(Descriptor$`agriculture [km²]`)])

##for population densities
unique(Descriptor$`Unique Sample Identifier`[is.na(Descriptor$`pop_dens`)])


##for temporal data: information is missing for 35 samples ==> to be added 
unique(Descriptor$`Unique Sample Identifier`[is.na(Descriptor$`tot_precip`)])
#Missing data: 
#"MIC_W_156" "MIC_W_158" "MIC_W_159" ==>ok, not all data available. 


##removing unnecessary columns
Descriptor_all<-Descriptor%>%
  select(c("Unique Sample Identifier","Sampling Location.x", "Matrix","Sediment type","Depth Sample (m)","Sampling Location specific", 33:86))%>%
  rename(`Sampling Location` = `Sampling Location.x`)%>%
  select(-c("sample ID","sampling location", "buffer", "year", "Sampling Location.y", "TIME"))

#export the dataset
write.csv(Descriptor_all, "Final analysis/Final dataset/Descriptor/Descriptor_all_2km_final.csv")





################################################
#Quality check for Descriptor_all dataset
################################################

Descriptor_all<-as.data.frame(read_csv("Final analysis/Final dataset/Descriptor/Descriptor_all_2km_final.csv"))


#1. Correct classes? 
str(Descriptor_all)
#changing land use to numeric
Descriptor_all$`agriculture [km²]`<-as.numeric(Descriptor_all$`agriculture [km²]`)
Descriptor_all$`industry  [km²]`<-as.numeric(Descriptor_all$`industry  [km²]`)
Descriptor_all$`transport  [km²]`<-as.numeric(Descriptor_all$`transport  [km²]`)
Descriptor_all$`urban  [km²]`<-as.numeric(Descriptor_all$`urban  [km²]`)
Descriptor_all$`water  [km²]`<-as.numeric(Descriptor_all$`water  [km²]`)
Descriptor_all$`nature  [km²]`<-as.numeric(Descriptor_all$`nature  [km²]`)
Descriptor_all$`recreation  [km²]`<-as.numeric(Descriptor_all$`recreation  [km²]`)
Descriptor_all$`waste  [km²]`<-as.numeric(Descriptor_all$`waste  [km²]`) 


#2. Missing values?

#general sampling information without missing values: 
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Sampling Location`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Matrix`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Sampling Location specific`)])

sum(is.na(Descriptor_all$`Unique Sample Identifier`))

#no missing values in the general sampling information

#Descriptors without missing values: 

unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Sediment type`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Depth Sample (m)`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Nearby vegetation`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`radius (degree)`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`radius (km)`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`area (km²)`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`RWZI [nr]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Waste facilities [nr]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`agriculture [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`industry  [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`meandering`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`transport  [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`urban  [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`waste  [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`water  [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`human foot print`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`ecotope`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`km² at flood risk`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`Pop_AVG `)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`pop_dens`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_temp_mean `)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`tidal data`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`nature  [km²]`)]) 
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`recreation  [km²]`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`NaturalBank`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_windspeed`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_winddirection`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_windspeed_mean`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_winddirection_mean`)]) 
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`width`)])
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`shortest distance from shore`)])


#To be checked:

#data from some samples missing (see before); but all not available so ok
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`tot_precip`)]) #"MIC_W_156" "MIC_W_158" "MIC_W_159"
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_temp`)]) #"MIC_W_156" "MIC_W_158" "MIC_W_159"

unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_pressure`)])#"MIC_W_156" "MIC_W_158" "MIC_W_159"
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`tot_precip_mean`)]) #"MIC_W_156" "MIC_W_158" "MIC_W_159"

unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_pressure_mean`)]) # "MIC_W_156" "MIC_W_158" "MIC_W_159"

#not always available but ok
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_cloudiness`)]) #lot of NA's
unique(Descriptor_all$`Unique Sample Identifier`[is.na(Descriptor_all$`mean_cloudiness_mean`)]) #lot of NA's



#3. outliers in any of the continuous values?
summary(Descriptor_all$`Depth Sample (m)`) #ranging 1 and 27 seems ok
summary(Descriptor_all$`Depth river`) #ranging between 70 and 2796 (cm) with a lot of NA's
summary(Descriptor_all$`Mean slope`) #ranging between 0.0 and 0.69 (unit?)
summary(Descriptor_all$`radius (degree)`)#all the same (0.025) correct
summary(Descriptor_all$`radius (km)`)#all the same (2.185) correct
summary(Descriptor_all$`area (km²)`)#all the same (15) correct
summary(Descriptor_all$`RWZI [nr]`)#ranging between 0 and 3, ok
ggplot(Descriptor_all, aes(x = `RWZI [nr]`)) + geom_histogram(bins = 20)
summary(Descriptor_all$`Waste facilities [nr]`) #ranging between 0 and 4
ggplot(Descriptor_all, aes(x = `Waste facilities [nr]`)) + geom_histogram(bins = 20)

summary(Descriptor_all$`agriculture [km²]`) #range between 0 and 7.248 ok
summary(Descriptor_all$`industry  [km²]`)#range between 0 and 7 ok
summary(Descriptor_all$`transport  [km²]`)#range between 0 and 1.91 ok
summary(Descriptor_all$`urban  [km²]`)#range between 0 and 6 ok
summary(Descriptor_all$`waste  [km²]`)#range between 0 and 1.08 ok
summary(Descriptor_all$`water  [km²]`)#range between 1 and 14 ok
summary(Descriptor_all$`nature  [km²]`) #range between 0 and 2.9 ok
summary(Descriptor_all$`recreation  [km²]`) #range between 0 and 1.41 ok

summary(Descriptor_all$`human foot print`) #ranging 0 to 59 
ggplot(Descriptor_all, aes(x = `human foot print`)) + geom_histogram(bins = 20)

summary(Descriptor_all$`km² at flood risk`) #ranging 0 to 7.02 

summary(Descriptor_all$`pop_dens`) #ranging 0-39023, large range but seems ok
ggplot(Descriptor_all, aes(x = `pop_dens`)) + geom_histogram(bins = 20)#large range, but various different locations, so ok


summary(Descriptor_all$`width`) #ranges between 12-13505  
ggplot(Descriptor_all, aes(x = `width`)) + geom_histogram(bins = 20)#one outlier
unique(Descriptor_all$`Sampling Location`[(Descriptor_all$`width` > 10000)]) 
#westerschelde (seems possible) 
unique(Descriptor_all$`Sampling Location`[(Descriptor_all$`width` < 100)])


summary(Descriptor_all$`shortest distance from shore`) #8.85 and 5282 
unique(Descriptor_all$`Sampling Location`[(Descriptor_all$`shortest distance from shore` > 5000)]) 
#again westerschelde


summary(Descriptor_all$`tot_precip`) #range between 0 and 32, 3 NA's
summary(Descriptor_all$`tot_precip_mean`) #range between 0 and 17, ok
summary(Descriptor_all$`tot_precip_TOT`) #range 0 to 52

summary(Descriptor_all$`mean_temp`) # 1.287 - 23.618==> acceptable range 
summary(Descriptor_all$`mean_temp_mean`) #range -0.6 to 26, acceptable range
unique(Descriptor_all$`Unique Sample Identifier`[(Descriptor_all$`mean_temp_mean` < 0)])
#check if all low temperatures are in winter months

summary(Descriptor_all$`mean_windspeed`) #range between 0.67 and 12.8, OK
summary(Descriptor_all$`mean_windspeed_mean`)#range between 1 and 10

summary(Descriptor_all$`mean_winddirection`) ##between 33 and 340 , ok
summary(Descriptor_all$`mean_winddirection_mean`) #range between 55 and 264, Ok 

summary(Descriptor_all$`mean_pressure`) #range bewteen 991 and 1038, OK
summary(Descriptor_all$`mean_pressure_mean`)#OK

summary(Descriptor_all$`mean_cloudiness`) #range between 0.91 and 8 (part of sky covered (okta)
summary(Descriptor_all$`mean_cloudiness_mean`)#ok


#4. all categorical values correct?
unique(Descriptor_all$`Sediment type`)
unique(Descriptor_all$`Nearby vegetation`)
unique(Descriptor_all$`meandering`)
unique(Descriptor_all$`ecotope`)
unique(Descriptor_all$`NaturalBank`)
#seem correct.

#5.all samples included
unique(Descriptor_all$`Unique Sample Identifier`) == unique(sample_metadata$`Unique Sample Identifier`)
#seems correct

length(unique(Descriptor_all$`Unique Sample Identifier`)) 
length(unique(sample_metadata$`Unique Sample Identifier`)) 

#yes, all samples included

#5. all dates correct
#Checked in Excel and all are fine!




################################################
## Fix some problems in colnames
################################################

colnames(Descriptor_all) <- make.names(colnames(Descriptor_all))




################################################
## Making subsets for water and sediment
################################################
Descriptor_2km_WAT<-Descriptor_all%>%
  filter(Matrix == "Water")

Descriptor_2km_SED<-Descriptor_all%>%
  filter(Matrix == "Sediment")






################################################
#Saving descriptor data
################################################

write.csv(Descriptor_2km_SED, "Final analysis/Final dataset/Descriptor/Descriptor_2km_SED.csv")

write.csv(Descriptor_2km_WAT, "Final analysis/Final dataset/Descriptor/Descriptor_2km_WAT.csv")











###################################################################################################
#########################adding additions to the descriptor data WAT
###################################################################################################



########################
#Datasets
########################

Descriptor_2km_WAT<-read.csv("Final analysis/Final dataset/Descriptor/Descriptor_2km_WAT.csv")
data_full<-read.csv("Final analysis/Final dataset/data_full_quality.csv")

########################
##Additions
########################

#adding proximity to the sea 
Descriptor_WAT_additions<-Descriptor_2km_WAT%>%
  mutate(Proximity_sea = ifelse(`Sampling.Location` == "Port of Oostende (B)"| 
                                  `Sampling.Location` == "Port of Oostende (A)"|
                                  `Sampling.Location` == "Port of Oostende (C)"|
                                  `Sampling.Location` == "Port of Nieuwpoort (A)"|
                                  `Sampling.Location` == "Port of Nieuwpoort (B)"|
                                  `Sampling.Location` == "Port of Nieuwpoort (C)"|
                                  `Sampling.Location` == "Zeebrugge (700)"|
                                  `Sampling.Location` == "River Scheldt (Westerschelde)"|
                                  `Sampling.Location` == "River Scheldt (Breskens)", "Coastal", "Inland"))




#adding port/river variable
Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  mutate(Port_river=ifelse(`Sampling.Location` == "River Scheldt (Antwerpen)"| 
                             `Sampling.Location` == "River Scheldt (Wintam)"| 
                             `Sampling.Location` == "River Scheldt (Doel)"| 
                             `Sampling.Location` == "River Scheldt (Temse)"| 
                             `Sampling.Location` == "River Scheldt (Breskens)"| 
                             `Sampling.Location` == "River Scheldt (Terneuzen)"| 
                             `Sampling.Location` == "River Scheldt (Westerschelde)"| 
                             `Sampling.Location` == "River Scheldt (Kruiningen)"| 
                             `Sampling.Location` == "River Scheldt (Bath)", "River", "Port"))







#Landuse
landuse<-Descriptor_WAT_additions%>%
  select(c("Unique.Sample.Identifier","Sampling.Location", "agriculture..km..",
           "transport...km..","urban...km..", "water...km..","nature...km..",
           "recreation...km..","waste...km.."))%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(Natural= sum(`water...km..`,`nature...km..`),
         NonNatural= sum (`agriculture..km..`,`transport...km..`, `urban...km..`, `recreation...km..`, `waste...km..`), 
         Total=sum(Natural, NonNatural), 
         NaturalPerc= Natural/Total*100,
         NonNaturalPerc= NonNatural/Total*100)%>%
  mutate(Natural_NonNatural= ifelse(NaturalPerc >= 70, "Highly Natural", 
                                    ifelse(NaturalPerc <70 & NaturalPerc >30, "Moderate natural", "Low natural")))
##summary
landuse%>%
  group_by(`Sampling.Location`)%>%
  summarise(unique(Natural_NonNatural))%>%
  print(n=26)

#join to the descriptor dataset
landuse<-landuse%>%
  select(c("Unique.Sample.Identifier","Natural_NonNatural"))  
Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  left_join(landuse, by="Unique.Sample.Identifier")







##Point discharge

"RWZI [nr]"
"Waste facilities [nr]"
"Active overflow"  

unique(Descriptor_WAT_additions$`RWZI..nr.`)
unique(Descriptor_WAT_additions$`Waste.facilities..nr.`)
unique(Descriptor_WAT_additions$`Active.overflow`)

Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(TotalPointDischarge=sum(`RWZI..nr.`,`Waste.facilities..nr.`,`Active.overflow`), 
         PointDischarge=ifelse(TotalPointDischarge>5, "High", 
                               ifelse(TotalPointDischarge>1 & TotalPointDischarge<=5, "Moderate", "Low")))






##PopulationDens
Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(PopulationDensity=ifelse(pop_dens>30000, "High", 
                                  ifelse(pop_dens>10000 & pop_dens<=30000, "Moderate", "Low")))







##HFP
Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(HFP=ifelse(`human.foot.print`>40, "High", 
                    ifelse(`human.foot.print`>20 & `human.foot.print`<=40, "Moderate", "Low")))







##Season
#based on full dataset
season<-data_full%>%
  select(`Unique.Sample.Identifier`,Date)%>%
  mutate(Season = ifelse(month(Date) %in% c(3,4,5), "Spring",
                         ifelse(month(Date) %in% c(6,7,8), "Summer",
                                ifelse(month(Date) %in% c(9,10,11), "Autumn",
                                       ifelse(month(Date) %in% c(12,1,2), "Winter", "Unknown")))))
season<-distinct(season ,`Unique.Sample.Identifier`,.keep_all=TRUE)

Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  left_join(season, by="Unique.Sample.Identifier")



########################
#save dataset
########################

write.csv(Descriptor_WAT_additions, "Final analysis/Final dataset/Descriptor/Descriptor_WAT_withadditions.csv")





###################################################################################################
#########################adding additions to the descriptor data SED
###################################################################################################



########################
#Datasets
########################

Descriptor_2km_SED<-read.csv("Final analysis/Final dataset/Descriptor/Descriptor_2km_SED.csv")
data_full<-read.csv("Final analysis/Final dataset/data_full_quality.csv")

########################
##Additions
########################

#adding proximity to the sea 
Descriptor_SED_additions<-Descriptor_2km_SED%>%
  mutate(Proximity_sea = ifelse(`Sampling.Location` == "Port of Oostende (B)"| 
                                  `Sampling.Location` == "Port of Oostende (A)"|
                                  `Sampling.Location` == "Port of Oostende (C)"|
                                  `Sampling.Location` == "Port of Nieuwpoort (A)"|
                                  `Sampling.Location` == "Port of Nieuwpoort (B)"|
                                  `Sampling.Location` == "Port of Nieuwpoort (C)"|
                                  `Sampling.Location` == "Zeebrugge (700)"|
                                  `Sampling.Location` == "River Scheldt (Westerschelde)"|
                                  `Sampling.Location` == "River Scheldt (Breskens)", "Coastal", "Inland"))




#adding port/river variable
Descriptor_SED_additions<-Descriptor_SED_additions%>%
  mutate(Port_river=ifelse(`Sampling.Location` == "River Scheldt (Antwerpen)"| 
                             `Sampling.Location` == "River Scheldt (Wintam)"| 
                             `Sampling.Location` == "River Scheldt (Doel)"| 
                             `Sampling.Location` == "River Scheldt (Temse)"| 
                             `Sampling.Location` == "River Scheldt (Breskens)"| 
                             `Sampling.Location` == "River Scheldt (Terneuzen)"| 
                             `Sampling.Location` == "River Scheldt (Westerschelde)"| 
                             `Sampling.Location` == "River Scheldt (Kruiningen)"| 
                             `Sampling.Location` == "River Scheldt (Bath)", "River", "Port"))







#Landuse
landuse<-Descriptor_SED_additions%>%
  select(c("Unique.Sample.Identifier","Sampling.Location", "agriculture..km..",
           "transport...km..","urban...km..", "water...km..","nature...km..",
           "recreation...km..","waste...km.."))%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(Natural= sum(`water...km..`,`nature...km..`),
         NonNatural= sum (`agriculture..km..`,`transport...km..`, `urban...km..`, `recreation...km..`, `waste...km..`), 
         Total=sum(Natural, NonNatural), 
         NaturalPerc= Natural/Total*100,
         NonNaturalPerc= NonNatural/Total*100)%>%
  mutate(Natural_NonNatural= ifelse(NaturalPerc >= 70, "Highly Natural", 
                                    ifelse(NaturalPerc <70 & NaturalPerc >30, "Moderate natural", "Low natural")))
##summary
landuse%>%
  group_by(`Sampling.Location`)%>%
  summarise(unique(Natural_NonNatural))%>%
  print(n=26)

#join to the descriptor dataset
landuse<-landuse%>%
  select(c("Unique.Sample.Identifier","Natural_NonNatural"))  
Descriptor_SED_additions<-Descriptor_SED_additions%>%
  left_join(landuse, by="Unique.Sample.Identifier")







##Point discharge

"RWZI [nr]"
"Waste facilities [nr]"
"Active overflow"  

unique(Descriptor_SED_additions$`RWZI..nr.`)
unique(Descriptor_SED_additions$`Waste.facilities..nr.`)
unique(Descriptor_SED_additions$`Active.overflow`)

Descriptor_SED_additions<-Descriptor_SED_additions%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(TotalPointDischarge=sum(`RWZI..nr.`,`Waste.facilities..nr.`,`Active.overflow`), 
         PointDischarge=ifelse(TotalPointDischarge>5, "High", 
                               ifelse(TotalPointDischarge>1 & TotalPointDischarge<=5, "Moderate", "Low")))






##PopulationDens
Descriptor_SED_additions<-Descriptor_SED_additions%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(PopulationDensity=ifelse(pop_dens>30000, "High", 
                                  ifelse(pop_dens>10000 & pop_dens<=30000, "Moderate", "Low")))







##HFP
Descriptor_SED_additions<-Descriptor_SED_additions%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(HFP=ifelse(`human.foot.print`>40, "High", 
                    ifelse(`human.foot.print`>20 & `human.foot.print`<=40, "Moderate", "Low")))







##Season
#based on full dataset
season<-data_full%>%
  select(`Unique.Sample.Identifier`,Date)%>%
  mutate(Season = ifelse(month(Date) %in% c(3,4,5), "Spring",
                         ifelse(month(Date) %in% c(6,7,8), "Summer",
                                ifelse(month(Date) %in% c(9,10,11), "Autumn",
                                       ifelse(month(Date) %in% c(12,1,2), "Winter", "Unknown")))))
season<-distinct(season ,`Unique.Sample.Identifier`,.keep_all=TRUE)

Descriptor_SED_additions<-Descriptor_SED_additions%>%
  left_join(season, by="Unique.Sample.Identifier")



########################
#save dataset
########################

write.csv(Descriptor_SED_additions, "Final analysis/Final dataset/Descriptor/Descriptor_SED_withadditions.csv")



