########################
#libraries
########################

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

##particle information #werken met de 'Particle data sheet' ipv 'Integrated dataset'
plastic<-as.data.frame(read_xlsx("PLUXIN-FinalAnalysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Particle data"))

##Sample metadata
sample_metadata<-as.data.frame(read_xlsx("PLUXIN-FinalAnalysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Sampling overview"))


#campaign metadata
campaign_metadata<-as.data.frame(read_xlsx("PLUXIN-FinalAnalysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Campaigns overview"))

##adding sampling area and alternative names of sampling locations 
SamplingArea<-as.data.frame(read_xlsx("PLUXIN-FinalAnalysis/Raw data/SamplingAreas.xlsx"))





########################
# Overview datasets
########################

##################### Information on the campaigns

#number of campaigns
length(unique(campaign_metadata$`Code Campaign`)) 
##51 campagnes

##number of tidal cycle measurements
sum(campaign_metadata$`Type of Campaign` == "Tidal cycle" ) 
##13 tidal campaigns

##number of spot sampling campaigns
sum(campaign_metadata$`Type of Campaign` == "Spotsampling")
#38 spot samplings

#sampling areas 
unique(campaign_metadata$`Sampling Area`)



##################### Information on the samples

#list of samples (can be used to check the data)
Summary_Samples<-sample_metadata%>%
  group_by(`Sample ID`, `Unique Sample Identifier`)%>%
  summarize("Aantal stalen" = n())
unique(Summary_Samples$`Aantal stalen`)
#all sample identifiers are unique

##Sampling locations 
sort(unique(sample_metadata$`Sampling Location`))

##number of samples
length(unique(sample_metadata$`Unique Sample Identifier`)) #355 samples

##number of sediment samples
sum(sample_metadata$Matrix == "Sediment") ##88
##number of water samples
sum(sample_metadata$Matrix == "Water") #267

##number of macroplastic samples
sum(sample_metadata$`Sample Type` == "Macroplastics (> 5 mm)") #191 samples

##number of microplastic samples  #164 samples
sum(sample_metadata$`Sample Type` == "Microplastics (> 100 µm)")

##number of sediment macro samples
sum(sample_metadata$Matrix == "Sediment" & sample_metadata$`Sample Type` == "Macroplastics (> 5 mm)") #61

##number of sediment micro samples
sum(sample_metadata$Matrix == "Sediment" & sample_metadata$`Sample Type` == "Microplastics (> 100 µm)") #27

##number of water macro samples
sum(sample_metadata$Matrix == "Water" & sample_metadata$`Sample Type` == "Macroplastics (> 5 mm)") #130

##number of water micro samples
sum(sample_metadata$Matrix == "Water" & sample_metadata$`Sample Type` == "Microplastics (> 100 µm)") #137


##check tidal campaigns and number of samples
sample_campaign<-sample_metadata%>%
  left_join(campaign_metadata, by = c("Code Campaign" = "Code Campaign"))

#number of macroplastic samples in spot samples
sum(sample_campaign$`Type of Campaign` == "Spotsampling" & sample_campaign$`Sample Type` == "Macroplastics (> 5 mm)") #133 samples

#number of macroplastic samples in spot samples in water
sum(sample_campaign$`Type of Campaign` == "Spotsampling" & sample_campaign$`Sample Type` == "Macroplastics (> 5 mm)" & sample_campaign$Matrix == "Water") #72 samples

#number of macroplastic samples in spot samples in sediment
sum(sample_campaign$`Type of Campaign` == "Spotsampling" & sample_campaign$`Sample Type` == "Macroplastics (> 5 mm)" & sample_campaign$Matrix == "Sediment") #61 samples

#number of microplastic samples in spot samples
sum(sample_campaign$`Type of Campaign` == "Spotsampling" & sample_campaign$`Sample Type` == "Microplastics (> 100 µm)") #91 samples

#number of microplastic samples in spot samples in water
sum(sample_campaign$`Type of Campaign` == "Spotsampling" & sample_campaign$`Sample Type` == "Microplastics (> 100 µm)" & sample_campaign$Matrix == "Water") #64 samples

#number of microplastic samples in spot samples in sediment
sum(sample_campaign$`Type of Campaign` == "Spotsampling" & sample_campaign$`Sample Type` == "Microplastics (> 100 µm)" & sample_campaign$Matrix == "Sediment") #27 samples


##################### Information from the Particle data (microplastic)

Summary_microplasticSamples<-plastic%>%
  group_by(`Sample ID`, `Unique Sample Identifier`)%>%
  summarize("Aantal particles" = n())

#total number of particles
length(plastic$`Unique Particle Identifier`)
#62116
length(unique(plastic$`Unique Particle Identifier`))
#62112 ==> aantal Unique particle identifiers are not unique!! ==> to be solved for final dataset


#number of samples #expected 355
length(unique(plastic$`Unique Sample Identifier`)) #only 235 samples are present ==> but if a sample doesn't contain any plastics ==>they will be removed ==> not necessarily correct here 
#235 of 355 samples contained plastic ==> 66% of the samples

##number of microplastic particles
sum(plastic$`Sample Type` == "Microplastics (> 100 µm)")
#61081 of the 62116 particles are microplastics ==> 98% of the particles are microplastics

##number of macroplastic particles
sum(plastic$`Sample Type` == "Macroplastics (> 5 mm)")
#only 1035 or 2% macroplastics



################################################
##Creating a full dataset
################################################
#use of "Unique Sample identifier" to merge with particle data
sample_metadata2<-sample_metadata%>% 
  select(-c("LongStart", "LatStart", "LatStop", "LongStop", "Sample Type", "Sample ID","Matrix" , 
            "High tide start (hh:mm)" ,"High tide end (hh:mm)", "Low tide start(hh:mm)", "Low tide end(hh:mm)", "Remarks", "Project"))


##preparing sampling area
SamplingArea<-SamplingArea%>%select(c("Unique Sample Identifier",  "Sampling Location specific", "Sampling Area general"))

#preparing plastic dataset
plastic2<-plastic%>%
  select(-c("Matrix","Project","Sample Type", "Sample ID"))



##Addition of all samples (based on Sample_metdata) to the particle_full dataset
data_full<-sample_metadata%>%
  left_join(plastic2, by = "Unique Sample Identifier")%>%
  left_join(campaign_metadata, by="Code Campaign")%>%
  left_join(SamplingArea, by="Unique Sample Identifier")



#expected to keep 62116 observations of the microplastic dataset + 120 rows of samples without plastics ==> 62236 observations
# have 47 columns (11 of microplastic +31 +5 + 3 ==> alles -3)


##microplastic ==> has 62116 obs 11 var
dim(plastic2)
##sample_metadata2==> 355 obs 31 variables
dim(sample_metadata)
##campaign_metadata==> 51 obs 5 variables
dim(campaign_metadata)
#sampling area ==> 355 rows, 3 columns
dim(SamplingArea)
##data_full ==> has 62236 obs 47 var
dim(data_full)


#check
#all samples 
length(unique(data_full$`Unique Sample Identifier`)) #==> 355 samples included in the data_full dataset


#how many samples in each category: 
#check how many samples in each category contains plastics: 
#macro_spot_wat
unique_macro_spot_wat<-data_full%>%
  filter(Matrix == "Water",`Sample Type` == "Macroplastics (> 5 mm)", `Type of Campaign` == "Spotsampling")%>%
  distinct(`Unique Sample Identifier`)
length(unique_macro_spot_wat$`Unique Sample Identifier`)
# 72 samples 

#macro_spot_sed
unique_macro_spot_sed<-data_full%>%
  filter(Matrix == "Sediment",`Sample Type` == "Macroplastics (> 5 mm)", `Type of Campaign` == "Spotsampling")%>%
  distinct(`Unique Sample Identifier`)
length(unique_macro_spot_sed$`Unique Sample Identifier`)
# 61 samples

#micro_spot_wat
unique_micro_spot_wat<-data_full%>%
  filter(Matrix == "Water",`Sample Type` == "Microplastics (> 100 µm)", `Type of Campaign` == "Spotsampling")%>%
  distinct(`Unique Sample Identifier`)
length(unique_micro_spot_wat$`Unique Sample Identifier`)
#64 samples

#micro_spot_sed
unique_micro_spot_sed<-data_full%>%
  filter(Matrix == "Sediment",`Sample Type` == "Microplastics (> 100 µm)", `Type of Campaign` == "Spotsampling")%>%
  distinct(`Unique Sample Identifier`)
length(unique_micro_spot_sed$`Unique Sample Identifier`)
# 27 samples 



######################################################
# Necessary additions or changes in the data 
######################################################



##checking classes of variables
data_full$Match<-as.numeric(data_full$Match)
data_full$`Lenght (mm)`<-as.numeric(data_full$`Lenght (mm)`)
data_full$`Width (mm)`<-as.numeric(data_full$`Width (mm)`)
data_full$`Maximal dimension (mm)`<-as.numeric(data_full$`Maximal dimension (mm)`)
data_full$`Weight (mg)`<-as.numeric(data_full$`Weight (mg)`)
data_full$`Volume (L)` <- as.numeric(data_full$`Volume (L)`, na.rm=T)
data_full$`DW  Sediment (kg)`<-as.numeric(data_full$`DW  Sediment (kg)`)
data_full$`Particle count`<-as.numeric(data_full$`Particle count`)
data_full$`Total plastic mass (kg)`<-as.numeric(data_full$`Total plastic mass (kg)`)
data_full$`Concentration (g plastic/kg dry weight sediment)`<-as.numeric(data_full$`Concentration (g plastic/kg dry weight sediment)`)
data_full$`Concentration (g plastic/L water)`<-as.numeric(data_full$`Concentration (g plastic/L water)`)
data_full$`Depth Sample (m)`<-as.numeric(data_full$`Depth Sample (m)`)



###toevoegen parameters plastics
data_full<-data_full%>%
  mutate(`LengthWidthratio`= `Lenght (mm)` / `Width (mm)`)%>%
  mutate(`estimated Hight (mm)` = (`Width (mm)`)^2 / (`Lenght (mm)`))%>% #koelmans
  mutate(CSF= `estimated Hight (mm)` / sqrt(`Width (mm)` * `Lenght (mm)`))%>% #koelmans ellipsoïde (CSF<1) of naar een sphere gaan (CSF=1): 
  mutate(Shape_est = ifelse(`LengthWidthratio` <0.375, "filaments", "Non-filaments")) ##cut-off 0.375 based on Koelmans publication



## adding discrete size classes
summary(data_full$`Maximal dimension (mm)`)
#sizes go between 0.1 to 327 mm

##size classes: 
##SC1=0.1mm-0.5mm (H0300011)
##SC2=0.5-1.5mm (H0300012)
##SC3= 1.5-5mm (H0300013)
##SC4=5-15mm (H0300014)
##SC5=15-50mm(H0300015)
##SC6=50-100mm (Not defined in NERC)
##SC7= >100mm (Not defined in NERC)

data_full<-data_full%>%
  mutate(SizeClass= ifelse(`Maximal dimension (mm)`<0.5, "SC1", 
                           ifelse(`Maximal dimension (mm)`>=0.5 & `Maximal dimension (mm)`<1.5, "SC2", 
                                  ifelse(`Maximal dimension (mm)`>=1.5 & `Maximal dimension (mm)`<4.99, "SC3", 
                                         ifelse(`Maximal dimension (mm)`>=5 & `Maximal dimension (mm)`<14.99, "SC4",
                                                ifelse(`Maximal dimension (mm)`>=15 & `Maximal dimension (mm)`<49.99, "SC5",
                                                       ifelse(`Maximal dimension (mm)`>=50 & `Maximal dimension (mm)`<99.99, "SC6",
                                                              ifelse(`Maximal dimension (mm)`>=100 , "SC7", NA))))))))

##if Shape_obs is pellet, the particle is classified as SC3 
data_full$SizeClass<-ifelse(data_full$Shape_obs=="Pellet"| data_full$Shape_obs=="PELLET", "SC3", data_full$SizeClass)
#Check
sum(is.na(data_full$SizeClass)) #122 partikels zonder size class
data_full$Shape_obs<-as.character(as.factor(data_full$Shape_obs))



##Transparency classes NERC: 
data_full$Transpar_NERC<-ifelse(data_full$Color == "TRANSPARENT" | 
                                  data_full$Color =="TRANSPARENT AND YELLOW"|
                                  data_full$Color =="TRANSPARENT AND BLUE"|
                                  data_full$Color =="TRANSPARENT AND GREEN"|
                                  data_full$Color =="TRANSP"|
                                  data_full$Color =="TRANSPARANT"|
                                  data_full$Color =="Transparant", 
                                    "transparent/translucent", 
                                    ifelse(!is.na(data_full$Color), "opaque", NA))

#if sample type is microplastic, the transpar_NERC should be opaque
data_full$Transpar_NERC<-ifelse(data_full$`Sample Type`=="Microplastics (> 100 µm)", "opaque", data_full$Transpar_NERC)
#some NA values as no color is mentioned


#Color NERC:
data_full$Color_NERC<-ifelse(data_full$Color== "TRANSPARENT" |
                                   data_full$Color== "TRANSP" |
                                   data_full$Color== "TRANSPARANT" |
                                   data_full$Color== "Transparant","COLOURLESS",
                                 ifelse(data_full$Color== "ORANGE","ORANGE",
                                        ifelse(data_full$Color== "WHITE" |
                                                 data_full$Color== "SILVER", "WHITE/CREAM",
                                               ifelse(data_full$Color== "BLUE" |
                                                        data_full$Color== "TRANSPARENT AND BLUE", "BLUE",
                                                      ifelse(data_full$Color== "GREY" |
                                                               data_full$Color== "BLACK" |
                                                               data_full$Color=="BLACK AND GREY", "BLACK/GREY",
                                                             ifelse(data_full$Color== "RED", "RED",
                                                                    ifelse(data_full$Color== "GREEN" |
                                                                             data_full$Color== "TRANSPARENT AND GREEN", "GREEN", 
                                                                           ifelse (data_full$Color == "YELLOW" |
                                                                                     data_full$Color == "TRANSPARENT AND YELLOW", "YELLOW",
                                                                                   ifelse(data_full$Color== "BROWN", "BROWN",
                                                                                          ifelse(data_full$Color== "PURPLE", "PURPLE",
                                                                                                 ifelse(data_full$Color== "PINK", "PINK",
                                                                                                        ifelse(data_full$Color== "NA"|
                                                                                                                 data_full$Color== "ND", NA, "MULTICOLOUR"))))))))))))

     
sum(is.na(data_full$Color_NERC))

## Color_NERC_Interreg
data_full$Color_NERC_Interreg<-ifelse(is.na(data_full$Color_NERC), "UNDEFINED",
                                          ifelse(data_full$Color_NERC == "BLUE" |data_full$Color_NERC == "GREEN", "BLUE/GREEN",
                                                 ifelse(data_full$Color_NERC == "ORANGE" | data_full$Color_NERC == "RED" | data_full$Color_NERC == "PURPLE" | data_full$Color_NERC == "PINK", "ORANGE/RED/PURPLE/PINK", data_full$Color_NERC)))
unique(data_full$Color_NERC_Interreg)
sum(is.na(data_full$Color_NERC_Interreg))


##Shape_NERC

data_full$Shape_NERC<-ifelse(data_full$Shape_obs == "Foil" |
                                   data_full$Shape_obs == "FOIL", "films",
                                 ifelse(data_full$Shape_obs == "Fiber" |
                                          data_full$Shape_obs == "FIBER", "filaments",
                                        ifelse(data_full$Shape_obs == "Pellet" |
                                                 data_full$Shape_obs == "PELLET", "pellets",
                                               ifelse(data_full$Shape_obs == "Foam", "foams",
                                                      ifelse(data_full$Shape_obs == "Styrofoam", "foams",
                                                             ifelse(data_full$Shape_obs == "Hard plastic" |
                                                                      data_full$Shape_obs == "HARD BLASTIC" |
                                                                      data_full$Shape_obs == "HARD PLASTIK" |
                                                                      data_full$Shape_obs == "HARD PLASTIC"|
                                                                      data_full$Shape_obs == "Fragment", "fragments",
                                                                    ifelse(data_full$Shape_obs == "Sphere", "granules", 
                                                                           ifelse(data_full$Shape_obs == "band-aid?", "J211",
                                                                                  ifelse(data_full$Shape_obs == "Bottle", "J8",
                                                                                         ifelse(data_full$Shape_obs == "Bubble wrap", "J67",
                                                                                                ifelse(data_full$Shape_obs == "Label Waterbottle", "films",
                                                                                                       ifelse(data_full$Shape_obs == "Sticker ('Bollo')", "films",
                                                                                                              ifelse(data_full$Shape_obs == "Packaging wash product ('Calgon')", "J9",
                                                                                                                     ifelse(data_full$Shape_obs == "Sticker ('Bollo')", "films",
                                                                                                                            ifelse(data_full$Shape_obs == "Clothing etiquette hanger", "fragments",
                                                                                                                                   ifelse(data_full$Shape_obs == "Bottle cap", "J23",
                                                                                                                                          ifelse(data_full$Shape_obs == "Ring", "J24","Undefined micro-litter items")))))))))))))))))



#Shape_tot
#if shape_obs is NA then use Shape_Est otherwise use shape_NERC
data_full$Shape_tot<- ifelse(is.na(data_full$Shape_obs), data_full$Shape_est, data_full$Shape_NERC)


#Polymer_NERC

data_full$Polymer_NERC<-ifelse(data_full$Polymer == "Polypropyleen", "polypropylene", 
                                   ifelse(data_full$Polymer == "Polystyreen"| data_full$Polymer == "Polystyrene", "polystyrene", 
                                          ifelse(data_full$Polymer == "Polyethyleen tereftalate" |data_full$Polymer == "Polyester", "polyester", 
                                                 ifelse(data_full$Polymer == "Polyamide", "polyamide (nylon)", 
                                                        ifelse(data_full$Polymer == "Polyurethaan", "polyurethane", 
                                                               ifelse(data_full$Polymer == "Polyvinyl Alcohol", "polyvinyl alcohol", 
                                                                      ifelse(data_full$Polymer == "Polyvinyl Chloride", "polychlorinated polymer", 
                                                                             ifelse(data_full$Polymer == "PBMA" | data_full$Polymer == "PMA", "polymethylacrylate", 
                                                                                    ifelse(data_full$Polymer == "Styrene-acrylonitrile copolymer", "acrylonitrile butadiene styrene", 
                                                                                           ifelse(data_full$Polymer == "Ethylene-vinyl-acetate", "ethylene-vinyl-acetate", 
                                                                                                  ifelse(data_full$Polymer == "Cellophane", "cellophane", 
                                                                                                         ifelse(data_full$Polymer == "Sodium sterate", "sodium sterate", 
                                                                                                                ifelse(data_full$Polymer == "Polybutadiene", "polybutadiene", 
                                                                                                                       ifelse(data_full$Polymer == "Polyacrilamide", "polyacrilamide", 
                                                                                                                              ifelse(data_full$Polymer == "Polyether urethane-polypropylene oxide - methylene", "polyether urethane-polypropylene oxide - methylene", 
                                                                                                                                     ifelse(data_full$Polymer == "Polyethyleen"| data_full$Polymer == "Ethyleen-Propyleen-Dieen-Monomeer" |data_full$Polymer == "Ethylene.propylene copolymer" |data_full$Polymer == "Monomer Ethylene-Propylene-Diene", "polyethylene", "undefined plastic"))))))))))))))))









################################################
## Clean-up of the data
################################################


# Removing columns
data_full<-data_full%>%
  select(-c("Sample ID", "Instrument", "Opening size (cm)",  "Mesh size" ,"Campaign Date start", "Campaign Date end","LongStart", "LatStart", "LatStop",
            "LongStop","High tide start (hh:mm)" ,"High tide end (hh:mm)", "Low tide start(hh:mm)", "Low tide end(hh:mm)", "Remarks", "Project"))

#T1 and T2 have the wrong date but correct time. ==> for calculation of Diff_time is no problem
##adding time difference for length of mantatrawl (T2-T1 in minutes)
data_full<- data_full%>%
  mutate(Diff_Time = as.numeric(difftime(`T2 (STOP)`, `T1 (START)`, units="mins")))





################################################
## Fix some problems in colnames
################################################

colnames(data_full) <- make.names(colnames(data_full))



################################################
## Save the data
################################################

write.csv(data_full, "PLUXIN-FinalAnalysis/Final dataset/data_full_cleaned.csv")

