########################
#libraries
########################
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggsci) #voor kleurpalet grafieken

options(scipen = 999, digits = 10)

##lay-out grafieken
gglayer_theme<-list(
  theme_classic(), 
  scale_fill_npg(),
  scale_color_npg(),
  theme(legend.position = "top", text= element_text(size=15)))


########################
#datasets
########################

data_full_SED_micro<- read.csv("Final analysis/Final dataset/Datasplits/data_full_SED_micro.csv")

data_full_WAT_micro_spot<-read.csv( "Final analysis/Final dataset/Datasplits/data_full_WAT_micro_spot.csv")

data_full_SED_macro<- read.csv("Final analysis/Final dataset/Datasplits/data_full_SED_macro.csv")

data_full_WAT_macro_spot<- read.csv("Final analysis/Final dataset/Datasplits/data_full_WAT_macro_spot.csv")


########################
#Merge all water and all sediment samples
########################


data_full_WAT_spot<-rbind(data_full_WAT_micro_spot, data_full_WAT_macro_spot)
data_full_SED_spot<-rbind(data_full_SED_micro, data_full_SED_macro)


#Datasets created: 
#1. all spot samples water: data_full_WAT_spot
#2. all spot samples sediment: data_full_SED_spot
#2. microplastics: data_full_WAT_micro_spot and data_full_SED_micro
#3. macroplastics: data_full_WAT_macro_spot and data_full_SED_macro


########################
#dataset with replicate and dry weight for SED_micro
########################

replicate_SED_metadata<-as.data.frame(read_xlsx("Final analysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Replicate data Sediment"))
## Fix some problems in colnames
colnames(replicate_SED_metadata) <- make.names(colnames(replicate_SED_metadata))

replicate_SED_metadata<-replicate_SED_metadata%>%
  select(c("Unique.Sample.Identifier","Replicate..sediment.", "Sample_replicate", "DW..Sediment..kg."))




############################################################
######Summary samples and number of particles
############################################################


#summary table for number of samples 
summary_spot_WAT<-data_full_WAT_spot%>%
  group_by( `Sample.Type`)%>%
  summarise(n=n_distinct(`Unique.Sample.Identifier`), n_loc=n_distinct(`Sampling.Location`), n_date=n_distinct(`Date`))

summary_spot_SED<-data_full_SED_spot%>%
  group_by( `Sample.Type`)%>%
  summarise(n=n_distinct(`Unique.Sample.Identifier`), n_loc=n_distinct(`Sampling.Location`), n_date=n_distinct(`Date`))


length(unique(data_full_WAT_spot$Sampling.Location.specific))

#number of particles per matrix
sum(!is.na(data_full_WAT_spot$Polymer))
sum(!is.na(data_full_SED_spot$Polymer))

total<-sum(!is.na(data_full_WAT_spot$Polymer)) + sum(!is.na(data_full_SED_spot$Polymer))
total

#number of particles per matrix and size
sum(!is.na(data_full_WAT_micro_spot$Polymer))
sum(!is.na(data_full_WAT_macro_spot$Polymer))
sum(!is.na(data_full_SED_macro$Polymer))
sum(!is.na(data_full_SED_micro$Polymer))

#Samples without polymers
num_samples_WAT_spot <- data_full_WAT_spot %>%
  filter(!is.na(Polymer)) %>%
  distinct(`Unique.Sample.Identifier`) %>%
  nrow()
print(num_samples_WAT_spot)
                                                   
num_samples_SED_spot <- data_full_SED_spot %>%
  filter(!is.na(Polymer)) %>%
  distinct(`Unique.Sample.Identifier`) %>%
  nrow()
print(num_samples_SED_spot)                                                
                                                   
num_samples_WAT_micro_spot <- data_full_WAT_micro_spot %>%
  filter(!is.na(Polymer)) %>%
  distinct(`Unique.Sample.Identifier`) %>%
  nrow()
print(num_samples_WAT_micro_spot)


num_samples_WAT_macro_spot <- data_full_WAT_macro_spot %>%
  filter(!is.na(Polymer)) %>%
  distinct(`Unique.Sample.Identifier`) %>%
  nrow()
print(num_samples_WAT_macro_spot)

num_samples_SED_macro <- data_full_SED_macro %>%
  filter(!is.na(Polymer)) %>%
  distinct(`Unique.Sample.Identifier`) %>%
  nrow()
print(num_samples_SED_macro)


num_samples_SED_micro <- data_full_SED_micro %>%
  filter(!is.na(Polymer)) %>%
  distinct(`Unique.Sample.Identifier`) %>%
  nrow()
print(num_samples_SED_micro)
                                       
                                                   
###########################################################################################################################################################################
##Water 
###########################################################################################################################################################################



#####################################
# Concentration
#####################################



#general plastic concentration calculation
conc_summary_WAT_spot<-data_full_WAT_spot%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", NA, 1))%>%
  group_by(`Unique.Sample.Identifier`, `Sampling.Area.general`,`Sampling.Location`, `Sample.Type`)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE), VolumeTot= mean(`Volume..L.`, na.rm=TRUE))%>%
  mutate(ConcMP= aantalMPs/VolumeTot, ConcMPm3=ConcMP*1000)

#average per micro/macro
conc_summary_WAT_micromacro<-conc_summary_WAT_spot%>%
  group_by(`Sample.Type`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPm3, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPm3, na.rm=T)),
            minMPconc=min(ConcMPm3, na.rm=T), maxMPconc=max(ConcMPm3, na.rm=T))

##Average concentration per Sampling area
conc_summary_WAT_area<-conc_summary_WAT_spot%>%
  group_by(`Sampling.Area.general`,`Sample.Type`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPm3, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPm3, na.rm=T)),
            minMPconc=min(ConcMPm3, na.rm=T), maxMPconc=max(ConcMPm3, na.rm=T), aantal=n())

##Average concentration per Sampling.Location
conc_summary_WAT_location<-conc_summary_WAT_spot%>%
  group_by(`Sampling.Location`,`Sample.Type`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPm3, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPm3, na.rm=T)),
            minMPconc=min(ConcMPm3, na.rm=T), maxMPconc=max(ConcMPm3, na.rm=T))

##graphs

#Sampling.Location
ggplot(conc_summary_WAT_spot, aes(`Sampling.Location`, ConcMPm3)) + geom_boxplot(aes(fill=factor(`Sample.Type`))) + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/m³)", x="Sampling.Location") + labs(fill="Sample.Type")


#sampling area
ggplot(conc_summary_WAT_spot, aes(`Sampling.Area.general`, ConcMPm3)) + geom_boxplot(aes(fill=factor(`Sample.Type`))) + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/m³)", x="Sampling.Location") + labs(fill="Sample.Type")    

##check water volume sampled for macro and micro
ggplot(conc_summary_WAT_spot, aes(x=`Sample.Type`, y= VolumeTot)) + geom_boxplot()+gglayer_theme
#macroplastic generally larger volume linked to larger mesh size and thus less clotting






#####################################
# PM composition
#####################################


#for this we first remove every row with NA for PM type

data_full_WAT_spot_adjusted<-data_full_WAT_spot%>%
  filter(!is.na(Polymer_NERC))


###General PM composition
SummaryPMtypes_WAT_General<-data_full_WAT_spot_adjusted%>%
  group_by(`Sample.Type`, Polymer_NERC)%>%
  summarise(aantalMPs=as.numeric(n()))%>%
  mutate(total=as.numeric(sum(aantalMPs)))%>%
  mutate(freq=as.numeric((aantalMPs/total)*100))
#pie chart
ggplot(SummaryPMtypes_WAT_General, aes("",  y=freq, fill=factor(Polymer_NERC))) + geom_bar(width=1, stat="identity") +
  coord_polar(theta = "y", start=0) + theme_bw()+
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(),legend.title=element_blank(), legend.text = element_text(size=9)) + 
  facet_grid(~`Sample.Type`)
#bar chart
ggplot(SummaryPMtypes_WAT_General, aes(x=`Sample.Type`, y=freq, fill=factor(Polymer_NERC))) + geom_bar(position="stack", stat="identity") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))


### per Sampling.Location
SummaryPMtypes_WAT<-data_full_WAT_spot_adjusted%>%
  group_by(`Sampling.Location`, `Sample.Type` ,Polymer_NERC)%>%
  summarise(aantalMPs=as.numeric(n()))%>%
  mutate(total=as.numeric(sum(aantalMPs)))%>%
  mutate(freq=as.numeric((aantalMPs/total)*100))
##no pie chart due to too many Sampling.Locations

###per sampling area
SummaryPMtypes_WAT2<-data_full_WAT_spot_adjusted%>%
  group_by(`Sampling.Area.general`, `Sample.Type` ,Polymer_NERC)%>%
  summarise(aantalMPs=as.numeric(n()))%>%
  mutate(total=as.numeric(sum(aantalMPs)))%>%
  mutate(freq=as.numeric((aantalMPs/total)*100))

#bar chart
ggplot(SummaryPMtypes_WAT2, aes(x=`Sample.Type`, y=freq, fill=factor(Polymer_NERC))) + geom_bar(position="stack", stat="identity") +
  facet_grid(~`Sampling.Area.general`) + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))

#pie chart
ggplot(SummaryPMtypes_WAT2, aes("",  y=freq, fill=factor(Polymer_NERC))) + geom_bar(width=1, stat="identity") +
  coord_polar(theta = "y", start=0) + theme_bw()+
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(),legend.title=element_blank(), legend.text = element_text(size=4))+
  facet_grid(`Sample.Type`~`Sampling.Area.general`) 




#####################################
# Size distribution
#####################################

##sizes go between 0.1 to 327 mm with very high concentration of microplastics and lower concentration of macroplastics
#With sizeClasses

#use also the data_full_WAT_adjusted dataset (removed NA's)

#General
ggplot(data_full_WAT_spot_adjusted, aes(SizeClass)) + geom_bar(aes(y=..prop..,group=1)) + 
  xlab("Size class")+ylab("Relative frequence ") + scale_y_continuous(labels = scales::percent)+gglayer_theme



##General with continuous variable 
ggplot(data_full_WAT_spot_adjusted, aes(`Lenght..mm.`)) + geom_histogram(aes(y = after_stat(count / sum(count)))) + 
  xlab("Length (mm)")+ylab("Relative frequence ") + scale_y_continuous(labels = scales::percent)+gglayer_theme
#clear over-representation of microplastics (0.1-5 mm) compared to macroplastics (>5 mm)
##General with continuous variable + color per Sample.Type
ggplot(data_full_WAT_spot_adjusted, aes(log(`Lenght..mm.`), fill=`Sample.Type`)) + geom_histogram(aes(y = after_stat(count / sum(count)))) + 
  xlab("Length (mm)")+ylab("Relative frequence ") + scale_y_continuous(labels = scales::percent)+gglayer_theme
#clear over-representation of microplastics (0.1-5 mm) compared to macroplastics (>5 mm)

#summary table per Sample.Type
length_summary<-data_full_WAT_spot_adjusted%>%
  group_by(`Sample.Type`)%>%
  summarise(mean=mean(`Lenght..mm.`, na.rm=T), sd=sd(`Lenght..mm.`, na.rm=T), 
            min=min(`Lenght..mm.`, na.rm=T), max=max(`Lenght..mm.`, na.rm=T))
#summary tabel per Sample.Type and sample area
length_summary_area<-data_full_WAT_spot_adjusted%>%
  group_by(`Sample.Type`, `Sampling.Area.general`)%>%
  summarise(mean=mean(`Lenght..mm.`, na.rm=T), sd=sd(`Lenght..mm.`, na.rm=T), 
            min=min(`Lenght..mm.`, na.rm=T), max=max(`Lenght..mm.`, na.rm=T))


### General with categorical variable (true microplastics (0.1-5 mm, according to definition of MP)) to allow for better visualization as macroplastics are underrepresented
WAT_size<-data.table(data_full_WAT_spot_adjusted)
WAT_size[, sizeGroup:=cut(`Lenght..mm.`, breaks=seq(0.1,5, 0.1), dig.lab=10)] 
WAT_size<-as.data.frame(WAT_size[, table(sizeGroup)])
WAT_size<-WAT_size%>%
  mutate(total=sum(Freq), relFreq=(Freq/total)*100)
#bar plot
ggplot(WAT_size, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity")+gglayer_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("Size class (mm)")+ylab("Relative frequence (%)")


###Size distribution per Sampling.Location

functie<-function(number){
  data_full_size<-data.table(data_full_WAT_spot_adjusted)
  data_full_size[, sizeGroup:=cut(`Lenght..mm.`, breaks=seq(0.1,5, 0.1), dig.lab=10)]  
  a<-subset(data_full_size, `Sampling.Location`  %in% number)
  Datasize<-as.data.frame(a[, table(sizeGroup)]) 
  Datasize<-Datasize%>%
    mutate(`Sampling.Location` = as.factor(number))%>%
    mutate(total=sum(Freq), relFreq=(Freq/total)*100)
  
}

y<-NULL

for(i in unique(data_full_WAT_spot_adjusted$`Sampling.Location`)){
  tmp<- functie(i)
  y<-rbind(y, tmp)
  y<-y%>%
    filter(!is.na(relFreq))
}
unique(y$`Sampling.Location`)


#Facet grid per locations, but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity",  aes(fill=`Sampling.Location`))+facet_grid(~`Sampling.Location`)+
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ylab("Relative frequency (%)")+xlab("Size classes (µm)")


###Size distribution per sampling Area

functie<-function(number){
  data_full_size<-data.table(data_full_WAT_spot_adjusted)
  data_full_size[, sizeGroup:=cut(`Lenght..mm.`, breaks=seq(0.1,3, 0.1), dig.lab=10)]  
  a<-subset(data_full_size, `Sampling.Area.general`  %in% number)
  Datasize<-as.data.frame(a[, table(sizeGroup)]) 
  Datasize<-Datasize%>%
    mutate(`Sampling.Area.general` = as.factor(number))%>%
    mutate(total=sum(Freq), relFreq=(Freq/total)*100)
  
}

y<-NULL
for(i in unique(data_full_WAT_spot_adjusted$`Sampling.Area.general`)){
  tmp<- functie(i)
  y<-rbind(y, tmp)
  y<-y%>%
    filter(!is.na(relFreq))
}
unique(y$`Sampling.Area.general`)

#One bar graph with different colors for locations but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity", position="dodge", aes(fill=`Sampling.Area.general`))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+ ylab("Relative frequency (%)")+xlab("Size classes (µm)")
#Facet grid per locations, but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity",  aes(fill=`Sampling.Area.general`))+facet_grid(~`Sampling.Area.general`)+ theme_bw()+
  theme(axis.text.x=element_text(angle=45, hjust=1, size=4), legend.position = "none") + ylab("Relative frequency (%)")+xlab("Size classes (µm)")






#####################################
# Size distribution - alfa values
#####################################

set.seed(123)
#Caculate general alfa value
## make into continuous powerlaw object
df.pl.full <- conpl$new(na.omit(data_full_WAT_spot_adjusted$`Lenght..mm.`))   
###########1 check power law distribution
bs.p.full <- bootstrap_p(df.pl.full, no_of_sims = 10, threads = 2, xmax = 2E12)
bs.p.full
#if p > 0.1 ==> power law distiribution is possible
########determine Xmin using KS statistics
estimate_xmin(df.pl.full)
######Determine alpha (mle method) with estimated xmin
bs.res.full <- bs.p.full$bootstraps## get results from the bootstrap
bs.res.full
mean(bs.res.full$pars)## gives individual fits (bs.res$pars), take the mean for your overall fit.
sd(bs.res.full$pars)
mean(bs.res.full$xmin)## gives individual fits (bs.res$pars), take the mean for your overall fit.
sd(bs.res.full$xmin)



##calculating the alfa values of the size distributions
##make a function for each Sampling.Location
calculate_alfa <- function(data, variable) {
  results <- data.frame(Location = character(),
                        Mean_alfa = numeric(),
                        SD_alfa = numeric(),
                        stringsAsFactors = FALSE)
  unique_locations <- unique(data$`Sampling.Location`)
  for (location in unique_locations) {
    subset_data <- data[data$`Sampling.Location` == location, ]
    df.pl <- conpl$new(na.omit(subset_data[[variable]]))
    bs <- bootstrap(df.pl, no_of_sims = 10, threads = 2, xmax = 2E12)
    bs.res <- bs$bootstraps
    mean_value <- mean(bs.res$pars)
    sd_value <- sd(bs.res$pars)
    result <- data.frame("Sampling.Location" = location, Mean_alfa = mean_value, SD_alfa = sd_value)
    results <- rbind(results, result)
  }
  return(results)
}

alfa_WAT <- calculate_alfa(data_full_WAT_spot_adjusted, "Lenght..mm.")

##calculating the alfa values of the size distributions per sampling area

##make a function for each sampling area
calculate_alfa_area <- function(data, variable) {
  results <- data.frame(Area = character(),
                        Mean_alfa = numeric(),
                        SD_alfa = numeric(),
                        stringsAsFactors = FALSE)
  unique_Area <- unique(data$`Sampling.Area.general`)
  for (Area in unique_Area) {
    subset_data <- data[data$`Sampling.Area.general` == Area, ]
    df.pl <- conpl$new(na.omit(subset_data[[variable]]))
    bs <- bootstrap(df.pl, no_of_sims = 10, threads = 2, xmax = 2E12)
    bs.res <- bs$bootstraps
    mean_value <- mean(bs.res$pars)
    sd_value <- sd(bs.res$pars)
    result <- data.frame("Sampling.Area.general" = Area, Mean_alfa = mean_value, SD_alfa = sd_value)
    results <- rbind(results, result)
  }
  return(results)
}

alfa_WAT_area <- calculate_alfa_area(data_full_WAT_spot_adjusted, "Lenght..mm.")








###########################################################################################################################################################################
##Sediment
###########################################################################################################################################################################





#####################################
# Concentration
#####################################



#1. Macroplastics
#general plastic concentration calculation
conc_summary_SED_macro<-data_full_SED_macro%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", NA, 1))%>%
  group_by(`Unique.Sample.Identifier`, `Sampling.Area.general`,`Sampling.Location`, `Sample.Type`)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE), DW=max(`DW..Sediment..kg.`, na.rm=T))%>%
  mutate(ConcMPkg= aantalMPs/DW)


#average per macro
conc_summary_SED_macro_total<-conc_summary_SED_macro%>%
  group_by(`Sample.Type`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPkg, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPkg, na.rm=T)),
            minMPconc=min(ConcMPkg, na.rm=T), maxMPconc=max(ConcMPkg, na.rm=T))

##Average concentration per Sampling area
conc_summary_SED_macro_area<-conc_summary_SED_macro%>%
  group_by(`Sampling.Area.general`,`Sample.Type`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPkg, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPkg, na.rm=T)),
            minMPconc=min(ConcMPkg, na.rm=T), maxMPconc=max(ConcMPkg, na.rm=T))

##Average concentration per Sampling.Location
conc_summary_SED_macro_location<-conc_summary_SED_macro%>%
  group_by(`Sampling.Location`,`Sample.Type`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPkg, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPkg, na.rm=T)),
            minMPconc=min(ConcMPkg, na.rm=T), maxMPconc=max(ConcMPkg, na.rm=T))

##graphs

#Sampling.Location
ggplot(conc_summary_SED_macro, aes(`Sampling.Location`, ConcMPkg)) + geom_boxplot(aes(fill=factor(`Sample.Type`))) + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/kg DW)", x="Sampling.Location") + labs(fill="Sample.Type")


#sampling area
ggplot(conc_summary_SED_macro, aes(`Sampling.Area.general`, ConcMPkg)) + geom_boxplot(aes(fill=factor(`Sample.Type`))) + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/m³)", x="Sampling.Location") + labs(fill="Sample.Type")    


#sampling area
ggplot(conc_summary_SED_macro, aes(`Sampling.Area.general`, ConcMPkg)) + geom_boxplot(aes(fill=factor(`Sample.Type`))) + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/m³)", x="Sampling.Location") + labs(fill="Sample.Type")+
  ylim(0,0.01)

ggplot(conc_summary_SED_macro_area, aes(`Sampling.Area.general`, avgMPconc))+ 
  geom_bar(stat="identity")+ 
  geom_errorbar(ymax=conc_summary_SED_macro_area$avgMPconc+conc_summary_SED_macro_area$sdMPconc, ymin=conc_summary_SED_macro_area$avgMPconc-conc_summary_SED_macro_area$sdMPconc, width=0.1) + 
  gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/m³)", x="Sampling.Location")





#2. microplastic



##Concentration per replicate
conc_summary_SED_micro<-data_full_SED_micro%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", NA, 1))%>%
  group_by(`Unique.Sample.Identifier`,`Replicate..sediment.`,`Sampling.Area.general`, `Sampling.Location`)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE))%>%
  left_join(replicate_SED_metadata, by=c("Unique.Sample.Identifier", "Replicate..sediment."))%>%#adding dry weigth based on both sample identifier and replicate
  mutate(concMPkg= aantalMPs/as.numeric(`DW..Sediment..kg.`))
length(unique(conc_summary_SED_micro$'Unique.Sample.Identifier')) ##correct, all samples are included in this new dataset


##concentration per sample
#decision on average or sum per replicate.  ==> average
conc_summary_SED_micro_sample<-conc_summary_SED_micro%>%
  group_by(`Unique.Sample.Identifier`, `Sampling.Area.general`, `Sampling.Location`)%>%
  summarise(ConcMPkg=mean(concMPkg, na.rm=T))

#average concentration microplastics
mean(as.numeric(conc_summary_SED_micro_sample$ConcMPkg), na.rm=T)
sd(as.numeric(conc_summary_SED_micro_sample$ConcMPkg), na.rm=T)
min(as.numeric(conc_summary_SED_micro_sample$ConcMPkg), na.rm=T)
max(as.numeric(conc_summary_SED_micro_sample$ConcMPkg), na.rm=T)


##Average concentration per Sampling area
conc_summary_SED_micro_area<-conc_summary_SED_micro_sample%>%
  group_by(`Sampling.Area.general`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPkg, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPkg, na.rm=T)),
            minMPconc=min(ConcMPkg, na.rm=T), maxMPconc=max(ConcMPkg, na.rm=T))
#plot
ggplot(conc_summary_SED_micro_sample, aes(`Sampling.Area.general`, ConcMPkg)) + geom_boxplot() + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/kg DW)", x="Sampling area") + labs(fill="Sample.Type")


##Average concentration per Sampling.Location
conc_summary_SED_micro_location<-conc_summary_SED_micro_sample%>%
  group_by(`Sampling.Location`)%>%
  summarise(avgMPconc=mean(as.numeric(ConcMPkg, na.rm=T)), sdMPconc=sd(as.numeric(ConcMPkg, na.rm=T)),
            minMPconc=min(ConcMPkg, na.rm=T), maxMPconc=max(ConcMPkg, na.rm=T))
#plot
ggplot(conc_summary_SED_micro_sample, aes(`Sampling.Location`, ConcMPkg)) + geom_boxplot() + gglayer_theme + theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(y="Plastic concentration (MP/kg DW)", x="Sampling.Location") + labs(fill="Sample.Type")






#####################################
# PM composition
#####################################
#macro en micro samen
#for this we first remove every row with NA for PM type
data_full_SED_adjusted<-data_full_SED_spot%>%
  filter(!is.na(Polymer_NERC))


###General PM composition
SummaryPMtypes_SED_General<-data_full_SED_adjusted%>%
  group_by(`Sample.Type`, Polymer_NERC)%>%
  summarise(aantalMPs=as.numeric(n()))%>%
  mutate(total=as.numeric(sum(aantalMPs)))%>%
  mutate(freq=as.numeric((aantalMPs/total)*100))
#pie chart
ggplot(SummaryPMtypes_SED_General, aes("",  y=freq, fill=factor(Polymer_NERC))) + geom_bar(width=1, stat="identity") +
  coord_polar(theta = "y", start=0) + theme_bw()+
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(),legend.title=element_blank(), legend.text = element_text(size=9)) + 
  facet_grid(~`Sample.Type`)
#bar chart
ggplot(SummaryPMtypes_SED_General, aes(x=`Sample.Type`, y=freq, fill=factor(Polymer_NERC))) + geom_bar(position="stack", stat="identity") +
  theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))


### per Sampling.Location
SummaryPMtypes_SED<-data_full_SED_adjusted%>%
  group_by(`Sampling.Location`, `Sample.Type` ,Polymer_NERC)%>%
  summarise(aantalMPs=as.numeric(n()))%>%
  mutate(total=as.numeric(sum(aantalMPs)))%>%
  mutate(freq=as.numeric((aantalMPs/total)*100))
##no pie chart due to too many Sampling.Locations

###per sampling area
SummaryPMtypes_SED2<-data_full_SED_adjusted%>%
  group_by(`Sampling.Area.general`, `Sample.Type` ,Polymer_NERC)%>%
  summarise(aantalMPs=as.numeric(n()))%>%
  mutate(total=as.numeric(sum(aantalMPs)))%>%
  mutate(freq=as.numeric((aantalMPs/total)*100))

#bar chart
ggplot(SummaryPMtypes_SED2, aes(x=`Sample.Type`, y=freq, fill=factor(Polymer_NERC))) + geom_bar(position="stack", stat="identity") +
  facet_grid(~`Sampling.Area.general`) + theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1))

#pie chart
ggplot(SummaryPMtypes_SED2, aes("",  y=freq, fill=factor(Polymer_NERC))) + geom_bar(width=1, stat="identity") +
  coord_polar(theta = "y", start=0) + theme_bw()+
  theme(axis.line = element_blank(), axis.text = element_blank(), axis.title = element_blank(),legend.title=element_blank(), legend.text = element_text(size=9))+
  facet_grid(`Sample.Type`~`Sampling.Area.general`) 



#####################################
# Size distribution
#####################################
#use also the data_full_SED_adjusted dataset (removed NA's)

#With sizeClasses
##size classes: 
##SC1=0.1mm-0.5mm
##SC2=0.5-1mm
##SC3= 1-5mm
##SC4=5-10mm
##SC5=10-50mm
##SC6=50-100mm
##SC7= >100mm 

summary(data_full_SED_adjusted$`Lenght..mm.`)
##sizes go between 0.1 to 118 mm with very high concentration of microplastics and lower concentration of macroplastics

#summary table per Sample.Type
length_summary_SED<-data_full_SED_adjusted%>%
  group_by(`Sample.Type`)%>%
  summarise(mean=mean(`Lenght..mm.`, na.rm=T), sd=sd(`Lenght..mm.`, na.rm=T), 
            min=min(`Lenght..mm.`, na.rm=T), max=max(`Lenght..mm.`, na.rm=T))
#summary tabel per Sample.Type and sample area
length_summary_SED_area<-data_full_SED_adjusted%>%
  group_by(`Sample.Type`, `Sampling.Area.general`)%>%
  summarise(mean=mean(`Lenght..mm.`, na.rm=T), sd=sd(`Lenght..mm.`, na.rm=T), 
            min=min(`Lenght..mm.`, na.rm=T), max=max(`Lenght..mm.`, na.rm=T))


#General
ggplot(data_full_SED_adjusted, aes(SizeClass)) + geom_bar(aes(y=..prop..,group=1)) + 
  xlab("Size class")+ylab("Relative frequence ") + scale_y_continuous(labels = scales::percent)+gglayer_theme

##General with continuous variable 
ggplot(data_full_SED_adjusted, aes(`Lenght..mm.`)) + geom_histogram(aes(y = after_stat(count / sum(count)))) + 
  xlab("Length (mm)")+ylab("Relative frequence ") + scale_y_continuous(labels = scales::percent)+gglayer_theme
#clear over-representation of microplastics (0.1-5 mm) compared to macroplastics (>5 mm)



### General with categorical variable (true microplastics (0.1-5 mm, according to definition of MP)) to allow for better visualization as macroplastics are underrepresented
SED_size<-data.table(data_full_SED_adjusted)
SED_size[, sizeGroup:=cut(`Lenght..mm.`, breaks=seq(0.1,5, 0.1), dig.lab=10)] 
SED_size<-as.data.frame(SED_size[, table(sizeGroup)])
SED_size<-SED_size%>%
  mutate(total=sum(Freq), relFreq=(Freq/total)*100)
#bar plot
ggplot(SED_size, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity")+gglayer_theme + 
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("Size class (mm)")+ylab("Relative frequence (%)")


###Size distribution per Sampling.Location

functie<-function(number){
  data_full_size<-data.table(data_full_SED_adjusted)
  data_full_size[, sizeGroup:=cut(`Lenght..mm.`, breaks=seq(0.1,5, 0.1), dig.lab=10)]  
  a<-subset(data_full_size, `Sampling.Location`  %in% number)
  Datasize<-as.data.frame(a[, table(sizeGroup)]) 
  Datasize<-Datasize%>%
    mutate(`Sampling.Location` = as.factor(number))%>%
    mutate(total=sum(Freq), relFreq=(Freq/total)*100)
  
}

y<-NULL
for(i in unique(data_full_SED_adjusted$`Sampling.Location`)){
  tmp<- functie(i)
  y<-rbind(y, tmp)
  y<-y%>%
    filter(!is.na(relFreq))
}
unique(y$`Sampling.Location`)

#One bar graph with different colors for locations but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity", position="dodge", aes(fill=`Sampling.Location`))+
  theme(axis.text.x=element_text(angle=45, hjust=1))+ ylab("Relative frequency (%)")+xlab("Size classes (µm)")+gglayer_theme
#Facet grid per locations, but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity",  aes(fill=`Sampling.Location`))+facet_grid(~`Sampling.Location`)+
  theme(axis.text.x=element_text(angle=45, hjust=1)) + ylab("Relative frequency (%)")+xlab("Size classes (µm)")+gglayer_theme


###Size distribution per sampling Area

functie<-function(number){
  data_full_size<-data.table(data_full_SED_adjusted)
  data_full_size[, sizeGroup:=cut(`Lenght..mm.`, breaks=seq(0.1,3, 0.1), dig.lab=10)]  
  a<-subset(data_full_size, `Sampling.Area.general`  %in% number)
  Datasize<-as.data.frame(a[, table(sizeGroup)]) 
  Datasize<-Datasize%>%
    mutate(`Sampling.Area.general` = as.factor(number))%>%
    mutate(total=sum(Freq), relFreq=(Freq/total)*100)
  
}

y<-NULL
for(i in unique(data_full_SED_adjusted$`Sampling.Area.general`)){
  tmp<- functie(i)
  y<-rbind(y, tmp)
  y<-y%>%
    filter(!is.na(relFreq))
}
unique(y$`Sampling.Area.general`)

#One bar graph with different colors for locations but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity", position="dodge", aes(fill=`Sampling.Area.general`))+
  ylab("Relative frequency (%)")+xlab("Size classes (µm)")+gglayer_theme+theme(axis.text.x=element_text(angle=45, hjust=1, size=4)) 
#Facet grid per locations, but to many locations for good overview
ggplot(y, aes(sizeGroup, relFreq)) + geom_bar(stat = "identity",  aes(fill=`Sampling.Area.general`))+facet_grid(~`Sampling.Area.general`)+ theme_bw()+
  ylab("Relative frequency (%)")+xlab("Size classes (µm)") +gglayer_theme +theme(axis.text.x=element_text(angle=45, hjust=1, size=4)) 





#####################################
# Size distribution - alfa values
#####################################
#only using microplastics
data_full_SED_micro_spot_adjusted<-data_full_SED_adjusted%>%
  filter(`Sample.Type`== "Microplastics (> 100 µm)")

set.seed(123)
#Caculate general alfa value
## make into continuous powerlaw object
df.pl <- conpl$new(na.omit(data_full_SED_micro_spot_adjusted$`Lenght..mm.`))   
###########1 check power law distribution
bs.p <- bootstrap_p(df.pl, no_of_sims = 10, threads = 2, xmax = 2E12)
bs.p
#if p > 0.1 ==> power law distiribution is possible
########determine Xmin using KS statistics
estimate_xmin(df.pl)
######Determine alpha (mle method) with estimated xmin
bs.res <- bs.p$bootstraps## get results from the bootstrap
bs.res
mean(bs.res$pars)## gives individual fits (bs.res$pars), take the mean for your overall fit.
sd(bs.res$pars)
mean(bs.res$xmin)## gives individual fits (bs.res$pars), take the mean for your overall fit.
sd(bs.res$xmin)

##Alfa value per Sampling.Location
alfa_SED <- calculate_alfa(data_full_SED_micro_spot_adjusted, "Lenght..mm.")
view(alfa_SED)

##Alfa value per sampling area
alfa_SED_area <- calculate_alfa_area(data_full_SED_micro_spot_adjusted, "Lenght..mm.")
view(alfa_SED_area)




