########################
#libraries
########################

library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(stringr)
library(lubridate)


########################
#datasets
########################


data <- read.csv("Final analysis/Final dataset/data_full_quality.csv")



########################
#dataset with replicate and dry weight for SED_micro
########################

replicate_SED_metadata<-as.data.frame(read_xlsx("Final analysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Replicate data Sediment"))
## Fix some problems in colnames
colnames(replicate_SED_metadata) <- make.names(colnames(replicate_SED_metadata))

replicate_SED_metadata<-replicate_SED_metadata%>%
  select(c("Unique.Sample.Identifier","Replicate..sediment.", "Sample_replicate", "DW..Sediment..kg."))





########################
#additions/changes
########################
#1. selecting spot samples DONE
#2. add own particle count and concentration per sample  DONE
#3. remove unnecessary columns
#4. check class of columns and change where needed
#5. do NA check
#6. check composition of data
#7. Get some typo's out + consistency in colnames
#8. change order of the columns to logic order

#check spelling of column names + change aantal MPs

########################################################################
#1. selecting spot samples 
########################################################################
data<-data%>%
  filter(`Type.of.Campaign` == "Spotsampling")


########################################################################
#2. add own particle count and concentration per sample + remove other columns
########################################################################
#calculate particle count and concentration for water samples
conc_summary_wat<-data%>%
  filter(Matrix == "Water")%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", 0, 1))%>%
  group_by(`Unique.Sample.Identifier`, `Sampling.Area.general`,`Sampling.Location`, `Sample.Type`)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE), Volume..L= mean(`Volume..L.`, na.rm=TRUE))%>%
  mutate(ConcMPL= aantalMPs/Volume..L, ConcMPm3=ConcMPL*1000)


#1. calculate particle count and concentration for  Sediment macroplastic samples
#general plastic concentration calculation
conc_summary_SED_macro<-data%>%
  filter(Matrix == "Sediment" & Sample.Type == "Macroplastics (> 5 mm)")%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", 0, 1))%>%
  group_by(`Unique.Sample.Identifier`, `Sampling.Area.general`,`Sampling.Location`, Sample.Type)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE), DW..Sediment..kg=max(`DW..Sediment..kg.`, na.rm=T))%>%
  mutate(concMPkg= aantalMPs/DW..Sediment..kg)%>%
  mutate(concMPkg = ifelse(is.infinite(DW..Sediment..kg), NA, concMPkg))


#sediment microplastics
conc_summary_SED_micro<-data%>%
  filter(Matrix == "Sediment"& Sample.Type == "Microplastics (> 100 µm)")%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", 0, 1))%>%
  group_by(`Unique.Sample.Identifier`,`Replicate..sediment.`,`Sampling.Area.general`, `Sampling.Location`, Sample.Type)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE))%>%
  left_join(replicate_SED_metadata, by=c("Unique.Sample.Identifier", "Replicate..sediment."))%>%#adding dry weigth based on both sample identifier and replicate
  mutate(concMPkg= aantalMPs/as.numeric(`DW..Sediment..kg.`))

conc_summary_SED_micro_pooled<-conc_summary_SED_micro%>%
  group_by(`Unique.Sample.Identifier`,`Sampling.Area.general`, `Sampling.Location`, Sample.Type)%>%
  summarise(aantalMPs=sum(aantalMPs, na.rm=T), DW..Sediment..kg = sum(as.numeric(DW..Sediment..kg.),na.rm=T), concMPkg = sum(concMPkg, na.rm=T))
  

#streamlinging columns so that we can create Rbind


colnames(conc_summary_SED_micro_pooled) == colnames(conc_summary_SED_macro)

conc_summary_SED<-rbind(conc_summary_SED_micro_pooled, conc_summary_SED_macro)

conc_summary_SED<-conc_summary_SED%>%
  mutate(Volume..L = NA, ConcMPL = NA, ConcMPm3 = NA)


conc_summary_wat<-conc_summary_wat%>%
  mutate(DW..Sediment..kg =NA, concMPkg = NA)%>%
  select("Unique.Sample.Identifier","Sampling.Area.general","Sampling.Location","Sample.Type",
         "aantalMPs","DW..Sediment..kg","concMPkg","Volume..L","ConcMPL","ConcMPm3")

colnames(conc_summary_wat) ==colnames(conc_summary_SED)

conc_summary<-rbind(conc_summary_wat, conc_summary_SED)
conc_summary<-conc_summary%>%
  ungroup()%>%
  select(-c("Sampling.Area.general","Sampling.Location","Sample.Type"))

#join with full dataset

data<-data%>%
  left_join(conc_summary, by = "Unique.Sample.Identifier")



########################################################################
#3. remove unnecessary columns 
########################################################################
#match is not necessary to report
#maximal dimension is mostly the same as length
#color is replaced by color_NERC
#T2 stop was used for calculation of diff_time and still present in the data
#Diff_time has an error and not usefull
#polymer is replaced by Polymer_NERC
#dw sediment, volume, particle count, mass based concentration: all replaced by new calculation

data<-data%>%
  select(-c("X", "X.1", "Match","Maximal.dimension..mm.", "Color","Polymer", "Diff_time",
            "DW..Sediment..kg.", "Volume..L.", "Particle.count", "Total.plastic.mass..kg.",
            "Concentration..g.plastic.L.water.", "Concentration..g.plastic.kg.dry.weight.sediment.",
            "T2..STOP..1"))
  


########################################################################
#4. check class of columns and change where needed
########################################################################

str(data)


#date to date format
data$Date <- as.Date(data$Date, origin = "1899-12-30")
data$Diff_Time<-as.numeric(data$Diff_Time)





########################################################################
#5. do NA check
########################################################################
generate_structure_checks <- function(data) {
  structure.checks <- data.frame(
    na.counts = sapply(data, function(x) sum(is.na(x))),
    na.percent = round(sapply(data, function(x) sum(is.na(x)) / nrow(data) * 100), digits = 1),
    n.levels = sapply(data, function(x) length(unique(x)))
  )
  return(structure.checks)
  
}

na.structure<-generate_structure_checks(data)
#all seems ok, water samples do not have info on sediment dry weight etc. and visa versa




########################################################################
#6. check composition of data
########################################################################
nrow(data)
#microplastic samples
data %>% filter(Sample.Type == "Microplastics (> 100 µm)") %>% nrow()
#microplastic samples sediment
data %>% filter(Sample.Type == "Microplastics (> 100 µm)" & Matrix == "Sediment") %>% nrow()
#microplastic samples water
data %>% filter(Sample.Type == "Microplastics (> 100 µm)" & Matrix == "Water") %>% nrow()

#macroplastic samples
data %>% filter(Sample.Type == "Macroplastics (> 5 mm)") %>% nrow()
#macroplastic samples sediment
data %>% filter(Sample.Type == "Macroplastics (> 5 mm)" & Matrix == "Sediment") %>% nrow()
#macroplastic samples water
data %>% filter(Sample.Type == "Macroplastics (> 5 mm)" & Matrix == "Water") %>% nrow()



########################################################################
#7. Get some typo's out + consistency in colnames
########################################################################


#typo's
colnames(data)[colnames(data) == "Lenght..mm."] <- "Length..mm."
colnames(data)[colnames(data) == "estimated.Hight..mm."] <- "estimated.Height..mm."
colnames(data)[colnames(data) == "CSF"] <- "Corey..Shape..factor"
colnames(data)[colnames(data) == "aantalMPs"] <- "Particle..Count"
colnames(data)[colnames(data) == "concMPkg"] <- "ConcMPkg"

#consistency, always use ..
colnames(data) <- gsub("_", "..", colnames(data))



########################################################################
#8. change order of the columns to logic order
########################################################################
#############campaign info
#"Code.Campaign" 
#"Date"
#"Type.of.Campaign" #removed as only spotsampling
 
#########Sampling info
#"Unique.Sample.Identifier"

#"Sampling.Area" 
#"Sampling.Area.general"
#"Sampling.Location"
#"Sampling.Location.specific"#maybe not necessary?

#"T1..START."                                       
#"T2..STOP." 
#"Diff..Time"
#"Tide" ==> remove as it only contains NA

#"Sample.Type"
#"Matrix" 
#"Depth.Sample..m."
#"Volume..L" 

#"Replicate..sediment." #not usefull as data is pooled for different replicates so removed
#"DW..Sediment..kg"
#"Sediment.type"

#########particle info
#"Unique.Particle.Identifier"  
#"Polymer_NERC"
#"Length..mm."
#"SizeClass"
#"Width..mm."
#"estimated.Height..mm."
#"Weight..mg."
#"LengthWidthratio"
#"Corey..Shape..factor"
#"Shape..obs"
#"Shape..est" 
#"Shape..NERC"
#"Shape..tot"
#"Transpar..NERC"                                    
#"Color..NERC"
#"Color..NERC..Interreg"                              
#"Particle.count"
#"ConcMPkg"
#"ConcMPL"
#"ConcMPm3"


data<-data%>%
  select("Code.Campaign","Date",
         "Unique.Sample.Identifier","Sampling.Area","Sampling.Area.general","Sampling.Location","Sampling.Location.specific", 
         "T1..START.","T2..STOP.","Diff..Time",
         "Sample.Type","Matrix","Depth.Sample..m.","Volume..L",
         "DW..Sediment..kg","Sediment.type",
         "Unique.Particle.Identifier","Polymer..NERC",
         "Length..mm.","SizeClass","Width..mm.","estimated.Height..mm.","Weight..mg.",
         "LengthWidthratio","Corey..Shape..factor","Shape..obs","Shape..est", "Shape..NERC","Shape..tot",
         "Transpar..NERC","Color..NERC","Color..NERC..Interreg",
         "Particle..Count","ConcMPkg","ConcMPL","ConcMPm3")
view(data)



########################################################################
## Exporting data
########################################################################   

write.csv(data, "Final analysis/Final dataset/data_full_publication.csv")
