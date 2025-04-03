########################
#libraries
########################


library(dplyr)
library(ggplot2)
library(tidyverse)
library(data.table)
library(stringr)
library(ggsci) #voor kleurpalet grafieken
library(patchwork)
library(lubridate)
library(factoextra)#for cluster analysis
library(corrplot)
library("ggpubr")
library(FactoMineR)
library(Factoshiny)
library(reshape2)
library(circlize)
library(ggthemes)#themes for maps


#setting encoding to UTF-8 to avoid some problems with 'µ'
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")

##lay-out grafieken
gglayer_theme<-list(
  theme_classic(), 
  scale_fill_npg(),
  scale_color_npg(),
  theme(legend.position = "top", text= element_text(size=15)))

set.seed(123) #for reproducibility

########################
#Datasets
########################

data_full_SED_micro<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Final dataset/Datasplits/data_full_SED_micro.csv"))

#descriptor
Descriptor_2km_SED<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Final dataset/Descriptor/Descriptor_SED_withadditions.csv"))


#dataset with replicate and dry weight for SED_micro
replicate_SED_metadata<-as.data.frame(read_xlsx("PLUXIN-FinalAnalysis/Raw data/OVAM PLUXIN merged dataset revAC_2 2024 manuscript_V3.xlsx", sheet = "Replicate data Sediment"))
## Fix some problems in colnames
colnames(replicate_SED_metadata) <- make.names(colnames(replicate_SED_metadata))

replicate_SED_metadata<-replicate_SED_metadata%>%
  select(c("Unique.Sample.Identifier","Replicate..sediment.", "Sample_replicate", "DW..Sediment..kg."))



#########################################
####    CLUSTER - preparation        ####
#########################################


#-------------------------------------Construction of the dataset 
#variables for cluster analysis: 
#SizeClass ==> concentration
#Polymer_NERC  ==> concentration
#Concentration==> #conc_summary_SED_micro2 ==> mean of concentrations of 3 replicates

#Not useful for micropl
#shape
#color_NERC(recalculated to fractions)
#Transpar_NERC(recalculated to fractions)


##cluster analysis doesn't work well with NA's 
sum(is.na(data_full_SED_micro$SizeClass)) #0 NA's
sum(is.na(data_full_SED_micro$Polymer_NERC))#no NA's

#average L/W ratio per location
LWratio_SED_micro<-data_full_SED_micro%>%
  filter(`Sample.Type` == "Microplastics (> 100 µm)")%>%
  group_by(`Unique.Sample.Identifier`)%>%
  summarise(avgLWratio=mean(LengthWidthratio, na.rm=T))


#concentration of microplastics in samples
conc_summary_SED_micro<-data_full_SED_micro%>%
  mutate(Aantal=ifelse(`Polymer` == "NA", NA, 1))%>%
  group_by(`Unique.Sample.Identifier`,`Replicate..sediment.`)%>%
  summarise(aantalMPs= sum(Aantal, na.rm=TRUE))%>%
  left_join(replicate_SED_metadata, by=c("Unique.Sample.Identifier", "Replicate..sediment."))%>%
  mutate(concMPkg= aantalMPs/as.numeric(`DW..Sediment..kg.`))
#calcuation of total dry weight of the 3 replicates per sample
conc_summary_SED_micro$`DW..Sediment..kg.`<-as.numeric(conc_summary_SED_micro$`DW..Sediment..kg.`)

Weight<-conc_summary_SED_micro%>%
  group_by(`Unique.Sample.Identifier`)%>%
  summarise(WeightTOT=sum(`DW..Sediment..kg.`, na.rm=T))


#preparing the dataset
data_full_SED_micro_cluster_sample<-data_full_SED_micro%>%
  group_by(`Unique.Sample.Identifier`)%>%
  mutate(PP_total = sum(Polymer_NERC == "polypropylene"),
         PE_total = sum(Polymer_NERC == "polyethylene"),
         PES_total = sum(Polymer_NERC == "polyester"),
         PS_total = sum(Polymer_NERC == "polystyrene"),
         PVC_total=sum (Polymer_NERC == "polychlorinated polymer"), 
         PAM_total=sum(Polymer_NERC == "polyacrilamide"),
         Others_total = sum(Polymer_NERC %in% c("ethylene-vinyl-acetate", "polyether urethane-polypropylene oxide - methylene", "polybutadiene" , "polyurethane", "acrylonitrile butadiene styrene", "polyvinyl alcohol","polymethylacrylate", "polyamide (nylon)", "cellophane", "sodium sterate")),
         Unknown_pm_total = sum(Polymer_NERC == "undefined plastic"),
         SC1_total = sum(SizeClass == "SC1"),
         SC2_total = sum(SizeClass == "SC2"),
         SC3_total = sum(SizeClass == "SC3"),
         SC4_total = sum(SizeClass == "SC4"),
         SC5_total = sum(SizeClass == "SC5"),
         SC6_total = sum(SizeClass == "SC6"),
         SC7_total = sum(SizeClass == "SC7"))%>%
  add_count(`Unique.Sample.Identifier`, name = "Loc_total")%>%
  left_join(Weight, by = "Unique.Sample.Identifier")%>%
  mutate_at(vars(PP_total:Loc_total),
            funs(./ WeightTOT)) %>%
  rename_with(~ paste0(., "_conc"), matches("_total"))%>%
  summarise(across(ends_with("_conc"), max, na.rm=T))%>%
  left_join(LWratio_SED_micro, by = "Unique.Sample.Identifier")


length(unique(data_full_SED_micro_cluster_sample$ `Unique.Sample.Identifier`))
#27 samples



#adding descriptor data (using both temporal and spatial)
SED_micro_cluster_sample<-data_full_SED_micro_cluster_sample%>%
  left_join(Descriptor_2km_SED, by ="Unique.Sample.Identifier")

#check that correct descriptor data is used
summary(SED_micro_cluster_sample$`radius..km.`)

summary(SED_micro_cluster_sample)
#based on first observation, a few fractions are not present: 
#PVC, Ohter, unknown
#SC5, 6 en 7

#remove variables that are not useful
#matrix, not useful
#outside_insideBend, not useful for first exploration
#"radius (degree)","radius (km)","area (km²)" chosen ourselves and similar for all data
#alternative measures for flood risk
#pop per year ==> to general, pop_dens gives specific information per locations
#...1

#extra changes: 
#removal of "Sampling Location specific" 
#removal of measured weather parameters (keep the 3dmean): "tot_precip", "mean_temp", "mean_windspeed", 
# "mean_winddirection" , "mean_pressure", "mean_cloudiness"
#removal of cloudiness in general: "mean_cloudiness_mean"
#removal of total precip.: "tot_precip_TOT"
#removal of mean pressure: "mean_pressure_mean"


SED_micro_cluster_sample<-SED_micro_cluster_sample%>%
  select(c("Unique.Sample.Identifier","PP_total_conc","PE_total_conc","PES_total_conc",
           "PS_total_conc","PAM_total_conc","SC1_total_conc","SC2_total_conc","SC3_total_conc",
           "SC4_total_conc","Loc_total_conc","avgLWratio","width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
           "agriculture..km..","industry...km..",
           "recreation...km..","transport...km..","urban...km..","nature...km..",
           "waste...km..","water...km..","human.foot.print",
           "pop_dens",
           "mean_winddirection_mean",
           "tot_precip_mean","mean_temp_mean","mean_windspeed_mean",
           "Active.overflow", "Mean.slope","Nearby.vegetation", "NaturalBank", "meandering", "ecotope", "Season", 
           "Sediment.type", "Depth.Sample..m."))


colnames(SED_micro_cluster_sample)

#######################
#Check collinearity
#######################
#change categorical in numeric
temporal<-SED_micro_cluster_sample%>%
  select(c("tot_precip_mean","mean_temp_mean", "mean_windspeed_mean" ,
           "mean_winddirection_mean", "pop_dens", "Season"))
Local<-SED_micro_cluster_sample%>%
  select(c("width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
           "agriculture..km..","industry...km..",
           "recreation...km..","transport...km..","urban...km..","nature...km..",
           "waste...km..","water...km..","human.foot.print",
           "pop_dens",
           "Active.overflow", "Mean.slope","Nearby.vegetation", "NaturalBank", "meandering", "ecotope"))



##change categorical variables to numeric
temporal_numeric <- temporal %>%
  mutate(across(where(is.character) , as.factor)) %>%  # Convert characters (except sample identifier) to factors
  mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))
Local_numeric <- Local %>%
  mutate(across(where(is.character) , as.factor)) %>%  # Convert characters (except sample identifier) to factors
  mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))


#We can visually look for correlations between variables:
heatmap(abs(cor(temporal_numeric)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 1, bty = "n",
       fill = rev(heat.colors(6)))

heatmap(abs(cor(Local_numeric)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 1, bty = "n",
       fill = rev(heat.colors(6)))


#Correlations between
# Depth and width ==> remove width
# width and distance from shore==> remove width and distance
# active overflow and waste facilities ==> remove active overflow
# active overlow and nature==> remove active overflow
# nearby vegetation and agriculutre
# depth and water
# water and shortest distance and width ==> remove width and distance
# ecotope and natural bank ==> remove natural bank
# transport and waste facilities


SED_micro_cluster_sample<-SED_micro_cluster_sample%>%
  select(-c(Active.overflow, shortest.distance.from.shore, width, NaturalBank, water...km..))

#########################################
####    CLUSTER - all parameters     ####
#########################################
#parameters for cluster analysis
#"PP_total_conc"                
#"PE_total_conc"                
#"PES_total_conc"              
##"PS_total_conc"                               
#"PAM_total_conc"               
#"SC1_total_conc"               
#"SC2_total_conc"              
#"SC3_total_conc"              
#"SC4_total_conc"               
#"Loc_total_conc"               
#"avgLWratio" 

#Supplementary variables
####quantitative (24)
#"width" 
#"shortest distance from shore" 
#"RWZI [nr]" 
#"Waste facilities [nr]"        
#"agriculture [km²]"            
#"industry  [km²]"              
#"transport  [km²]"            
#"urban  [km²]"                 
#"water  [km²]"                 
#"nature  [km²]"                
#"recreation  [km²]"           
#"waste  [km²]"                
#"human foot print"
#"km² at flood risk" ==> remove   
#"pop_dens"                     
#"tot_precip_mean"             
#"mean_temp_mean"               
#"mean_windspeed_mean"          
#"mean_winddirection_mean"      
#"Active overflow" 
#"Mean slope"  ==> remove              
#"Depth river"
#"Depth Sample (m)"

####qualitative  (5)
#"Nearby vegetation"            
#"NaturalBank"                 
#"meandering"                                     
#"ecotope" 
#Sediment type"


#change to dataframe if necessary         
class(SED_micro_cluster_sample)
SED_micro_cluster_sample<-as.data.frame(SED_micro_cluster_sample)

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(SED_micro_cluster_sample) <- SED_micro_cluster_sample[,"Unique.Sample.Identifier"]
#verwijderen kolom sample names
SED_micro_cluster_sample<-SED_micro_cluster_sample%>%
  select(-"Unique.Sample.Identifier")
colnames(SED_micro_cluster_sample)

#-------------------------------------PCA analyse
#test voor nodig aantal dimensies om mee te nemen naar clustering, default = 5
#voor consistentie kiezen we altijd min. 95% van de variantie beschreven.
quanti.sup<-c( "width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
               "agriculture..km..","industry...km..",
               "recreation...km..","transport...km..","urban...km..","nature...km..",
               "waste...km..","water...km..","human.foot.print",
               "pop_dens",
               "mean_winddirection_mean",
               "tot_precip_mean","mean_temp_mean","mean_windspeed_mean",
               "Active.overflow", "Depth.river","Depth.Sample..m.", "Mean.slope")


quali.sup<-c("Nearby.vegetation", "NaturalBank", "meandering", "ecotope","Sediment.type", "Season")

set.seed(123)
#PCA uitvoeren maar de juiste variablen als supplementary ingeven (alle descriptors)
res<-PCA(SED_micro_cluster_sample,ncp=Inf,quanti.sup=quanti.sup, quali.sup=quali.sup,graph=FALSE)
res$eig #6 dimensions explain 95.76% of variance 


#nieuwe analyse met bepaald aantal dimensies
res.pca.micro.sed<-PCA(SED_micro_cluster_sample,ncp=6,quanti.sup=quanti.sup, quali.sup=quali.sup,graph=FALSE)
res.pca.micro.sed$eig

#verkennende plots
plot.PCA(res.pca.micro.sed,choix='var')
plot.PCA(res.pca.micro.sed,invisible=c('ind.sup'),label =c('ind'), loadings=TRUE)

#ophalen gegevens PCA
#res.pca.micro.sed$quali.sup$v.test
#res.pca.micro.sed$quanti.sup
#dimdesc(res.pca.micro.wat, axes=c(1,2)) #doesn't work yet





#-------------------------------------hierarchical clustering based on PCA
set.seed(123) #for reproducibility
##cluster analysis 
res.hcpc<-HCPC(res.pca.micro.sed, nb.clust=-1, kk=Inf) #let the code determine own number of clusters

#visualize the clusters (dendrogram)
#fviz_dend(res.hcpc,cex = 0.7,                     # Label size
#         palette = "jco",               # Color palette see ?ggpubr::ggpar
#         rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
#         rect_border = "jco",           # Rectangle color
#          labels_track_height = 1.6      # Augment the room for labels
#)
#visualize the clusters (factor map)
#fviz_cluster(res.hcpc,
#             geom="point",
#             repel = TRUE,            # Avoid label overlapping
#             show.clust.cent = TRUE, # Show cluster centers
#             palette = "jco",         # Color palette see ?ggpubr::ggpar
#             ggtheme = theme_minimal(),
#             main = "Factor map"
#)


## 6 main clusters would end up in a more relevant analysis? 
res.hcpc<-HCPC(res.pca.micro.sed, nb.clust=6, kk=Inf) 


#visualize the clusters (dendrogram)
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 1.6      # Augment the room for labels
)
#visualize the clusters (factor map)
fviz_cluster(res.hcpc,
             geom="point",
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco",         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)


#details of the clusters
res.hcpc$desc.var
res.hcpc$desc.axes
#res.hcpc$desc.ind
res.hcpc$call$t #details of clustering
res.hcpc$data.clust

res.hcpc$desc.var$test.chi2
res.hcpc$desc.var$quali
res.hcpc$desc.var$quanti
res.hcpc$call$var
table(res.hcpc$data.clust$clust)
#Investigate(res.pca.micro.wat) ##small report on the analysis


#visualize the clusters (black and white dendrogram)
fviz_dend(res.hcpc, 
          cex = 0.9,                     # Label size
          palette = "black",               # Color palette see ?ggpubr::ggpar
          rect = TRUE,                   # Add rectangle around groups
          labels_track_height = 1.6,      # Augment the room for labels
          main=""
)

##################################
##exporting results
#################################


result_HCPC_micro_SED<- res.hcpc$data.clust
#add sample ID as column
result_HCPC_micro_SED$Unique.Sample.Identifier<-rownames(result_HCPC_micro_SED)

write.csv(result_HCPC_micro_SED, "PLUXIN-FinalAnalysis/Results/result_HCPC_micro_SED.csv")
