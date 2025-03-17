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

data_full_WAT_micro_spot<-as.data.frame(read_csv("Final analysis/Final dataset/Datasplits/data_full_WAT_micro_spot.csv"))

#descriptor
Descriptor_2km_WAT<-as.data.frame(read_csv("Final analysis/Final dataset/Descriptor/Descriptor_2km_WAT.csv"))




#########################################
####    CLUSTER - preparation        ####
#########################################


#-------------------------------------Construction of the dataset 
#variables for cluster analysis: 
#SizeClass (recalculated to fractions)
#Polymer_NERC (recalculated to fractions) 
#Concentration (continuous)
#shape (L/W ratio;continuous)


#using concentrations instead of fractions to avoid dependence of the different groups of variables (eg. different polymer types) 

##using only spot sample data
length(unique(data_full_WAT_micro_spot$ `Unique.Sample.Identifier`))
#64 samples

#check for NA's
sum(is.na(data_full_WAT_micro_spot$`Polymer_NERC`))
sum(is.na(data_full_WAT_micro_spot$`SizeClass`))#1 NA ==> remove
data_full_WAT_micro_spot<-data_full_WAT_micro_spot%>%
  filter(SizeClass != "NA")
#removing samples with no plastics which might cause confounding? but were looking at plastic fingerprint... 

##removing or adding 'na.rm=T' has same effect, but with this more control
data_full_WAT_micro_cluster_sample <- data_full_WAT_micro_spot %>%
  group_by(`Unique.Sample.Identifier`) %>%
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
         SC5_total = sum(SizeClass == "SC5")) %>%
  add_count(`Unique.Sample.Identifier`, name = "Loc_total")%>%
  mutate(VolumeTOT= sum(unique(`Volume..L.`)))%>%
  mutate_at(vars(PP_total:Loc_total),
            funs((. / VolumeTOT)*1000))%>%
  rename_with(~ paste0(., "_conc"), matches("_total"))%>%
  summarise(across(ends_with("_conc"), max, na.rm=T))

#average L/W ratio per location
LWratio_WAT_micro<-data_full_WAT_micro_spot%>%
  filter(`Sample.Type` == "Microplastics (> 100 µm)")%>%
  group_by(`Unique.Sample.Identifier`)%>%
  summarise(avgLWratio=mean(LengthWidthratio, na.rm=T))

#average Microplastic concentration per sample ##not necessary as it is the same as Loc_total_conc
#making full dataset
WAT_micro_cluster_sample<-data_full_WAT_micro_cluster_sample%>%
  left_join(LWratio_WAT_micro, by= "Unique.Sample.Identifier")

length(unique(WAT_micro_cluster_sample$ `Unique.Sample.Identifier`))
#62 samples


summary(WAT_micro_cluster_sample)
#adding descriptor data (using both temporal and spatial)
WAT_micro_cluster_sample_2km<-WAT_micro_cluster_sample%>%
  left_join(Descriptor_2km_WAT, by ="Unique.Sample.Identifier")
#check that correct descriptor data is used
summary(WAT_micro_cluster_sample_2km$`radius..km.`)


#verwijderen data die we niet nodig hebben 
#verwijderen PVC aangezien geen enkele meting
#matrix, not useful 
#depth of sampling, sediment type not useful for water samples
#outside_insideBend, not useful for first exploration
#"radius (degree)","radius (km)","area (km²)" chosen ourselves and similar for all data
#alternative measures for flood risk
#pop per year ==> to general, pop_dens gives specific information per locations

#extra changes: 
#removal of "Sampling Location specific" 
#removal of measured weather parameters (keep the 3dmean): "tot_precip", "mean_temp", "mean_windspeed", 
# "mean_winddirection" , "mean_pressure", "mean_cloudiness"
#removal of cloudiness in general: "mean_cloudiness_mean"
#removal of total precip.: "tot_precip_TOT"
#removal of mean pressure: "mean_pressure_mean"

WAT_micro_cluster_sample_2km<-WAT_micro_cluster_sample_2km%>%
  select(-c( "PVC_total_conc","Matrix", "Sampling.Location", "Sediment.type",  "Depth.Sample..m.",
             "Outside_insideBend","radius..degree.","radius..km.","area..km..",
             "pixels.at.flood.risk", "X..buffer.at.flood.risk", "total.flooding.depth.in.buffer",
             "average.flooding.depth.in.buffer", "pop20", "pop22" , "pop21","Pop_AVG", 
             "tot_precip", "mean_temp", "mean_windspeed", 
             "mean_winddirection" , "mean_pressure", "mean_cloudiness", 
             "mean_cloudiness_mean", "tot_precip_TOT","mean_pressure_mean"))                  

colnames(WAT_micro_cluster_sample_2km)
summary(WAT_micro_cluster_sample_2km)




#########################################
####    CLUSTER - correlation        ####
#########################################


##---------------------------------------------test correlation between temporal variables and local continuous variables
temporal<-WAT_micro_cluster_sample_2km%>%
  select(c("tot_precip_mean","mean_temp_mean", "mean_windspeed_mean" ,"pop_dens", "mean_winddirection_mean"))
Local<-WAT_micro_cluster_sample_2km%>%
  select(c("Mean.slope","Depth.river","Active.overflow",
           "width",
           "shortest.distance.from.shore","RWZI..nr.","Waste.facilities..nr.",
           "agriculture..km..","industry...km.." ,
           "transport...km..","urban...km..", "water...km.." ,"nature...km..","recreation...km..",
           "waste...km..","human.foot.print","km..at.flood.risk"))


#Checking correlation of the variables correlation matrix
data.cor.temporal<-cor(temporal, method="spearman")
corrplot(data.cor.temporal)

#Checking correlation of the local variables correlation matrix
data.cor.local<-cor(Local, method="spearman")
corrplot(data.cor.local)
#recreation and km² at flood risk are highly negatively correlated





#########################################
####    CLUSTER - all parameters     ####
#########################################
#parameters for cluster analysis
#"PP_total_conc"                
#"PE_total_conc"                
#"PES_total_conc"              
##"PS_total_conc"                               
#"PAM_total_conc"               
#"Others_total_conc"        
#"Unknown_pm_total_conc"        
#"SC1_total_conc"               
#"SC2_total_conc"              
#"SC3_total_conc"              
#"SC4_total_conc"               
#"SC5_total_conc"               
#"Loc_total_conc"               
#"avgLWratio" 

#Supplementary variables
####quantitative
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
#"km² at flood risk" ==> Remove
#"pop_dens"                     
#"tot_precip_mean"             
#"mean_temp_mean"               
#"mean_windspeed_mean"          
#"mean_winddirection_mean"      
#"Active overflow" 
#"Mean slope"               
#"Depth river"

####qualitative 
#"Nearby vegetation"            
#"NaturalBank"                 
#"meandering"                                     
#"ecotope"                                 

#samplig location specific  as we have descriptors for locations

WAT_micro_cluster_sample_2km<-WAT_micro_cluster_sample_2km%>%
  select(-c( "Sampling.Location.specific","km..at.flood.risk"))            

#change to dataframe if necessary         
class(WAT_micro_cluster_sample_2km)
WAT_micro_cluster_sample_2km<-as.data.frame(WAT_micro_cluster_sample_2km)

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(WAT_micro_cluster_sample_2km) <- WAT_micro_cluster_sample_2km[,1]
#verwijderen kolom sample names
WAT_micro_cluster_sample_2km<-WAT_micro_cluster_sample_2km[,-1]

colnames(WAT_micro_cluster_sample_2km)
summary(WAT_micro_cluster_sample_2km)



#-------------------------------------PCA analyse
set.seed(123) #for reproducibility
quanti.sup<-c( "width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
               "agriculture...km..","industry...km..",
               "recreation...km..","transport...km..","urban...km..","nature...km..",
               "waste...km..","water...km..","human.foot.print",
               "pop_dens",
               "mean_winddirection_mean",
               "tot_precip_mean","mean_temp_mean","mean_windspeed_mean",
               "Active.overflow", "Mean.slope","Depth.river")

quali.sup<-c("Nearby.vegetation", "NaturalBank", "meandering", "ecotope")


#PCA uitvoeren maar de juiste variablen als supplementary ingeven (alle descriptors)
res<-PCA(WAT_micro_cluster_sample_2km,ncp=Inf,quanti.sup=quanti.sup, quali.sup=quali.sup,graph=FALSE)
res$eig #9 dimensions explain 95% of variance 

##only 3 samples that do not have a temp and precipitation data ==> should be fine

#nieuwe analyse met bepaald aantal dimensies
res.pca.micro.wat<-PCA(WAT_micro_cluster_sample_2km,ncp=9,quanti.sup=quanti.sup, quali.sup=quali.sup,graph=FALSE)
res.pca.micro.wat$eig

#verkennende plots
plot.PCA(res.pca.micro.wat,choix='var')
plot.PCA(res.pca.micro.wat,invisible=c('ind.sup'),label =c('ind'), loadings=TRUE)



#ophalen gegevens PCA
#res.pca.micro.wat$quali.sup$v.test
#res.pca.micro.wat$quanti.sup
#dimdesc(res.pca.micro.wat, axes=c(1,2)) #doens't work yet


#-------------------------------------hierarchical clustering based on PCA
set.seed(123) #for reproducibility
##cluster analysis 
res.hcpc<-HCPC(res.pca.micro.wat, nb.clust=3, kk=Inf) #fixed at 3 clusters

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

str(res.hcpc$desc.var)


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


result_HCPC_micro_WAT<- res.hcpc$data.clust
#add sample ID as column
result_HCPC_micro_WAT$Unique.Sample.Identifier<-rownames(result_HCPC_micro_WAT)

write.csv(result_HCPC_micro_WAT, "Final analysis/Results/result_HCPC_micro_WAT.csv")
