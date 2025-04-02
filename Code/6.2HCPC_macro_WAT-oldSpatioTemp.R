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

data_full_WAT_macro_spot<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Final dataset/Datasplits/data_full_WAT_macro_spot.csv"))

#descriptor
Descriptor_2km_WAT<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Final dataset/Descriptor/Descriptor_2km_WAT.csv"))




#########################################
####    CLUSTER - preparation        ####
#########################################

#variables for cluster analysis: 
#SizeClass (recalculated to fractions)
#Polymer_NERC (recalculated to fractions) 
#shape
#transparency and color
#Concentration (continuous)


#using concentrations instead of fractions to avoid dependence of the different groups of variables (eg. different polymer types) 
##using only spot sample data


length(unique(data_full_WAT_macro_spot$ `Unique.Sample.Identifier`))
#72 samples

#check for NA's
sum(is.na(data_full_WAT_macro_spot$`Polymer_NERC`))
sum(is.na(data_full_WAT_macro_spot$`SizeClass`))
sum(is.na(data_full_WAT_macro_spot$`Shape_tot`))
#34 samples with NA (no plastic found)
#remove samples without plastics
data_full_WAT_macro_spot<-data_full_WAT_macro_spot%>%
  filter(SizeClass != "NA")


#reforming the dataset
data_full_WAT_macro_cluster_sample <- data_full_WAT_macro_spot %>%
  group_by(`Unique.Sample.Identifier`) %>%
  mutate(PP_total = sum(Polymer_NERC == "polypropylene"),
         PE_total = sum(Polymer_NERC == "polyethylene"),
         PES_total = sum(Polymer_NERC == "polyester"),
         PS_total = sum(Polymer_NERC == "polystyrene"),
         PVC_total=sum (Polymer_NERC == "polychlorinated polymer"), 
         PAM_total=sum(Polymer_NERC == "polyacrilamide"),
         Others_total = sum(Polymer_NERC %in% c("ethylene-vinyl-acetate", "polyether urethane-polypropylene oxide - methylene", "polybutadiene" , "polyurethane", "acrylonitrile butadiene styrene", "polyvinyl alcohol","polymethylacrylate", "polyamide (nylon)", "cellophane", "sodium sterate")),
         Unknown_pm_total = sum(Polymer_NERC == "undefined plastic"),
         Unknown_shape_total = sum(Shape_tot == "Undefined micro-litter items", na.rm=T),
         fragments_total = sum(Shape_tot == "fragments", na.rm=T),
         filaments_total = sum(Shape_tot == "filaments", na.rm=T),
         pellets_total = sum(Shape_tot == "pellets", na.rm=T),
         films_total = sum(Shape_tot == "films", na.rm=T),
         foams_total = sum(Shape_tot == "foams", na.rm=T),
         granules_total = sum(Shape_tot == "granules", na.rm=T),
         others_shape_total = sum(Shape_tot %in% c("J211","J8", "J67", "J9", "J23", "J24")),
         SC1_total = sum(SizeClass == "SC1"),
         SC2_total = sum(SizeClass == "SC2"),
         SC3_total = sum(SizeClass == "SC3"),
         SC4_total = sum(SizeClass == "SC4"),
         SC5_total = sum(SizeClass == "SC5"),
         SC6_total = sum(SizeClass == "SC6"),
         SC7_total = sum(SizeClass == "SC7"))%>%
  add_count(`Unique.Sample.Identifier`, name = "Loc_total") %>%
  mutate(VolumeTOT= sum(unique(`Volume..L.`)))%>%
  mutate_at(vars(PP_total:Loc_total),
            funs((. / VolumeTOT)*1000)) %>%
  rename_with(~ paste0(., "_conc"), matches("_total"))%>%
  summarise(across(ends_with("_conc"), max, na.rm=T))

length(unique(data_full_WAT_macro_cluster_sample$ `Unique.Sample.Identifier`))
#38 samples, correct

#adding descriptor data
WAT_macro_cluster_sample_2km<-data_full_WAT_macro_cluster_sample%>%
  left_join(Descriptor_2km_WAT, by ="Unique.Sample.Identifier")
colnames(WAT_macro_cluster_sample_2km)

#check that correct descriptor data is used
summary(WAT_macro_cluster_sample_2km$`radius..km.`)
#correct, 2.185 km radius

summary(WAT_macro_cluster_sample_2km)
#Removal of data that isn't necessary
#based on visual observation, a few fractions are not present: 
#PES_total
#PVC_total
#PAM_total
#granules_total
#SC1_total
#SC2_total

#samplig location as we have descriptors for locations
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

#remove ...1 and ...2 columns


WAT_macro_cluster_sample_2km<-WAT_macro_cluster_sample_2km%>%
  select(-c("PES_total_conc", "PAM_total_conc","PVC_total_conc", "granules_total_conc","SC1_total_conc","SC2_total_conc",
            "Sampling.Location","Matrix","Sediment.type",  "Depth.Sample..m.",
            "Outside_insideBend","radius..degree.","radius..km.","area..km..","Sampling.Location.specific",
            "pixels.at.flood.risk", "X..buffer.at.flood.risk", "total.flooding.depth.in.buffer",
            "average.flooding.depth.in.buffer", "pop20", "pop22" , "pop21","Pop_AVG", 
            "tot_precip", "mean_temp", "mean_windspeed", 
            "mean_winddirection" , "mean_pressure", "mean_cloudiness", 
            "mean_cloudiness_mean", "tot_precip_TOT","mean_pressure_mean",
            "...1", "...2"))


colnames(WAT_macro_cluster_sample_2km)
summary(WAT_macro_cluster_sample_2km)



#########################################
####    CLUSTER - correlation        ####
#########################################


##---------------------------------------------test correlation between temporal variables and local continuous variables
temporal<-WAT_macro_cluster_sample_2km%>%
  select(c("tot_precip_mean","mean_temp_mean", "mean_windspeed_mean" ,
           "mean_winddirection_mean", "pop_dens"))
Local<-WAT_macro_cluster_sample_2km%>%
  select(c("Mean.slope","Depth.river","Active.overflow",
           "width",
           "shortest.distance.from.shore","RWZI..nr.","Waste.facilities..nr.",
           "agriculture..km..","industry...km.." ,
           "transport...km..","urban...km..", "water...km.." ,"nature...km..","recreation...km..",
           "waste...km..","human.foot.print","km..at.flood.risk"))


#Checking correlation of the variables correlation matrix
data.cor.temporal<-cor(temporal, method="spearman")
corrplot(data.cor.temporal)
#no correlations

#Checking correlation of the local variables correlation matrix
data.cor.local<-cor(Local, method="spearman")
corrplot(data.cor.local)
#recreation and km² at flood risk are highly negatively correlated
#industry is neg correlated with mean slope

#km² at flood risk removed
#mean slope removed

WAT_macro_cluster_sample_2km<-WAT_macro_cluster_sample_2km%>%
  select(-c("km..at.flood.risk", "Mean.slope"))




#########################################
####    CLUSTER - all parameters     ####
#########################################



#parameters for cluster analysis
#"PP_total_conc"                
#"PE_total_conc"                
##"PS_total_conc"                               
#"Others_total_conc"        
#"Unknown_pm_total_conc"        
#"SC3_total_conc"              
#"SC4_total_conc"               
#"SC5_total_conc"   
#"SC6_total_conc"
#"SC7_total_conc"
#"Loc_total_conc"               
#"Unknown_shape_total_conc"     
#"fragments_total_conc"        
#"filaments_total_conc"         
#"pellets_total_conc"           
#"films_total_conc"             
#"foams_total_conc"            
#"others_shape_total_conc"

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


#change to dataframe if necessary         
class(WAT_macro_cluster_sample_2km)
WAT_macro_cluster_sample_2km<-as.data.frame(WAT_macro_cluster_sample_2km)
#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(WAT_macro_cluster_sample_2km) <- WAT_macro_cluster_sample_2km[,1]
#verwijderen kolom sample names
WAT_macro_cluster_sample_2km<-WAT_macro_cluster_sample_2km[,-1]
view(WAT_macro_cluster_sample_2km)
colnames(WAT_macro_cluster_sample_2km)


#-------------------------------------PCA analyse


quanti.sup<-c( "width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
               "agriculture..km..","industry...km..",
               "recreation...km..","transport...km..","urban...km..","nature...km..",
               "waste...km..","water...km..","human.foot.print",
               "pop_dens",
               "mean_winddirection_mean",
               "tot_precip_mean","mean_temp_mean","mean_windspeed_mean",
               "Active.overflow","Depth.river")

quali.sup<-c("Nearby.vegetation", "NaturalBank", "meandering", "ecotope")



#PCA uitvoeren maar de juiste variablen als supplementary ingeven (alle descriptors)
res<-PCA(WAT_macro_cluster_sample_2km,ncp=Inf,quanti.sup=quanti.sup, quali.sup=quali.sup,graph=FALSE)
res$eig #3 dimensions explain 95.72% of variance 


#nieuwe analyse met bepaald aantal dimensies
res.pca.macro.wat<-PCA(WAT_macro_cluster_sample_2km,ncp=3,quanti.sup=quanti.sup, quali.sup=quali.sup,graph=FALSE)
res.pca.macro.wat$eig

#verkennende plots
plot.PCA(res.pca.macro.wat,choix='var')
plot.PCA(res.pca.macro.wat,invisible=c('ind.sup'),label =c('ind'), loadings=TRUE)

#ophalen gegevens PCA
#res.pca.macro.wat$quali.sup$v.test
#res.pca.macro.wat$quanti.sup
#dimdesc(res.pca.macro.wat, axes=c(1,2)) #doens't work yet


#-------------------------------------hierarchical clustering based on PCA
##cluster analysis 
res.hcpc<-HCPC(res.pca.macro.wat, nb.clust=-1, kk=Inf) #let the code determine own number of 

##############3 clusters
#visualize the clusters (dendrogram)
fviz_dend(res.hcpc, 
          cex = 0.7,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 1.6    # Augment the room for labels
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


result_HCPC_macro_WAT<- res.hcpc$data.clust
#add sample ID as column
result_HCPC_macro_WAT$Unique.Sample.Identifier<-rownames(result_HCPC_macro_WAT)

write.csv(result_HCPC_macro_WAT, "PLUXIN-FinalAnalysis/Results/result_HCPC_macro_WAT.csv")

