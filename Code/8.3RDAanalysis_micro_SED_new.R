########################
#libraries
########################

library(dplyr)
library(ggplot2)
library(tidyverse)
library(vegan)


set.seed(123) #for reproducibility

########################
#Datasets
########################


#microplastic data
data_full_SED_micro_spot<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Final dataset/Datasplits/data_full_SED_micro.csv"))


#cluster
cluster<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Final dataset/Cluster_mic_sed.csv"))
cluster<-cluster%>%
  select(-`...1`)

###Descriptor data
#based on full Descriptor dataset with additions
Descriptor_SED_additions<- read.csv("PLUXIN-FinalAnalysis/Final dataset/Descriptor/Descriptor_SED_withadditions.csv")




#############################################
#Necessary changes in plastic dataset
#############################################
#concentration of microplastics in samples
conc_summary_SED_micro<-data_full_SED_micro_spot%>%
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


data_test<- data_full_SED_micro_spot%>%
  mutate(Particle_identity=str_c(Polymer_NERC, SizeClass, sep="-"))%>%
  group_by(Unique.Sample.Identifier,Particle_identity)%>%
  summarise(aantalMPs= n())%>%
  left_join(Weight, by="Unique.Sample.Identifier")%>%
  mutate(ConcMP= aantalMPs/WeightTOT)

data_conc_particle_id_pivot<-data_test%>%
  select(Unique.Sample.Identifier, Particle_identity, ConcMP)%>%
  pivot_wider(names_from = Particle_identity, values_from = ConcMP)%>%
  mutate_all(~replace(., is.na(.), 0))%>%
  ungroup()
data_conc_particle_id_pivot<-as.data.frame(data_conc_particle_id_pivot)




# SED_micro_selectionSamples<-SED_micro_RDA%>%
#   select(`Unique.Sample.Identifier`)

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(data_conc_particle_id_pivot) <- data_conc_particle_id_pivot[,"Unique.Sample.Identifier"]


data_conc_particle_id_pivot<-data_conc_particle_id_pivot%>%
  select(-`Unique.Sample.Identifier`)
# Hellinger transform the community data
SED_micro_RDA_hel <- decostand(data_conc_particle_id_pivot, method = "hellinger")


#check colinearity between plastics
heatmap(abs(cor(SED_micro_RDA_hel)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
#no high collinearity 

#############################################
#Necessary changes in descriptor dataset
#############################################

#Environmental variables
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
#"Depth Sample (m)"

####qualitative  (5)
#"Nearby vegetation"            
#"NaturalBank"                 
#"meandering"                                     
#"ecotope" 
#Sediment type"


#NAs in diepte (redelijk veel) en regenval en temperatuur
#Diepte weglaten 
#extra kolommen weglaten 

Descriptor_SED_additions<-Descriptor_SED_additions%>%
  select(c("Unique.Sample.Identifier","width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
           "agriculture..km..","industry...km..",
           "recreation...km..","transport...km..","urban...km..","nature...km..",
           "waste...km..","water...km..","human.foot.print",
           "pop_dens",
           "mean_winddirection_mean",
           "tot_precip_mean","mean_temp_mean","mean_windspeed_mean",
           "Active.overflow", "Mean.slope","Nearby.vegetation", "NaturalBank", "meandering", "ecotope", "Season", 
           "Sediment.type", "Depth.Sample..m."))

#temp en neerslag: aanvullen met gemiddelde in het seizoen
mean_autumn_temp <- Descriptor_SED_additions %>%
  filter(Season == "Autumn" & !is.na(mean_temp_mean)) %>%
  summarize(mean_temp = mean(mean_temp_mean)) %>%
  pull(mean_temp)

Descriptor_SED_additions<- Descriptor_SED_additions %>%
  group_by(Unique.Sample.Identifier)%>%
  mutate(mean_temp_mean = ifelse(is.na(mean_temp_mean) & Season == "Autumn", 
                                 mean_autumn_temp, 
                                 mean_temp_mean))

mean_autumn_precip <- Descriptor_SED_additions %>%
  filter(Season == "Autumn" & !is.na(tot_precip_mean)) %>%
  summarize(mean_precip = mean(tot_precip_mean)) %>%
  pull(mean_precip)

Descriptor_SED_additions <- Descriptor_SED_additions %>%
  group_by(Unique.Sample.Identifier)%>%
  mutate(tot_precip_mean = ifelse(is.na(tot_precip_mean) & Season == "Autumn", 
                                  mean_autumn_precip, 
                                  tot_precip_mean))

Descriptor_SED_additions<-SED_micro_selectionSamples%>%
  left_join(Descriptor_SED_additions, by = "Unique.Sample.Identifier")

#NA's in sediment type replacen
Descriptor_SED_additions$`Sediment.type`[is.na(Descriptor_SED_additions$`Sediment.type`)] <- "Unknown"



# Standardize quantitative environmental data
Descriptor_SED_additions$width <- decostand(Descriptor_SED_additions$width, method = "standardize")
Descriptor_SED_additions$width <- round(Descriptor_SED_additions$width, 6)

Descriptor_SED_additions$`shortest.distance.from.shore` <- decostand(Descriptor_SED_additions$`shortest.distance.from.shore`, method = "standardize")
Descriptor_SED_additions$shortest.distance.from.shore <- round(Descriptor_SED_additions$shortest.distance.from.shore, 6)

Descriptor_SED_additions$`RWZI..nr.` <- decostand(Descriptor_SED_additions$`RWZI..nr.`, method = "standardize")
Descriptor_SED_additions$RWZI..nr. <- round(Descriptor_SED_additions$RWZI..nr., 6)

Descriptor_SED_additions$`Waste.facilities..nr.` <- decostand(Descriptor_SED_additions$`Waste.facilities..nr.`, method = "standardize")
Descriptor_SED_additions$Waste.facilities..nr. <- round(Descriptor_SED_additions$Waste.facilities..nr., 6)

Descriptor_SED_additions$`agriculture..km..` <- decostand(Descriptor_SED_additions$`agriculture..km..`, method = "standardize")
Descriptor_SED_additions$agriculture..km.. <- round(Descriptor_SED_additions$agriculture..km.., 6)

Descriptor_SED_additions$`industry...km..` <- decostand(Descriptor_SED_additions$`industry...km..`, method = "standardize")
Descriptor_SED_additions$industry...km.. <- round(Descriptor_SED_additions$industry...km.., 6)

Descriptor_SED_additions$`transport...km..` <- decostand(Descriptor_SED_additions$`transport...km..`, method = "standardize")
Descriptor_SED_additions$transport...km.. <- round(Descriptor_SED_additions$transport...km.., 6)

Descriptor_SED_additions$`urban...km..` <- decostand(Descriptor_SED_additions$`urban...km..`, method = "standardize")
Descriptor_SED_additions$urban...km.. <- round(Descriptor_SED_additions$urban...km.., 6)

Descriptor_SED_additions$`nature...km..` <- decostand(Descriptor_SED_additions$`nature...km..`, method = "standardize")
Descriptor_SED_additions$nature...km.. <- round(Descriptor_SED_additions$nature...km.., 6)

Descriptor_SED_additions$`recreation...km..` <- decostand(Descriptor_SED_additions$`recreation...km..`, method = "standardize")
Descriptor_SED_additions$recreation...km.. <- round(Descriptor_SED_additions$recreation...km.., 6)

Descriptor_SED_additions$`waste...km..` <- decostand(Descriptor_SED_additions$`waste...km..`, method = "standardize")
Descriptor_SED_additions$waste...km.. <- round(Descriptor_SED_additions$waste...km.., 6)

Descriptor_SED_additions$`human.foot.print` <- decostand(Descriptor_SED_additions$`human.foot.print`, method = "standardize")
Descriptor_SED_additions$human.foot.print <- round(Descriptor_SED_additions$human.foot.print, 6)

Descriptor_SED_additions$`pop_dens` <- decostand(Descriptor_SED_additions$`pop_dens`, method = "standardize")
Descriptor_SED_additions$pop_dens <- round(Descriptor_SED_additions$pop_dens, 6)

Descriptor_SED_additions$`tot_precip_mean` <- decostand(Descriptor_SED_additions$`tot_precip_mean`, method = "standardize")
Descriptor_SED_additions$tot_precip_mean <- round(Descriptor_SED_additions$tot_precip_mean, 6)

Descriptor_SED_additions$`mean_temp_mean` <- decostand(Descriptor_SED_additions$`mean_temp_mean`, method = "standardize")
Descriptor_SED_additions$mean_temp_mean <- round(Descriptor_SED_additions$mean_temp_mean, 6)

Descriptor_SED_additions$`mean_windspeed_mean` <- decostand(Descriptor_SED_additions$`mean_windspeed_mean`, method = "standardize")
Descriptor_SED_additions$mean_windspeed_mean <- round(Descriptor_SED_additions$mean_windspeed_mean, 6)

Descriptor_SED_additions$`mean_winddirection_mean` <- decostand(Descriptor_SED_additions$`mean_winddirection_mean` , method = "standardize")
Descriptor_SED_additions$mean_winddirection_mean <- round(Descriptor_SED_additions$mean_winddirection_mean, 6)

Descriptor_SED_additions$`Active.overflow` <- decostand(Descriptor_SED_additions$`Active.overflow`, method = "standardize")
Descriptor_SED_additions$Active.overflow <- round(Descriptor_SED_additions$Active.overflow, 6)

Descriptor_SED_additions$`Mean.slope` <- decostand(Descriptor_SED_additions$`Mean.slope`, method = "standardize")
Descriptor_SED_additions$Mean.slope <- round(Descriptor_SED_additions$Mean.slope, 6)

Descriptor_SED_additions$`Depth.Sample..m.` <- decostand(Descriptor_SED_additions$`Depth.Sample..m.`, method = "standardize")
Descriptor_SED_additions$Depth.Sample..m. <- round(Descriptor_SED_additions$Depth.Sample..m., 6)



class(Descriptor_SED_additions)

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(Descriptor_SED_additions)<- Descriptor_SED_additions[,"Unique.Sample.Identifier"]
#verwijderen kolom sample names
Descriptor_SED_additions<-Descriptor_SED_additions%>%
  select(-Unique.Sample.Identifier)



################################################
##Redundancy analysis
################################################
#Redundancy Analysis (RDA) is a direct extension of multiple regression, as it models the effect 
#of an explanatory matrix X on a response matrix  Y
# https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html




########################
##Check collinearity
########################
#change categorical in numeric
Descriptor_numeric <- Descriptor_SED_additions %>%
  mutate(across(where(is.character), as.factor)) %>%  # Convert characters to factors
  mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))

# We can visually look for correlations between variables:
heatmap(abs(cor(Descriptor_numeric)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
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


Descriptor_SED_additions<-Descriptor_SED_additions%>%
  select(-c(Active.overflow, shortest.distance.from.shore, width, NaturalBank, water...km..))


####################################
##initial RDA: full model
####################################
set.seed(123) #for reproducibility
# Initial RDA with ALL of the environmental data
micro.SED.rda <- rda(SED_micro_RDA_hel ~ ., data = Descriptor_SED_additions)
summary(micro.SED.rda)

# Find the adjusted R2 of the model with the retained env
# variables
RsquareAdj(micro.SED.rda)$adj.r.squared
#calculating the proportion of the variation of  Y explained by the variables in  X

#test model significance
anova.cca(micro.SED.rda, step = 1000)

#You can also test the significance of each variable
anova.cca(micro.SED.rda, step = 1000, by = "term")

# RDA plot
ordiplot(micro.SED.rda, scaling = 2)




####################################
##model optimization 
####################################

micro.SED.rda.test1 <- rda(SED_micro_RDA_hel ~  RWZI..nr., data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test1)$adj.r.squared


micro.SED.rda.test2 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens, data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test2)$adj.r.squared

micro.SED.rda.test3 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season, data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test3)$adj.r.squared


micro.SED.rda.test4 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season + ecotope , data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test4)$adj.r.squared


micro.SED.rda.test5 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season + ecotope +meandering, data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test5)$adj.r.squared

micro.SED.rda.test6 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season + ecotope +meandering + 
                             Waste.facilities..nr., data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test6)$adj.r.squared

micro.SED.rda.test7 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season + ecotope +meandering + 
                             Waste.facilities..nr. + industry...km.. , data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test7)$adj.r.squared

micro.SED.rda.test8 <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season + ecotope +meandering + 
                             Waste.facilities..nr. + industry...km.. +mean_winddirection_mean, data = Descriptor_SED_additions)
RsquareAdj(micro.SED.rda.test8)$adj.r.squared



#adding transport...km../recreation...km../agriculture..km../urban...km../nature/waste/human.foot.print/tot_precip_mean
#tot_precip_mean/mean_temp_mean/ mean_windspeed_mean/Mean.slope/Nearby.vegetation/Sediment.type/Depth.Sample..m. didn't improve the model




######################
#Final model
########################

set.seed(123) #for reproducibility

# Initial RDA with ALL of the environmental data
micro.SED.rda <- rda(SED_micro_RDA_hel ~  RWZI..nr. + pop_dens + Season + ecotope +meandering + 
                       Waste.facilities..nr. + industry...km.. +mean_winddirection_mean, data = Descriptor_SED_additions)
summary(micro.SED.rda)

# Find the adjusted R2 of the model with the retained env
# variables
RsquareAdj(micro.SED.rda)$adj.r.squared
#calculating the proportion of the variation of  Y explained by the variables in  X

#test model significance
anova.cca(micro.SED.rda, step = 1000, permutation=999)

#test model significance
anova.cca(micro.SED.rda, by= "axis", step = 1000, permutation=999)

#You can also test the significance of each variable
anova.cca(micro.SED.rda, step = 1000, by = "term", permutation=999)

ordiplot(micro.SED.rda, scaling=1)
ordiplot(micro.SED.rda, scaling=2)


#plot Scaling 2
env_scores <- scores(micro.SED.rda, display = "bp", scaling = 2)
env_df <- as.data.frame(env_scores)
env_df$Variable <- rownames(env_df)

#with clusters
ggplot() +
  geom_point(data = as.data.frame(scores(micro.SED.rda, display = "sites", scaling = 2)), 
             aes(x = RDA1, y = RDA2, color = as.factor(cluster$clust)), size=3) +
  geom_segment(data = env_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = env_df, 
            aes(x = RDA1, y = RDA2, label = Variable), 
            vjust = -1, color = "black") +
  theme_bw() +geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "RDA1 (32.32%)", y = "RDA2 (14.82%)")+ labs(color='Clusters') 


#plot Scaling 1
env_scores <- scores(micro.SED.rda, display = "bp", scaling = 1)
env_df <- as.data.frame(env_scores)
env_df$Variable <- rownames(env_df)

#with clusters
ggplot() +
  geom_point(data = as.data.frame(scores(micro.SED.rda, display = "sites", scaling = 1)), 
             aes(x = RDA1, y = RDA2, color = as.factor(cluster$clust)), size=3) +
  geom_segment(data = env_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = env_df, 
            aes(x = RDA1, y = RDA2, label = Variable), 
            vjust = -1, color = "black") +
  theme_bw() +geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "RDA1 (32.32%)", y = "RDA2 (14.82%)")+ labs(color='Clusters') 


















# #problems: 
# #overfitting
# #remove redundant variables
# #imbalance between predictor and response
# #check shared variance via variance partitioning.
# 
# ##########################################################################################
# ##RDA with reduced number of predictor variables to increase fitting of the model.
# ##########################################################################################
# set.seed(123) #for reproducibility
# ###Descriptor data
# #based on full Descriptor dataset with additions
# Descriptor_SED_additions<- read.csv("PLUXIN-FinalAnalysis/Final dataset/Descriptor/Descriptor_SED_withadditions.csv")
# 
# 
# #############################################
# #Necessary changes in descriptor dataset
# #############################################
# 
# #Environmental variables
# ####quantitative
# #"width" ==> expected to be less important for sediment 
# #"shortest distance from shore"==> expected to be less important for sediment 
# #"Nearby vegetation"   ==> expected to be less important for sediment          
# #"NaturalBank" ==> expected to be less important for sediment                
# #"meandering"  ==> keep                                   
# #"ecotope" ==> keep
# 
# 
# #"RWZI [nr]" 
# #"Waste facilities [nr]" 
# #"Active overflow" 
# #replace by "TotalPointDischarge"
# 
# 
# #"agriculture [km²]"            
# #"industry  [km²]"              
# #"transport  [km²]"            
# #"urban  [km²]"                 
# #"water  [km²]"                 
# #"nature  [km²]"                
# #"recreation  [km²]"           
# #"waste  [km²]"  
# #landuse==> can be captured by Natural_non natural
# 
# #"human foot print" ==> repetitive, remove
# #"pop_dens"       ==> keep              
# 
# #"tot_precip_mean"             
# #"mean_temp_mean"               
# #"mean_windspeed_mean"          
# #"mean_winddirection_mean" 
# #replace by season
# 
# 
# #"Depth Sample (m)"
# #Sediment type"
# 
# 
# #NAs in diepte (redelijk veel) en regenval en temperatuur
# #Diepte waterkolom weglaten 
# #extra kolommen weglaten 
# 
# Descriptor_SED_reduced<-Descriptor_SED_additions%>%
#   select(c("Unique.Sample.Identifier", "meandering", "ecotope","TotalPointDischarge", 
#            "Natural_NonNatural","pop_dens", "Season","Depth.Sample..m.", "Sediment.type"))
# 
# #remove macroplastic samples
# Descriptor_SED_reduced <- Descriptor_SED_reduced %>%
#   filter(!grepl("^MAC", Unique.Sample.Identifier))
# 
# 
# #NA's in sediment type replacen
# Descriptor_SED_reduced$`Sediment.type`[is.na(Descriptor_SED_reduced$`Sediment.type`)] <- "Unknown"
# 
# 
# # Standardize quantitative environmental data
# Descriptor_SED_reduced$`pop_dens` <- decostand(Descriptor_SED_reduced$`pop_dens`, method = "standardize")
# Descriptor_SED_reduced$`Depth.Sample..m.` <- decostand(Descriptor_SED_reduced$`Depth.Sample..m.`, method = "standardize")
# Descriptor_SED_reduced$TotalPointDischarge <- decostand(Descriptor_SED_reduced$TotalPointDischarge, method = "standardize")
# 
# 
# 
# class(Descriptor_SED_reduced)
# 
# #veranderen rijnamen (nodig voor PCA/Cluster)
# rownames(Descriptor_SED_reduced)<- Descriptor_SED_reduced[,"Unique.Sample.Identifier"]
# #verwijderen kolom sample names
# Descriptor_SED_reduced<-Descriptor_SED_reduced%>%
#   select(-Unique.Sample.Identifier)
# 
# 
# 
# ################################################
# ##Redundancy analysis
# 
# #Redundancy Analysis (RDA) is a direct extension of multiple regression, as it models the effect 
# #of an explanatory matrix X on a response matrix  Y
# # https://r.qcbs.ca/workshop10/book-en/redundancy-analysis.html
# 
# 
# 
# 
# ########################
# ##Check collinearity
# #change categorical in numeric
# Descriptor_numeric <- Descriptor_SED_reduced %>%
#   mutate(across(where(is.character), as.factor)) %>%  # Convert characters to factors
#   mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))
# 
# # We can visually look for correlations between variables:
# heatmap(abs(cor(Descriptor_numeric)), 
#         # Compute pearson correlation (note they are absolute values)
#         col = rev(heat.colors(6)), 
#         Colv = NA, Rowv = NA)
# legend("topright", 
#        title = "Absolute Pearson R",
#        legend =  round(seq(0,1, length.out = 6),1),
#        y.intersp = 0.7, bty = "n",
#        fill = rev(heat.colors(6)))
# 
# 
# set.seed(123) #for reproducibility
# # Initial RDA with ALL of the environmental data
# micro.SED.rda <- rda(SED_micro_RDA_hel ~  ., data = Descriptor_SED_reduced)
# summary(micro.SED.rda)
# 
# # Find the adjusted R2 of the model with the retained env
# # variables
# RsquareAdj(micro.SED.rda)$adj.r.squared
# #calculating the proportion of the variation of  Y explained by the variables in  X
# 
# #test model significance
# anova.cca(micro.SED.rda, step = 1000)
# 
# #You can also test the significance of each variable
# anova.cca(micro.SED.rda, step = 1000, by = "term")
# 
# ordiplot(micro.SED.rda, scaling=1)
# 
# 
# 
# ####################################
# ##model optimization based on highest R² adj
# ####################################
# 
# micro.SED.rda.test1 <- rda(SED_micro_RDA_hel ~ ecotope + Season, data = Descriptor_SED_reduced)
# RsquareAdj(micro.SED.rda.test1)$adj.r.squared
# 
# 
# micro.SED.rda.test2 <- rda(SED_micro_RDA_hel ~ ecotope + Season + meandering, data = Descriptor_SED_reduced)
# RsquareAdj(micro.SED.rda.test2)$adj.r.squared
# 
# micro.SED.rda.test3 <- rda(SED_micro_RDA_hel ~ ecotope + Season + meandering + pop_dens, data = Descriptor_SED_reduced)
# RsquareAdj(micro.SED.rda.test3)$adj.r.squared
# 
# micro.SED.rda.test4 <- rda(SED_micro_RDA_hel ~ ecotope + Season + meandering + pop_dens + Depth.Sample..m., data = Descriptor_SED_reduced)
# RsquareAdj(micro.SED.rda.test4)$adj.r.squared
# 
# #adding sediment type/TotalPointDischarge/Natural_NonNatural/Depth.Sample..m. didn't improve the model
# 
# 
# #test3 is the best model and the final model
#       

