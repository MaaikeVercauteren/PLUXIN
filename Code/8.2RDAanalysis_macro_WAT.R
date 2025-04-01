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


#Plastic characteristics
#based on dataset with results of the clusters
WAT_macro_RDA<-as.data.frame(read_csv("Final analysis/Results/result_HCPC_macro_WAT.csv"))


###Descriptor data
#based on full Descriptor dataset with additions
Descriptor_WAT_additions<- read.csv("Final analysis/Final dataset/Descriptor/Descriptor_WAT_withadditions.csv")




#############################################
#Necessary changes in plastic dataset
#############################################
WAT_macro_selectionSamples<-WAT_macro_RDA%>%
  select(`Unique.Sample.Identifier`)

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(WAT_macro_RDA) <- WAT_macro_RDA[,"Unique.Sample.Identifier"]

#Plastic parameters 
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

WAT_macro_RDA<-WAT_macro_RDA%>%
  select("PP_total_conc","PE_total_conc","PS_total_conc", "Others_total_conc","Unknown_pm_total_conc",
         "Unknown_shape_total_conc","fragments_total_conc","filaments_total_conc","pellets_total_conc","films_total_conc","foams_total_conc",
         "others_shape_total_conc",
         "SC3_total_conc","SC4_total_conc","SC5_total_conc","SC6_total_conc" ,"SC7_total_conc",
         "Loc_total_conc")


# Hellinger transform the community data
WAT_macro_RDA_hel <- decostand(WAT_macro_RDA, method = "hellinger")




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
#"Mean slope"               


####qualitative 
#"Nearby vegetation"            
#"NaturalBank"                 
#"meandering"                                     
#"ecotope"  

#NAs in diepte (redelijk veel) en regenval en temperatuur
#Diepte weglaten 

Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  select(c("Unique.Sample.Identifier","width", "shortest.distance.from.shore", "RWZI..nr.", "Waste.facilities..nr.",
           "agriculture..km..","industry...km..",
           "recreation...km..","transport...km..","urban...km..","nature...km..",
           "waste...km..","water...km..","human.foot.print",
           "pop_dens",
           "mean_winddirection_mean",
           "tot_precip_mean","mean_temp_mean","mean_windspeed_mean",
           "Active.overflow", "Mean.slope","Nearby.vegetation", "NaturalBank", "meandering", "ecotope", "Season"))


#temp en neerslag: aanvullen met gemiddelde in het seizoen

mean_autumn_temp <- Descriptor_WAT_additions %>%
  filter(Season == "Autumn" & !is.na(mean_temp_mean)) %>%
  summarize(mean_temp = mean(mean_temp_mean)) %>%
  pull(mean_temp)

Descriptor_WAT_additions<- Descriptor_WAT_additions %>%
  group_by(Unique.Sample.Identifier)%>%
  mutate(mean_temp_mean = ifelse(is.na(mean_temp_mean) & Season == "Autumn", 
                                 mean_autumn_temp, 
                                 mean_temp_mean))

mean_autumn_precip <- Descriptor_WAT_additions %>%
  filter(Season == "Autumn" & !is.na(tot_precip_mean)) %>%
  summarize(mean_precip = mean(tot_precip_mean)) %>%
  pull(mean_precip)

Descriptor_WAT_additions <- Descriptor_WAT_additions %>%
  group_by(Unique.Sample.Identifier)%>%
  mutate(tot_precip_mean = ifelse(is.na(tot_precip_mean) & Season == "Autumn", 
                                  mean_autumn_precip, 
                                  tot_precip_mean))
#selection of only macrooplastic samples
Descriptor_WAT_additions<-WAT_macro_selectionSamples%>%
  left_join(Descriptor_WAT_additions, by = "Unique.Sample.Identifier")
 




# Standardize quantitative environmental data
Descriptor_WAT_additions$width <- decostand(Descriptor_WAT_additions$width, method = "standardize")
Descriptor_WAT_additions$width <- round(Descriptor_WAT_additions$width, 6)

Descriptor_WAT_additions$`shortest.distance.from.shore` <- decostand(Descriptor_WAT_additions$`shortest.distance.from.shore`, method = "standardize")
Descriptor_WAT_additions$shortest.distance.from.shore <- round(Descriptor_WAT_additions$shortest.distance.from.shore, 6)

Descriptor_WAT_additions$`RWZI..nr.` <- decostand(Descriptor_WAT_additions$`RWZI..nr.`, method = "standardize")
Descriptor_WAT_additions$RWZI..nr. <- round(Descriptor_WAT_additions$RWZI..nr., 6)

Descriptor_WAT_additions$`Waste.facilities..nr.` <- decostand(Descriptor_WAT_additions$`Waste.facilities..nr.`, method = "standardize")
Descriptor_WAT_additions$Waste.facilities..nr. <- round(Descriptor_WAT_additions$Waste.facilities..nr., 6)

Descriptor_WAT_additions$`agriculture..km..` <- decostand(Descriptor_WAT_additions$`agriculture..km..`, method = "standardize")
Descriptor_WAT_additions$agriculture..km.. <- round(Descriptor_WAT_additions$agriculture..km.., 6)

Descriptor_WAT_additions$`industry...km..` <- decostand(Descriptor_WAT_additions$`industry...km..`, method = "standardize")
Descriptor_WAT_additions$industry...km.. <- round(Descriptor_WAT_additions$industry...km.., 6)

Descriptor_WAT_additions$`transport...km..` <- decostand(Descriptor_WAT_additions$`transport...km..`, method = "standardize")
Descriptor_WAT_additions$transport...km.. <- round(Descriptor_WAT_additions$transport...km.., 6)

Descriptor_WAT_additions$`urban...km..` <- decostand(Descriptor_WAT_additions$`urban...km..`, method = "standardize")
Descriptor_WAT_additions$urban...km.. <- round(Descriptor_WAT_additions$urban...km.., 6)

Descriptor_WAT_additions$`nature...km..` <- decostand(Descriptor_WAT_additions$`nature...km..`, method = "standardize")
Descriptor_WAT_additions$nature...km.. <- round(Descriptor_WAT_additions$nature...km.., 6)

Descriptor_WAT_additions$`recreation...km..` <- decostand(Descriptor_WAT_additions$`recreation...km..`, method = "standardize")
Descriptor_WAT_additions$recreation...km.. <- round(Descriptor_WAT_additions$recreation...km.., 6)

Descriptor_WAT_additions$`waste...km..` <- decostand(Descriptor_WAT_additions$`waste...km..`, method = "standardize")
Descriptor_WAT_additions$waste...km.. <- round(Descriptor_WAT_additions$waste...km.., 6)

Descriptor_WAT_additions$`human.foot.print` <- decostand(Descriptor_WAT_additions$`human.foot.print`, method = "standardize")
Descriptor_WAT_additions$human.foot.print <- round(Descriptor_WAT_additions$human.foot.print, 6)

Descriptor_WAT_additions$`pop_dens` <- decostand(Descriptor_WAT_additions$`pop_dens`, method = "standardize")
Descriptor_WAT_additions$pop_dens <- round(Descriptor_WAT_additions$pop_dens, 6)

Descriptor_WAT_additions$`tot_precip_mean` <- decostand(Descriptor_WAT_additions$`tot_precip_mean`, method = "standardize")
Descriptor_WAT_additions$tot_precip_mean <- round(Descriptor_WAT_additions$tot_precip_mean, 6)

Descriptor_WAT_additions$`mean_temp_mean` <- decostand(Descriptor_WAT_additions$`mean_temp_mean`, method = "standardize")
Descriptor_WAT_additions$mean_temp_mean <- round(Descriptor_WAT_additions$mean_temp_mean, 6)

Descriptor_WAT_additions$`mean_windspeed_mean` <- decostand(Descriptor_WAT_additions$`mean_windspeed_mean`, method = "standardize")
Descriptor_WAT_additions$mean_windspeed_mean <- round(Descriptor_WAT_additions$mean_windspeed_mean, 6)

Descriptor_WAT_additions$`mean_winddirection_mean` <- decostand(Descriptor_WAT_additions$`mean_winddirection_mean` , method = "standardize")
Descriptor_WAT_additions$mean_winddirection_mean <- round(Descriptor_WAT_additions$mean_winddirection_mean, 6)

Descriptor_WAT_additions$`Active.overflow` <- decostand(Descriptor_WAT_additions$`Active.overflow`, method = "standardize")
Descriptor_WAT_additions$Active.overflow <- round(Descriptor_WAT_additions$Active.overflow, 6)

Descriptor_WAT_additions$`Mean.slope` <- decostand(Descriptor_WAT_additions$`Mean.slope`, method = "standardize")
Descriptor_WAT_additions$Mean.slope <- round(Descriptor_WAT_additions$Mean.slope, 6)


class(Descriptor_WAT_additions)

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(Descriptor_WAT_additions)<- Descriptor_WAT_additions[,"Unique.Sample.Identifier"]
#verwijderen kolom sample names
Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
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

##change categorical variables to numeric
Descriptor_numeric <- Descriptor_WAT_additions %>%
  mutate(across(where(is.character), as.factor)) %>%  # Convert characters (except sample identifier) to factors
  mutate(across(where(is.factor), ~ as.numeric(as.factor(.))))


# We can visually look for correlations between variables:
heatmap(abs(cor(Descriptor_numeric)), 
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)), 
        Colv = NA, Rowv = NA)
legend("topright", 
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 1, bty = "n",
       fill = rev(heat.colors(6)))
#correlations between 
#width and urban
#water and distance from shore ==> remove distance
#pop dens and transport 
#active overflow and nature 
#nature and meandering
#mean slope and industry
#ecotope and natural bank==> remove natural bank
#natural bank and nearby vegetation==> remove natural bank
#population density and human foot print ==> remove human foot print

#remove collinear variables
Descriptor_WAT_additions<-Descriptor_WAT_additions%>%
  select(-c(NaturalBank, human.foot.print, shortest.distance.from.shore))

colnames(Descriptor_WAT_additions)

#included descriptor variables:

# "width"                   "RWZI..nr."               "Waste.facilities..nr."  
#"agriculture..km.."       "industry...km.."         "recreation...km.."      
#"transport...km.."        "urban...km.."            "nature...km.."          
#"waste...km.."            "water...km.."            "pop_dens"               
#"mean_winddirection_mean" "tot_precip_mean"         "mean_temp_mean"         
#"mean_windspeed_mean"     "Active.overflow"         "Mean.slope"             
#"Nearby.vegetation"       "meandering"              "ecotope"                
#"Season"                 

####################################
##initial RDA: full model
####################################
set.seed(123) #for reproducibility
#remove NA columns to avoid errors in handling NA
Descriptor_WAT_additions <- na.omit(Descriptor_WAT_additions)
# Initial RDA with ALL of the environmental data
macro.wat.rda <- rda(WAT_macro_RDA_hel ~ ., data = Descriptor_WAT_additions)
summary(macro.wat.rda)

# Find the adjusted R2 of the model with the retained env
# variables
RsquareAdj(macro.wat.rda)$adj.r.squared
#calculating the proportion of the variation of  Y explained by the variables in  X

#test model significance
anova.cca(macro.wat.rda, step = 1000, permutations = 999)

#You can also test the significance of each variable
anova.cca(macro.wat.rda, step = 1000, by = "term", permutations = 999)

# RDA plot
ordiplot(macro.wat.rda, scaling = 2)




####################################
##initial RDA: model with significant variables
####################################
set.seed(123) #for reproducibility
# Initial RDA with ALL of the environmental data
macro.wat.rda <- rda(WAT_macro_RDA_hel ~ meandering+ agriculture..km.. + industry...km..+ 
                       recreation...km..+urban...km..  +mean_windspeed_mean  + mean_winddirection_mean + Season, data = Descriptor_WAT_additions)
summary(macro.wat.rda)

# Find the adjusted R2 of the model with the retained env
# variables
RsquareAdj(macro.wat.rda)$adj.r.squared
#calculating the proportion of the variation of  Y explained by the variables in  X

#test model significance
anova.cca(macro.wat.rda, step = 1000, permutations = 999)

#You can also test the significance of each variable
anova.cca(macro.wat.rda, step = 1000, by = "term", permutations = 999)
#urban no longer significant

ordiplot(macro.wat.rda, scaling=1)

#plot scaling 2
env_scores <- scores(macro.wat.rda, display = "bp", scaling = 2)
env_df <- as.data.frame(env_scores)
env_df$Variable <- rownames(env_df)

#with clusters
ggplot() +
  geom_point(data = as.data.frame(scores(macro.wat.rda, display = "sites", scaling = 1)), 
             aes(x = RDA1, y = RDA2, color = as.factor(result_HCPC_macro_WAT$clust)), size=3) +
  geom_segment(data = env_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = env_df, 
            aes(x = RDA1, y = RDA2, label = Variable), 
            vjust = -1, color = "black") +
  theme_bw() +geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "RDA1 (XX %)", y = "RDA2 (XX %)")+ labs(color='Clusters') 

#with indication of subclusters
ggplot() +
  geom_point(data = as.data.frame(scores(macro.wat.rda, display = "sites", scaling = 1)), 
             aes(x = RDA1, y = RDA2, color = result_HCPC_macro_WAT$cluster), size=3) +
  geom_segment(data = env_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = env_df, 
            aes(x = RDA1, y = RDA2, label = Variable), 
            vjust = -1, color = "black") +
  theme_bw() +geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "RDA1 (XX %)", y = "RDA2 (XX %)")+ labs(color='Clusters and subclusters') 

#plot scaling 1
env_scores <- scores(macro.wat.rda, display = "bp", scaling = 1)
env_df <- as.data.frame(env_scores)
env_df$Variable <- rownames(env_df)

#with clusters
ggplot() +
  geom_point(data = as.data.frame(scores(macro.wat.rda, display = "sites", scaling = 1)), 
             aes(x = RDA1, y = RDA2, color = as.factor(result_HCPC_macro_WAT$clust)), size=3) +
  geom_segment(data = env_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = env_df, 
            aes(x = RDA1, y = RDA2, label = Variable), 
            vjust = -1, color = "black") +
  theme_bw() +geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "RDA1 (XX %)", y = "RDA2 (XX %)")+ labs(color='Clusters') 

#with indication of subclusters
ggplot() +
  geom_point(data = as.data.frame(scores(macro.wat.rda, display = "sites", scaling = 1)), 
             aes(x = RDA1, y = RDA2, color = result_HCPC_macro_WAT$cluster), size=3) +
  geom_segment(data = env_df, 
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black") +
  geom_text(data = env_df, 
            aes(x = RDA1, y = RDA2, label = Variable), 
            vjust = -1, color = "black") +
  theme_bw() +geom_hline(yintercept=0) + geom_vline(xintercept=0) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(x = "RDA1 (XX %)", y = "RDA2 (XX %)")+ labs(color='Clusters and subclusters') 

