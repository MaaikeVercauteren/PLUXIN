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
library(patchwork)
library(lubridate)
library(corrplot)
library("ggpubr")
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
#resutls cluster
result_HCPC_micro_SED<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Results/result_HCPC_micro_SED.csv"))

#adding sampling areas
SamplingArea<-as.data.frame(read_xlsx("PLUXIN-FinalAnalysis/Raw data/SamplingAreas.xlsx"))

## Fix some problems in colnames
colnames(SamplingArea) <- make.names(colnames(SamplingArea))

#needed for RDA
cluster<-result_HCPC_micro_SED%>%
  select(Unique.Sample.Identifier, clust)
write.csv(cluster, "PLUXIN-FinalAnalysis/Final dataset/Cluster_mic_sed.csv")
######################################################################################################
####Visualization clusters
######################################################################################################


############################################
####    CLUSTER - Bargraphs             ####
############################################
colnames(result_HCPC_micro_SED)

result_bar<-result_HCPC_micro_SED%>%
  select(c("PP_total_conc","PE_total_conc","PES_total_conc","PS_total_conc","PAM_total_conc",
           "SC1_total_conc","SC2_total_conc" ,"SC3_total_conc","SC4_total_conc",
           "Loc_total_conc","avgLWratio","clust"))


##change the column names according to the variable's category 

# Define the columns to prefix
columns_PM_ <- c("PP_total_conc","PE_total_conc","PES_total_conc","PS_total_conc","PAM_total_conc")

columns_SC_ <- c("SC1_total_conc", "SC2_total_conc","SC3_total_conc", "SC4_total_conc") 

columns_Conc_ <- c("Loc_total_conc")

columns_shape_ <- c("avgLWratio")

# Function to add the prefix
new_column_names <- sapply(names(result_bar), function(c) {
  if (c %in% columns_PM_) {
    paste0("PM_", c)
  } else if (c %in% columns_SC_) {
    paste0("SC_", c)
  } else if (c %in% columns_Conc_) {
    paste0("Conc_", c)
  } else if (c %in% columns_shape_) {
    paste0("Shape_", c)
  } else {
    c
  }
})

# Assign the new column names to the dataframe
names(result_bar) <- new_column_names

# Convert data to long format for ggplot2
data_long <- melt(result_bar, id.vars = "clust")

data_long<-data_long%>%
  mutate(group=ifelse(grepl("PM_", variable), "PM",
                      ifelse(grepl("SC_", variable), "SC",
                             ifelse(grepl("Conc_", variable), "Conc",
                                    ifelse(grepl("Shape_", variable), "Shape", NA)))))




data_long$clust <- factor(data_long$clust, levels = c("6", "1","2","5","3","4"))

options(scipen=999)

data_mean<-data_long%>%
  group_by(clust,group, variable)%>%
  summarise(value_mean=mean(value))


## bargraph all variables
#not relative but absolute
ggplot(data_mean, aes(x = variable, y = value_mean, fill = clust)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clust, nrow=1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

#bargraph relative all variables
data_mean<-data_mean%>%
  group_by(clust,group)%>%
  mutate(value_rel= (value_mean / sum(value_mean)) * 100)%>%
  ungroup()

ggplot(data_mean, aes(x = variable, y = value_rel, fill = clust)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clust, nrow=1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Relative Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")
#doesn't make sense for conc and shape


## bargraph polymertypes
data_mean_PM<-data_mean%>%
  filter(group=="PM")%>%
  group_by(clust)%>%
  mutate(value_rel= (value_mean / sum(value_mean)) * 100)%>%
  ungroup()

ggplot(data_mean_PM, aes(x = variable, y = value_rel, fill = clust)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clust, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Relative Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

## bargraph size groups
data_mean_SC<-data_mean%>%
  filter(group=="SC")%>%
  group_by(clust)%>%
  mutate(value_rel= (value_mean / sum(value_mean)) * 100)%>%
  ungroup()

ggplot(data_mean_SC, aes(x = variable, y = value_rel, fill = clust)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clust, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Relative Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

## bargraph shape
data_mean_Shape<-data_mean%>%
  filter(group=="Shape")

ggplot(data_mean_Shape, aes(x = variable, y = value_mean, fill = clust)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clust, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")

## bargraph concentration 
data_mean_Conc<-data_mean%>%
  filter(group=="Conc")

ggplot(data_mean_Conc, aes(x = variable, y = value_mean, fill = clust)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~clust, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")





# ###########################################
# #visualization PCA plot
# ###########################################
result_HCPC_micro_SED_PCA<-as.data.frame(read_csv("PLUXIN-FinalAnalysis/Results/result_HCPC_micro_SED.csv"))


#merging with dataset containing sample areas
result_HCPC_micro_SED_PCA<-result_HCPC_micro_SED_PCA%>%
  left_join(SamplingArea, by ="Unique.Sample.Identifier")

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(result_HCPC_micro_SED_PCA) <- result_HCPC_micro_SED_PCA[,"Unique.Sample.Identifier"]


colnames(result_HCPC_micro_SED_PCA)

##selection of variables that were significantly associated with the clusters
result_HCPC_micro_SED_PCA<-result_HCPC_micro_SED_PCA%>%
  select (c("PP_total_conc","PE_total_conc","PES_total_conc","PS_total_conc","PAM_total_conc",
            "SC1_total_conc","SC2_total_conc" ,"SC3_total_conc","SC4_total_conc",
            "Loc_total_conc","avgLWratio",
            "ecotope","Sediment.type",
            "recreation...km..","urban...km..","RWZI..nr.","Depth.Sample..m.","mean_windspeed_mean",
            "width","agriculture..km..","industry...km..","waste...km..",
            "pop_dens","human.foot.print","mean_temp_mean", "clust",  "Sampling.Area.general"))


#defining quantitative and qualitative variables
quanti.sup.final<-c("recreation...km..","urban...km..","RWZI..nr.","Depth.Sample..m.","mean_windspeed_mean",
                    "width","agriculture..km..","industry...km..","waste...km..",
                    "pop_dens","human.foot.print","mean_temp_mean")
quali.sup.final<-c("ecotope","Sediment.type"? "clust",  "Sampling.Area.general")




res.pca<-PCA(result_HCPC_micro_SED_PCA, ncp=6,quanti.sup=quanti.sup.final, quali.sup=c(quali.sup.final, "clust", "Sampling.Area.general"),graph=FALSE)
res.pca$eig

plot.PCA(res.pca,choix='var')
plot.PCA(res.pca,invisible=c('ind.sup'),label =c('ind'), loadings=TRUE)



#pca
fviz_pca_biplot(res.pca,
                obs.scale=1, var.scale=1, alpha=0.5,col.quanti.sup="darkgrey",
                label = c("quali.sup.final", "quanti.sup.final"),
                invisible = "var") +
  geom_point(size=4,aes(color=factor(result_HCPC_micro_SED_PCA$`Sampling.Area.general`),
                        shape=factor(result_HCPC_micro_SED_PCA$clust)))+
  guides(shape = guide_legend(title = "Cluster"),
         colour = guide_legend(title = "Sampling location"))+ labs(title="PCA plot of clusters")


#pca
fviz_pca_biplot(res.pca,
                obs.scale=1, var.scale=1, alpha=0.5,col.quanti.sup="darkgrey",
                label = c("quali.sup.final", "quanti.sup.final"),
                invisible = "var") +
  geom_point(size=4,aes(shape=factor(result_HCPC_micro_SED_PCA$clust)))+
  guides(shape = guide_legend(title = "Cluster"))+ labs(title="PCA plot of clusters")




# ######################################################################################
# ##visualisation contribution of quali variables to PCA
# ######################################################################################

# Extract contributions of variables to the first two principal components
qualisup_contrib <- as.data.frame(res.pca$quali.sup$v.test)

# Keep only PC1 and PC2 (or choose more components if needed)
qualisup_contrib <- qualisup_contrib[, c("Dim.1", "Dim.2")]

# Add variable names as a new column
qualisup_contrib$Variable <- rownames(qualisup_contrib)


#Select variables
qualisup_contrib<-qualisup_contrib%>%
  filter(Variable %in% c("Antropogeen","Diep Subtidaal","Hoogdynamisch sublittoraal","Matig Diep Subtidaal", "Sand", "Sludge", "Sludge+sand"))
                         
                 

# Convert data to long format for plotting
qualisup_contrib_long <- melt(qualisup_contrib, id.vars = "Variable")


ggplot(qualisup_contrib_long, aes(x = Variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~variable, scales = "free_x", nrow=1) +  # Separate plots for PC1 and PC2
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#E69F00")) + # Custom colors
  labs(x = "", y = "Contribution (%)", title = "Variable Contributions to PCA") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

# ######################################################################################
# ##visualisation contribution of quanti to PCA
# ######################################################################################


# Extract contributions of variables to the first two principal components
quantisup_contrib <- as.data.frame(res.pca$quanti.sup$cos2) #same result using cos
#represents the quality of representation for variables on the factor map. It’s calculated as the squared coordinates: var.cos2 = var.coord * var.coord.


# Keep only PC1 and PC2 (or choose more components if needed)
quantisup_contrib <- quantisup_contrib[, c("Dim.1", "Dim.2")]

# Add variable names as a new column
quantisup_contrib$Variable <- rownames(quantisup_contrib)

# Convert data to long format for plotting
quantisup_contrib_long <- melt(quantisup_contrib, id.vars = "Variable")


ggplot(quantisup_contrib_long, aes(x = Variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~variable, scales = "free_x", nrow=1) +  # Separate plots for PC1 and PC2
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#E69F00")) + # Custom colors
  labs(x = "", y = "Contribution (%)", title = "Variable Contributions to PCA") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")



