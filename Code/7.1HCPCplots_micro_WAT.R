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
result_HCPC_micro_WAT<-as.data.frame(read_csv("Final analysis/Results/result_HCPC_micro_WAT.csv"))

#adding sampling areas
SamplingArea<-as.data.frame(read_xlsx("Final analysis/Raw data/SamplingAreas.xlsx"))

## Fix some problems in colnames
colnames(SamplingArea) <- make.names(colnames(SamplingArea))


######################################################################################################
####Visualization main clusters
######################################################################################################


############################################
####    CLUSTER - Bargraphs             ####
############################################
result_bar<-result_HCPC_micro_WAT%>%
  select(c("PES_total_conc", "PS_total_conc", "PAM_total_conc","Others_total_conc","PE_total_conc",
           "PP_total_conc", "Unknown_pm_total_conc",
           "SC1_total_conc", "SC2_total_conc","SC3_total_conc", "SC4_total_conc", "SC5_total_conc","Loc_total_conc", 
           "avgLWratio","clust"))


##change the column names according to the variable's category 

# Define the columns to prefix
columns_PM_ <- c("PES_total_conc", "PS_total_conc", "PAM_total_conc","Others_total_conc","PE_total_conc",
                 "PP_total_conc", "Unknown_pm_total_conc")

columns_SC_ <- c("SC1_total_conc", "SC2_total_conc","SC3_total_conc", "SC4_total_conc", "SC5_total_conc") 

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

data_long$clust <- factor(data_long$clust, levels = c("3","1", "2"))

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
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) + # Customize colors
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
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) + # Customize colors
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
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) + # Customize colors
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
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) + # Customize colors
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
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) + # Customize colors
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
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")





###########################################
#visualization PCA plot
###########################################
result_HCPC_micro_WAT_PCA<-as.data.frame(read_csv("Final analysis/Results/result_HCPC_micro_WAT.csv"))


#merging with dataset containing sample areas
result_HCPC_micro_WAT_PCA<-result_HCPC_micro_WAT_PCA%>%
  left_join(SamplingArea, by ="Unique.Sample.Identifier")

#veranderen rijnamen (nodig voor PCA/Cluster)
rownames(result_HCPC_micro_WAT_PCA) <- result_HCPC_micro_WAT_PCA[,"Unique.Sample.Identifier"]


colnames(result_HCPC_micro_WAT_PCA)

##selection of variables that were significantly associated with the clusters
result_HCPC_micro_WAT_PCA<-result_HCPC_micro_WAT_PCA%>%
  select (c("PES_total_conc", "PS_total_conc", "PAM_total_conc","Others_total_conc","PE_total_conc",
            "PP_total_conc", "Unknown_pm_total_conc",
            "SC1_total_conc", "SC2_total_conc","SC3_total_conc", "SC4_total_conc", "SC5_total_conc",
            "avgLWratio","Loc_total_conc", 
            "NaturalBank", "meandering", "ecotope", 
            "width", "recreation...km..",  "Waste.facilities..nr.", 
            "nature...km..", "Active.overflow", "mean_temp_mean","pop_dens","transport...km..", "clust", "Sampling.Area.general"))

#defining quantitative and qualitative variables
quanti.sup.final<-c("width", "recreation...km..",  "Waste.facilities..nr.", 
                    "nature...km..", "Active.overflow", "mean_temp_mean","pop_dens","transport...km..")
quali.sup.final<-c("NaturalBank", "meandering", "ecotope")




res.pca<-PCA(result_HCPC_micro_WAT_PCA, ncp=9,quanti.sup=quanti.sup.final, quali.sup=c(quali.sup.final, "clust", "Sampling.Area.general"),graph=FALSE)
res.pca$eig

plot.PCA(res.pca,choix='var')
plot.PCA(res.pca,invisible=c('ind.sup'),label =c('ind'), loadings=TRUE)



#pca
fviz_pca_biplot(res.pca, 
                obs.scale=1, var.scale=1, alpha=0.5,col.quanti.sup="darkgrey",
                label = c("quali.sup", "quanti.sub"), 
                invisible = "var") +
  geom_point(size=4,aes(color=factor(result_HCPC_micro_WAT_PCA$`Sampling.Area.general`),
                        shape=factor(result_HCPC_micro_WAT_PCA$clust)))+
  guides(shape = guide_legend(title = "Cluster"),
         colour = guide_legend(title = "Sampling location"))+ labs(title="PCA plot of clusters")





######################################################################################
##visualisation contribution of quali variables to PCA
######################################################################################

# Extract contributions of variables to the first two principal components
qualisup_contrib <- as.data.frame(res.pca$quali.sup$v.test)

# Keep only PC1 and PC2 (or choose more components if needed)
qualisup_contrib <- qualisup_contrib[, c("Dim.1", "Dim.2")]

# Add variable names as a new column
qualisup_contrib$Variable <- rownames(qualisup_contrib)


#Select variables
qualisup_contrib<-qualisup_contrib%>%
  filter(Variable %in% c("NaturalBank_Yes", "meandering_Yes", "Antropogeen","Diep Subtidaal","Hoogdynamisch sublittoraal","Matig Diep Subtidaal"))

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

######################################################################################
##visualisation contribution of quanti to PCA
######################################################################################


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







######################################################################################################
####Visualization main clusters with subclusters
######################################################################################################



############################################
####    CLUSTER - Defining subclusters  ####
############################################


sum(result_HCPC_micro_WAT$clust==1)#48 samples
sum(result_HCPC_micro_WAT$clust==2)#3 samples
sum(result_HCPC_micro_WAT$clust==3)#11 samples
subcluster1<-c("MIC_W_017", "MIC_W_042", "MIC_W_104", 
               "MIC_W_055",  "MIC_W_019","MIC_W_035",
               "MIC_W_020","MIC_W_044", "MIC_W_077",
               "MIC_W_041","MIC_W_067", "MIC_W_078",
               "MIC_W_040","MIC_W_016", "MIC_W_015")
               
 

subcluster2<-c("MIC_W_046","MIC_W_087", "MIC_W_054",
               "MIC_W_009","MIC_W_027", "MIC_W_028")
  

subcluster3<- c("MIC_W_148", "MIC_W_156", "MIC_W_144","MIC_W_043","MIC_W_026", "MIC_W_018","MIC_W_013",  "MIC_W_145")

subcluster4<-c("MIC_W_146", "MIC_W_011","MIC_W_065", "MIC_W_069", "MIC_W_023","MIC_W_012","MIC_W_064",
               "MIC_W_066","MIC_W_039", "MIC_W_037","MIC_W_072","MIC_W_038","MIC_W_076", "MIC_W_086", "MIC_W_024", 
               "MIC_W_036","MIC_W_088","MIC_W_045","MIC_W_025")


result_HCPC_micro_WAT<-result_HCPC_micro_WAT%>%
  mutate(cluster = ifelse(clust== 2 |clust == 3, as.factor(result_HCPC_micro_WAT$clust),
                          ifelse(result_HCPC_micro_WAT$Unique.Sample.Identifier %in% subcluster1, "S1", 
                                 ifelse(result_HCPC_micro_WAT$Unique.Sample.Identifier %in% subcluster2, "S2", 
                                        ifelse(result_HCPC_micro_WAT$Unique.Sample.Identifier %in% subcluster3, "S3", "S4")))))
                 



############################################
####    CLUSTER - Bargraphs             ####
############################################
result_bar<-result_HCPC_micro_WAT%>%
  select(c("PES_total_conc", "PS_total_conc", "PAM_total_conc","Others_total_conc","PE_total_conc",
           "PP_total_conc", "Unknown_pm_total_conc",
           "SC1_total_conc", "SC2_total_conc","SC3_total_conc", "SC4_total_conc", "SC5_total_conc","Loc_total_conc", 
           "avgLWratio","cluster"))


##change the column names according to the variable's category 

# Define the columns to prefix
columns_PM_ <- c("PES_total_conc", "PS_total_conc", "PAM_total_conc","Others_total_conc","PE_total_conc",
                 "PP_total_conc", "Unknown_pm_total_conc")

columns_SC_ <- c("SC1_total_conc", "SC2_total_conc","SC3_total_conc", "SC4_total_conc", "SC5_total_conc") 

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
data_long <- melt(result_bar, id.vars = "cluster")

data_long<-data_long%>%
  mutate(group=ifelse(grepl("PM_", variable), "PM",
                      ifelse(grepl("SC_", variable), "SC",
                             ifelse(grepl("Conc_", variable), "Conc",
                                    ifelse(grepl("Shape_", variable), "Shape", NA)))))

data_long$cluster <- factor(data_long$cluster, levels = c("3","S1", "S2", "S3", "S4", "2"))

options(scipen=999)

data_mean<-data_long%>%
  group_by(cluster,group, variable)%>%
  summarise(value_mean=mean(value))


## bargraph all variables
#not relative but absolute
ggplot(data_mean, aes(x = variable, y = value_mean, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cluster, nrow=1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00","yellow", "red", "pink")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

#bargraph relative all variables
data_mean<-data_mean%>%
  group_by(cluster,group)%>%
  mutate(value_rel= (value_mean / sum(value_mean)) * 100)%>%
  ungroup()

ggplot(data_mean, aes(x = variable, y = value_rel, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cluster, nrow=1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Relative Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")
#doesn't make sence for conc and shape


## bargraph polymertypes
data_mean_PM<-data_mean%>%
  filter(group=="PM")%>%
  group_by(cluster)%>%
  mutate(value_rel= (value_mean / sum(value_mean)) * 100)%>%
  ungroup()

ggplot(data_mean_PM, aes(x = variable, y = value_rel, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cluster, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00","yellow", "red", "pink")) + # Customize colors
  labs(x = "", y = "Relative Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

## bargraph size groups
data_mean_SC<-data_mean%>%
  filter(group=="SC")%>%
  group_by(cluster)%>%
  mutate(value_rel= (value_mean / sum(value_mean)) * 100)%>%
  ungroup()

ggplot(data_mean_SC, aes(x = variable, y = value_rel, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cluster, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Relative Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14),legend.position = "none")

## bargraph shape
data_mean_Shape<-data_mean%>%
  filter(group=="Shape")

ggplot(data_mean_Shape, aes(x = variable, y = value_mean, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cluster, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00","yellow", "red", "pink")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")

## bargraph concentration 
data_mean_Conc<-data_mean%>%
  filter(group=="Conc")

ggplot(data_mean_Conc, aes(x = variable, y = value_mean, fill = cluster)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~cluster, nrow = 1)+
  coord_flip() +  # Flip for horizontal bars
  scale_fill_manual(values = c("#009E73", "#56B4E9", "#E69F00", "yellow","red", "pink")) + # Customize colors
  labs(x = "", y = "Absolute Concentration") +
  theme_minimal() +
  theme(text = element_text(size = 14), legend.position = "none")



