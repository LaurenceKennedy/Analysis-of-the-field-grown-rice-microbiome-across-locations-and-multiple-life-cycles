--------------------------------------------------
#R Library preparation + data loading and cleaning
--------------------------------------------------
  
#load in the relevant packages for dataset analysis
library(vegan)
library(vioplot)
library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)

#Load in the OTU dataset and metadata
OTU_data<-read.delim("C:/Users/laure/Documents/R files/Rice_Microbiome_Project/data/raw_data/lc_study_otu_table/lc_study_otu_table.tsv",header=T)
n=unlist(OTU_data)
Meta<-read.delim("C:/Users/laure/Documents/R files/Rice_Microbiome_Project/data/raw_data/lc_study_mapping_files/lc_study_mapping_file.tsv",header=T)

#logical matrix to classify the OTU_data dataset
OTU_presence <- OTU_data > 0
prevalence <- rowSums(OTU_presence) / ncol(OTU_data)

# Removing OTUs with >= 5% prevalence
filtered_otus <- which(prevalence >= 0.05)
OTU_data_filtered <- OTU_data[filtered_otus, ]

#Defining the core OTUs present in 80% of samples (302/30916)
presence_threshold.8 <- 0.8
core_OTUs <- which(prevalence >= presence_threshold.8)
core_otu_table <- OTU_data_filtered[core_OTUs, ]



---------------------------------------
#Analysis of cross season data for M206
---------------------------------------
  
#Creating data subset for M206 cultivar
Meta_M206 <- Meta[Meta$Cultivar == "M206", ]
OTU_data_M206 <- OTU_data_filtered[, Meta_M206$SampleID]
CultivarM<-which(Meta$Cultivar=="M206")
#Creating data subset for M206 cultivar
Meta_2016 <- Meta[Meta$Season == "2016", ]
OTU_data_2016 <- OTU_data_filtered[, Meta_2016$SampleID]
Season2016=which(Meta$Season=='2016')

#Transforming dataset into binary data
bin_M206_OTU_data=OTU_data_filtered[,CultivarM]
bin_M206_OTU_data[bin_M206_OTU_data>1]=1
core_bin=core_otu_table[,CultivarM]
core_bin[core_bin>1]=1


#Multidimensional scaling analysis of dataset
bray_dist <- vegdist(t(bin_M206_OTU_data), method = "bray")
mds_result <- cmdscale(bray_dist, k = 3)
# Create a data frame with coordinates
mds_df <- as.data.frame(mds_result)
mds_df$SampleID <- rownames(mds_result)  
# Merge with metadata using the correct column
mds_merged <- merge(mds_df, Meta_M206, by = "SampleID")



----------------------------------------
#Cross-Cultivar analysis of 2016 samples
----------------------------------------
  
#Creating data subset for 2016 season
Meta_2016 <- Meta[Meta$Season == "2016", ]
OTU_data_2016 <- OTU_data_filtered[, Meta_2016$SampleID]
Season2016=which(Meta$Season=='2016')

#Transforming dataset into binary data
bin_2016_OTU_data=OTU_data_filtered[,Season2016]
bin_2016_OTU_data[bin_2016_OTU_data>1]=1
core_bin2=core_otu_table[,Season2016]
core_bin2[core_bin2>1]=1


#Multidimensional scaling analysis of dataset
bray_dist2 <- vegdist(t(bin_2016_OTU_data), method = "bray")
mds_result2 <- cmdscale(bray_dist2, k = 3)
# Create a data frame with coordinates
mds_df2 <- as.data.frame(mds_result2)
mds_df2$SampleID <- rownames(mds_result2)  
# Merge with metadata using the correct column
mds_merged2 <- merge(mds_df2, Meta_2016, by = "SampleID")



---------------------------------------------
#Fisher's alpha diversity metric calculations
---------------------------------------------  
  
#M206 diversity
  
alphaFM206<-array()
for(i in 1:(ncol(OTU_data_M206))){
  n=OTU_data_M206[,i]
  n=n[n>0]
  alphaFM206[i]=fisher.alpha(n)
}
names(alphaFM206)=colnames(OTU_data_M206)[1:841]
NF=names(alphaFM206)
rownames(Meta)=Meta$SampleID

# Ensure sample names match between alphaFM206 and Meta
alpha_df <- data.frame(SampleID = names(alphaFM206),
                       FisherAlpha = alphaFM206)

# Merge with metadata and group age values into 'bins'
alpha_merged <- merge(alpha_df, Meta_M206, by = "SampleID")
alpha_merged$AgeBin <- cut(alpha_merged$Age,
                           breaks = seq(0, 140, by = 20),
                           include.lowest = TRUE,
                           right = FALSE,
                           labels = c("0-19", "20-39", "40-59", "60-79", "80-99", "100-119", "120-139"))


#2016 diversity

alphaF2016<-array()
for(i in 1:(ncol(OTU_data_2016))){
  n=OTU_data_2016[,i]
  n=n[n>0]
  alphaF2016[i]=fisher.alpha(n)
}
names(alphaF2016)=colnames(OTU_data_2016)[1:788]
NF2=names(alphaF2016)
rownames(Meta)=Meta$SampleID

# Ensure sample names match between alphaFM206 and Meta
alpha_df2 <- data.frame(SampleID = names(alphaF2016),
                        FisherAlpha = alphaF2016)

# Merge with metadata and group age values into 'bins'
alpha_merged2 <- merge(alpha_df2, Meta_2016, by = "SampleID")
alpha_merged2$AgeBin <- cut(alpha_merged2$Age,
                            breaks = seq(0, 140, by = 20),
                            include.lowest = TRUE,
                            right = FALSE,
                            labels = c("0-19", "20-39", "40-59", "60-79", "80-99", "100-119", "120-139"))



-----------------------------------------------------------------------------------------------
#Statistical analysis of PCoA and Diversity results (Linear mixed effects regression modelling)
-----------------------------------------------------------------------------------------------

#M206 lmer modeling across all 3 axis
  
PCoAlmerk1=lmer(mds_result[,1]~ Age + Season + FinalDepth + (1|Compartment),data=Meta_M206)
summary(PCoAlmerk1)
ranef(PCoAlmerk1)
rand(PCoAlmerk1)
r.squaredGLMM(PCoAlmerk1)

PCoAlmerk2=lmer(mds_result[,2]~ Age + Season + FinalDepth + (1|Compartment),data=Meta_M206)
summary(PCoAlmerk2)
ranef(PCoAlmerk2)
rand(PCoAlmerk2)
r.squaredGLMM(PCoAlmerk2)

PCoAlmerk3=lmer(mds_result[,3]~ Age + Season + FinalDepth + (1|Compartment),data=Meta_M206)
summary(PCoAlmerk3)
ranef(PCoAlmerk3)
rand(PCoAlmerk3)
r.squaredGLMM(PCoAlmerk3)

#Test for whether assumptions are met. Vary to test k1, k2 and k3
plot(PCoAlmerk3)            # Residuals vs fitted
qqnorm(resid(PCoAlmerk3))   # Q-Q plot for normality
qqline(resid(PCoAlmer))



#2016 lmer modeling across all 3 axis

PCoAlmer2k1=lmer(mds_result2[,1]~ Age + Site + Cultivar + FinalDepth + (1|Compartment),data=Meta_2016)
summary(PCoAlmer2k1)
ranef(PCoAlmer2k1)
rand(PCoAlmer2k1)
r.squaredGLMM(PCoAlmer2k1)

PCoAlmer2k2=lmer(mds_result2[,2]~ Age + Site + Cultivar + FinalDepth + (1|Compartment),data=Meta_2016)
summary(PCoAlmer2k2)
ranef(PCoAlmer2k2)
rand(PCoAlmer2k2)
r.squaredGLMM(PCoAlmer2k2)

PCoAlmer2k3=lmer(mds_result2[,3]~ Age + Site + Cultivar + FinalDepth + (1|Compartment),data=Meta_2016)
summary(PCoAlmer2k3)
ranef(PCoAlmer2k3)
rand(PCoAlmer2k3)
r.squaredGLMM(PCoAlmer2k3)

#Test for whether assumptions are met. Vary to test k1, k2 and k3
plot(PCoAlmer2k1)            # Residuals vs fitted
qqnorm(resid(PCoAlmer2k1))   # Q-Q plot for normality
qqline(resid(PCoAlmer2k1))



#Bayesian lmer model in case dataset is not behaving.
fdlmer=blmer(mds_merged2[,"V2"]~ Age + (1|Compartment) + (1|Cultivar) + (1|Site) ,data=Meta_2016)
summary(fdlmer)
ranef(fdlmer)
rand(fdlmer)
r.squaredGLMM(fdlmer)


#save a fdlmer for both axis and use the fitted function on each object. Can then plot the fitted values against each other
#this should produce a similar graph to the MDS of 2016 Samples by Cultivar graph with 4 clumped groups
fitted(fdlmer)
plot(fitted(fdlmer),fitted(fdlmerforsecondaxis))

fitted(PCoAlmerk1)
fitted(PCoAlmerk2)
fitted(PCoAlmerk3)
plot(fitted(PCoAlmerk1),fitted(PCoAlmerk3))








-------------
#WIP Analyses
-------------


#Core microbiome analysis yet to do (for M206 cultivar)

# Transpose OTU table to long format with sample metadata
otu_long <- as.data.frame(t(OTU_data_M206))
otu_long$SampleID <- rownames(otu_long)
otu_meta <- left_join(otu_long, Meta_M206, by = "SampleID")

# Pivot longer to calculate presence per OTU per season
otu_longer <- pivot_longer(otu_meta, 
                           cols = 1:ncol(OTU_data_M206),  # Adjust if needed
                           names_to = "OTU",
                           values_to = "abundance"
)
otu_longer <- mutate(otu_longer, presence = abundance > 0)

# Calculate prevalence per season
core_by_season <- otu_longer %>%
  group_by(year, OTU) %>%
  summarize(prevalence = sum(presence) / n(), .groups = 'drop') %>%
  filter(prevalence >= 0.8)  # Core threshold (80%)

# Example: list of core OTUs in each season
core_2014 <- core_by_season$OTU[core_by_season$year == 2014]
core_2015 <- core_by_season$OTU[core_by_season$year == 2015]
core_2016 <- core_by_season$OTU[core_by_season$year == 2016]






