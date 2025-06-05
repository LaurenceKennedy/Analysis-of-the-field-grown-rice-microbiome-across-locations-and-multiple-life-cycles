#load in the relevant packages for dataset analysis
library(vegan)
library(vioplot)
library(lme4)
library(lmerTest)
library(MuMIn)

#Load in the OTU dataset and metadata
OTU_data<-read.delim("C:/Users/laure/Documents/R files/Rice_Microbiome_Project/data/raw_data/lc_study_otu_table/lc_study_otu_table.tsv",header=T)
n=unlist(OTU_data)
Meta<-read.delim("C:/Users/laure/Documents/R files/Rice_Microbiome_Project/data/raw_data/lc_study_mapping_files/lc_study_mapping_file.tsv",header=T)

#Tables are useful to give scope of dataset and indicate avenues of interest. For example:

table(Meta[,"Cultivar"],Meta[,"Season"])
table(Meta[,"Age"],Meta[,"Cultivar"])
table(Meta[,"Cultivar"],Meta[,"Compartment"],Meta[,"Season"])
# the above tables indicate where the bulk of data that is relevant to our hypotheses are held:
# across cultivar within the 2016 season    OR
# across season within the M206 cultivar 

#Store objects that will allow analysis of these reduced facets of the dataset
CultivarM=which(Meta$Cultivar=="M206")
Season2016=which(Meta$Season=='2016')



#OTU dataset is very large, explore ways that it could be reduced
sum(rowSums(OTU_data)==0)
sum(colSums(OTU_data)==0)

binaryOTU_data=OTU_data[,CultivarM]
binaryOTU_data[binaryOTU_data>1]=1

cultivar_OTU_data_filtered=OTU_data_filtered[,CultivarM]
cultivar_OTU_data_filtered[cultivar_OTU_data_filtered>1]=1
table(rowSums(binaryOTU_data))
table(colSums(binaryOTU_data))




#setup for PCA on reduced OTU data only using rows with large reads
w=which(rowSums(binaryOTU_data)>420)
k=cmdscale(vegdist(t(binaryOTU_data[w,])))


compartment=as.integer(as.factor(Meta$Compartment[CultivarM]))
table(compartment)
jpeg('PCA of bacterial OTUs based on compartment.jpg') #code for saving plot
plot(k,xlab="Dim.1",ylab="Dim.2",main="PCA of bacterial OTU presence in M206 rice microbiomes",pch=compartment,col=hsv(h=(compartment-1)/5))
legend(-0.285,0.17,legend=c('Bulk Soil','Endosphere','Rhizoplane','Rhizosphere'),box.lwd=0,bty='n',cex=0.6,pch=c(1,2,3,4),col=hsv(h=(1:4-1)/5))
dev.off()

summary(lm(k[,1]~as.factor(compartment)))
fdlmer=lmer(k[,1]~FinalDepth + (1|Compartment) + (1|Age) + (1|Season),data=Meta_M206)
summary(fdlmer)
ranef(fdlmer)
rand(fdlmer)
r.squaredGLMM(fdlmer)


summary(lm(k[,1]~as.factor(compartment)))
fdlmer=lmer(k[,2]~ FinalDepth + Age + I(Age^2) + (1|Compartment) + (1|Season),data=Meta_M206)
summary(fdlmer)
ranef(fdlmer)
rand(fdlmer)
r.squaredGLMM(fdlmer)
#ranef results are comparable to each other even across factor being analysed
#mixed effects modelling can produce the same results as an ANOVA but also giving fixed effects (random effects akin to ANOVA results) 
#important to compare the results of both dimensions K1 and K2, explain 
#can look for interactions by using * instead of + when looking at 
summary(lm(k[,1]~as.factor(compartment)))
slmer=lmer(k[,2]~Season + (1|FinalDepth) + (1|Age) + (1|compartment),data=Meta_M206)
summary(slmer)
ranef(slmer)
rand(slmer)
r.squaredGLMM(slmer)
culseason=as.numeric(unlist(ranef(slmer)$Age))
plot(culseason)












#Defining the core OTUs present in 80% of samples (302/30916)
presence_threshold <- 0.8
otu_presence <- OTU_data > 0  # Logical matrix: TRUE where abundance > 0
prevalence <- rowSums(otu_presence) / ncol(OTU_data)

core_OTUs <- which(prevalence >= presence_threshold)
core_otu_table <- OTU_data[core_OTUs, ]

# Step 2: Filter OTUs with >= 5% prevalence
filtered_otus <- which(prevalence >= 0.05)
OTU_data_filtered <- OTU_data[filtered_otus, ]



# Bray-Curtis dissimilarity
bray <- vegdist(t(OTU_data), method = "bray")

# PCoA
pcoa <- cmdscale(bray, eig = TRUE, k = 2)
pcoa_df <- data.frame(SampleID = rownames(pcoa$points), 
                      PC1 = pcoa$points[,1], 
                      PC2 = pcoa$points[,2])
pcoa_df <- left_join(pcoa_df, metadata, by = "SampleID")

# Plot
library(ggplot2)
ggplot(pcoa_df, aes(x = PC1, y = PC2, color = age)) +
  geom_point() +
  facet_wrap(~year) +
  theme_minimal()













#Shannons H Diversity index for alphaF comparison
shannonH<-array()
for(i in 1:(ncol(OTU_data_filtered)-1)){
  n=OTU_data[,i]
  n=n[n>0]
  shannonH[i]= diversity(n,index="shannon")
}
names(shannonH)=colnames(OTU_data_filtered)[1:1511]
NH=names(shannonH)
rownames(Meta)=Meta$SampleID

#Simpsons D 
simpsonD<-array()
for(i in 1:(ncol(OTU_data_filtered)-1)){
  n=OTU_data[,i]
  n=n[n>0]
  simpsonD[i]= diversity(n,index="invsimpson")
}
names(simpsonD)=colnames(OTU_data_filtered)[1:1511]
ND=names(simpsonD)
rownames(Meta)=Meta$SampleID


is.data.frame(Meta)
summary(lm(log(alphaF)~Cultivar + Season + Compartment,data=Meta))
summary(lm(log(alphaF)~Cultivar + Season + Compartment + Age + Lane + Plot,data=Meta))

#Diversity indexs looking at only M206 cultivar
MetaM206=Meta[Meta$Cultivar=="M206",]
alphaFM206=alphaF[CultivarM]
shannonHM206=exp(shannonH[CultivarM])
simpsonDM206=simpsonD[CultivarM]

AFlmer=lmer(log(alphaFM206)~FinalDepth + (1|Compartment) + (1|Age) + (1|Season),data=Meta_M206)
SHlmer=lmer(log(shannonHM206)~FinalDepth + (1|Compartment) + (1|Age) + (1|Season),data=Meta_M206)
SDlmer=lmer(log(simpsonDM206)~FinalDepth + (1|Compartment) + (1|Age) + (1|Season),data=Meta_M206)

summary(AFlmer)
ranef(AFlmer)
rand(AFlmer)
r.squaredGLMM(AFlmer)
plot(alphaFM206,exp(fitted(AFlmer)),log="xy")
abline(0,1)

plot(alphaFM206,exp(fitted(AFlmer)),log="xy")
abline(0,1)


summary(SHlmer)
ranef(SHlmer)
rand(SHlmer)
r.squaredGLMM(SHlmer)
plot(alphaFM206,fitted(AFlmer))

summary(SDlmer)
ranef(SDlmer)
rand(SDlmer)
r.squaredGLMM(SDlmer)
plot(alphaFM206,fitted(AFlmer))

culage=as.numeric(unlist(ranef(AFlmer)$Age))
cor.test(1:32,culage,method="s")
#spearmans correlation test that proves that alphaF increases with plant age
plot(culage)
#must be cleaned up ^)

culage2=as.numeric(unlist(ranef(SHlmer)$Age))
cor.test(1:32,culage2,method="s")
#spearmans correlation test that proves that alphaF increases with plant age
plot(culage2)
#must be cleaned up ^)

#tested other two variables and found no interesting pattern
culseason=as.numeric(unlist(ranef(AFlmer)$Season))
cor.test(1:3,culseason,method="s")
plot(culseason)
culcompart=as.numeric(unlist(ranef(AFlmer)$Compartment))
cor.test(1:4,culcompart,method="s")
plot(culcompart)










#checking variables that were throwing off the linear mixed effects modeling
sum(is.na(MetaM206$FinalDepth))
sum(is.na(MetaM206$Compartment))
sum(is.na(MetaM206$Age))
sum(is.na(MetaM206$Plot))



#spearmans correlation
cor.test(alphaF,Meta[NF,"FinalDepth"],method="s")
cor.test(alphaF,Meta[NF,"Age"],method="s")
#data is too skewed when not logged (asymettric distribution with outliers, log transformation should be used)
hist(Meta[NF,"FinalDepth"])
hist(log(Meta[NF,"FinalDepth"]))
hist(Meta[NF,"Age"])
hist(log(Meta[NF,"Age"]))






#Testing diversity across categories - Vary 'h' to alter the vioplot structure
vioplot(alphaF~Meta[NF,"Compartment"],h=5)
vioplot(alphaF~Meta[NF,"Cultivar"],h=10)
vioplot(alphaF~Meta[NF,"Season"],h=5)


#visualising the diversity analysis outputs

plot(alphaF,Meta[NF,"FinalDepth"])
plot(shannonH,Meta[NH,"FinalDepth"])
plot(simpsonD,Meta[ND,"FinalDepth"])
#however, 'depth' relates to depth of sequencing reads, not wholly useful

plot(alphaF,Meta[NF,"Age"])
plot(shannonH,Meta[NH,"Age"])
plot(simpsonD,Meta[ND,"Age"])

plot(alphaF,Meta[NF,"Season"])
plot(shannonH,Meta[NH,"Season"])
plot(simpsonD,Meta[ND,"Season"])






