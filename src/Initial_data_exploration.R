#download/use R package 'vegan'
install.packages('vegan')
library(vegan)


#note that objects have not yet been properly stored in this script. 


OTU_data<-read.delim("C:/Users/laure/Documents/R data sets/Dryad download/lc_study_otu_table.tsv/lc_study_otu_table.tsv",header=T)
n=unlist(OTU_data)

alphaF<-array()
for(i in 1:ncol(OTU_data)){
  n=OTU_data[,i]
  n=n[n>0]
  alphaF[i]=fisher.alpha(n)
}

PMetaData<-read.delim("C:Users/laure/Documents/R data sets/Dryad download/lc_study_mapping_file.tsv",header=T)

names(alphaF)=colnames(OTU_data)
N=names(alphaF)
rownames(PMetaData)=PMetaData$SampleID

plot(alphaF,PMetaData[N,"FinalDepth"])

#spearmans correlation
cor.test(alphaF,PMetaData[N,"FinalDepth"],method="s")
