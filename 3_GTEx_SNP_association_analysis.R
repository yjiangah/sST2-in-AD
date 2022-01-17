#####
##### Code for association anlayses between SNP and ST2 isoform expressions in multiple tissues in GTEx datasets
#####-------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------
##### Input merged files
# -------------------------------------------------------------------------------------
Merge_ST2_transcript = read.table('./GTEx_RNAseq_ST2L_sST2_level.txt', header = T)
Tissue_types = read.table('./GTEx_tissue_types.txt', header = T)
Merge_Phenotype = read.table('./GTEx_ID_Phenotype_Genotype.txt', header = T)

# -------------------------------------------------------------------------------------
##### Seperate expression file by tissue types
# -------------------------------------------------------------------------------------
setwd('./Tissues')
for (i in 1:length(Tissue_types[,1]))
{
  tissue_ST2_transcript<-subset(Merge_ST2_transcript,Merge_ST2_transcript$Note==toString(Tissue_types[i,1]))
  write.table(tissue_ST2_transcript, file = paste(toString(Tissue_types[i,1]),"_GTEx_RNAseq_ST2L_sST2_level",".txt", sep = ""),row.names = F)
}

# -------------------------------------------------------------------------------------
##### Seperate phenotype file by tissue-specific ID list
# -------------------------------------------------------------------------------------
for (i in 1:length(Tissue_types[,1]))
{
  m=1
  tissue_ST2_transcript<-subset(Merge_ST2_transcript,Merge_ST2_transcript$Note==toString(Tissue_types[i,1]))
  tissue_ID<-matrix(,ncol=15,nrow=length(tissue_ST2_transcript[,1]))
  colnames(tissue_ID)<-c("SUBJID","COHORT","SEX","AGE","RACE","ETHNCTY","BMI","PC1","PC2","PC3","PC4","TRISCH","TRCHSTIN","ID","SNP")
  for (j in 1:length(tissue_ST2_transcript[,1]))
  {
    for (k in 1:length(Merge_Phenotype[,1]))
    {
      if (Merge_Phenotype[k,1]==tissue_ST2_transcript[j,8])
      {
        tissue_ID[m,1:15]=as.vector(Merge_Phenotype[k,1:15])
        m=m+1
      }
    }
  }
  tissue_ID=na.omit(tissue_ID)
  write.table(tissue_ID, file = paste(toString(Tissue_types[i,1]),"IDs_GTEx_Phenotype.txt", sep = "_"),row.names = F)
}

# -------------------------------------------------------------------------------------
##### Expression abundance of ST2 isoforms in multiple tissues
# -------------------------------------------------------------------------------------
ST2_tissue_abundance<-matrix(,ncol=5,nrow=length(Tissue_types[,1]))
colnames(ST2_tissue_abundance)<-c("ST2L_mean","ST2L_SE","sST2_mean","sST2_SE","Sample_size")
for (i in 1:length(Tissue_types[,1]))
{
  tissue_ST2_transcript<-subset(Merge_ST2_transcript,Merge_ST2_transcript$Note==toString(Tissue_types[i,1]))
  ST2_tissue_abundance[i,1]=apply(tissue_ST2_transcript[c(2)],2,mean)
  ST2_tissue_abundance[i,2]=apply(tissue_ST2_transcript[c(2)],2,sd) / sqrt(length(tissue_ST2_transcript[,1]))
  ST2_tissue_abundance[i,3]=apply(tissue_ST2_transcript[c(3)],2,mean)
  ST2_tissue_abundance[i,4]=apply(tissue_ST2_transcript[c(3)],2,sd) / sqrt(length(tissue_ST2_transcript[,1]))
  ST2_tissue_abundance[i,5]=length(tissue_ST2_transcript[,1])
}
ST2_isoform_tissue_abundance=cbind(Tissue_types,data.frame(ST2_tissue_abundance))

write.table(ST2_isoform_tissue_abundance, './Summary_GTEx_ST2_isoform_tissue_abundance.txt',row.names = F)

# -------------------------------------------------------------------------------------
##### SNP ~ ST2 isoforms level in multiple tissues
# -------------------------------------------------------------------------------------
setwd('./Tissues')

library(MASS)
library(GenABEL.data)
library(GenABEL)
library(ggplot2)

Tissue_sample_size=read.table( './Summary_GTEx_ST2_isoform_tissue_abundance.txt', header = T)
# Include tissues with sample size > 50
Tissue_sample_size_clean<-subset(Tissue_sample_size,Tissue_sample_size$Sample_size>50)
Tissue_types_clean=Tissue_sample_size_clean[c(1)]

# -------------------------------------------------------------------------------------
# SNP effects on sST2
SNP_effects_on_sST2=data.frame(matrix(,ncol=9,nrow=length(Tissue_types_clean[,1])))
colnames(SNP_effects_on_sST2)<-c("SNP_on_sST2_Estimate","SNP_on_sST2_SE","SNP_on_sST2_T_value","SNP_on_sST2_P_value","ST2L_mean","ST2L_SE","sST2_mean","sST2_SE","Sample_size")

for (i in 1:length(Tissue_types_clean[,1]))
{
  tissue_ST2_transcript=read.table(paste(toString(Tissue_types_clean[i,1]),"_GTEx_RNAseq_ST2L_sST2_level",".txt", sep = ""), header = T)
  tissue_ID=read.table(paste(toString(Tissue_types_clean[i,1]),"IDs_GTEx_Phenotype.txt", sep = "_"), header = T)
  Merge_tissue_file=cbind(tissue_ST2_transcript,tissue_ID)
  
  tissue_sST2_SNP_test<-lmrob(rntransform(ENST00000311734)~SNP+SEX+AGE+RIN+PC1+PC2+PC3+PC4,data=Merge_tissue_file,k.max=900000)
  
  SNP_effects_on_sST2[i,1:4]=as.vector(t(summary(tissue_sST2_SNP_test)$coefficients[2,1:4]))
  SNP_effects_on_sST2[i,5]=apply(tissue_ST2_transcript[c(2)],2,mean)
  SNP_effects_on_sST2[i,6]=apply(tissue_ST2_transcript[c(2)],2,sd) / sqrt(length(tissue_ST2_transcript[,1]))
  SNP_effects_on_sST2[i,7]=apply(tissue_ST2_transcript[c(3)],2,mean)
  SNP_effects_on_sST2[i,8]=apply(tissue_ST2_transcript[c(3)],2,sd) / sqrt(length(tissue_ST2_transcript[,1]))
  SNP_effects_on_sST2[i,9]=length(tissue_ST2_transcript[,1])
}
Summary_SNP_effects_on_sST2=cbind(Tissue_types_clean,SNP_effects_on_sST2)
Summary_SNP_effects_on_sST2$SNP_on_sST2_FDR=p.adjust(Summary_SNP_effects_on_sST2$SNP_on_sST2_P_value,method='fdr')

write.table(Summary_SNP_effects_on_sST2, './Summary_GTEx_SNP_effects_on_sST2.txt',row.names = F)

# -------------------------------------------------------------------------------------
# SNP effects on ST2L
SNP_effects_on_ST2L=data.frame(matrix(,ncol=6,nrow=length(Tissue_types_clean[,1])))
colnames(SNP_effects_on_ST2L)<-c("SNP_on_ST2L_Estimate","SNP_on_ST2L_SE","SNP_on_ST2L_T_value","SNP_on_ST2L_P_value","ST2L_expressing_percentage","sST2_expressing_percentage")

for (i in 1:length(Tissue_types_clean[,1]))
{
  tryCatch({
    tissue_ST2_transcript=read.table(paste(toString(Tissue_types_clean[i,1]),"_GTEx_RNAseq_ST2L_sST2_level",".txt", sep = ""), header = T)
    tissue_ID=read.table(paste(toString(Tissue_types_clean[i,1]),"IDs_GTEx_Phenotype.txt", sep = "_"), header = T)
    Merge_tissue_file=cbind(tissue_ST2_transcript,tissue_ID)
    
    tissue_ST2L_pos=Merge_tissue_file[c(2)]
    tissue_ST2L_pos[tissue_ST2L_pos==0]<-NA
    tissue_ST2L_pos=na.omit(tissue_ST2L_pos)
  
    tissue_sST2_pos=Merge_tissue_file[c(3)]
    tissue_sST2_pos[tissue_sST2_pos==0]<-NA
    tissue_sST2_pos=na.omit(tissue_sST2_pos)
    
    SNP_effects_on_ST2L[i,5]= length(tissue_ST2L_pos[,1]) / length(Merge_tissue_file[,1])
    SNP_effects_on_ST2L[i,6]= length(tissue_sST2_pos[,1]) / length(Merge_tissue_file[,1])
    
    tissue_ST2L_SNP_test<-lmrob(rntransform(ENST00000233954)~SNP+SEX+AGE+RIN+PC1+PC2+PC3+PC4,data=Merge_tissue_file,k.max=900000)
    
    SNP_effects_on_ST2L[i,1:4]=as.vector(t(summary(tissue_ST2L_SNP_test)$coefficients[2,1:4]))
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

SNP_effects_on_ST2_isoforms=cbind(SNP_effects_on_ST2L,SNP_effects_on_sST2)
Summary_SNP_effects_on_ST2_isoforms=cbind(Tissue_types_clean,SNP_effects_on_ST2_isoforms)

write.table(Summary_SNP_effects_on_ST2_isoforms, './Summary_GTEx_SNP_effects_on_ST2_isoforms.txt',row.names = F)








