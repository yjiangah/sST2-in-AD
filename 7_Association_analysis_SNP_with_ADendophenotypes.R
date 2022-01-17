#####
##### Code for association analysis between candidate SNPs and AD-related endophenotypes
#####-------------------------------------------------------------------------------------------------

#####---------------------------------------------------
## SNP ~ Onset age of dementia
#####---------------------------------------------------

library(survival)

# Data input
Data_input=read.table('./Input_file_Onset_age_Phenotype_Genotype.txt',header=T)

# Data stratification by Sex / Apoe4 genotypes
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")

#####----------------------------------------------
### In overall: SNP ~ Onset age of dementia
test_all<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+Gender+PC_1+PC_2+PC_3+PC_4+PC_5,Data_input)
# Test the proportional hazard assumption
test.ph_all=data.frame(as.matrix(cox.zph(test_all)$table))

#####-------------
### In female and male: SNP ~ Onset age of dementia
test_F<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+PC_1+PC_2+PC_3+PC_4+PC_5,Data_female)
test_M<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+PC_1+PC_2+PC_3+PC_4+PC_5,Data_male)

# Test the proportional hazard assumption
test.ph_F=data.frame(as.matrix(cox.zph(test_F)$table))
test.ph_M=data.frame(as.matrix(cox.zph(test_M)$table))

#####-------------
## In female and male APOE4 carriers: SNP ~ Onset age of dementia
test_F_APOE4<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+PC_1+PC_2+PC_3+PC_4+PC_5,Data_female_APOE4)
test_M_APOE4<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+PC_1+PC_2+PC_3+PC_4+PC_5,Data_male_APOE4)

# Test the proportional hazard assumption
test.ph_F_APOE4=data.frame(as.matrix(cox.zph(test_F_APOE4)$table))
test.ph_M_APOE4=data.frame(as.matrix(cox.zph(test_M_APOE4)$table))

#####-------------
## In female and male APOE4 noncarriers: SNP ~ Onset age of dementia
test_F_nonAPOE4<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+PC_1+PC_2+PC_3+PC_4+PC_5,Data_female_nonAPOE4)
test_M_nonAPOE4<-coxph(Surv(as.numeric(AgeDem),as.numeric(Status)) ~ SNP+PC_1+PC_2+PC_3+PC_4+PC_5,Data_male_nonAPOE4)

# Test the proportional hazard assumption
test.ph_F_nonAPOE4=data.frame(as.matrix(cox.zph(test_F_nonAPOE4)$table))
test.ph_M_nonAPOE4=data.frame(as.matrix(cox.zph(test_M_nonAPOE4)$table))

#####----------------------------------------------
## Merge summary and FDR correction
Summary_statistics=data.frame(matrix(,ncol=5,nrow=7))
colnames(Summary_statistics)=c("beta","HR","SE","Z_score","P_value")
rownames(Summary_statistics)=c("overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Summary_statistics[1,1:5]=as.vector(t(summary(test_all)$coefficients[1,1:5]))
Summary_statistics[2,1:5]=as.vector(t(summary(test_M)$coefficients[1,1:5]))
Summary_statistics[3,1:5]=as.vector(t(summary(test_F)$coefficients[1,1:5]))
Summary_statistics[4,1:5]=as.vector(t(summary(test_M_APOE4)$coefficients[1,1:5]))
Summary_statistics[5,1:5]=as.vector(t(summary(test_M_nonAPOE4)$coefficients[1,1:5]))
Summary_statistics[6,1:5]=as.vector(t(summary(test_F_APOE4)$coefficients[1,1:5]))
Summary_statistics[7,1:5]=as.vector(t(summary(test_F_nonAPOE4)$coefficients[1,1:5]))

Summary_statistics$FDR=p.adjust(Summary_statistics$P_value,method='fdr')


#####---------------------------------------------------
## SNP ~ MMSE score
#####---------------------------------------------------
# Data input
Data_input=read.table('./Input_file_MMSE_Phenotype_Genotype.txt',header=T)

# Data stratification by Sex / Apoe4 genotypes
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")

#####---------------------------------
Test_all<-lm(MMSE~SNP+Age+Gender+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_input)
Test_F<-lm(MMSE~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_female)
Test_M<-lm(MMSE~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_male)
Test_F_APOE4<-lm(MMSE~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_female_APOE4)
Test_M_APOE4<-lm(MMSE~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_male_APOE4)
Test_F_nonAPOE4<-lm(MMSE~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_female_nonAPOE4)
Test_M_nonAPOE4<-lm(MMSE~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5, data = Data_male_nonAPOE4)

#####----------------------------------------------
## Merge summary and FDR correction
Merge_summary_test=data.frame(matrix(,ncol=4,nrow=7))
rownames(Merge_summary_test)=c("overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")
colnames(Merge_summary_test)=c("Estimate","SE","T","P_value")

Merge_summary_test[1,1:4]=as.vector(summary(Test_all)$coefficient[2,1:4])
Merge_summary_test[2,1:4]=as.vector(summary(Test_M)$coefficient[2,1:4])
Merge_summary_test[3,1:4]=as.vector(summary(Test_F)$coefficient[2,1:4])
Merge_summary_test[4,1:4]=as.vector(summary(Test_M_APOE4)$coefficient[2,1:4])
Merge_summary_test[5,1:4]=as.vector(summary(Test_M_nonAPOE4)$coefficient[2,1:4])
Merge_summary_test[6,1:4]=as.vector(summary(Test_F_APOE4)$coefficient[2,1:4])
Merge_summary_test[7,1:4]=as.vector(summary(Test_F_nonAPOE4)$coefficient[2,1:4])

Merge_summary_test$FDR=p.adjust(Merge_summary_test$P_value,method='fdr')


#####---------------------------------------------------
## SNP ~ Brain region volumes
#####---------------------------------------------------
library(MASS)
library(GenABEL.data)
library(GenABEL)

# Data input
Data_input=read.table('./Input_file_Brain_region_volumes_Phenotype_Genotype.txt',header=T)

# Data stratification by Sex / Apoe4 genotypes
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")

#####---------------------------------
## Statistics in overall
Test_all=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_all)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_all)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_input[,i])~SNP+Age+Gender+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_input)
  Test_all[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_all$FDR=p.adjust(Test_all$P_value,method='fdr')

#####---------------------------------
## Statistics in female
Test_F=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_F)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_F)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_female[,i])~SNP+Age+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_female)
  Test_F[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F$FDR=p.adjust(Test_F$P_value,method='fdr')

#####---------------------------------
## Statistics in male
Test_M=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_M)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_M)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_male[,i])~SNP+Age+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_male)
  Test_M[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_M$FDR=p.adjust(Test_M$P_value,method='fdr')

#####---------------------------------
## Statistics in female APOE4 carriers
Test_F_APOE4=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_F_APOE4)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_F_APOE4)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_female_APOE4[,i])~SNP+Age+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_female_APOE4)
  Test_F_APOE4[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F_APOE4$FDR=p.adjust(Test_F_APOE4$P_value,method='fdr')

#####---------------------------------
## Statistics in male APOE4 carriers
Test_M_APOE4=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_M_APOE4)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_M_APOE4)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_male_APOE4[,i])~SNP+Age+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_male_APOE4)
  Test_M_APOE4[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_M_APOE4$FDR=p.adjust(Test_M_APOE4$P_value,method='fdr')

#####---------------------------------
## Statistics in female APOE4 noncarriers
Test_F_nonAPOE4=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_F_nonAPOE4)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_F_nonAPOE4)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_female_nonAPOE4[,i])~SNP+Age+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_female_nonAPOE4)
  Test_F_nonAPOE4[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F_nonAPOE4$FDR=p.adjust(Test_F_nonAPOE4$P_value,method='fdr')

#####---------------------------------
## Statistics in male APOE4 noncarriers
Test_M_nonAPOE4=data.frame(matrix(,ncol = 4,nrow = 6))
rownames(Test_M_nonAPOE4)=c("Ventricles","Hippocampus","Entorhinal_cortex","Fusiform","MidTemp","Whole_brain")
colnames(Test_M_nonAPOE4)=c("Estimate","SE","T","P_value")

for (i in 1:6){
  test<-lm(rntransform(Data_male_nonAPOE4[,i])~SNP+Age+ICV+MRI_platform+dementia_stages+PC1+PC2+PC3+PC4+PC5,data=Data_male_nonAPOE4)
  Test_M_nonAPOE4[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_M_nonAPOE4$FDR=p.adjust(Test_M_nonAPOE4$P_value,method='fdr')


#####---------------------------------------------------
## SNP ~ Plasma p-Tau181 and Nfl levels
#####---------------------------------------------------
# Data input
Data_input=read.table('./Input_file_plasma_AD_biomarkers_Phenotype_Genotype.txt',header=T)

# Data stratification by Sex / Apoe4 genotypes
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")

#####---------------------------------
## Statistics in overall
Test_pTau_all<-lm(Ptau181_pg.mL~SNP+Age+Gender+PC1+PC2+PC3+PC4+PC5,data = Data_input)
Test_NfL_all<-lm(NfL_pg.mL~SNP+Age+Gender+PC1+PC2+PC3+PC4+PC5,data = Data_input)

#####---------------------------------
## Statistics in male
Test_pTau_M<-lm(Ptau181_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_male)
Test_NfL_M<-lm(NfL_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_male)

#####---------------------------------
## Statistics in female
Test_pTau_F<-lm(Ptau181_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_female)
Test_NfL_F<-lm(NfL_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_female)

#####---------------------------------
## Statistics in male APOE4 carriers
Test_pTau_M_APOE4<-lm(Ptau181_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_male_APOE4)
Test_NfL_M_APOE4<-lm(NfL_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_male_APOE4)

#####---------------------------------
## Statistics in female APOE4 carriers
Test_pTau_F_APOE4<-lm(Ptau181_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_female_APOE4)
Test_NfL_F_APOE4<-lm(NfL_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_female_APOE4)

#####---------------------------------
## Statistics in male APOE4 noncarriers
Test_pTau_M_APOE4<-lm(Ptau181_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_male_APOE4)
Test_NfL_M_APOE4<-lm(NfL_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_male_APOE4)

#####---------------------------------
## Statistics in female APOE4 noncarriers
Test_pTau_F_nonAPOE4<-lm(Ptau181_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_female_nonAPOE4)
Test_NfL_F_nonAPOE4<-lm(NfL_pg.mL~SNP+Age+PC1+PC2+PC3+PC4+PC5,data = Data_female_nonAPOE4)

#####----------------------------------------------
## Merge summary and FDR correction
Merge_pTau181_summary=data.frame(matrix(,ncol = 4,nrow = 7))
colnames(Merge_pTau181_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_pTau181_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Merge_pTau181_summary[1,1:4]=as.vector(summary(Test_pTau_all)$coefficient[2,1:4])
Merge_pTau181_summary[2,1:4]=as.vector(summary(Test_pTau_M)$coefficient[2,1:4])
Merge_pTau181_summary[3,1:4]=as.vector(summary(Test_pTau_F)$coefficient[2,1:4])
Merge_pTau181_summary[4,1:4]=as.vector(summary(Test_pTau_M_APOE4)$coefficient[2,1:4])
Merge_pTau181_summary[5,1:4]=as.vector(summary(Test_pTau_M_nonAPOE4)$coefficient[2,1:4])
Merge_pTau181_summary[6,1:4]=as.vector(summary(Test_pTau_F_APOE4)$coefficient[2,1:4])
Merge_pTau181_summary[7,1:4]=as.vector(summary(Test_pTau_F_nonAPOE4)$coefficient[2,1:4])

Merge_pTau181_summary$FDR=p.adjust(Merge_pTau181_summary$P_value,method = "fdr")

#------
Merge_NfL_summary=data.frame(matrix(,ncol = 4,nrow = 7))
colnames(Merge_NfL_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_NfL_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Merge_NfL_summary[1,1:4]=as.vector(summary(Test_NfL_all)$coefficient[2,1:4])
Merge_NfL_summary[2,1:4]=as.vector(summary(Test_NfL_M)$coefficient[2,1:4])
Merge_NfL_summary[3,1:4]=as.vector(summary(Test_NfL_F)$coefficient[2,1:4])
Merge_NfL_summary[4,1:4]=as.vector(summary(Test_NfL_M_APOE4)$coefficient[2,1:4])
Merge_NfL_summary[5,1:4]=as.vector(summary(Test_NfL_M_nonAPOE4)$coefficient[2,1:4])
Merge_NfL_summary[6,1:4]=as.vector(summary(Test_NfL_F_APOE4)$coefficient[2,1:4])
Merge_NfL_summary[7,1:4]=as.vector(summary(Test_NfL_F_nonAPOE4)$coefficient[2,1:4])

Merge_NfL_summary$FDR=p.adjust(Merge_NfL_summary$P_value,method = "fdr")


#####---------------------------------------------------
## SNP ~ Abeta plaque load / Microglial coverage of Abeta
#####---------------------------------------------------
# Data input
Data_input=read.table('./Input_file_AD_Brain_Abeta_microglia_Phenotype_Genotype.txt',header=T)

# Data stratification by Sex / Apoe4 genotypes
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")
Data_female_SNP<-subset(Data_female,Data_female$SNP!=0)
Data_female_nonSNP<-subset(Data_female,Data_female$SNP==0)
Data_male_SNP<-subset(Data_male,Data_male$SNP!=0)
Data_male_nonSNP<-subset(Data_male,Data_male$SNP==0)

#####---------------------------------
## Statistics of Abeta load: male AD vs. female AD
Test_in_all<-lm(Abeta_area~SEX+AGE+PMD,data = Data_input)
Test_in_APOE4<-lm(Abeta_area~SEX+AGE+PMD,data = Data_APOE4)
Test_in_nonAPOE4<-lm(Abeta_area~SEX+AGE+PMD,data = Data_nonAPOE4)

Merge_sex_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_sex_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_sex_summary)=c("Overall","APOE4","nonAPOE4")
Merge_sex_summary[1,1:4]=as.vector(summary(Test_in_all)$coefficient[2,1:4])
Merge_sex_summary[2,1:4]=as.vector(summary(Test_in_APOE4)$coefficient[2,1:4])
Merge_sex_summary[3,1:4]=as.vector(summary(Test_in_nonAPOE4)$coefficient[2,1:4])

Merge_sex_summary$FDR=p.adjust(Merge_sex_summary$P_value,method = "fdr")

#####---------------------------------
## Statistics of Abeta load / microglia coverage of Abeta: APOE4 vs. nonAPOE4
# In overall AD
Test_all=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_all)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_all)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_input[,i]~APOE4+AGE+SEX+PMD,data=Data_input)
  Test_all[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_all$FDR=p.adjust(Test_all$P_value,method='fdr')

# In male AD
Test_M=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_M)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_M)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_male[,i]~APOE4+AGE+PMD,data=Data_male)
  Test_M[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_M$FDR=p.adjust(Test_M$P_value,method='fdr')

# In female AD
Test_F=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_F)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_F)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_female[,i]~APOE4+AGE+PMD,data=Data_female)
  Test_F[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F$FDR=p.adjust(Test_F$P_value,method='fdr')

# In female AD SNP carriers
Test_F_SNP=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_F_SNP)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_F_SNP)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_female_SNP[,i]~APOE4+AGE+PMD,data=Data_female_SNP)
  Test_F_SNP[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F_SNP$FDR=p.adjust(Test_F_SNP$P_value,method='fdr')

# In female AD SNP noncarriers
Test_F_nonSNP=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_F_nonSNP)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_F_nonSNP)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_female_nonSNP[,i]~APOE4+AGE+PMD,data=Data_female_nonSNP)
  Test_F_nonSNP[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F_nonSNP$FDR=p.adjust(Test_F_nonSNP$P_value,method='fdr')

#####---------------------------------
## Statistics: SNP ~ Abeta load / microglia coverage of Abeta in female APOE4 carriers and noncarriers
# In female AD APOE4 carriers
Test_F_APOE4=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_F_APOE4)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_F_APOE4)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_female_APOE4[,i]~SNP+AGE+PMD,data=Data_female_APOE4)
  Test_F_APOE4[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F_APOE4$FDR=p.adjust(Test_F_APOE4$P_value,method='fdr')

# In female AD APOE4 noncarriers
Test_F_nonAPOE4=data.frame(matrix(,ncol = 4,nrow = 4))
rownames(Test_F_nonAPOE4)=c("Abeta_area","Abeta_number","Abeta_median_size","Microglial_coverage_Abeta")
colnames(Test_F_nonAPOE4)=c("Estimate","SE","T","P_value")

for (i in 1:4){
  test<-lm(Data_female_nonAPOE4[,i]~SNP+AGE+PMD,data=Data_female_nonAPOE4)
  Test_F_nonAPOE4[i,1:4]=as.vector(summary(test)$coefficient[2,1:4])
}
Test_F_nonAPOE4$FDR=p.adjust(Test_F_nonAPOE4$P_value,method='fdr')



