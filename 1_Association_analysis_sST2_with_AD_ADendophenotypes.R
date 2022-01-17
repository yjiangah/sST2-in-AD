#####
##### Code for association analysis between plasma and CSF sST2 levels and AD/AD-related endophenotypes
#####-------------------------------------------------------------------------------------------------

#####----------------------------------------------
#### Plasma sST2 ~ AD / AD-related endophenotypes
#####----------------------------------------------

# Data input
Data_input=read.table('./Input_file_plasma_sST2_AD_Graymatter_pTau181_NfL.txt',header=T)
Data_input$Phenotype<-relevel(Data_input$Phenotype,"NC")

# Data stratification by AD diagnosis / Sex / Apoe4 genotypes
Data_AD<-subset(Data_input,Data_input$Phenotype=="AD")
Data_NC<-subset(Data_input,Data_input$Phenotype=="NC")
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_female_AD<-subset(Data_AD,Data_AD$Gender=="2")
Data_female_NC<-subset(Data_NC,Data_NC$Gender=="2")
Data_male_AD<-subset(Data_AD,Data_AD$Gender=="1")
Data_male_NC<-subset(Data_NC,Data_NC$Gender=="1")
Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")

#####----------------------------------------------
## Statistics on sST2: male vs female
Test_in_NC<-lm(sST2_ng.mL~Gender+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_NC)
Test_in_AD<-lm(sST2_ng.mL~Gender+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_AD)
Test_in_all<-lm(sST2_ng.mL~Gender+Phenotype+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_input)

Merge_sex_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_sex_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_sex_summary)=c("Overall","NC","AD")
Merge_sex_summary[1,1:4]=as.vector(summary(Test_in_all)$coefficient[2,1:4])
Merge_sex_summary[2,1:4]=as.vector(summary(Test_in_NC)$coefficient[2,1:4])
Merge_sex_summary[3,1:4]=as.vector(summary(Test_in_AD)$coefficient[2,1:4])

Merge_sex_summary$FDR=p.adjust(Merge_sex_summary$P_value,method = "fdr")

#####----------------------------------------------
## Statistics on sST2: APOE4 vs nonAPOE4
Test_in_female_NC<-lm(sST2_ng.mL~APOE4+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_NC)
Test_in_female_AD<-lm(sST2_ng.mL~APOE4+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_AD)
Test_in_male_NC<-lm(sST2_ng.mL~APOE4+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_NC)
Test_in_male_AD<-lm(sST2_ng.mL~APOE4+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_AD)
Test_in_NC<-lm(sST2_ng.mL~APOE4+Gender+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_NC)
Test_in_AD<-lm(sST2_ng.mL~APOE4+Gender+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_AD)
Test_in_all<-lm(sST2_ng.mL~APOE4+Phenotype+Gender+Age+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_input)

Merge_APOE4_summary=data.frame(matrix(,ncol = 4,nrow = 7))
colnames(Merge_APOE4_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_APOE4_summary)=c("Overall","NC","AD","F_NC","F_AD","M_NC","M_AD")
Merge_APOE4_summary[1,1:4]=as.vector(summary(Test_in_all)$coefficient[2,1:4])
Merge_APOE4_summary[2,1:4]=as.vector(summary(Test_in_NC)$coefficient[2,1:4])
Merge_APOE4_summary[3,1:4]=as.vector(summary(Test_in_AD)$coefficient[2,1:4])
Merge_APOE4_summary[4,1:4]=as.vector(summary(Test_in_female_NC)$coefficient[2,1:4])
Merge_APOE4_summary[5,1:4]=as.vector(summary(Test_in_female_AD)$coefficient[2,1:4])
Merge_APOE4_summary[6,1:4]=as.vector(summary(Test_in_male_NC)$coefficient[2,1:4])
Merge_APOE4_summary[7,1:4]=as.vector(summary(Test_in_male_AD)$coefficient[2,1:4])

Merge_APOE4_summary$FDR=p.adjust(Merge_APOE4_summary$P_value,method = "fdr")

#####----------------------------------------------
## Statistics in overall: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_all<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_input)
Test_Graymatter_all<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_input)
Test_pTau_all<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_input)
Test_NfL_all<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_input)

## Statistics in male: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_M<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male)
Test_pTau_M<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male)
Test_NfL_M<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male)
Test_Graymatter_M<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male)

## Statistics in female: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_F<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female)
Test_pTau_F<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female)
Test_NfL_F<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female)
Test_Graymatter_F<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female)

## Statistics in female APOE4 carriers: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_F_APOE4<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_APOE4)
Test_pTau_F_APOE4<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_APOE4)
Test_NfL_F_APOE4<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_APOE4)
Test_Graymatter_F_APOE4<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_APOE4)

## Statistics in male APOE4 carriers: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_M_APOE4<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_APOE4)
Test_pTau_M_APOE4<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_APOE4)
Test_NfL_M_APOE4<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_APOE4)
Test_Graymatter_M_APOE4<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_APOE4)

## Statistics in female APOE4 noncarriers: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_F_nonAPOE4<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_nonAPOE4)
Test_pTau_F_nonAPOE4<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_nonAPOE4)
Test_NfL_F_nonAPOE4<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_nonAPOE4)
Test_Graymatter_F_nonAPOE4<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_female_nonAPOE4)

## Statistics in male APOE4 noncarriers: sST2 ~ AD / Gray matter volume / p-Tau181 / NfL
Test_AD_M_nonAPOE4<-lm(sST2_ng.mL~Phenotype+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_nonAPOE4)
Test_pTau_M_nonAPOE4<-lm(sST2_ng.mL~Ptau181_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_nonAPOE4)
Test_NfL_M_nonAPOE4<-lm(sST2_ng.mL~NfL_pg.mL+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_nonAPOE4)
Test_Graymatter_M_nonAPOE4<-lm(sST2_ng.mL~Gray_matter_.ICV+Age+Gender+Heart_Disease+Hypertension+Diabete_Mellitus+Hyperlipidaemia+Edu_yr+BMI,data = Data_male_nonAPOE4)

#####----------------------------------------------
## Merge summary and FDR correction
# sST2 ~ AD
Merge_AD_summary=data.frame(matrix(,ncol = 4,nrow = 7))
colnames(Merge_AD_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_AD_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Merge_AD_summary[1,1:4]=as.vector(summary(Test_AD_all)$coefficient[2,1:4])
Merge_AD_summary[2,1:4]=as.vector(summary(Test_AD_M)$coefficient[2,1:4])
Merge_AD_summary[3,1:4]=as.vector(summary(Test_AD_F)$coefficient[2,1:4])
Merge_AD_summary[4,1:4]=as.vector(summary(Test_AD_M_APOE4)$coefficient[2,1:4])
Merge_AD_summary[5,1:4]=as.vector(summary(Test_AD_M_nonAPOE4)$coefficient[2,1:4])
Merge_AD_summary[6,1:4]=as.vector(summary(Test_AD_F_APOE4)$coefficient[2,1:4])
Merge_AD_summary[7,1:4]=as.vector(summary(Test_AD_F_nonAPOE4)$coefficient[2,1:4])

Merge_AD_summary$FDR=p.adjust(Merge_AD_summary$P_value,method = "fdr")

# # sST2 ~ Gray matter volume
Merge_Graymatter_summary=data.frame(matrix(,ncol = 4,nrow = 7))
colnames(Merge_Graymatter_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_Graymatter_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Merge_Graymatter_summary[1,1:4]=as.vector(summary(Test_Graymatter_all)$coefficient[2,1:4])
Merge_Graymatter_summary[2,1:4]=as.vector(summary(Test_Graymatter_M)$coefficient[2,1:4])
Merge_Graymatter_summary[3,1:4]=as.vector(summary(Test_Graymatter_F)$coefficient[2,1:4])
Merge_Graymatter_summary[4,1:4]=as.vector(summary(Test_Graymatter_M_APOE4)$coefficient[2,1:4])
Merge_Graymatter_summary[5,1:4]=as.vector(summary(Test_Graymatter_M_nonAPOE4)$coefficient[2,1:4])
Merge_Graymatter_summary[6,1:4]=as.vector(summary(Test_Graymatter_F_APOE4)$coefficient[2,1:4])
Merge_Graymatter_summary[7,1:4]=as.vector(summary(Test_Graymatter_F_nonAPOE4)$coefficient[2,1:4])

Merge_Graymatter_summary$FDR=p.adjust(Merge_Graymatter_summary$P_value,method = "fdr")

# sST2 ~ Plasma p-Tau181 levels
Merge_pTau_summary=data.frame(matrix(,ncol = 4,nrow = 7))
colnames(Merge_pTau_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_pTau_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Merge_pTau_summary[1,1:4]=as.vector(summary(Test_pTau_all)$coefficient[2,1:4])
Merge_pTau_summary[2,1:4]=as.vector(summary(Test_pTau_M)$coefficient[2,1:4])
Merge_pTau_summary[3,1:4]=as.vector(summary(Test_pTau_F)$coefficient[2,1:4])
Merge_pTau_summary[4,1:4]=as.vector(summary(Test_pTau_M_APOE4)$coefficient[2,1:4])
Merge_pTau_summary[5,1:4]=as.vector(summary(Test_pTau_M_nonAPOE4)$coefficient[2,1:4])
Merge_pTau_summary[6,1:4]=as.vector(summary(Test_pTau_F_APOE4)$coefficient[2,1:4])
Merge_pTau_summary[7,1:4]=as.vector(summary(Test_pTau_F_nonAPOE4)$coefficient[2,1:4])

Merge_pTau_summary$FDR=p.adjust(Merge_pTau_summary$P_value,method = "fdr")

# sST2 ~ Plasma NfL levels
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

#####----------------------------------------------
#### CSF sST2 ~ AD / AD-related endophenotypes
#####----------------------------------------------

# Data input
Data_input=read.table('./Input_file_CSF_sST2_AD_Abeta.txt',header=T)
Data_input$DIAG<-relevel(Data_input$DIAG,"CONTROL")

# Data stratification by AD diagnosis / Sex
Data_AD<-subset(Data_input,Data_input$DIAG=="AD")
Data_NC<-subset(Data_input,Data_input$DIAG=="CONTROL")
Data_male<-subset(Data_input,Data_input$SEX=="1")
Data_female<-subset(Data_input,Data_input$SEX=="2")

Data_female_AD<-subset(Data_AD,Data_AD$SEX=="2")
Data_female_NC<-subset(Data_NC,Data_NC$SEX=="2")
Data_male_AD<-subset(Data_AD,Data_AD$SEX=="1")
Data_male_NC<-subset(Data_NC,Data_NC$SEX=="1")

#####----------------------------------------------
## Statistics on sST2: male vs female
Test_in_NC<-lm(sST2_ng.mL~SEX+AGE+PMD,data = Data_NC)
Test_in_AD<-lm(sST2_ng.mL~SEX+AGE+PMD,data = Data_AD)
Test_in_all<-lm(sST2_ng.mL~SEX+AGE+PMD,data = Data_input)

Merge_sex_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_sex_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_sex_summary)=c("Overall","NC","AD")
Merge_sex_summary[1,1:4]=as.vector(summary(Test_in_all)$coefficient[2,1:4])
Merge_sex_summary[2,1:4]=as.vector(summary(Test_in_NC)$coefficient[2,1:4])
Merge_sex_summary[3,1:4]=as.vector(summary(Test_in_AD)$coefficient[2,1:4])

Merge_sex_summary$FDR=p.adjust(Merge_sex_summary$P_value,method = "fdr")

#####----------------------------------------------
## Statistics in overall / male / female: sST2 ~ AD
Test_AD_all<-lm(sST2_ng.mL~DIAG+AGE+SEX+PMD,data = Data_input)
Test_AD_M<-lm(sST2_ng.mL~DIAG+AGE+SEX+PMD,data = Data_male)
Test_AD_F<-lm(sST2_ng.mL~DIAG+AGE+SEX+PMD,data = Data_female)

Merge_AD_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_AD_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_AD_summary)=c("Overall","Male","Female")

Merge_AD_summary[1,1:4]=as.vector(summary(Test_AD_all)$coefficient[2,1:4])
Merge_AD_summary[2,1:4]=as.vector(summary(Test_AD_M)$coefficient[2,1:4])
Merge_AD_summary[3,1:4]=as.vector(summary(Test_AD_F)$coefficient[2,1:4])

Merge_AD_summary$FDR=p.adjust(Merge_AD_summary$P_value,method = "fdr")

#---------------------------------------------------
#### Cutoff point of CSF sST2 to define Low and High CSF sST2 group using Youden's method
library(OptimalCutpoints)
test<-optimal.cutpoints(X="sST2_ng.mL",status = "DIAG", tag.healthy = "CONTROL", methods="Youden",data=Data_input, direction="<")
summary(test)

#---------------------------------------------------
####  CSF sST2 ~ Abeta deposition in AD

# In overall: Abeta %area / number / median size ~ sST2
sST2_AD_Abeta_test <-lm(Abeta_area~sST2_ng.mL+AGE+SEX+PMD,data=Data_AD)
sST2_AD_Abeta_number_test <-lm(Abeta_number_OLD~sST2_ng.mL+AGE+SEX+PMD,data=Data_AD)
sST2_AD_Abeta_size_test <-lm(Abeta_median_size~sST2_ng.mL+AGE+SEX+PMD,data=Data_AD)

# In male: Abeta %area / number / median size ~ sST2
sST2_Male_AD_Abeta_test <-lm(Abeta_area~sST2_ng.mL+AGE+SEX+PMD,data=Data_male_AD)
sST2_Male_AD_Abeta_number_test <-lm(Abeta_number_OLD~sST2_ng.mL+AGE+SEX+PMD,data=Data_male_AD)
sST2_Male_AD_Abeta_size_test <-lm(Abeta_median_size~sST2_ng.mL+AGE+SEX+PMD,data=Data_male_AD)

# In female: Abeta %area / number / median size ~ sST2
sST2_Female_AD_Abeta_test <-lm(Abeta_area~sST2_ng.mL+AGE+SEX+PMD,data=Data_female_AD)
sST2_Female_AD_Abeta_number_test <-lm(Abeta_number_OLD~sST2_ng.mL+AGE+SEX+PMD,data=Data_female_AD)
sST2_Female_AD_Abeta_size_test <-lm(Abeta_median_size~sST2_ng.mL+AGE+SEX+PMD,data=Data_female_AD)

#####----------------------------------------------
## Merge summary and FDR correction
# Percentage area of abeta
Merge_AD_Abeta_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_AD_Abeta_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_AD_Abeta_summary)=c("Overall","Male","Female")

Merge_AD_Abeta_summary[1,1:4]=as.vector(summary(sST2_AD_Abeta_test)$coefficient[2,1:4])
Merge_AD_Abeta_summary[2,1:4]=as.vector(summary(sST2_Male_AD_Abeta_test)$coefficient[2,1:4])
Merge_AD_Abeta_summary[3,1:4]=as.vector(summary(sST2_Female_AD_Abeta_test)$coefficient[2,1:4])

Merge_AD_Abeta_summary$FDR=p.adjust(Merge_AD_Abeta_summary$P_value,method = "fdr")

# Number of abeta
Merge_AD_Abeta_number_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_AD_Abeta_number_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_AD_Abeta_number_summary)=c("Overall","Male","Female")

Merge_AD_Abeta_number_summary[1,1:4]=as.vector(summary(sST2_AD_Abeta_number_test)$coefficient[2,1:4])
Merge_AD_Abeta_number_summary[2,1:4]=as.vector(summary(sST2_Male_AD_Abeta_number_test)$coefficient[2,1:4])
Merge_AD_Abeta_number_summary[3,1:4]=as.vector(summary(sST2_Female_AD_Abeta_number_test)$coefficient[2,1:4])

Merge_AD_Abeta_number_summary$FDR=p.adjust(Merge_AD_Abeta_number_summary$P_value,method = "fdr")

# Median size of abeta
Merge_AD_Abeta_size_summary=data.frame(matrix(,ncol = 4,nrow = 3))
colnames(Merge_AD_Abeta_size_summary)=c("Estimate","SE","T","P_value")
rownames(Merge_AD_Abeta_size_summary)=c("Overall","Male","Female")

Merge_AD_Abeta_size_summary[1,1:4]=as.vector(summary(sST2_AD_Abeta_size_test)$coefficient[2,1:4])
Merge_AD_Abeta_size_summary[2,1:4]=as.vector(summary(sST2_Male_AD_Abeta_size_test)$coefficient[2,1:4])
Merge_AD_Abeta_size_summary[3,1:4]=as.vector(summary(sST2_Female_AD_Abeta_size_test)$coefficient[2,1:4])

Merge_AD_Abeta_size_summary$FDR=p.adjust(Merge_AD_Abeta_size_summary$P_value,method = "fdr")


