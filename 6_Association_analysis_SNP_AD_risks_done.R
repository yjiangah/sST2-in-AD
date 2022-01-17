#####
##### Code for association analysis between candidate SNPs and AD risks
#####-------------------------------------------------------------------------------------------------

library(questionr)

# Data input
Data_input=read.table('./Input_file_AD_Phenotype_Genotype.txt',header=T)

# Data stratification by AD diagnosis / Sex / Apoe4 genotypes
Data_AD<-subset(Data_input,Data_input$Phenotype=="AD")
Data_NC<-subset(Data_input,Data_input$Phenotype=="NC")
Data_APOE4<-subset(Data_input,Data_input$APOE4!=0)
Data_nonAPOE4<-subset(Data_input,Data_input$APOE4==0)
Data_male<-subset(Data_input,Data_input$Gender=="1")
Data_female<-subset(Data_input,Data_input$Gender=="2")

Data_male_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="1")
Data_female_APOE4<-subset(Data_APOE4,Data_APOE4$Gender=="2")
Data_male_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="1")
Data_female_nonAPOE4<-subset(Data_nonAPOE4,Data_nonAPOE4$Gender=="2")

#####----------------------------------------------
## In overall: SNP ~ AD
All_AD_SNP_test <- glm(Data_input$Phenotype==2~SNP+Age+Gender+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_input,family=binomial("logit"))
OR_All<-odds.ratio(All_AD_SNP_test)

## In female and male: SNP ~ AD
Female_AD_SNP_test <- glm(Data_female$Phenotype==2~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_female,family=binomial("logit"))
Male_AD_SNP_test <- glm(Data_male$Phenotype==2~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_male,family=binomial("logit"))

OR_F<-odds.ratio(Female_AD_SNP_test)
OR_M<-odds.ratio(Male_AD_SNP_test)

## In female and male APOE4 carriers: SNP ~ AD
Female_APOE4_AD_SNP_test <- glm(Data_female_APOE4$Phenotype==2~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_female_APOE4,family=binomial("logit"))
Male_APOE4_AD_SNP_test <- glm(Data_male_APOE4$Phenotype==2~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_male_APOE4,family=binomial("logit"))

OR_F_APOE4<-odds.ratio(Female_APOE4_AD_SNP_test)
OR_M_APOE4<-odds.ratio(Male_APOE4_AD_SNP_test)

## In female and male APOE4 carriers: SNP ~ AD
Female_nonAPOE4_AD_SNP_test <- glm(Data_female_nonAPOE4$Phenotype==2~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_female_nonAPOE4,family=binomial("logit"))
Male_nonAPOE4_AD_SNP_test <- glm(Data_male_nonAPOE4$Phenotype==2~SNP+Age+PC_1+PC_2+PC_3+PC_4+PC_5,data=Data_male_nonAPOE4,family=binomial("logit"))

OR_F_nonAPOE4<-odds.ratio(Female_nonAPOE4_AD_SNP_test)
OR_M_nonAPOE4<-odds.ratio(Male_nonAPOE4_AD_SNP_test)

#####----------------------------------------------
## Merge summary and FDR correction
Merge_AD_summary=data.frame(matrix(,ncol = 5,nrow = 7))
colnames(Merge_AD_summary)=c("Estimate","SE","T","P_value","OR")
rownames(Merge_AD_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")

Merge_AD_summary[1,1:4]=as.vector(summary(All_AD_SNP_test)$coefficient[2,1:4])
Merge_AD_summary[2,1:4]=as.vector(summary(Male_AD_SNP_test)$coefficient[2,1:4])
Merge_AD_summary[3,1:4]=as.vector(summary(Female_AD_SNP_test)$coefficient[2,1:4])
Merge_AD_summary[4,1:4]=as.vector(summary(Male_APOE4_AD_SNP_test)$coefficient[2,1:4])
Merge_AD_summary[5,1:4]=as.vector(summary(Male_nonAPOE4_AD_SNP_test)$coefficient[2,1:4])
Merge_AD_summary[6,1:4]=as.vector(summary(Female_APOE4_AD_SNP_test)$coefficient[2,1:4])
Merge_AD_summary[7,1:4]=as.vector(summary(Female_nonAPOE4_AD_SNP_test)$coefficient[2,1:4])

Merge_AD_summary[1,5]=OR_All[2,1]
Merge_AD_summary[2,5]=OR_M[2,1]
Merge_AD_summary[3,5]=OR_F[2,1]
Merge_AD_summary[4,5]=OR_M_APOE4[2,1]
Merge_AD_summary[5,5]=OR_M_nonAPOE4[2,1]
Merge_AD_summary[6,5]=OR_F_APOE4[2,1]
Merge_AD_summary[7,5]=OR_F_nonAPOE4[2,1]

Merge_AD_summary$FDR=p.adjust(Merge_AD_summary$P_value,method = "fdr")








