#####
##### Code for two-sample Mendelian randomization analysis of sST2 in AD
#####-------------------------------------------------------------------------------------------------
library(TwoSampleMR)

setwd("./Mendelian_Randomization")
File_groups=c("Overall","Male","Female","Male_APOE4","Male_nonAPOE4","Female_APOE4","Female_nonAPOE4")

MR_test_summary=data.frame(matrix(,ncol = 5,nrow = 7))
colnames(MR_test_summary)=c("Method","Number_of_SNP","Estimate","SE","P_value")
rownames(MR_test_summary)=c("Overall","Male","Female","M_APOE4","M_nonAPOE4","F_APOE4","F_nonAPOE4")
k=1

for (i in File_groups){
  Data_input=read.table(paste("MR_input_file_",toString(i),".txt",sep = ""),header = T)
  MR_test<-mr(Data_input)
  MR_test_summary[k,1:5]=as.vector(MR_test[3,5:9])
  k=k+1
}
MR_test_summary$FDR=p.adjust(MR_test_summary$P_value,method = "fdr")





