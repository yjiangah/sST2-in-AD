#####
##### Code for data visualization of sST2 GWAS and proporition variance explained
#####-------------------------------------------------------------------------------------------------

#####---------------------------------------------------
## Plots for sST2 GWAS
#####---------------------------------------------------
## Manhattan plot for sST2 GWAS
library(qqman)

gwasResults=read.table('./Input_file_sST2_GWAS_results.txt',header=T)
gwasResults=na.omit(gwasResults)

png(filename="sST2_GWAS_Manhattan_Plot.png", width = 900, height = 800)
manhattan(gwasResults)
dev.off()

#####-------------------
## Q-Q plot for sST2 GWAS
library(ggplot2)
library(qqman)

png(filename="sST2_GWAS_QQ_Plot.png", width = 900, height = 800)
qq(gwasResults$P, main = "Q-Q plot of sST2 p-values", pch=18, cex=1.5, las=1)
dev.off()

#lambda
p_value=gwasResults$P
z=qnorm(p_value/2)
lambda=round(median(z^2, na.rm = TRUE) / 0.454, 3)
lambda

#####-------------------
## Bubble plot at IL1RL1 locus for fine mapping analysis
library(plotly)

Fine_mapped_IL1RL1_locus=read.table('./sST2_GWAS_fine_mapped_IL1RL1_locus.txt',header = T)

png(filename="sST2_GWAS_Fined_mapped_IL1RL1_locus_bubble_Plot.png", width = 900, height = 800)

test<-plot_ly(Fine_mapped_IL1RL1_locus,x = ~BP, y= ~BETA, text= ~SNP, type='scatter',mode='markers',color= ~-BETA, colors='RdBu',sizes=c(50,1000),size=~Fine_mapped_P,marker=list(opacity=0.5)) %>%
  layout(title='IL1RL1 locus fine mapping',
         xaxis= list(title='Position on chr2 (Mb)',showline=T, linewidth=3,zeroline=FALSE,ticklen=5,tickwidth=3,showgrid=FALSE,tickfont=list(size=20),titlefont=list(size=30)),
         yaxis= list(title='Estimate (Effect size)',showline=T, linewidth=3,zerolinewidth=3,ticklen=5,tickwidth=3,showgrid=FALSE,tickfont=list(size=20),titlefont=list(size=30),side="left"),
         paper_bgcolor='rgb(256,256,256)',
         plot_bgcolor='rgb(256,256,256)')
test
dev.off()

#####---------------------------------------------------
## Proportion of variance of sST2 levels explained by genetic and nongenetic factors
#####---------------------------------------------------
library(relaimpo)

Data_input=read.table('./Input_file_sST2_Phenotype_Genotype.txt',header=T)


linear_test<-lm(sST2_level~Age+Gender+SNP, Data_input)
Proportion_variance_explained<-calc.relimp(linear_test,type = c("lmg"),rela=F)

png(filename="Proportion_variance_explained_by_candidate_factors.png", width = 900, height = 800)
plot(Proportion_variance_explained)
dev.off()



