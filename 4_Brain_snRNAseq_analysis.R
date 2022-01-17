#####
##### Code for brain snRNA-seq analysis
#####-------------------------------------------------------------------------------------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SoupX)

# -------------------------------------------------------------------------------------
# Input individual files
# -------------------------------------------------------------------------------------
ID1 <- load10X("./ID1_snRNAseq/outs/")
ID1 = autoEstCont(ID1)
ID1 = adjustCounts(ID1)
ID1 <- CreateSeuratObject(counts = ID1, project = "SNP_noncarrier", min.cells = 5, min.features = 200)

ID2 <- load10X("./ID2_snRNAseq/outs/")
ID2 = autoEstCont(ID2)
ID2 = adjustCounts(ID2)
ID2 <- CreateSeuratObject(counts = ID2, project = "SNP_carrier", min.cells = 5, min.features = 200)

ID3 <- load10X("./ID3_snRNAseq/outs/")
ID3 = autoEstCont(ID3)
ID3 = adjustCounts(ID3)
ID3 <- CreateSeuratObject(counts = ID3, project = "SNP_noncarrier", min.cells = 5, min.features = 200)

# -------------------------------------------------------------------------------------
# Merge datasets
# -------------------------------------------------------------------------------------
data <- merge(ID1, y= c(ID2,ID3), 
              add.cell.ids = c("ID1", "ID2", "ID3"), 
              project = "Merged_snRNAseq")

head(colnames(data))
table(data$orig.ident)
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data <- subset(data, subset = nFeature_RNA > 200 & percent.mt < 10)
table(data$orig.ident)

# -------------------------------------------------------------------------------------
# Normalizing data
# -------------------------------------------------------------------------------------
data <- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
data <- RunPCA(data, features = VariableFeatures(object = data))
saveRDS(data, file = "./RDS_brain_snRNAseq_beforecluster.rds")

# -------------------------------------------------------------------------------------
# Determine cell type
# -------------------------------------------------------------------------------------
data <- FindNeighbors(data, dims = 1:20)
data <- FindClusters(data, resolution = 0.1)
data <- RunUMAP(data, dims = 1:20)

DefaultAssay(data) <- "RNA"

png(filename="UMAP_all_cells.png", width = 900, height = 800)
DimPlot(data, label = T, label.size = 8, pt.size = 1.5)
dev.off()

png(filename="UMAP_Endothelial_marker_gene.png", width = 1000, height = 1000)
FeaturePlot(data, features = c("PRKCH", "CLDN5", "PDGFRB", "SLC2A1"), label = T, order = T, pt.size = 0.5)
dev.off()

png(filename="UMAP_Oligodendrocyte_marker_gene.png", width = 1000, height = 1000)
FeaturePlot(data, features = c("CLDN11", "OLIG2", "PHLPP1", "TF"), label = T, order = T, pt.size = 0.5)
dev.off()

png(filename="UMAP_OPC_marker_gene.png", width = 1000, height = 1000)
FeaturePlot(data, features = c("VCAN", "CSPG4", "XYLT1", "PTPRZ1"), label = T, order = T, pt.size = 0.5)
dev.off()

png(filename="UMAP_Astrocyte_marker_gene.png", width = 1000, height = 1000)
FeaturePlot(data, features = c("ATP1A2", "PREX2", "CTNNA2", "AQP4"), label = T, order = T, pt.size = 0.5)
dev.off()

png(filename="UMAP_Microglia_marker_gene.png", width = 1000, height = 1000)
FeaturePlot(data, features = c("NAV3", "MEF2A", "ST6GAL1", "PTPRC"), label = T, order = T, pt.size = 0.5)
dev.off()

png(filename="UMAP_Excitatory_Inhibitory_neuron_marker_gene.png", width = 1000, height = 1000)
FeaturePlot(data, features = c("SLC17A7", "SATB2", "GAD1", "GAD2", "GABRA1", "SLC6A1"), label = T, order = T, pt.size = 0.5)
dev.off()

# -------------------------------------------------------------------------------------
# Assign cell type
# -------------------------------------------------------------------------------------
data <- RenameIdents(data,
                        `1`="Excit",`6`="Excit",`7`="Excit",`9`="Excit",`11`="Excit",`13`="Excit",
                        `3`="Inhib",`5`="Inhib",`10`="Inhib",
                        `2`="Astro",
                        `8`="Micro",`12`="Endo",
                        `0`="Oligo",`4`="OPC")

saveRDS(data, file = "./RDS_brain_snRNAseq_aftercluster_celltypes.rds")

# UMAP plot of cell types
png(filename="UMAP_brain_snRNAseq_cell_type.png", width = 900, height = 800)
DimPlot(data,reduction="umap", label = T, label.size = 6, pt.size = 1.5)
dev.off()

# UMAP plot of candidate genes
png(filename="UMAP_brain_snRNAseq_sST2.png", width = 900, height = 800)
FeaturePlot(data, features = c("SOLUBLEST2"), order = T, pt.size = 1.5, cols=c("grey","red"),split.by="orig.ident")
dev.off()

# -------------------------------------------------------------------------------------
# Examine SNP effect on sST2 expression in brain endothelial cells
# -------------------------------------------------------------------------------------
Endo_data<-subset(data,idents="Endo")

my_genes = c("CLDN5", "SOLUBLEST2")
Endo_data@assays$RNA[my_genes]
write.csv(Endo_data@assays$RNA[my_genes], "./Brain_snRNAseq_Endothelial_cells_sST2_expression_singlecells.csv")

# Statistical analysis on the effects of SNP on sST2 expression
Data = read.table('./Brain_snRNAseq_ECs_sST2_CLDN5_expression_proportion_per_ID.txt',header = T)

test_sST2_expression<-lm(Average_sST2_expression_in_POS_cells~SNP+AD_diagnosis+AGE+SEX+PMD, data= Data)
test_sST2_proportion<-lm(sST2EC_proportion~SNP+AD_diagnosis+AGE+SEX+PMD, data= Data)
test_CLDN5_expression<-lm(Average_CLDN5_expression_in_POS_cells~SNP+AD_diagnosis+AGE+SEX+PMD, data= Data)
test_CLDN5_proportion<-lm(CLDN5EC_proportion~SNP+AD_diagnosis+AGE+SEX+PMD, data= Data)

# Bubble plot for the effects of SNP on sST2 expression
library(ggplot2)

Data = read.table('./Brain_snRNAseq_ECs_sST2_CLDN5_for_bubble_plot.txt',header = T)

png(filename="Brain_snRNAseq_ECs_sST2_CLDN5_SNP_bubble_plot.png", width = 900, height = 800)
test<-ggplot(Data,aes(y=factor(Genes),x=factor(Genotype_order)))+
  geom_point(aes(colour = Expression,size = Proportion))+
  scale_color_gradient(low='yellow',high='red')+
  scale_size(range=c(1, 25))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.ticks.x=element_blank())
test
dev.off()

# -------------------------------------------------------------------------------------
# Examine SNP effect on microglial transcriptome
# -------------------------------------------------------------------------------------
Micro_data<-subset(my_data,idents="Micro")
# Extract all microglial genes (expressed by > 1% microglia)
All_microglia_genes<-as.matrix(GetAssayData(object = Micro_data,slot="counts"))
All_microglia_genes_pos=All_microglia_genes
All_microglia_genes_pos[All_microglia_genes_pos==0]<-NA

Number=data.frame(rowSums(!is.na(All_microglia_genes_pos)))
colnames(Number)<-c("Expressing_cells")
Merge=cbind(Number,All_microglia_genes_pos)
Merge_filtered_microglia_onepercent_cell_genes<-subset(Merge,Merge$Expressing_cells>length(All_microglia_genes[1,])*0.01)

write.table(t(Merge_filtered_microglia_onepercent_cell_genes),"./Brain_snRNAseq_Microglia_filtered_gene_expression_singlecells.txt",row.names = T)

# Statistical analysis on the effects of SNP on microglial transcriptome
Data = read.table('./Brain_snRNAseq_Microglial_transcriptome_expression_proportion_per_ID.txt',header = T)

Test_microglial_transcriptome<-data.frame(matrix(,ncol=4,nrow=(length(Data[,1])-6)))
colnames(Test_microglial_transcriptome)<-c("Estimate","SE","T_value","P_value")
rownames(Test_microglial_transcriptome)<-names(Data[c(1:(length(Data[1,])-6))])

for (i in 1:(length(Data[1,])-6))
{
  All_test<-lm(Data[,i]~SNP+AD_diagnosis+AGE+SEX+PMD,data=Data)
  Test_microglial_transcriptome[i,1:4]=as.vector(t(summary(All_test)$coefficients[2,1:4]))
}

# Volcano plot for the effects of SNP on microglial transcriptome
library(MASS)
library(calibrate)

Data=read.table('./Brain_snRNAseq_SNP_effects_on_Microglia_transcriptome_for_volcano_plot.txt',header=T)

png(filename="Brain_snRNAseq_Microglial_transcriptome_SNP_volcano_plot.png", width = 900, height = 800)

with(Data,plot(Estimate,log_FDR,pch=20,main="SNP effects on microglial genes",ylim=c(0,40),xlim=c(-1.5,1.5)))
with(subset(Data,Estimate>0),points(Estimate,log_FDR,pch=21,col="black",bg="tomato2",cex=log_FDR/10+0.5))
with(subset(Data,Estimate<0),points(Estimate,log_FDR,pch=21,col="black",bg="steelblue2",cex=log_FDR/10+0.5))
with(subset(Data,log_FDR>10),textxy(Estimate,log_FDR,labs=Genes,cex=0.4))
dev.off()

# Bubble plot for the effects of SNP on candidate microglial activation genes
library(ggplot2)

Data = read.table('./Brain_snRNAseq_Microglial_activation_genes_for_bubble_plot.txt',header = T)

png(filename="Brain_snRNAseq_Microglial_activation_genes_SNP_bubble_plot.png", width = 900, height = 800)
test<-ggplot(Data,aes(y=factor(Genes),x=factor(Genotype_order)))+
  geom_point(aes(colour = Expression,size = Proportion))+
  scale_color_gradient(low='yellow',high='red')+
  scale_size(range=c(1, 25))+
  theme_bw()+ theme(panel.border = element_rect(colour = "black", fill=NA, size=1), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),axis.ticks.x=element_blank())
test
dev.off()


