

library(fishpond)
library(tximport)
library(devtools)
library(ggploT2)
library(Seurat)

library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(corrplot)
library(MatrixGenerics)
library(dplyr)
library(FateID)


library(slingshot)

library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(Matrix)


library(scater)

library(cowplot)


library("ComplexHeatmap")
library("zoo")
library(circlize)
library("Seurat")


### demultiplex and process all the gene counts output by Alevin

filesD1_2 <- file.path("lacbD1_2_alevin_out_iso_trinity/alevin/quants_mat.gz")
file.exists(filesD1_2)
txi_la_D1_2<-tximport(filesD1_2, type = "alevin")
library(Seurat)
library(Matrix)
mat_D1_2<-Matrix(txi_la_D1_2$counts, sparse = TRUE)


D1_2_psums<-colSums(mat_D1_2)



# knee plot for D1_2
sortedD1_2 <- sort(D1_2_psums,decreasing = TRUE)
plot(sortedD1_2, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot D1_2",
     col = ifelse(sortedD1_2 > 1000,'red','black'),
     xlim = c(1,length(sortedD1_2)), 
     ylim = c(1,sortedD1_2[[1]]),
     log = "xy"
)

keepD1_2 <- colSums(mat_D1_2) > 1000
countsD1_2_filter <- mat_D1_2[,keepD1_2]

sumD1_2 <- colSums(countsD1_2_filter)
range(sumD1_2)
length(sumD1_2)



UMIspercellD1_2_sorted <- sort(sumD1_2,decreasing = TRUE)

pD1_2.a <- plot(UMIspercellD1_2_sorted, 
                xlab = "Cells",
                ylab = "UMI counts",
                main = "Frequency distribution of UMIs per cell D1_2",
                xlim = c(1,length(UMIspercellD1_2_sorted)), 
                ylim = c(1,UMIspercellD1_2_sorted[[1]]),
)

pD1_2.b <- hist(log10(UMIspercellD1_2_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellD1_2 <- apply(countsD1_2_filter > 0, 2, sum)
hist(log10(genespercellD1_2), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes D1_2")


plot(UMIspercellD1_2_sorted, genespercellD1_2, log='xy')
title('counts vs genes per cell D1_2')

write.table(countsD1_2_filter, file = "D1_2_iso.mtx", append = FALSE, quote = TRUE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, 
            qmethod = c("escape", "double"), fileEncoding = "")

libD2_2 <- CreateSeuratObject(counts = countsD2_2_filter, min.cells = 3)

CellsMetaD1_2 = libD1_2_iso@meta.data
Library_source <- rep("libD1_2",dim(libD1_2_iso@meta.data)[1])
CellsMetaD1_2 <- cbind(CellsMetaD1_2, Library_source)
libD1_2_iso <- AddMetaData(libD1_2_iso, CellsMetaD1_2)




#make D1_1

files_la_D1_1 <- file.path("lacbD1_1_alevin_out_iso_trinity/alevin/quants_mat.gz")
file.exists(files_la_D1_1)
txi_la_D1_1<-tximport(files_la_D1_1, type = "alevin")


mat_D1_1<-Matrix(txi_la_D1_1$counts, sparse = TRUE)
D1_1_psums<-colSums(mat_D1_1)

sortedD1_1 <- sort(D1_1_psums,decreasing = TRUE)
plot(sortedD1_1, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot D1_1",
     col = ifelse(sortedD1_1 > 1000,'red','black'),
     xlim = c(1,length(sortedD1_1)), 
     ylim = c(1,sortedD1_1[[1]]),
     log = "xy"
)

keepD1_1 <- colSums(mat_D1_1) > 1000
countsD1_1_filter <- mat_D1_1[,keepD1_1]

sumD1_1 <- colSums(countsD1_1_filter)
range(sumD1_1)
length(sumD1_1)
UMIspercellD1_1_sorted <- sort(sumD1_1,decreasing = TRUE)

pD1_1.1 <- plot(UMIspercellD1_1_sorted, 
                xlab = "Cells",
                ylab = "UMI counts",
                main = "Frequency distribution of UMIs per cell D1_1",
                xlim = c(1,length(UMIspercellD1_1_sorted)), 
                ylim = c(1,UMIspercellD1_1_sorted[[1]]),
)

pD1_1.b <- hist(log10(UMIspercellD1_1_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellD1_1 <- apply(countsD1_1_filter > 0, 2, sum)
hist(log10(genespercellD1_1), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes D1_1")


plot(UMIspercellD1_1_sorted, genespercellD1_1, log='xy')
title('counts vs genes per cell D1_1')
#important filter step
libD2_1_iso <- CreateSeuratObject(counts = countsD2_1_filter, min.cells = 3)

CellsMetaD1_1 = libD1_1_iso@meta.data
Library_source <- rep("libD1_1",dim(libD1_1_iso@meta.data)[1])
CellsMetaD1_1 <- cbind(CellsMetaD1_1, Library_source)
libD1_1_iso <- AddMetaData(libD1_1_iso, CellsMetaD1_1)

D1_x_iso_AI <- merge(libD1_2_iso, y = libD1_1_iso, add.cell.ids = c("D1_2", "D1_1"), project = "Dx_Pleuro")





mito.genes <- read.table(file = "mitochondrial_genes_updated.txt")
mito.genes <- as.character(mito.genes[,1])
mito.genesD1_D1_iso <- mito.genes[mito.genes %in% rownames(D1_D1_iso)]
D1_D1_iso_AI[["percent.mt"]] <- PercentageFeatureSet(D1_D1_iso_AI, features = mito.genesD1_D1_iso)

D1_D1_iso_AI<-subset(x=D1_D1_iso_AI, subset = percent.mt< 15)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
D1_D1_iso_AI <- CellCycleScoring(D1_D1_iso_AI, s.features = s.genes, g2m.features = g2m.genes)


library(tximport)
filesD2_2 <- file.path("lacbD2_2_alevin_out_iso_trinity/alevin/quants_mat.gz")
file.exists(filesD2_2)
txi_la_D2_2<-tximport(filesD2_2, type = "alevin")
library(Seurat)
library(Matrix)
mat_D2_2<-Matrix(txi_la_D2_2$counts, sparse = TRUE)


D2_2_psums<-colSums(mat_D2_2)



# knee plot for D2_2
sortedD2_2 <- sort(D2_2_psums,decreasing = TRUE)
plot(sortedD2_2, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot D2_2",
     col = ifelse(sortedD2_2 > 1000,'red','black'),
     xlim = c(1,length(sortedD2_2)), 
     ylim = c(1,sortedD2_2[[1]]),
     log = "xy"
)

keepD2_2 <- colSums(mat_D2_2) > 1000
countsD2_2_filter <- mat_D2_2[,keepD2_2]

sumD2_2 <- colSums(countsD2_2_filter)
range(sumD2_2)
length(sumD2_2)



UMIspercellD2_2_sorted <- sort(sumD2_2,decreasing = TRUE)

pD2_2.a <- plot(UMIspercellD2_2_sorted, 
                xlab = "Cells",
                ylab = "UMI counts",
                main = "Frequency distribution of UMIs per cell D2_2",
                xlim = c(1,length(UMIspercellD2_2_sorted)), 
                ylim = c(1,UMIspercellD2_2_sorted[[1]]),
)

pD2_2.b <- hist(log10(UMIspercellD2_2_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellD2_2 <- apply(countsD2_2_filter > 0, 2, sum)
hist(log10(genespercellD2_2), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes D2_2")


plot(UMIspercellD2_2_sorted, genespercellD2_2, log='xy')
title('counts vs genes per cell D2_2')

write.table(countsD2_2_filter, file = "D2_2_iso.mtx", append = FALSE, quote = TRUE, sep = " ", 
            eol = "\n", na = "NA", dec = ".", row.names = TRUE, col.names = TRUE, 
            qmethod = c("escape", "double"), fileEncoding = "")

libD2_2 <- CreateSeuratObject(counts = countsD2_2_filter, min.cells = 3)

CellsMetaD2_2 = libD2_2_iso@meta.data
Library_source <- rep("libD2_2",dim(libD2_2_iso@meta.data)[1])
CellsMetaD2_2 <- cbind(CellsMetaD2_2, Library_source)
libD2_2_iso <- AddMetaData(libD2_2_iso, CellsMetaD2_2)




#make D2_1

files_la_D2_1 <- file.path("lacbD2_1_alevin_out_iso_trinity/alevin/quants_mat.gz")
file.exists(files_la_D2_1)
txi_la_D2_1<-tximport(files_la_D2_1, type = "alevin")


mat_D2_1<-Matrix(txi_la_D2_1$counts, sparse = TRUE)
D2_1_psums<-colSums(mat_D2_1)

sortedD2_1 <- sort(D2_1_psums,decreasing = TRUE)
plot(sortedD2_1, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot D2_1",
     col = ifelse(sortedD2_1 > 1000,'red','black'),
     xlim = c(1,length(sortedD2_1)), 
     ylim = c(1,sortedD2_1[[1]]),
     log = "xy"
)

keepD2_1 <- colSums(mat_D2_1) > 1000
countsD2_1_filter <- mat_D2_1[,keepD2_1]

sumD2_1 <- colSums(countsD2_1_filter)
range(sumD2_1)
length(sumD2_1)
UMIspercellD2_1_sorted <- sort(sumD2_1,decreasing = TRUE)

pD2_1.1 <- plot(UMIspercellD2_1_sorted, 
                xlab = "Cells",
                ylab = "UMI counts",
                main = "Frequency distribution of UMIs per cell D2_1",
                xlim = c(1,length(UMIspercellD2_1_sorted)), 
                ylim = c(1,UMIspercellD2_1_sorted[[1]]),
)

pD2_1.b <- hist(log10(UMIspercellD2_1_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellD2_1 <- apply(countsD2_1_filter > 0, 2, sum)
hist(log10(genespercellD2_1), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes D2_1")


plot(UMIspercellD2_1_sorted, genespercellD2_1, log='xy')
title('counts vs genes per cell D2_1')


#important filter step
libD2_1_iso <- CreateSeuratObject(counts = countsD2_1_filter, min.cells = 3)

CellsMetaD2_1 = libD2_1_iso@meta.data
Library_source <- rep("libD2_1",dim(libD2_1_iso@meta.data)[1])
CellsMetaD2_1 <- cbind(CellsMetaD2_1, Library_source)
libD2_1_iso <- AddMetaData(libD2_1_iso, CellsMetaD2_1)



D2_x_iso_AI <- merge(libD2_2_iso, y = libD2_1_iso, add.cell.ids = c("D2_2", "D2_1"), project = "Dx_Pleuro")


CellsMetaD2_x = D2_x_iso_AI@meta.data
Stage <- rep("D2_x",dim(D2_x_iso_AI@meta.data)[1])
CellsMetaD2_x <- cbind(CellsMetaD2_x, Stage)
D2_x_iso_AI <- AddMetaData(D2_x_iso_AI, CellsMetaD2_x)


CellsMetaD1_x = D1_x_iso_AI@meta.data
Stage <- rep("D1_x",dim(D1_x_iso_AI@meta.data)[1])
CellsMetaD1_x <- cbind(CellsMetaD1_x, Stage)
D1_x_iso_AI <- AddMetaData(D1_x_iso_AI, CellsMetaD1_x)


D1_D2_iso_AI <- merge(D1_x_iso_AI, y =D2_x_iso_AI, add.cell.ids = c("D1_x", "D2_x"), project = "D1_D2_Pleuro")

#prepare to run classifieer on D1_D2
pleuro.2<-D1_D2_iso_AI
# Create a data frame containing metadata
data.info <- pleuro.2@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")]
cell.names <- colnames(pleuro.2)

# Here you add all the idents from bad clusters
bad.names <- colnames(subset(x = pleuro.2, idents = c(0, 28, 7, 10)))
good.names <- cell.names[!cell.names %in% bad.names] # The good clusters must be the rest of the cells

train.bad.names<-sample(bad.names,110) # selects randomly X bad cells for training set. You will want to train your set with ~10% of the total number of "bad cells" that you have (look at length(bad.names))
train.good.names<-sample(good.names,110) # selects randomly X good cells for training set. Keep the same number of cells as above.

train.bad<-data.info[train.bad.names,]
train.good<-data.info[train.good.names,]
train<-rbind(train.good,train.bad)
ytrain <- matrix(c(rep(1,110),rep(-1,110))) # Limits are the same as your number of trainig cells

train.final = data.frame(train, ytrain=as.factor(ytrain))

library(e1071)
svmfit=svm(ytrain~., data=train.final, kernel="linear", cost=0.1, scale=F)

set.seed (1)
tune.out=tune(svm,ytrain~.,data=train.final,kernel="linear", ranges=list(cost=c(0.1,1,10,100,500)))
summary(tune.out)


bestmod=tune.out$best.model
summary(bestmod)

test.bad.names<-bad.names[!bad.names %in% train.bad.names] 
test.good.names<-good.names[!good.names %in% train.good.names] 
test.bad<-data.info[test.bad.names,] 
test.good<-data.info[test.good.names,] 
test<-rbind(test.good,test.bad) 

test<-data.info[!rownames(data.info) %in% rownames(train),] 

test<-as.matrix(test)
ytest<-matrix(c(rep(1,length(test.good.names)),rep(-1,length(test.bad.names)))) 

test.final = data.frame(test, ytest=as.factor(ytest))

ypred=predict(bestmod,test.final)
table(predict=ypred, truth=test.final$ytest)

ypred<-as.data.frame(ypred)
pleuro.2$ypred <- ypred

pred<-as.data.frame(pleuro.2$ypred)
rownames(pred)<-rownames(pleuro.2[[]])
pred[train.bad.names,]<-"-1"
pred[train.good.names,]<-"1"
colnames(pred)<-c("svm_class")
pleuro.2$svm_class <- pred

DimPlot(pleuro.2, split.by="svm_class", pt.size=0.7, label = T, group.by = "Stage") #+ NoLegend()


pleuro.3 <- subset(pleuro.2, cells=rownames(pleuro.2@meta.data[pleuro.2@meta.data$svm_class=="1",]))

pleuro.2 <- pleuro.3
rm(pleuro.3)

pleuro.2 <- subset(pleuro.2, cells=rownames(pleuro.2@meta.data[pleuro.2@meta.data$nFeature_RNA>800,]))


pleuro.2 <- SCTransform(pleuro.2, vars.to.regress = c("percent.mt","nCount_RNA","Stage"))

pleuro.2 <- RunPCA(pleuro.2, npcs=100)

ElbowPlot(pleuro.2, ndims=100)
ndims= 100

pleuro.2 <- RunUMAP(pleuro.2, dims=1:30)
pleuro.2 <- FindNeighbors(pleuro.2, dims=1:30)

pleuro.2 <- FindClusters(pleuro.2, resolution = 1)

DimPlot(pleuro.2, reduction="umap", label=T, pt.size = 0.3) + NoLegend()


D1_D2<-pleuro.2



library(tximport)

st36_T1.files <- file.path("ST36_T1_alevin_out_cat_MPX/alevin/quants_mat.gz")

st36_M1.files <- file.path("multiplex_st36_M1_cmo_out/alevin/quants_mat.gz")
file.exists(c(st36_T1.files, st36_M1.files))

st36_T1.txi <- tximport(files = st36_T1.files, type = "alevin")

st36_M1.txi <- tximport(files = st36_M1.files, type = "alevin")

# the below file can be found at https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
map <- read.table("3M-february-2018.txt.gz")
st36_T1ToADT <- map$V1
names(st36_T1ToADT) <- map$V2

st36_M1 <- tximport(files = st36_M1.files, type = "alevin")$counts

st36_T1<-st36_T1.txi$counts





colnames(st36_T1) <- st36_T1ToADT[colnames(st36_T1)]
length(intersect(colnames(st36_T1), colnames(st36_M1)))



common.cells<-(intersect(colnames(st36_T1),colnames(st36_M1)))

length(common.cells)




st_36_1_O <- CreateSeuratObject(st36_T1[, common.cells])


st36_M1<-st36_M1[1:4,]

st_36_1_O[["HTO"]] <- CreateAssayObject(counts = st36_M1[, common.cells])





VlnPlot(st_36_1_O, features = c("nCount_HTO"), log = TRUE )
st_36_1_O <- subset(st_36_1_O, subset = nCount_HTO < 2e4)
# st36_M1 Normalization
DefaultAssay(st_36_1_O) <- "HTO"
st_36_1_O <- NormalizeData(st_36_1_O,   normalization.method = "CLR", margin = 2, verbose = F)
VariableFeatures(st_36_1_O) <- rownames(st_36_1_O[["HTO"]]@counts)
st_36_1_O <- ScaleData(st_36_1_O, assay = "HTO", verbose = F)




st_36_1_O <- HTODemux(st_36_1_O, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st36_T1_class.pdf"), width=10, height=10)
Idents(st_36_1_O) <- "HTO_classification.global"
VlnPlot(st_36_1_O, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_36_1_O, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


### determine filter threshold
st36_RNA_psums<-colSums(st36_T1)

pdf(paste0("st36_norm2.pdf"), width=10, height=10)

# knee plot for E7_RNA
sortedSt_36_1_RNA <- sort(st36_RNA_psums, decreasing = TRUE)
plot(sortedSt_36_1_RNA, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot St_36_1_RNA",
     col = ifelse(sortedSt_36_1_RNA > 700,'red','black'),
     xlim = c(1,length(sortedSt_36_1_RNA)), 
     ylim = c(1,sortedSt_36_1_RNA[[1]]),
     log = "xy"
)

keepSt_36_1_RNA <- colSums(st36_T1) > 700
countsSt_36_1_RNA_filter <- st36_T1[,keepSt_36_1_RNA]

sumSt_36_1_RNA <- colSums(countsSt_36_1_RNA_filter)
range(sumSt_36_1_RNA)
length(sumSt_36_1_RNA)



UMIspercellSt_36_1_RNA_sorted <- sort(sumSt_36_1_RNA,decreasing = TRUE)

pSt_36_1_RNA.a <- plot(UMIspercellSt_36_1_RNA_sorted, 
                       xlab = "Cells",
                       ylab = "UMI counts",
                       main = "Frequency distribution of UMIs per cell St_36_1_RNA",
                       xlim = c(1,length(UMIspercellSt_36_1_RNA_sorted)), 
                       ylim = c(1,UMIspercellSt_36_1_RNA_sorted[[1]]),
)

pSt_36_1_RNA.b <- hist(log10(UMIspercellSt_36_1_RNA_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellSt_36_1_RNA <- apply(countsSt_36_1_RNA_filter > 0, 2, sum)
hist(log10(genespercellSt_36_1_RNA), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes St_36_1_RNA")


plot(UMIspercellSt_36_1_RNA_sorted, genespercellSt_36_1_RNA, log='xy')
title('counts vs genes per cell St_36_1_RNA')






dev.off()


#apply filter threshold to St36 and demultiplex data 



st_36_1_OC<-subset(x=st_36_1_O, subset = nCount_RNA > 700)


st_36_1_OC <- HTODemux(st_36_1_OC, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st36_T1_class_clean.pdf"), width=10, height=10)
Idents(st_36_1_OC) <- "HTO_classification.global"
VlnPlot(st_36_1_OC, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_36_1_OC, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()






saveRDS(st_36_1_OC, file = "St36_T1_mpx_clean.rds")
saveRDS(st_36_1_O, file="St_36_T1_mpx_full.rds")

doublets<-subset(x=st_36_1_OC, idents = "Doublet")

St_36_1_singlets<-subset(x=st_36_1_OC, idents = "Singlet")



mito.genes <- read.table(file = "mitochondrial_genes_updated.txt")
mito.genes <- as.character(mito.genes[,1])
mito.genesSt_36_1_singlets <- mito.genes[mito.genes %in% rownames(St_36_1_singlets)]
St_36_1_singlets[["percent.mt"]] <- PercentageFeatureSet(St_36_1_singlets, features = mito.genesSt_36_1_singlets)

St_36_1_singlets<-subset(x=St_36_1_singlets, subset= percent.mt<10)

St_36_1_singlets <- SCTransform(St_36_1_singlets, vars.to.regress = c("percent.mt", "nCount_RNA"))
St_36_1_singlets <- RunPCA(St_36_1_singlets, features = VariableFeatures(object = St_36_1_singlets))
St_36_1_singlets <- FindNeighbors(St_36_1_singlets, dims = 1:35)
St_36_1_singlets <- FindClusters(St_36_1_singlets, resolution = 1.5)
St_36_1_singlets <- RunUMAP(St_36_1_singlets, dims = 1:35)
DimPlot(St_36_1_singlets, reduction = "umap", label = T)

FeaturePlot(St_36_1_singlets, features=c("nFeature_RNA","nCount_RNA","percent.mt"))

FeaturePlot(St_36_1_singlets, features=c("RELN","TP73","NEUROD6","TBR1"))

saveRDS(St_36_1_singlets, file = "St_36_1_singlets_clean.rds")







###clean ST 36 singlets data using classifier

pleuro.4<-St_36_1_singlets

# Create a data frame containing metadata
data.info <- pleuro.4@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")]
cell.names <- colnames(pleuro.4)

# Here you add all the idents from bad clusters

VlnPlot(pleuro.4, features = c("nCount_RNA","nFeature_RNA", "percent.mt"))
bad.names <- colnames(subset(x = pleuro.4, idents = c(8, 21, 24)))
good.names <- cell.names[!cell.names %in% bad.names] # The good clusters must be the rest of the cells

train.bad.names<-sample(bad.names,50) # selects randomly X bad cells for training set. You will want to train your set with ~10% of the total number of "bad cells" that you have (look at length(bad.names))
train.good.names<-sample(good.names,50) # selects randomly X good cells for training set. Keep the same number of cells as above.

train.bad<-data.info[train.bad.names,]
train.good<-data.info[train.good.names,]
train<-rbind(train.good,train.bad)
ytrain <- matrix(c(rep(1,50),rep(-1,50))) # Limits are the same as your number of trainig cells

train.final = data.frame(train, ytrain=as.factor(ytrain))

#Run the clasifier 
library(e1071)
svmfit=svm(ytrain~., data=train.final, kernel="linear", cost=0.1, scale=F)

set.seed (1)
tune.out=tune(svm,ytrain~.,data=train.final,kernel="linear", ranges=list(cost=c(0.1,1,10,100,500)))
summary(tune.out)


bestmod=tune.out$best.model
summary(bestmod)

test.bad.names<-bad.names[!bad.names %in% train.bad.names] 
test.good.names<-good.names[!good.names %in% train.good.names] 
test.bad<-data.info[test.bad.names,] 
test.good<-data.info[test.good.names,] 
test<-rbind(test.good,test.bad) 

test<-data.info[!rownames(data.info) %in% rownames(train),] 

test<-as.matrix(test)
ytest<-matrix(c(rep(1,length(test.good.names)),rep(-1,length(test.bad.names)))) 

test.final = data.frame(test, ytest=as.factor(ytest))

ypred=predict(bestmod,test.final)
table(predict=ypred, truth=test.final$ytest)

ypred<-as.data.frame(ypred)


pleuro.4$ypred <- ypred

pred<-as.data.frame(pleuro.4$ypred)
rownames(pred)<-rownames(pleuro.4[[]])
pred[train.bad.names,]<-"-1"
pred[train.good.names,]<-"1"
colnames(pred)<-c("svm_class")
pleuro.4$svm_class <- pred


table(pleuro.4@meta.data$svm_class)



pleuro.5 <- subset(pleuro.4, cells=rownames(pleuro.4@meta.data[pleuro.4@meta.data$svm_class=="1",]))

pleuro.4 <- pleuro.5
rm(pleuro.5)



### check cleaned object
pleuro.4 <- SCTransform(pleuro.4, vars.to.regress = c("percent.mt","nCount_RNA","percent.mt"))

pleuro.4 <- RunPCA(pleuro.4, npcs=100)

ElbowPlot(pleuro.4, ndims=100)
ndims= 100

pleuro.4 <- RunUMAP(pleuro.4, dims=1:35)
pleuro.4 <- FindNeighbors(pleuro.4, dims=1:35)

pleuro.4 <- FindClusters(pleuro.4, resolution = 1.5)

DimPlot(pleuro.4, reduction="umap", label=T, pt.size = 0.3) + NoLegend()


St_36_clean_singles_doublets <- RunUMAP(St_36_clean_singles_doublets, dims=1:35)
St_36_clean_singles_doublets <- FindNeighbors(St_36_clean_singles_doublets, dims=1:35)

St_36_clean_singles_doublets <- FindClusters(St_36_clean_singles_doublets, resolution = 1.5)

DimPlot(St_36_clean_singles_doublets, reduction="umap", label=T, pt.size = 0.3) + NoLegend()


DimPlot(St_36_clean_singles_doublets, group.by = "HTO_classification.global")

# going to proceed with all labeled singlets

pleuro.4$Stage<-"36"

D1_D2<-readRDS("D1_D2_clean.rds")

DimPlot(D1_D2, reduction="umap", label=T, pt.size = 0.3) + NoLegend()


###merge stage 36 with D1_D2(41 and 50)
all.dev.old<-merge(pleuro.4, y=D1_D2)

saveRDS(all.dev, file="all.dev.4-6-22.rds")


st36_T2.files <- file.path("ST36_T2_alevin_out_cat_MPX/alevin/quants_mat.gz")

st36_M2.files <- file.path("multiplex_st36_M2_cmo_out/alevin/quants_mat.gz")
file.exists(c(st36_T2.files, st36_M2.files))

st36_T2.txi <- tximport(files = st36_T2.files, type = "alevin")

st36_M2.txi <- tximport(files = st36_M2.files, type = "alevin")


map <- read.table("3M-february-2018.txt.gz")
st36_T2ToADT <- map$V1
names(st36_T2ToADT) <- map$V2

st36_M2 <- tximport(files = st36_M2.files, type = "alevin")$counts

st36_T2<-st36_T2.txi$counts





colnames(st36_T2) <- st36_T2ToADT[colnames(st36_T2)]
length(intersect(colnames(st36_T2), colnames(st36_M2)))



common.cells<-(intersect(colnames(st36_T2),colnames(st36_M2)))

length(common.cells)




st_36_2_O <- CreateSeuratObject(st36_T2[, common.cells])


st36_M2<-st36_M2[1:4,]

st_36_2_O[["HTO"]] <- CreateAssayObject(counts = st36_M2[, common.cells])





VlnPlot(st_36_2_O, features = c("nCount_HTO"), log = TRUE )
st_36_2_O <- subset(st_36_2_O, subset = nCount_HTO < 2e4)
# st36_M2 Normalization
DefaultAssay(st_36_2_O) <- "HTO"
st_36_2_O <- NormalizeData(st_36_2_O,   normalization.method = "CLR", margin = 2, verbose = F)
VariableFeatures(st_36_2_O) <- rownames(st_36_2_O[["HTO"]]@counts)
st_36_2_O <- ScaleData(st_36_2_O, assay = "HTO", verbose = F)




st_36_2_O <- HTODemux(st_36_2_O, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st36_T2_class.pdf"), width=10, height=10)
Idents(st_36_2_O) <- "HTO_classification.global"
VlnPlot(st_36_2_O, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_36_2_O, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


### determine filter threshold
st36_RNA_psums<-colSums(st36_T2)



# knee plot for E7_RNA
sortedSt_36_2_RNA <- sort(st36_RNA_psums, decreasing = TRUE)
plot(sortedSt_36_2_RNA, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot St_36_2_RNA",
     col = ifelse(sortedSt_36_2_RNA > 700,'red','black'),
     xlim = c(1,length(sortedSt_36_2_RNA)), 
     ylim = c(1,sortedSt_36_2_RNA[[1]]),
     log = "xy"
)

keepSt_36_2_RNA <- colSums(st36_T2) > 700
countsSt_36_2_RNA_filter <- st36_T2[,keepSt_36_2_RNA]

sumSt_36_2_RNA <- colSums(countsSt_36_2_RNA_filter)
range(sumSt_36_2_RNA)
length(sumSt_36_2_RNA)



UMIspercellSt_36_2_RNA_sorted <- sort(sumSt_36_2_RNA,decreasing = TRUE)

pSt_36_2_RNA.a <- plot(UMIspercellSt_36_2_RNA_sorted, 
                       xlab = "Cells",
                       ylab = "UMI counts",
                       main = "Frequency distribution of UMIs per cell St_36_2_RNA",
                       xlim = c(1,length(UMIspercellSt_36_2_RNA_sorted)), 
                       ylim = c(1,UMIspercellSt_36_2_RNA_sorted[[1]]),
)

pSt_36_2_RNA.b <- hist(log10(UMIspercellSt_36_2_RNA_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellSt_36_2_RNA <- apply(countsSt_36_2_RNA_filter > 0, 2, sum)
hist(log10(genespercellSt_36_2_RNA), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes St_36_2_RNA")


plot(UMIspercellSt_36_2_RNA_sorted, genespercellSt_36_2_RNA, log='xy')
title('counts vs genes per cell St_36_2_RNA')






dev.off()


#apply filter threshold to St36 and demultiplex data 



st_36_2_OC<-subset(x=st_36_2_O, subset = nCount_RNA > 700)


st_36_2_OC <- HTODemux(st_36_2_OC, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st36_T2_class_clean.pdf"), width=10, height=10)
Idents(st_36_2_OC) <- "HTO_classification.global"
VlnPlot(st_36_2_OC, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_36_2_OC, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()






saveRDS(st_36_2_OC, file = "St36_T2_mpx_clean.rds")
saveRDS(st_36_2_O, file="St_36_T2_mpx_full.rds")

doublets<-subset(x=st_36_2_OC, idents = "Doublet")

St_36_2_singlets<-subset(x=st_36_2_OC, idents = "Singlet")



mito.genes <- read.table(file = "mitochondrial_genes_updated.txt")
mito.genes <- as.character(mito.genes[,1])
mito.genesSt_36_2_singlets <- mito.genes[mito.genes %in% rownames(St_36_2_singlets)]
St_36_2_singlets[["percent.mt"]] <- PercentageFeatureSet(St_36_2_singlets, features = mito.genesSt_36_2_singlets)

St_36_2_singlets<-subset(x=St_36_2_singlets, subset= percent.mt<15)

St_36_2_singlets <- SCTransform(St_36_2_singlets, vars.to.regress = c("percent.mt", "nCount_RNA"))
St_36_2_singlets <- RunPCA(St_36_2_singlets, features = VariableFeatures(object = St_36_2_singlets))
St_36_2_singlets <- FindNeighbors(St_36_2_singlets, dims = 1:35)
St_36_2_singlets <- FindClusters(St_36_2_singlets, resolution = 1.5)
St_36_2_singlets <- RunUMAP(St_36_2_singlets, dims = 1:35)
DimPlot(St_36_2_singlets, reduction = "umap", label = T)

FeaturePlot(St_36_2_singlets, features=c("nFeature_RNA","nCount_RNA","percent.mt"))

FeaturePlot(St_36_2_singlets, features=c("RELN","TP73","NEUROD6","TBR1"))

saveRDS(St_36_2_singlets, file = "St_36_2_singlets_clean.rds")




#Cell cleaning

pleuro.4<-St_36_2_singlets

# Create a data frame containing metadata
data.info <- pleuro.4@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")]
cell.names <- colnames(pleuro.4)

# Here you add all the idents from bad clusters

VlnPlot(pleuro.4, features = c("nCount_RNA","nFeature_RNA", "percent.mt"))
bad.names <- colnames(subset(x = pleuro.4, idents = c(19,24)))
good.names <- cell.names[!cell.names %in% bad.names] # The good clusters must be the rest of the cells

train.bad.names<-sample(bad.names,40) # selects randomly X bad cells for training set. You will want to train your set with ~10% of the total number of "bad cells" that you have (look at length(bad.names))
train.good.names<-sample(good.names,40) # selects randomly X good cells for training set. Keep the same number of cells as above.

train.bad<-data.info[train.bad.names,]
train.good<-data.info[train.good.names,]
train<-rbind(train.good,train.bad)
ytrain <- matrix(c(rep(1,40),rep(-1,40))) # Limits are the same as your number of trainig cells

train.final = data.frame(train, ytrain=as.factor(ytrain))

#Run the clasifier 
library(e1071)
svmfit=svm(ytrain~., data=train.final, kernel="linear", cost=0.1, scale=F)

set.seed (1)
tune.out=tune(svm,ytrain~.,data=train.final,kernel="linear", ranges=list(cost=c(0.1,1,10,100,500)))
summary(tune.out)


bestmod=tune.out$best.model
summary(bestmod)

test.bad.names<-bad.names[!bad.names %in% train.bad.names] 
test.good.names<-good.names[!good.names %in% train.good.names] 
test.bad<-data.info[test.bad.names,] 
test.good<-data.info[test.good.names,] 
test<-rbind(test.good,test.bad) 

test<-data.info[!rownames(data.info) %in% rownames(train),] 

test<-as.matrix(test)
ytest<-matrix(c(rep(1,length(test.good.names)),rep(-1,length(test.bad.names)))) 

test.final = data.frame(test, ytest=as.factor(ytest))

ypred=predict(bestmod,test.final)
table(predict=ypred, truth=test.final$ytest)

ypred<-as.data.frame(ypred)


pleuro.4$ypred <- ypred

pred<-as.data.frame(pleuro.4$ypred)
rownames(pred)<-rownames(pleuro.4[[]])
pred[train.bad.names,]<-"-1"
pred[train.good.names,]<-"1"
colnames(pred)<-c("svm_class")
pleuro.4$svm_class <- pred


table(pleuro.4@meta.data$svm_class)



pleuro.5 <- subset(pleuro.4, cells=rownames(pleuro.4@meta.data[pleuro.4@meta.data$svm_class=="1",]))

pleuro.4 <- pleuro.5
rm(pleuro.5)



### check cleaned object
pleuro.4 <- SCTransform(pleuro.4, vars.to.regress = c("percent.mt","nCount_RNA","percent.mt"))

pleuro.4 <- RunPCA(pleuro.4, npcs=100)

ElbowPlot(pleuro.4, ndims=100)
ndims= 100

pleuro.4 <- RunUMAP(pleuro.4, dims=1:35)
pleuro.4 <- FindNeighbors(pleuro.4, dims=1:35)

pleuro.4 <- FindClusters(pleuro.4, resolution = 1.5)

DimPlot(pleuro.4, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

ST_36_2<-pleuro.4

ST_36_2$Stage<-"36"
all.dev<-merge(all.dev, y=ST_36_2)

all.dev <- SCTransform(all.dev, vars.to.regress = c("nCount_RNA","Stage","percent.mt"))
all.dev <- RunPCA(all.dev, npcs=100)
ElbowPlot(all.dev, ndims=100)
ndims= 100
all.dev <- RunUMAP(all.dev, dims=1:30)
all.dev <- FindNeighbors(all.dev, dims=1:30)
all.dev <- FindClusters(all.dev, resolution = 1.5)
DimPlot(all.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()
all.dev <- SCTransform(all.dev, vars.to.regress = c("nCount_RNA","Stage","percent.mt"))
saveRDS(all.dev, file="all.dev.4-6-22.rds")




st41_T1.files <- file.path("ST41_T1_alevin_out_cat_MPX/alevin/quants_mat.gz")

st41_M1.files <- file.path("multiplex_st41_M1_cmo_out/alevin/quants_mat.gz")
file.exists(c(st41_T1.files, st41_M1.files))

st41_T1.txi <- tximport(files = st41_T1.files, type = "alevin")

st41_M1.txi <- tximport(files = st41_M1.files, type = "alevin")


map <- read.table("3M-february-2018.txt.gz")
st41_T1ToADT <- map$V1
names(st41_T1ToADT) <- map$V2

st41_M1 <- tximport(files = st41_M1.files, type = "alevin")$counts

st41_T1<-st41_T1.txi$counts





colnames(st41_T1) <- st41_T1ToADT[colnames(st41_T1)]
length(intersect(colnames(st41_T1), colnames(st41_M1)))



common.cells<-(intersect(colnames(st41_T1),colnames(st41_M1)))

length(common.cells)




st_41_1_O <- CreateSeuratObject(st41_T1[, common.cells])


st41_M1<-st41_M1[5:8,]

st_41_1_O[["HTO"]] <- CreateAssayObject(counts = st41_M1[, common.cells])





VlnPlot(st_41_1_O, features = c("nCount_HTO"), log = TRUE )
st_41_1_O <- subset(st_41_1_O, subset = nCount_HTO < 2e4)
# st41_M1 Normalization
DefaultAssay(st_41_1_O) <- "HTO"
st_41_1_O <- NormalizeData(st_41_1_O,   normalization.method = "CLR", margin = 2, verbose = F)
VariableFeatures(st_41_1_O) <- rownames(st_41_1_O[["HTO"]]@counts)
st_41_1_O <- ScaleData(st_41_1_O, assay = "HTO", verbose = F)




st_41_1_O <- HTODemux(st_41_1_O, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st41_T1_class.pdf"), width=10, height=10)
Idents(st_41_1_O) <- "HTO_classification.global"
VlnPlot(st_41_1_O, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_41_1_O, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


### determine filter threshold
st41_RNA_psums<-colSums(st41_T1)



# knee plot for E7_RNA
sortedSt_41_1_RNA <- sort(st41_RNA_psums, decreasing = TRUE)
plot(sortedSt_41_1_RNA, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot St_41_1_RNA",
     col = ifelse(sortedSt_41_1_RNA > 700,'red','black'),
     xlim = c(1,length(sortedSt_41_1_RNA)), 
     ylim = c(1,sortedSt_41_1_RNA[[1]]),
     log = "xy"
)

keepSt_41_1_RNA <- colSums(st41_T1) > 700
countsSt_41_1_RNA_filter <- st41_T1[,keepSt_41_1_RNA]

sumSt_41_1_RNA <- colSums(countsSt_41_1_RNA_filter)
range(sumSt_41_1_RNA)
length(sumSt_41_1_RNA)



UMIspercellSt_41_1_RNA_sorted <- sort(sumSt_41_1_RNA,decreasing = TRUE)

pSt_41_1_RNA.a <- plot(UMIspercellSt_41_1_RNA_sorted, 
                       xlab = "Cells",
                       ylab = "UMI counts",
                       main = "Frequency distribution of UMIs per cell St_41_1_RNA",
                       xlim = c(1,length(UMIspercellSt_41_1_RNA_sorted)), 
                       ylim = c(1,UMIspercellSt_41_1_RNA_sorted[[1]]),
)

pSt_41_1_RNA.b <- hist(log10(UMIspercellSt_41_1_RNA_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellSt_41_1_RNA <- apply(countsSt_41_1_RNA_filter > 0, 2, sum)
hist(log10(genespercellSt_41_1_RNA), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes St_41_1_RNA")


plot(UMIspercellSt_41_1_RNA_sorted, genespercellSt_41_1_RNA, log='xy')
title('counts vs genes per cell St_41_1_RNA')






dev.off()


#apply filter threshold to St41 and demultiplex data 



st_41_1_OC<-subset(x=st_41_1_O, subset = nCount_RNA > 700)


st_41_1_OC <- HTODemux(st_41_1_OC, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st41_T1_class_clean.pdf"), width=10, height=10)
Idents(st_41_1_OC) <- "HTO_classification.global"
VlnPlot(st_41_1_OC, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_41_1_OC, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()






saveRDS(st_41_1_OC, file = "St41_T1_mpx_clean.rds")
saveRDS(st_41_1_O, file="St_41_T1_mpx_full.rds")

doublets<-subset(x=st_41_1_OC, idents = "Doublet")

St_41_1_singlets<-subset(x=st_41_1_OC, idents = "Singlet")

St_41_1_singlets_telen.cells<-subset(St_41_1_singlets@meta.data, St_41_1_singlets@meta.data$HTO_classification%in%c("CMO305","CMO306"))
St_41_1_singlets_telen.SO<-subset(St_41_1_singlets, cells = rownames(St_41_1_singlets_telen.cells))










mito.genes <- read.table(file = "mitochondrial_genes_updated.txt")
mito.genes <- as.character(mito.genes[,1])
mito.genesSt_41_1_singlets_telen.SO <- mito.genes[mito.genes %in% rownames(St_41_1_singlets_telen.SO)]
St_41_1_singlets_telen.SO[["percent.mt"]] <- PercentageFeatureSet(St_41_1_singlets_telen.SO, features = mito.genesSt_41_1_singlets_telen.SO)

St_41_1_singlets_telen.SO<-subset(x=St_41_1_singlets_telen.SO, subset= percent.mt<15)

St_41_1_singlets_telen.SO <- SCTransform(St_41_1_singlets_telen.SO, vars.to.regress = c("percent.mt", "nCount_RNA"))
St_41_1_singlets_telen.SO <- RunPCA(St_41_1_singlets_telen.SO, features = VariableFeatures(object = St_41_1_singlets_telen.SO))
St_41_1_singlets_telen.SO <- FindNeighbors(St_41_1_singlets_telen.SO, dims = 1:35)
St_41_1_singlets_telen.SO <- FindClusters(St_41_1_singlets_telen.SO, resolution = 1.5)
St_41_1_singlets_telen.SO <- RunUMAP(St_41_1_singlets_telen.SO, dims = 1:35)
DimPlot(St_41_1_singlets_telen.SO, reduction = "umap", label = T)

FeaturePlot(St_41_1_singlets_telen.SO, features=c("nFeature_RNA","nCount_RNA","percent.mt"))

FeaturePlot(St_41_1_singlets_telen.SO, features=c("RELN","TP73","NEUROD6","TBR1"))

saveRDS(St_41_1_singlets_telen.SO, file = "St_41_1_singlets_telen.SO_clean.rds")


st46_T1.files <- file.path("ST46_T1_alevin_out_cat_MPX/alevin/quants_mat.gz")

st46_M1.files <- file.path("multiplex_st46_M1_cmo_out/alevin/quants_mat.gz")
file.exists(c(st46_T1.files, st46_M1.files))

st46_T1.txi <- tximport(files = st46_T1.files, type = "alevin")

st46_M1.txi <- tximport(files = st46_M1.files, type = "alevin")

# the below file can be found at https://github.com/10XGenomics/cellranger/blob/master/lib/python/cellranger/barcodes/translation/3M-february-2018.txt.gz
map <- read.table("3M-february-2018.txt.gz")
st46_T1ToADT <- map$V1
names(st46_T1ToADT) <- map$V2

st46_M1 <- tximport(files = st46_M1.files, type = "alevin")$counts

st46_T1<-st46_T1.txi$counts





colnames(st46_T1) <- st46_T1ToADT[colnames(st46_T1)]
length(intersect(colnames(st46_T1), colnames(st46_M1)))



common.cells<-(intersect(colnames(st46_T1),colnames(st46_M1)))

length(common.cells)




st_46_1_O <- CreateSeuratObject(st46_T1[, common.cells])


st46_M1<-st46_M1[5:8,]

st_46_1_O[["HTO"]] <- CreateAssayObject(counts = st46_M1[, common.cells])





VlnPlot(st_46_1_O, features = c("nCount_HTO"), log = TRUE )
st_46_1_O <- subset(st_46_1_O, subset = nCount_HTO < 2e4)
# st46_M1 Normalization
DefaultAssay(st_46_1_O) <- "HTO"
st_46_1_O <- NormalizeData(st_46_1_O,   normalization.method = "CLR", margin = 2, verbose = F)
VariableFeatures(st_46_1_O) <- rownames(st_46_1_O[["HTO"]]@counts)
st_46_1_O <- ScaleData(st_46_1_O, assay = "HTO", verbose = F)




st_46_1_O <- HTODemux(st_46_1_O, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st46_T1_class.pdf"), width=10, height=10)
Idents(st_46_1_O) <- "HTO_classification.global"
VlnPlot(st_46_1_O, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_46_1_O, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


### determine filter threshold
st46_RNA_psums<-colSums(st46_T1)



# knee plot for E7_RNA
sortedSt_46_1_RNA <- sort(st46_RNA_psums, decreasing = TRUE)
plot(sortedSt_46_1_RNA, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot St_46_1_RNA",
     col = ifelse(sortedSt_46_1_RNA > 700,'red','black'),
     xlim = c(1,length(sortedSt_46_1_RNA)), 
     ylim = c(1,sortedSt_46_1_RNA[[1]]),
     log = "xy"
)

keepSt_46_1_RNA <- colSums(st46_T1) > 700
countsSt_46_1_RNA_filter <- st46_T1[,keepSt_46_1_RNA]

sumSt_46_1_RNA <- colSums(countsSt_46_1_RNA_filter)
range(sumSt_46_1_RNA)
length(sumSt_46_1_RNA)



UMIspercellSt_46_1_RNA_sorted <- sort(sumSt_46_1_RNA,decreasing = TRUE)

pSt_46_1_RNA.a <- plot(UMIspercellSt_46_1_RNA_sorted, 
                       xlab = "Cells",
                       ylab = "UMI counts",
                       main = "Frequency distribution of UMIs per cell St_46_1_RNA",
                       xlim = c(1,length(UMIspercellSt_46_1_RNA_sorted)), 
                       ylim = c(1,UMIspercellSt_46_1_RNA_sorted[[1]]),
)

pSt_46_1_RNA.b <- hist(log10(UMIspercellSt_46_1_RNA_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellSt_46_1_RNA <- apply(countsSt_46_1_RNA_filter > 0, 2, sum)
hist(log10(genespercellSt_46_1_RNA), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes St_46_1_RNA")


plot(UMIspercellSt_46_1_RNA_sorted, genespercellSt_46_1_RNA, log='xy')
title('counts vs genes per cell St_46_1_RNA')






dev.off()


#apply filter threshold to St46 and demultiplex data 



st_46_1_OC<-subset(x=st_46_1_O, subset = nCount_RNA > 700)


st_46_1_OC <- HTODemux(st_46_1_OC, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st46_T1_class_clean.pdf"), width=10, height=10)
Idents(st_46_1_OC) <- "HTO_classification.global"
VlnPlot(st_46_1_OC, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_46_1_OC, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()






saveRDS(st_46_1_OC, file = "St46_T1_mpx_clean.rds")
saveRDS(st_46_1_O, file="St_46_T1_mpx_full.rds")

doublets<-subset(x=st_46_1_OC, idents = "Doublet")

St_46_1_singlets<-subset(x=st_46_1_OC, idents = "Singlet")

St_46_1_singlets_telen.cells<-subset(St_46_1_singlets@meta.data, St_46_1_singlets@meta.data$HTO_classification%in%c("CMO305","CMO306"))
St_46_1_singlets_telen.SO<-subset(St_46_1_singlets, cells = rownames(St_46_1_singlets_telen.cells))




st46_T2.files <- file.path("ST46_T2_alevin_out_cat_MPX/alevin/quants_mat.gz")

st46_M2.files <- file.path("multiplex_st46_M2_cmo_out/alevin/quants_mat.gz")
file.exists(c(st46_T2.files, st46_M2.files))

st46_T2.txi <- tximport(files = st46_T2.files, type = "alevin")

st46_M2.txi <- tximport(files = st46_M2.files, type = "alevin")


map <- read.table("3M-february-2018.txt.gz")
st46_T2ToADT <- map$V1
names(st46_T2ToADT) <- map$V2

st46_M2 <- tximport(files = st46_M2.files, type = "alevin")$counts

st46_T2<-st46_T2.txi$counts





colnames(st46_T2) <- st46_T2ToADT[colnames(st46_T2)]
length(intersect(colnames(st46_T2), colnames(st46_M2)))



common.cells<-(intersect(colnames(st46_T2),colnames(st46_M2)))

length(common.cells)




st_46_2_O <- CreateSeuratObject(st46_T2[, common.cells])


st46_M2<-st46_M2[5:8,]

st_46_2_O[["HTO"]] <- CreateAssayObject(counts = st46_M2[, common.cells])





VlnPlot(st_46_2_O, features = c("nCount_HTO"), log = TRUE )
st_46_2_O <- subset(st_46_2_O, subset = nCount_HTO < 2e4)
# st46_M2 Normalization
DefaultAssay(st_46_2_O) <- "HTO"
st_46_2_O <- NormalizeData(st_46_2_O,   normalization.method = "CLR", margin = 2, verbose = F)
VariableFeatures(st_46_2_O) <- rownames(st_46_2_O[["HTO"]]@counts)
st_46_2_O <- ScaleData(st_46_2_O, assay = "HTO", verbose = F)




st_46_2_O <- HTODemux(st_46_2_O, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st46_T2_class.pdf"), width=10, height=10)
Idents(st_46_2_O) <- "HTO_classification.global"
VlnPlot(st_46_2_O, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_46_2_O, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()


### determine filter threshold
st46_RNA_psums<-colSums(st46_T2)



# knee plot for E7_RNA
sortedSt_46_2_RNA <- sort(st46_RNA_psums, decreasing = TRUE)
plot(sortedSt_46_2_RNA, 
     xlab = "Barcodes",
     ylab = "UMI counts",
     main = "Barcode Rank Plot St_46_2_RNA",
     col = ifelse(sortedSt_46_2_RNA > 700,'red','black'),
     xlim = c(1,length(sortedSt_46_2_RNA)), 
     ylim = c(1,sortedSt_46_2_RNA[[1]]),
     log = "xy"
)

keepSt_46_2_RNA <- colSums(st46_T2) > 700
countsSt_46_2_RNA_filter <- st46_T2[,keepSt_46_2_RNA]

sumSt_46_2_RNA <- colSums(countsSt_46_2_RNA_filter)
range(sumSt_46_2_RNA)
length(sumSt_46_2_RNA)



UMIspercellSt_46_2_RNA_sorted <- sort(sumSt_46_2_RNA,decreasing = TRUE)

pSt_46_2_RNA.a <- plot(UMIspercellSt_46_2_RNA_sorted, 
                       xlab = "Cells",
                       ylab = "UMI counts",
                       main = "Frequency distribution of UMIs per cell St_46_2_RNA",
                       xlim = c(1,length(UMIspercellSt_46_2_RNA_sorted)), 
                       ylim = c(1,UMIspercellSt_46_2_RNA_sorted[[1]]),
)

pSt_46_2_RNA.b <- hist(log10(UMIspercellSt_46_2_RNA_sorted), breaks = 200, xlab = "UMIs expressed", main = "Distribution of UMIs")

genespercellSt_46_2_RNA <- apply(countsSt_46_2_RNA_filter > 0, 2, sum)
hist(log10(genespercellSt_46_2_RNA), breaks = 200, xlab = "Numer of Genes Expressed", main = "Distribution of Expressed Genes St_46_2_RNA")


plot(UMIspercellSt_46_2_RNA_sorted, genespercellSt_46_2_RNA, log='xy')
title('counts vs genes per cell St_46_2_RNA')






dev.off()


#apply filter threshold to St46 and demultiplex data 



st_46_2_OC<-subset(x=st_46_2_O, subset = nCount_RNA > 700)


st_46_2_OC <- HTODemux(st_46_2_OC, assay = "HTO", positive.quantile = 0.99, verbose = F)
pdf(paste0("st46_T2_class_clean.pdf"), width=10, height=10)
Idents(st_46_2_OC) <- "HTO_classification.global"
VlnPlot(st_46_2_OC, features = "nCount_HTO", pt.size = 0.1, log = TRUE)

VlnPlot(st_46_2_OC, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()






saveRDS(st_46_2_OC, file = "St46_T2_mpx_clean.rds")
saveRDS(st_46_2_O, file="St_46_T2_mpx_full.rds")

doublets<-subset(x=st_46_2_OC, idents = "Doublet")

St_46_2_singlets<-subset(x=st_46_2_OC, idents = "Singlet")

St_46_2_singlets_telen.cells<-subset(St_46_2_singlets@meta.data, St_46_2_singlets@meta.data$HTO_classification%in%c("CMO305","CMO306"))
St_46_2_singlets_telen.SO<-subset(St_46_2_singlets, cells = rownames(St_46_2_singlets_telen.cells))




St_46_1_singlets_telen.SO@meta.data$Stage<-"46"
St_46_2_singlets_telen.SO@meta.data$Stage<-"46"
St_41_1_singlets_telen.SO@meta.data$Stage<-"41"

merge_data <- merge(x = St_46_1_singlets_telen.SO, y = list(St_46_2_singlets_telen.SO, St_41_1_singlets_telen.SO), add.cell.ids = c("46_1", "46_2", "41_1"))

saveRDS(merge_data, file = "merge_new_data_6-9-22.rds")

mito.genes <- read.table(file = "mitochondrial_genes_updated.txt")
mito.genes <- as.character(mito.genes[,1])
mito.genesmerge_data <- mito.genes[mito.genes %in% rownames(merge_data)]
merge_data[["percent.mt"]] <- PercentageFeatureSet(merge_data, features = mito.genesmerge_data)

merge_data<-subset(x=merge_data, subset= percent.mt<15)
merge_data@meta.data$percent.mt

merge_data <- SCTransform(merge_data, vars.to.regress = c("percent.mt", "nCount_RNA"))
merge_data <- RunPCA(merge_data, features = VariableFeatures(object = merge_data))
merge_data <- FindNeighbors(merge_data, dims = 1:35)
merge_data <- FindClusters(merge_data, resolution = 1.5)
merge_data <- RunUMAP(merge_data, dims = 1:35)
DimPlot(merge_data, reduction = "umap", label = T)

FeaturePlot(merge_data, features=c("nFeature_RNA","nCount_RNA","percent.mt"))

FeaturePlot(merge_data, features=c("RELN","TP73","NEUROD6","TBR1"))

saveRDS(merge_data, file = "merge_data_clean.6-10-22.rds")



#Cell cleaning

pleuro.4<-merge_data

rm(merge_data)

# Create a data frame containing metadata
data.info <- pleuro.4@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt")]
cell.names <- colnames(pleuro.4)

# Here you add all the idents from bad clusters

VlnPlot(pleuro.4, features = c("nCount_RNA","nFeature_RNA", "percent.mt"))
bad.names <- colnames(subset(x = pleuro.4, idents = c(4,15,26,30)))
good.names <- cell.names[!cell.names %in% bad.names] # The good clusters must be the rest of the cells

train.bad.names<-sample(bad.names,100) # selects randomly X bad cells for training set. You will want to train your set with ~10% of the total number of "bad cells" that you have (look at length(bad.names))
train.good.names<-sample(good.names,100) # selects randomly X good cells for training set. Keep the same number of cells as above.

train.bad<-data.info[train.bad.names,]
train.good<-data.info[train.good.names,]
train<-rbind(train.good,train.bad)
ytrain <- matrix(c(rep(1,100),rep(-1,100))) # Limits are the same as your number of trainig cells

train.final = data.frame(train, ytrain=as.factor(ytrain))

#Run the clasifier 
library(e1071)
svmfit=svm(ytrain~., data=train.final, kernel="linear", cost=0.1, scale=F)

set.seed (1)
tune.out=tune(svm,ytrain~.,data=train.final,kernel="linear", ranges=list(cost=c(0.1,1,10,100,500)))
summary(tune.out)


bestmod=tune.out$best.model
summary(bestmod)

test.bad.names<-bad.names[!bad.names %in% train.bad.names] 
test.good.names<-good.names[!good.names %in% train.good.names] 
test.bad<-data.info[test.bad.names,] 
test.good<-data.info[test.good.names,] 
test<-rbind(test.good,test.bad) 

test<-data.info[!rownames(data.info) %in% rownames(train),] 

test<-as.matrix(test)
ytest<-matrix(c(rep(1,length(test.good.names)),rep(-1,length(test.bad.names)))) 

test.final = data.frame(test, ytest=as.factor(ytest))

ypred=predict(bestmod,test.final)
table(predict=ypred, truth=test.final$ytest)

ypred<-as.data.frame(ypred)


pleuro.4$ypred <- ypred

pred<-as.data.frame(pleuro.4$ypred)
rownames(pred)<-rownames(pleuro.4[[]])
pred[train.bad.names,]<-"-1"
pred[train.good.names,]<-"1"
colnames(pred)<-c("svm_class")
pleuro.4$svm_class <- pred


table(pleuro.4@meta.data$svm_class)



pleuro.5 <- subset(pleuro.4, cells=rownames(pleuro.4@meta.data[pleuro.4@meta.data$svm_class=="1",]))

pleuro.4 <- pleuro.5
rm(pleuro.5)



### check cleaned object
pleuro.4 <- SCTransform(pleuro.4, vars.to.regress = c("percent.mt","nCount_RNA","percent.mt"))

pleuro.4 <- RunPCA(pleuro.4, npcs=100)

ElbowPlot(pleuro.4, ndims=100)
ndims= 100

pleuro.4 <- RunUMAP(pleuro.4, dims=1:35)
pleuro.4 <- FindNeighbors(pleuro.4, dims=1:35)

pleuro.4 <- FindClusters(pleuro.4, resolution = 1.5)

DimPlot(pleuro.4, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

FeaturePlot(pleuro.4, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(pleuro.4, features = c("SOX6","TP73","RELN","FEZF2"))


saveRDS(pleuro.4, file="SVM_clean_6-10_newdata_merge.rds")

all.dev.new<-pleuro.4



all.dev<-merge(all.dev.old, all.dev.new)



rm(all.dev.new)
rm(all.dev.old)

#############
#############
###begin clustering and analyzing data

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
all.dev <- CellCycleScoring(all.dev, s.features = s.genes, g2m.features = g2m.genes)



all.dev <- SCTransform(all.dev, vars.to.regress = c("nCount_RNA","Stage","percent.mt","G2M.Score","S.Score"), variable.features.n = 3000 )

all.dev <- RunPCA(all.dev, npcs=100)

ElbowPlot(all.dev, ndims=100)
ndims= 100

all.dev <- RunUMAP(all.dev, dims=1:30)
all.dev <- FindNeighbors(all.dev, dims=1:30)

all.dev <- FindClusters(all.dev, resolution = 1.5)

DimPlot(all.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(all.dev, group.by  ="Stage", label=T, pt.size = 0.3)

FeaturePlot(all.dev, features = c("GAD1","NEUROD6","TP73","TBR1"))

FeaturePlot(all.dev, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(all.dev, features = c("FOXG1","NEUROD6","GAD1","TBX21"))
FeaturePlot(all.dev, features = c("LHX1","SLC17A6"))
FeaturePlot(all.dev, features = c("GFAP","PCNA", "SNAP25"), label = T)




#note, it is possible that the actual cluster numbers will change below depending on which cells were filtered during SVM. If re-running analysis must check that clusters you pull out still make sense
markers_33<-FindMarkers(all.dev,"33") #CPXM1= endothelial (remove)
markers_37<-FindMarkers(all.dev,"37")  #SPTB retina (keep)
markers_36<-FindMarkers(all.dev,"36") #sox10 UGT8 OPC remove
markers_38<-FindMarkers(all.dev,"38") # ARHGAP25 (remove)
markers_39<-FindMarkers(all.dev,"39")   #FLI1 SLC2A2 adipose/epithelia 

neur.dev<-subset(all.dev, ident=c("33","36","38","39"),invert=T)



neur.dev <- SCTransform(neur.dev, vars.to.regress = c("nCount_RNA","Stage","percent.mt", "S.Score","G2M.Score"), variable.features.n = 3000)

neur.dev <- RunPCA(neur.dev, npcs=100)

ElbowPlot(neur.dev, ndims=100)
ndims= 100
#orig 25 dim
neur.dev <- RunUMAP(neur.dev, dims=1:35)
neur.dev <- FindNeighbors(neur.dev, dims=1:35)

neur.dev <- FindClusters(neur.dev, resolution = 3)

DimPlot(neur.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(neur.dev, group.by  ="Stage", label=T, pt.size = 0.3)
DimPlot(neur.dev, group.by  ="Phase", label=T, pt.size = 0.3)




FeaturePlot(neur.dev, features = c("GAD1","NEUROD6","TP73","TBR1"))

FeaturePlot(neur.dev, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(neur.dev, features = c("FOXG1","NEUROD6","GAD1","TBX21"))
FeaturePlot(neur.dev, features = c("LHX1","SLC17A6"))
FeaturePlot(neur.dev, features = c("GFAP","PCNA", "NEUROG2", "SNAP25"))

## This vln plot as well as below label transfer used to identify telenephalic clusters
VlnPlot(neur.dev, features = "FOXG1")

saveRDS(neur.dev, file="neur.dev.big.data.6-11.rds")



#### below code is run on cluser. It is label transfer of adult to development
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(corrplot)
library(MatrixGenerics)

query<-readRDS("neur.dev.big.data.6-11.rds")
reference<-readRDS("EG_pleuro_high_annot.rds")


Idents(reference) <- reference$high.level.annot.tosches
query<-AddMetaData(query, (paste('D', as.character(Idents(query), sep="-"))),  col.name = "cluster_label")

reference<-AddMetaData(reference, (paste(as.character(Idents(reference)))), col.name ="cluster_label")

reference <- RunUMAP(reference, dims = 1:180, reduction = "pca", return.model = TRUE )

neur.anchors <- FindTransferAnchors(reference = reference, query = query,
                                    dims = 1:50, reference.reduction = "pca", k.anchor=10)
# k.weight default is 50; setting it to 20 in this case because otherwise it throws an error
predictions <- TransferData(anchorset = neur.anchors, refdata = reference$cluster_label,
                            dims = 1:50, k.weight=20,  weight.reduction = "pcaproject")
query <- AddMetaData(query, metadata = predictions)



cells <- rownames(query@meta.data[query@meta.data$prediction.score.max>0.75,])



pdf(paste0("6-11LT_predictedID_neur.dev.pleur.pdf"), width=20, height=10)
DimPlot(query, cells= cells, group.by="predicted.id", label=T)
DimPlot(query, label = T) 
dev.off()

pdf(paste0("vlnplots_6-15.pdf"))
VlnPlot(neur.dev, "nFeature_RNA", pt.size = 0)
VlnPlot(neur.dev, "nCount_RNA", pt.size = 0)
VlnPlot(neur.dev, "percent.mt", pt.size = 0)
dev.off()


pdf(paste0("neur.dev.plot.pdf"), width=8, height=6)
DimPlot(neur.dev, label = T)+NoLegend()
FeaturePlot(neur.dev, "FOXG1")
dev.off()


#neur.dev<-readRDS("neur.dev.big.data.6-11.rds")


## return to running loccally 

telen.dev<-subset(neur.dev, ident=c("49","46","45","39","35","38","47","30","53","44","22","37"), invert=T)



rm(neur.dev)



telen.dev <- SCTransform(telen.dev, vars.to.regress = c("nCount_RNA","percent.mt", "Stage", "S.Score","G2M.Score"), variable.features.n = 3000)

telen.dev <- RunPCA(telen.dev, npcs=100)

ElbowPlot(telen.dev, ndims=100)
ndims= 100

telen.dev <- RunUMAP(telen.dev, dims=1:35)
telen.dev <- FindNeighbors(telen.dev, dims=1:35)

telen.dev <- FindClusters(telen.dev, resolution = 1.5)

DimPlot(telen.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(telen.dev, group.by  ="Stage", label=T, pt.size = 0.3)
DimPlot(telen.dev, group.by  ="Phase", label=T, pt.size = 0.3)
FeaturePlot(telen.dev, features = c("GAD1","NEUROD6","TP73","TBR1"))

FeaturePlot(telen.dev, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(telen.dev, features = c("FOXG1","NEUROD6","GAD1","TBX21"))
FeaturePlot(telen.dev, features = c("LHX1","SLC17A6"))
FeaturePlot(telen.dev, features = c("EMX1","DLX1","NEUROG2"), label = T)

FeaturePlot(telen.dev, features = c("LHX9","NFIX","PROX1", "ETV1"))

FeaturePlot(telen.dev, features = c("SIM1"))

VlnPlot(telen.dev, "NEUROG2")


saveRDS(telen.dev, file="6-13.telen.dev.big.data.rds")


DimPlot(telen.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()
#do LT for telen dev against all adult cells


query<-readRDS("6-13.telen.dev.big.data.rds")
reference<-readRDS("EG_pleuro_high_annot.2.rds")


Idents(reference) <- reference$high.level.annot.tosches

query<-AddMetaData(query, (paste('D', as.character(Idents(query), sep="-"))),  col.name = "cluster_label")

reference<-AddMetaData(reference, (paste(as.character(Idents(reference)))), col.name ="cluster_label")

reference <- RunUMAP(reference, dims = 1:180, reduction = "pca", return.model = TRUE )

neur.anchors <- FindTransferAnchors(reference = reference, query = query,
                                    dims = 1:50, reference.reduction = "pca", k.anchor=10)
# k.weight default is 50; setting it to 20 in this case because otherwise it throws an error
predictions <- TransferData(anchorset = neur.anchors, refdata = reference$cluster_label,
                            dims = 1:50, k.weight=20,  weight.reduction = "pcaproject")
query <- AddMetaData(query, metadata = predictions)



cells <- rownames(query@meta.data[query@meta.data$prediction.score.max>0.75,])



pdf(paste0("6-13LT_telen.dev.35.pleur.pdf"), width=20, height=10)
DimPlot(query, cells= cells, group.by="predicted.id", label=T)
DimPlot(query, label = T) 
dev.off()

telen.dev<-readRDS("6-13LT_telen.dev.35.pleur.rds")



pal.dev<-subset(telen.dev, ident=c("16","5","24","7","15","14","12","0","27","30","29","8"))


pal.dev <- SCTransform(pal.dev, vars.to.regress = c("nCount_RNA","percent.mt", "Stage", "S.Score","G2M.Score"), variable.features.n = 3000)

pal.dev <- RunPCA(pal.dev, npcs=100)

ElbowPlot(pal.dev, ndims=100)
ndims= 100



pal.dev <- RunUMAP(pal.dev, dims=1:53)
pal.dev <- FindNeighbors(pal.dev, dims=1:53)

pal.dev <- FindClusters(pal.dev, resolution = 1)

DimPlot(pal.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

saveRDS(pal.dev, file="pal.dev.50.neurog2.6-13.rds")

DimPlot(pal.dev, group.by  ="Stage", label=T, pt.size = 0.3)

FeaturePlot(pal.dev, features = c("SIM1","PCNA","EOMES"), label=T)



###LT of adult neurons onto pallial development (run on HPC Cluster)


query<-readRDS("pal.dev.50.neurog2.6-13.rds")
##neuronal object from adult analysis
reference<-readRDS("EG_neurons_final6.2.rds")

Idents(reference) <- reference$high_anno

query<-AddMetaData(query, (paste('D', as.character(Idents(query), sep="-"))),  col.name = "cluster_label")

reference<-AddMetaData(reference, (paste(as.character(Idents(reference)))), col.name ="cluster_label")

reference <- RunUMAP(reference, dims = 1:180, reduction = "pca", return.model = TRUE )

neur.anchors <- FindTransferAnchors(reference = reference, query = query,
                                    dims = 1:50, reference.reduction = "pca", k.anchor=10)
# k.weight default is 50; setting it to 20 in this case because otherwise it throws an error
predictions <- TransferData(anchorset = neur.anchors, refdata = reference$cluster_label,
                            dims = 1:50, k.weight=20,  weight.reduction = "pcaproject")
query <- AddMetaData(query, metadata = predictions)



cells <- rownames(query@meta.data[query@meta.data$prediction.score.max>0.75,])



pdf(paste0("6-13LT_pal.dev.neurog2.50.neurons.pdf"), width=20, height=10)
DimPlot(query, cells= cells, group.by="predicted.id", label=T)
DimPlot(query, label = T) 
dev.off()


saveRDS(query, "6-13_LT.pal.dev.neurog2.50.neur.rds" )



####back locally





pal.dev<-readRDS("6-13_LT.pal.dev.neurog2.50.neur.rds")


#Select cells that belong to terminal points for trajectory analysis 

VP_cells <- rownames(pal.dev@meta.data[pal.dev@meta.data$prediction.score.VP>0.75,])
MP_cells <- rownames(pal.dev@meta.data[pal.dev@meta.data$prediction.score.MP>0.75,])
DP_cells <- rownames(pal.dev@meta.data[pal.dev@meta.data$prediction.score.DP>0.75,])
LP_cells <- rownames(pal.dev@meta.data[pal.dev@meta.data$prediction.score.LP>0.75,])

AMY_cells <- rownames(pal.dev@meta.data[pal.dev@meta.data$prediction.score.AMY>0.75,])

# MT_cells <- rownames(pal.dev@meta.data[pal.dev@meta.data$prediction.score.MT>0.75,])



length(VP_cells)
length(MP_cells)
length(DP_cells)
length(LP_cells)
length(AMY_cells)



pal.dev@meta.data$new_ident<-"NA"
#condMT <- rownames(pal.dev@meta.data)%in%MT_cells
condAMY <- rownames(pal.dev@meta.data)%in%AMY_cells
condVP <- rownames(pal.dev@meta.data)%in%VP_cells
condLP <- rownames(pal.dev@meta.data)%in%LP_cells
condMP <- rownames(pal.dev@meta.data)%in%MP_cells
condDP <- rownames(pal.dev@meta.data)%in%DP_cells



#Assign
#pal.dev@meta.data$new_ident[condMT]<-'MT'
pal.dev@meta.data$new_ident[condAMY]<-'AMY'
pal.dev@meta.data$new_ident[condVP]<-'VP'
pal.dev@meta.data$new_ident[condLP]<-'LP'
pal.dev@meta.data$new_ident[condMP]<-'MP'
pal.dev@meta.data$new_ident[condDP]<-'DP'


pal.dev@meta.data$new_ident<-na_if(pal.dev@meta.data$new_ident, "NA")

pal.dev@meta.data<-(pal.dev@meta.data %>%
                      mutate(new_ident = coalesce(new_ident, seurat_clusters)))

Idents(pal.dev)<-pal.dev@meta.data$new_ident


DimPlot(pal.dev, group.by = "new_ident", label=T, pt.size = 0.3)

#remove cells from terminal cluster that were not above .75 predicted id for a pallial region

pal.dev.clean<-subset(pal.dev, ident=c("0","14","8","15","9","18","16"), invert=T)


DimPlot(pal.dev.clean, label=T, pt.size = 0.3)






pal.dev.clean.sce <- as.SingleCellExperiment(pal.dev.clean)

seu<-pal.dev.clean.sce

umap_embed<-pal.dev.clean@reductions$umap@cell.embeddings
reducedDims(pal.dev.clean.sce)<-list(UMAP = umap_embed)

reducedDim(pal.dev.clean.sce)

seu<-pal.dev.clean.sce



sds <- slingshot(seu, reducedDim= "UMAP", clusterLabels = seu$ident,
                 start.clus = c("5"), end.clus = c("1","MP","DP","LP","VP","AMY"), reweight=T, reassign=T, extend='n',maxit=50 )

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}


cell_colors_clust <- cell_pal(seu$ident, hue_pal())




nc <- 3
pt <- slingPseudotime(SlingshotDataSet(sds))
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(SlingshotDataSet(sds), lwd=2,  col='black')
}



cw<-slingCurveWeights(sds)
cw.df<-data.frame(cw)
pt<-slingPseudotime(sds)
pt.df<-data.frame(pt)


pal.dev.clean<-AddMetaData(pal.dev.clean, cw.df)
FeaturePlot(pal.dev.clean, features = c("Lineage1","Lineage2","Lineage4","Lineage5"))


telen.dev.pal<-subset(telen.dev, cells= Cells(pal.dev.clean))

DimPlot(telen.dev.pal)


telen.dev.pal.cw<-AddMetaData(telen.dev.pal, cw.df)
telen.dev.pal.pt<-AddMetaData(telen.dev.pal, pt.df)



#plot trajectories for figure
pdf(paste0("spectral_traj_6-15.telen_UMAP.pdf"),width = 8, height=6)
FeaturePlot(telen.dev.pal.cw, features = c("Lineage1")) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(telen.dev.pal.cw, features = c("Lineage2")) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(telen.dev.pal.cw, features = c("Lineage4")) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
FeaturePlot(telen.dev.pal.cw, features = c("Lineage5")) + scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "Spectral")))
dev.off()




pdf(paste0("6-15pseudotime_map_lines.pdf"), width=8, height=6)
FeaturePlot(telen.dev.pal.pt, features = "Lineage1", label = F, cols = pal)
FeaturePlot(telen.dev.pal.pt, features = "Lineage2", label = F, cols = pal)
FeaturePlot(telen.dev.pal.pt, features = "Lineage4", label = F, cols = pal)
FeaturePlot(telen.dev.pal.pt, features = "Lineage5", label = F, cols = pal)
dev.off()





#rm(cw.df)


cw.df$max<-colnames(cw.df)[max.col(cw.df,ties.method="random")]

#saveRDS(cw.df, file="cw.df.6-16.rds")

DP<-subset(cw.df, cw.df$max=="Lineage2" )
VP<-subset(cw.df, cw.df$max=="Lineage4" )

DP.95<-subset(DP ,(DP$Lineage2>.95))
VP.95<-subset(VP ,(VP$Lineage4>.95))


no.prog.DP<-subset(pt.df, (pt.df$Lineage2>7))
no.prog.VP<-subset(pt.df, (pt.df$Lineage4>7))






DP.95.11<-subset(DP.95, rownames(DP.95)%in%rownames(no.prog.DP))
VP.95.11<-subset(VP.95, rownames(VP.95)%in%rownames(no.prog.VP))

share_DV<-intersect(rownames(DP.95.11), rownames(VP.95.11))




pal.dev.clean@meta.data$traj<-"NA"

condDP <- rownames(pal.dev.clean@meta.data)%in%rownames(DP.95.11)
condVP <- rownames(pal.dev.clean@meta.data)%in%rownames(VP.95.11)


#Assign
pal.dev.clean@meta.data$traj[condDP]<-'DP.traj'
pal.dev.clean@meta.data$traj[condVP]<-'VP.traj'



pal.dev.clean@meta.data$traj<-na_if(pal.dev.clean@meta.data$traj, "NA")

#pal.dev.clean@meta.data<-(pal.dev.clean@meta.data %>%
#  mutate(traj = coalesce(traj, seurat_clusters)))

Idents(pal.dev.clean)<-pal.dev.clean@meta.data$traj


DimPlot(pal.dev.clean, group.by = "traj", label=T, pt.size = 0.3)





DP.VP.genes<-FindMarkers(pal.dev.clean, "DP.traj","VP.traj")



DP.VP.genes$pct.dif<-abs(DP.VP.genes$pct.1- DP.VP.genes$pct.2)
DP.VP.genes$pct.ratio<-(DP.VP.genes$pct.1/DP.VP.genes$pct.2)
DP.VP.genes$gene<-rownames(DP.VP.genes)
top_DP.VP.genes<-subset(DP.VP.genes, (DP.VP.genes$pct.dif>.1 &(DP.VP.genes$pct.ratio>2 | DP.VP.genes$pct.ratio<(1/2))))

#View(top_DP.VP.genes)



#View(pt.df)



pt.df.DP<-subset(pt.df, rownames(pt.df)%in%rownames(DP.95))
pt.df.VP<-subset(pt.df, rownames(pt.df)%in%rownames(VP.95))

pt.df.DP.only<-pt.df.DP[2]
pt.df.VP.only<-pt.df.VP[4]


pt.df.DP.only.order <- pt.df.DP.only[order(pt.df.DP.only[, "Lineage2"]), , drop = FALSE]

pt.df.VP.only.order <- pt.df.VP.only[order(pt.df.VP.only[, "Lineage4"]), , drop = FALSE]


pal.dev.clean.mx<-as.matrix(pal.dev.clean@assays$RNA@counts)

DP.dev.clean.mx<-pal.dev.clean.mx[,rownames(pt.df.DP.only.order)]
VP.dev.clean.mx<-pal.dev.clean.mx[,rownames(pt.df.VP.only.order)]



#pt.df.DP<-pt.df[2]
#pt.df.VP<-pt.df[4]


#pt.df.DP.clean<-subset(pt.df.DP, rownames(pt.df.DP)%in%rownames(DP.95.11))
#pt.df.VP.clean<-subset(pt.df.VP, rownames(pt.df.VP)%in%rownames(VP.95.11))


#pt.df.DP.clean <- pt.df.DP.clean[order(pt.df.DP.clean[, "Lineage2"]), , drop = FALSE]

#pt.df.VP.clean <- pt.df.VP.clean[order(pt.df.VP.clean[, "Lineage4"]), , drop = FALSE]

#dev.DP.traj<-subset(pal.dev.clean, idents="DP.traj")
#dev.VP.traj<-subset(pal.dev.clean, idents="VP.traj")

dev.DP.traj.df<-data.frame(dev.DP.traj@assays$RNA@counts)
dev.VP.traj.df<-data.frame(dev.VP.traj@assays$RNA@counts)

top_DP.genes<-subset(top_DP.VP.genes, top_DP.VP.genes$avg_log2FC>0)
top_VP.genes<-subset(top_DP.VP.genes, top_DP.VP.genes$avg_log2FC<0)


DP.colnames<-colnames(dev.DP.traj.df)
DP.colnames.fix<-gsub("X", "", DP.colnames)
colnames(dev.DP.traj.df)<-DP.colnames.fix

VP.colnames<-colnames(dev.VP.traj.df)
VP.colnames.fix<-gsub("X", "", VP.colnames)
colnames(dev.VP.traj.df)<-VP.colnames.fix


dev.DP.traj.df.sub<-subset(dev.DP.traj.df, rownames(dev.DP.traj.df)%in%top_DP.genes$gene)
dev.VP.traj.df.sub<-subset(dev.VP.traj.df, rownames(dev.VP.traj.df)%in%top_VP.genes$gene)

View(dev.DP.traj.df.sub)


dev.DP.traj.df.sub <- dev.DP.traj.df.sub[, rownames(pt.df.DP.clean)]
dev.VP.traj.df.sub <- dev.VP.traj.df.sub[, rownames(pt.df.VP.clean)]

library(zoo)


#max_cell_DP<-data.frame(rownames(dev.DP.traj.df.sub))
#max_cell_DP$max_cell<-max.col(dev.DP.traj.df.sub)


#max_cell_VP<-data.frame(rownames(dev.VP.traj.df.sub))
#max_cell_VP$max_cell<-max.col(dev.VP.traj.df.sub)

avg_roll_DP<-t(rollmean(t(dev.DP.traj.df.sub), 5))
avg_roll_VP<-t(rollmean(t(dev.VP.traj.df.sub), 5))




max_cell_DP<-data.frame(rownames(avg_roll_DP))
max_cell_DP$max_cell<-max.col(avg_roll_DP)


max_cell_VP<-data.frame(rownames(avg_roll_VP))
max_cell_VP$max_cell<-max.col(avg_roll_VP)


max_cell_DP_order <- max_cell_DP[order(max_cell_DP[, "max_cell"]), , drop = FALSE]
max_cell_VP_order <- max_cell_VP[order(max_cell_VP[, "max_cell"]), , drop = FALSE]








pt.df.DP.only.order.neg<-pt.df.DP.only.order*(-1)

curve.DP<-data.frame("cell"=rownames(pt.df.DP.only.order.neg), "pt"=pt.df.DP.only.order.neg$Lineage2)

curve.VP<-data.frame("cell"=rownames(pt.df.VP.only.order), "pt"=pt.df.VP.only.order$Lineage4)

rownames(curve.DP)<-curve.DP$cell
rownames(curve.VP)<-curve.VP$cell

curve.merge<-rbind(curve.DP, curve.VP)



curve.merge.order <- curve.merge[order(curve.merge[, "pt"]), , drop = FALSE]

rownames(curve.merge.order)<-curve.merge.order$cell

merge.t<-curve.merge.order$pt


#dev_curve_merge<-subset(pal.dev.clean, cells= curve.merge.order$cell)
pal.dev.clean.mx.order<-pal.dev.clean.mx[,curve.merge.order$cell]

#dev_curve_merge


TFs_mouse <- read.table("~/Toshces_lab/Tosches_r/TFs_mouse.txt", quote="\"", comment.char="")


order_genes_DP<-max_cell_DP_order$rownames.avg_roll_DP.
order_genes_DP_label<-data.frame("gene"=order_genes_DP, "order"=c(1:76))

top_TFs_DP_order<-subset(order_genes_DP_label, order_genes_DP_label$gene%in%TFs_mouse$V1)



order_genes_VP<-max_cell_VP_order$rownames.avg_roll_VP.
order_genes_VP_label<-data.frame("gene"=order_genes_VP, "order"=c(1:97))

top_TFs_VP_order<-subset(order_genes_VP_label, order_genes_VP_label$gene%in%TFs_mouse$V1)








col_fun = colorRamp2(c(-2,0, 6), c("blue", "white", "red"))
col_fun(seq(-3, 3))



heatdata_merge <- as.matrix(pal.dev.clean.mx.order[rownames(pal.dev.clean.mx.order) %in% top_DP.genes$gene, order(merge.t, na.last = NA)])

heatdata_merge_order<- heatdata_merge[max_cell_DP_order$rownames.avg_roll_DP., ]

heatdata_merge_order_scale = t(scale(t(heatdata_merge_order)))


#library("ComplexHeatmap")
#library("zoo")


#note: some values below will need to change each time due to random assignment of ties at the start

ha = rowAnnotation(foo = anno_mark(at = c(1,2,5,8,10,13,38,58,69), labels = top_TFs_DP_order$gene))


pdf(paste0("6-16.heatmap_pseudotime_merge_DP.95.7.2xb.test.smooth.final.pdf"), width=10, height=10)
Heatmap(
  heatdata_merge_order_scale,
  cluster_rows=FALSE, cluster_columns=FALSE,  show_row_names = FALSE,
  right_annotation = ha,
  show_column_names=FALSE, heatmap_height = unit(20, "cm"),
  name = "heatmap", column_split =  rep(c("DP Trajectory","VP Trajectory"),c(737,810)))
dev.off()





heatdata_merge <- as.matrix(pal.dev.clean.mx.order[rownames(pal.dev.clean.mx.order) %in% top_VP.genes$gene, order(merge.t, na.last = NA)])

heatdata_merge_order<- heatdata_merge[max_cell_VP_order$rownames.avg_roll_VP., ]

heatdata_merge_order_scale = t(scale(t(heatdata_merge_order)))



hv = rowAnnotation(foo = anno_mark(at = c(4,8,10,12,23,35), labels = top_TFs_VP_order$gene))


pdf(paste0("6-16.heatmap_pseudotime_merge_VP.95.7.2xb.pdf"), width=10, height=10)
Heatmap(
  heatdata_merge_order_scale,
  cluster_rows=FALSE, cluster_columns=FALSE, show_row_names = FALSE,
  right_annotation = hv, show_column_names=FALSE, heatmap_height = unit(20, "cm"),
  name = "heatmap", column_split =  rep(c("DP Trajectory","VP Trajectory"),c(737,810)))
dev.off()


# make feature plots for trajecory genes 

pdf(paste0("6-16VP_fp.pdf"), width=8, height=6)
FeaturePlot(telen.dev.pal, features = "PBX3")
FeaturePlot(telen.dev.pal, features = "SOX6")
FeaturePlot(telen.dev.pal, features = "EBF4")
FeaturePlot(telen.dev.pal, features = "TSHZ2")
FeaturePlot(telen.dev.pal, features = "TSHZ3")
FeaturePlot(telen.dev.pal, features = "PCP4")
dev.off()

pdf(paste0("6-16DP_fp.pdf"), width=8, height=6)
FeaturePlot(telen.dev.pal, features = "NHLH1")
FeaturePlot(telen.dev.pal, features = "FEZF2")
FeaturePlot(telen.dev.pal, features = "LHX2")
FeaturePlot(telen.dev.pal, features = "SOX8")
FeaturePlot(telen.dev.pal, features = "FOXO6")
FeaturePlot(telen.dev.pal, features = "NFIX")
FeaturePlot(telen.dev.pal, features = "ST18")
FeaturePlot(telen.dev.pal, features = "ETV1")
FeaturePlot(telen.dev.pal, features = "FOXP1")
dev.off()


#### make transfer plot for telendev


reference<-subset(reference, idents=c("Inhibitory telencephalic", "NO-TE OTHER","Non-telencephalic", "Olfactory ensheating cells","TE OTHER","Vascular leptomeningeal cells"), invert=T)

query<-readRDS("6-13.telen.dev.big.data.rds")



query<-AddMetaData(query, (paste('D', as.character(Idents(query), sep="-"))),  col.name = "cluster_label")

reference<-AddMetaData(reference, (paste(as.character(Idents(reference)))), col.name ="cluster_label")

reference <- RunUMAP(reference, dims = 1:180, reduction = "pca", return.model = TRUE )

neur.anchors <- FindTransferAnchors(reference = reference, query = query,
                                    dims = 1:50, reference.reduction = "pca", k.anchor=10)
# k.weight default is 50; setting it to 20 in this case because otherwise it throws an error
predictions <- TransferData(anchorset = neur.anchors, refdata = reference$cluster_label,
                            dims = 1:50, k.weight=20,  weight.reduction = "pcaproject")
query <- AddMetaData(query, metadata = predictions)



cells <- rownames(query@meta.data[query@meta.data$prediction.score.max>0.75,])



pdf(paste0("6-15LT.telen.dev.clean.cells.pdf"), width=20, height=10)
DimPlot(query, cells= cells, group.by="predicted.id", label=T)
DimPlot(query, label = T) 
dev.off()



Idents(telen.dev)<-"predicted.id"

levels(telen.dev)
#cols=c("#17B5A9","#EA5F5B","#FF6699","#DE4AF0","MT","#33CC99","#CCBB44","IT","#8666A9","OBG","#6773B5","#458AC9", "TE_other","#2DA736")

cols=c("#ABBCBB","#EA5F5B","#FF6699","#33CC99","#DE4AF0","#3F270A","#CCBB44","#8666A9","#934790","#6773B5","#458AC9","#2DA736")

DimPlot(telen.dev, label = T, cols=cols)



Idents(telen.dev)<-"Stage"
telen.dev<-RenameIdents(telen.dev, `36`="36", `D1_x`="41", `D2_x`="50", `41`="41")

pdf(paste0("Stage_UMAP_TELEN.pdf"), width=8, height = 6)
DimPlot(telen.dev, shuffle = T)
dev.off()


telen.dev<-readRDS("6-13.telen.dev.big.data.rds")

telen.dev.46.cells<-rownames(subset(telen.dev@meta.data, telen.dev@meta.data$Stage=="46"))
telen.dev.46<-subset(telen.dev, cells=telen.dev.46.cells)
telen.dev.46.mx<-as.matrix(telen.dev.46@assays$RNA@counts)
View(telen.dev.46.mx)
saveRDS(telen.dev.46.mx, file="telen.dev.46.mx")


telen.dev.36.cells<-rownames(subset(telen.dev@meta.data, telen.dev@meta.data$Stage=="36"))
telen.dev.36<-subset(telen.dev, cells=telen.dev.36.cells)
telen.dev.36.mx<-as.matrix(telen.dev.36@assays$RNA@counts)

saveRDS(telen.dev.36.mx, file="telen.dev.36.mx.rds")



telen.dev.41.cells<-rownames(subset(telen.dev@meta.data, (telen.dev@meta.data$Stage=="41"| telen.dev@meta.data$Stage=="D1_x")))
telen.dev.41<-subset(telen.dev, cells=telen.dev.41.cells)
telen.dev.41.mx<-as.matrix(telen.dev.41@assays$RNA@counts)

saveRDS(telen.dev.41.mx, file="telen.dev.41.mx.rds")



