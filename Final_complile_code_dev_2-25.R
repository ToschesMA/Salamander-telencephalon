

library(tximport)
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



### bring back in doublets to identify potential doublets. 

#Probabaly going to skip this aspect as seems like very minor population of doublets
mito.genes <- read.table(file = "mitochondrial_genes_updated.txt")
mito.genes <- as.character(mito.genes[,1])

mito.genesSt_36_1_doublets <- mito.genes[mito.genes %in% rownames(doublets)]
doublets[["percent.mt"]] <- PercentageFeatureSet(doublets, features = mito.genesSt_36_1_doublets)



St_36_clean_singles_doublets <- merge(pleuro.4, y =doublets)


St_36_clean_singles_doublets <- SCTransform(St_36_clean_singles_doublets, vars.to.regress = c("nFeature_RNA","nCount_RNA","percent.mt"))

St_36_clean_singles_doublets <- RunPCA(St_36_clean_singles_doublets, npcs=100)

ElbowPlot(St_36_clean_singles_doublets, ndims=100)
ndims= 100

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
all.dev<-merge(pleuro.4, y=D1_D2)


all.dev <- SCTransform(all.dev, vars.to.regress = c("nCount_RNA","Stage","percent.mt"))

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
FeaturePlot(all.dev, features = c("GFAP","PCNA"))


#note, it is possible that the actual cluster numbers will change below depending on which cells were filtered during SVM. If re-running analysis must check that clusters you pull out still make sense
markers_27<-FindMarkers(all.dev,"27")
markers_28<-FindMarkers(all.dev,"28") #CPXM1 so likely endothelial cells, FOXG1 negative 
markers_29<-FindMarkers(all.dev,"29")# MYH11= smooth muscle cells. FOXG1 negative
markers_30<-FindMarkers(all.dev,"30") #FL1 =blood. FOXG1 negative
markers_0<-FindMarkers(all.dev,"0")   #LHX1= hypothalmus? FOXG1 negative
markers_25<-FindMarkers(all.dev,"25")




neur.dev<-subset(x=all.dev, idents = c("28", "29", "30"), invert=T)
 

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
neur.dev <- CellCycleScoring(neur.dev, s.features = s.genes, g2m.features = g2m.genes)


neur.dev <- SCTransform(neur.dev, vars.to.regress = c("nCount_RNA","Stage","percent.mt", "S.Score","G2M.Score"))

neur.dev <- RunPCA(neur.dev, npcs=100)

ElbowPlot(neur.dev, ndims=100)
ndims= 100
#orig 25 dim
neur.dev <- RunUMAP(neur.dev, dims=1:25)
neur.dev <- FindNeighbors(neur.dev, dims=1:25)

neur.dev <- FindClusters(neur.dev, resolution = 2)

DimPlot(neur.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(neur.dev, group.by  ="Stage", label=T, pt.size = 0.3)
DimPlot(neur.dev, group.by  ="Phase", label=T, pt.size = 0.3)


FeaturePlot(neur.dev, features = c("GAD1","NEUROD6","TP73","TBR1"))

FeaturePlot(neur.dev, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(neur.dev, features = c("FOXG1","NEUROD6","GAD1","TBX21"))
FeaturePlot(neur.dev, features = c("LHX1","SLC17A6"))
FeaturePlot(neur.dev, features = c("GFAP","PCNA", "NEUROG2", "NEUROD6"))

saveRDS(neur.dev, file = "neur.dev.rds")
saveRDS(all.dev, file ="all.dev.rds")


### run label transfer from adult to development
library(Seurat)
library(ggplot2)
library(SeuratDisk)
library(corrplot)
library(MatrixGenerics)

query<-neur.dev
reference<-readRDS("grouped_neurons_final.rds")

query<-AddMetaData(query, (paste('D', as.character(Idents(query), sep="-"))),  col.name = "cluster_label")

reference<-AddMetaData(reference, (paste(as.character(Idents(reference)))), col.name ="cluster_label")

reference <- RunUMAP(reference, dims = 1:180, reduction = "pca", return.model = TRUE )

neur.anchors <- FindTransferAnchors(reference = reference, query = query,
                                    dims = 1:50, reference.reduction = "pca", k.anchor=10)
# k.weight default is 50; setting it to 20 in this case because otherwise it throws an error
predictions <- TransferData(anchorset = neur.anchors, refdata = reference$cluster_label,
                            dims = 1:50, k.weight=20,  weight.reduction = "pcaproject")
query <- AddMetaData(query, metadata = predictions)



high_score_cells <- rownames(query@meta.data[query@meta.data$prediction.score.max>0.90,])
p1 <- DimPlot(query, cells=high_score_cells, group.by = "predicted.id", label = T)+NoLegend()
p2 <- DimPlot(query, cells=high_score_cells, group.by="cluster_label", label = T) +NoLegend()

pdf(paste0("2-17-Neur_dev_LT_cluster_pred_score_90.pdf"), width=20, height=10)
p1+p2
dev.off()



reference <- RunUMAP(reference, dims = 1:180, reduction = "pca", return.model = TRUE)
query@tools$TransferData <- NULL
query <- MapQuery(anchorset = neur.anchors, reference = reference, query = query,
                  refdata = list(celltype = "cluster_label"), reference.reduction = "pca", reduction.model = "umap")


xmin = min(reference[["umap"]]@cell.embeddings[,1],query[["ref.umap"]]@cell.embeddings[,1])
xmax = max(reference[["umap"]]@cell.embeddings[,1],query[["ref.umap"]]@cell.embeddings[,1])
ymin = min(reference[["umap"]]@cell.embeddings[,2],query[["ref.umap"]]@cell.embeddings[,2])
ymax = max(reference[["umap"]]@cell.embeddings[,2],query[["ref.umap"]]@cell.embeddings[,2])


pal <- c("#F8766D", "#0CB702", "#00A9FF", "#C77CFF", "#ED68ED", "#ABA300", "#8494FF", "#00B8E7", "#CD9600", "#00BFC4", "#7CAE00", "#E68613", "#FF68A1", "#00C19A", "#00BE67", "#FF61CC")



p1 <- DimPlot(reference, reduction = "umap", label = TRUE, label.size = 2,
              repel = F, cols = pal) + ggtitle("Reference annotations") + xlim(xmin, xmax) + ylim(ymin, ymax)

cells <- rownames(query@meta.data[query@meta.data$prediction.score.max>0.90,])

p2 <- DimPlot(query, cells=cells, reduction = "ref.umap", group.by = "cluster_label", label = TRUE,
              label.size = 2, repel = F)  + ggtitle("Query transferred labels") + xlim(xmin, xmax) + ylim(ymin, ymax)

pdf(paste0("2-17-full_adult_LT_pred_score_90p_color_TELE_legend.pdf"), width=20, height=10)
p1
p2
dev.off()


TransferPlot<-function(reduction.rot,prediction){
        transfer.score <- prediction
        pl <- ggplot(reduction.rot, aes(x=UMAP_1, y=UMAP_2, color=transfer.score)) + geom_point(size=.75)
        med <- median(transfer.score)
        pl + scale_color_gradient(low="grey80",high="magenta", space ="Lab" , limits=c(0,1.05))	
}


#must change VP to MP, DP etc where ever it apears below to get each plot 
pdf(paste0("t_score_3VCP.pdf"), width=11.5, height=10)
cluster_to_plot = "VCP"
TransferPlot(as.data.frame(Embeddings(query, reduction="umap")), query@meta.data$prediction.score.VCP) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + ggtitle(paste("transfer score", cluster_to_plot))

dev.off()





#dev_tip<-subset(neur.dev, idents = c("10","17","19","20","14","9","22","16","11","12","24")




telen.dev<-subset(x=neur.dev, idents = c("18", "10", "27"), invert=T)



telen.dev <- SCTransform(telen.dev, vars.to.regress = c("nCount_RNA","percent.mt", "Stage", "S.Score","G2M.Score"))

telen.dev <- RunPCA(telen.dev, npcs=100)

ElbowPlot(telen.dev, ndims=100)
ndims= 100

telen.dev <- RunUMAP(telen.dev, dims=1:25)
telen.dev <- FindNeighbors(telen.dev, dims=1:25)

telen.dev <- FindClusters(telen.dev, resolution = 2)

DimPlot(telen.dev, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(telen.dev, group.by  ="Stage", label=T, pt.size = 0.3)

FeaturePlot(telen.dev, features = c("GAD1","NEUROD6","TP73","TBR1"))

FeaturePlot(telen.dev, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(telen.dev, features = c("FOXG1","NEUROD6","GAD1","TBX21"))
FeaturePlot(telen.dev, features = c("LHX1","SLC17A6"))
FeaturePlot(telen.dev, features = c("EMX1","DLX1"))

FeaturePlot(telen.dev, features = c("LHX9","NFIX","PROX1", "ETV1"))

FeaturePlot(telen.dev, features = c("RORB","FOXO6","FEZF2", "ZBTB20"))



saveRDS(telen.dev, file="telen.dev.rds")

telen.dev<-readRDS("telen.dev.rds")
telen.dev.2<-subset(x= telen.dev, idents="31", invert=T )



telen.dev.2 <- SCTransform(telen.dev.2, vars.to.regress = c("nCount_RNA","percent.mt", "Stage", "S.Score","G2M.Score"))

telen.dev.2 <- RunPCA(telen.dev.2, npcs=100)

ElbowPlot(telen.dev.2, ndims=100)
ndims= 100

telen.dev.2 <- RunUMAP(telen.dev.2, dims=1:25)
telen.dev.2 <- FindNeighbors(telen.dev.2, dims=1:25)

telen.dev.2 <- FindClusters(telen.dev.2, resolution = 2.5)

DimPlot(telen.dev.2, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(telen.dev.2, group.by  ="Stage", label=T, pt.size = 0.3)

FeaturePlot(telen.dev.2, features = c("GAD1","NEUROD6","RELN","SLC17A7"))

FeaturePlot(telen.dev.2, features = c("percent.mt","nCount_RNA","nFeature_RNA"))

FeaturePlot(telen.dev.2, features = c("FOXG1","NEUROD6","GAD1","TBX21"))
FeaturePlot(neur.dev, features = c("SIM1"))
FeaturePlot(telen.dev.2, features = c("EMX1","SLC1A3","NEUROG2","ASCL1"), label=T)

FeaturePlot(telen.dev.2, features = c("OTP","SLC17A6","SIM1", "LHX9"))
FeaturePlot(telen.dev.2, features =c("NKX2-1"))

#Take NEUROG2 postive progenitors
### subset out progenitors 

saveRDS(telen.dev.2, file="telen.dev.2.rds")

dev.prog<-subset(x=telen.dev.2, idents=c("5","7","3","0"))



dev.prog <- RunPCA(dev.prog,features=VariableFeatures(dev.prog), npcs=100)

ElbowPlot(dev.prog, ndims=100)
ndims= 100

dev.prog <- RunUMAP(dev.prog, dims=1:25)
dev.prog <- FindNeighbors(dev.prog, dims=1:25)

dev.prog <- FindClusters(dev.prog, resolution = 1.5)

DimPlot(dev.prog, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(dev.prog, group.by  ="Stage", label=T, pt.size = 0.3)

FeaturePlot(dev.prog, features = c("EMX1","POU3F1","NEUROG2", "FEZF2"), label = T)


dev.prog.pal<-subset(dev.prog, idents=c("14","1","12"))


#dev.pal<-subset(telen.dev.2, idents = c("2","17","1","16","12","14","20","8","18"))
dev.pal<-subset(telen.dev.2, idents = c("20","21","13","3","31","16","2","27","14","9","8","4","17"))



dev.traj<-merge(dev.pal, y= dev.prog.pal)





### check merged object


dev.pal <- FindNeighbors(dev.pal, dims=1:25)

dev.pal <- FindClusters(dev.pal, resolution = .5)

DimPlot(dev.pal, reduction="umap", label=T, pt.size = 0.3) + NoLegend()

DimPlot(dev.pal, group.by  ="Stage", label=T, pt.size = 0.3)



library(slingshot
)

library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)


library(scater)

library(cowplot)







#make trajectories using origianl projection pattern from telen.dev.2











plot <- DimPlot(object = dev.pal)
cells.located <- CellSelector(plot = plot)
cells.located



dev.pal.2<-subset(dev.pal, cells = cells.located ,invert= T)



DimPlot(dev.pal.2, reduction="umap", label=T, pt.size = 0.3) + NoLegend()



plot <- DimPlot(object = dev.pal.2)
cells.located <- CellSelector(plot = plot)

dev.pal.2<-subset(dev.pal.2, cells = cells.located ,invert= T)



dev.pal.2$seurat_clusters<-Idents(dev.pal.2)

dev.pal.2.sce <- as.SingleCellExperiment(dev.pal.2)

seu<-dev.pal.2.sce

umap_embed<-dev.pal.2@reductions$umap@cell.embeddings
reducedDims(dev.pal.2.sce)<-list(UMAP = umap_embed)

reducedDim(dev.pal.2.sce)




seu<-dev.pal.2.sce



sds <- slingshot(seu, reducedDim= "UMAP", clusterLabels = seu$seurat_clusters, 
                 start.clus = c("9"), end.clus = c("7","6","3"))

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


cell_colors_clust <- cell_pal(seu$seurat_clusters, hue_pal())




pdf(paste0("slingshot_lines_2-17.pdf"), width=10, height=10)
plot(reducedDim(sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(SlingshotDataSet(sds), lwd=2,  col='black')
dev.off()



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



dev.pal.2$curve3<-pt.df$curve3
dev.pal.2$curve4<-pt.df$curve4
pdf(paste0("2-18pseudotime_map_lines.pdf"), width=10, height=10)
FeaturePlot(dev.pal.2, features = "curve3", label = F, cols = pal)
FeaturePlot(dev.pal.2, features = "curve4", label = F, cols = pal)
dev.off()

dev.pal.2$curve1<-pt.df$curve1
dev.pal.2$curve2<-pt.df$curve2
pdf(paste0("2-18pseudotime_map_MT.pdf"), width=10, height=10)
FeaturePlot(dev.pal.2, features = "curve1", label = F, cols = pal)
FeaturePlot(dev.pal.2, features = "curve2", label = F, cols = pal)
dev.off()





FeaturePlot(dev.pal.2, features = c("SOX6","POU3F1","PROX1", "FEZF2"), label = T)


### start identifying differentially expressed genes for the heatmap 



pt.df<-as.data.frame(pt)


pt.df_curve3<-subset(pt.df, !pt.df$curve3=="NA")

pt.df_curve4<-subset(pt.df, !pt.df$curve4=="NA")




install.packages("gam")
library(gam)

# Only look at the 1,000 most variable genes when identifying temporally expressesd genes.
# Identify the variable genes by ranking all genes by their variance.
Y <- log2(counts(seu) + 1)
var1K <- names(sort(apply(Y, 1, var),decreasing = TRUE))[1:4000]
Y <- Y[var1K, ]  # only counts for variable genes
Y1<-Y[,rownames(pt.df_curve3)]
Y2<-Y[,rownames(pt.df_curve4)]

FeaturePlot(dev.pal.2, features = c("PROX1","SOX6","RELN","FEZF2"))


markers6v7<-FindMarkers(dev.pal.2, ident.1= "6", ident.2 ="7")
markers7v6<-FindMarkers(dev.pal.2, ident.1= "7", ident.2 ="6")
markers0v5<-FindMarkers(dev.pal.2, ident.1= "0", ident.2 ="5")
markers5v0<-FindMarkers(dev.pal.2, ident.1= "5", ident.2 ="0")


markers6v7$pct.dif<-(markers6v7$pct.1- markers6v7$pct.2)
markers6v7$pct.ratio<-(markers6v7$pct.1/markers6v7$pct.2)
markers6v7$gene<-rownames(markers6v7)
top_markers6v7<-subset(markers6v7, (markers6v7$pct.dif>.21 & markers6v7$pct.ratio>3))

markers7v6$pct.dif<-(markers7v6$pct.1 - markers7v6$pct.2)
markers7v6$pct.ratio<-(markers7v6$pct.1/markers7v6$pct.2)
markers7v6$gene<-rownames(markers7v6)
top_markers7v6<-subset(markers7v6, (markers7v6$pct.dif>.21 & markers7v6$pct.ratio>3))


markers0v5$pct.dif<-(markers0v5$pct.1- markers0v5$pct.2)
markers0v5$pct.ratio<-(markers0v5$pct.1/markers0v5$pct.2)
markers0v5$gene<-rownames(markers0v5)
top_markers0v5<-subset(markers0v5, (markers0v5$pct.dif>.21 & markers0v5$pct.ratio>3))

markers5v0$pct.dif<-(markers5v0$pct.1- markers5v0$pct.2)
markers5v0$pct.ratio<-(markers5v0$pct.1/markers5v0$pct.2)
markers5v0$gene<-rownames(markers5v0)
top_markers5v0<-subset(markers5v0, (markers5v0$pct.dif>.21 & markers5v0$pct.ratio>3))


TFs_mouse <- read.table("~/Toshces_lab/Tosches_r/TFs_mouse.txt", quote="\"", comment.char="")

top_markers_venlat <- data.frame(gene = c(top_markers6v7[,"gene"], top_markers0v5[,"gene"]))

top_markers_dormed <- data.frame(gene = c(top_markers7v6[,"gene"], top_markers5v0[,"gene"]))

top_markers_dormed_unique<-unique(top_markers_dormed$gene)

top_markers_venlat_unique<-unique(top_markers_venlat$gene)

top_markers_dormed_unique_TF<-subset(top_markers_dormed_unique, top_markers_dormed_unique%in%TFs_mouse$V1)
top_markers_venlat_unique_TF<-subset(top_markers_venlat_unique, top_markers_venlat_unique%in%TFs_mouse$V1)


t1 <- pt.df_curve3$curve3
t2 <- pt.df_curve4$curve4




pt.df_curve3_unique<-subset(pt.df_curve3, !rownames(pt.df_curve3) %in% rownames(pt.df_curve4))

pt.df_curve4_unique<-subset(pt.df_curve4, !rownames(pt.df_curve4) %in% rownames(pt.df_curve3))




pt.df_curve3 <- pt.df_curve3[order(pt.df_curve3[, "curve3"]), , drop = FALSE]
curve3_order<-rownames(pt.df_curve3)


D1_D2_curve3<-subset(dev.pal.2, cells= curve3_order)

df_curve3<-as.data.frame(D1_D2_curve3@assays$RNA@data)

df_curve3 <- df_curve3[, curve3_order]

max_cell_curve3<-data.frame(rownames(df_curve3))

max_cell_curve3$max_cell<-max.col(df_curve3)

max_cell_curve3_sub<-subset(max_cell_curve3, max_cell_curve3$rownames.df_curve3. %in% top_markers_dormed_unique)

max_cell_curve3_sub_order <- max_cell_curve3_sub[order(max_cell_curve3_sub[, "max_cell"]), , drop = FALSE]





curve4_cells<-rownames(pt.df_curve4)

dev_curve4<-subset(dev.pal.2, cells= curve4_cells)

pt.df_curve4 <- pt.df_curve4[order(pt.df_curve4[, "curve4"]), , drop = FALSE]
curve4_order<-rownames(pt.df_curve4)

df_curve4<-as.data.frame(dev_curve4@assays$RNA@data)

df_curve4 <- df_curve4[, curve4_order]

max_cell_curve4<-data.frame(rownames(df_curve4))

max_cell_curve4$max_cell<-max.col(df_curve4)

max_cell_curve4_sub<-subset(max_cell_curve4, max_cell_curve4$rownames.df_curve4. %in% top_markers_venlat_unique)

max_cell_curve4_sub_order <- max_cell_curve4_sub[order(max_cell_curve4_sub[, "max_cell"]), , drop = FALSE]






View(pt.df)
curve3_curve4_intesect<-subset(pt.df, ((!pt.df$curve3=="NA")&(!pt.df$curve4=="NA")))

set.seed(42)
rand <- sample(nrow(curve3_curve4_intesect))

rand_curve3_curve4_intersect<-curve3_curve4_intesect[rand,]

View(rand_curve3_curve4_intersect)
rand_intesect_curve3<-rand_curve3_curve4_intersect[1:423,]
rand_intesect_curve4<-rand_curve3_curve4_intersect[424:846,]

rand_intesect_curve3<-data.frame("cell"=rownames(rand_intesect_curve3), "pt"= rand_intesect_curve3$curve3)

rand_intesect_curve4<-data.frame("cell"=rownames(rand_intesect_curve4), "pt"= rand_intesect_curve4$curve4)

rownames(rand_intesect_curve3)<-rand_intesect_curve3$cell


rand_intersect_neg_curve3<-rand_intesect_curve3$pt*(-1)
rand_intesect_curve3$pt<-rand_intersect_neg_curve3


pt.df_curve3_unique_neg<-pt.df_curve3_unique*(-1)
curve3<- data.frame("cell"=rownames(pt.df_curve3_unique_neg), "pt"=pt.df_curve3_unique_neg$curve3)

curve4<- data.frame("cell"=rownames(pt.df_curve4_unique), "pt"=pt.df_curve4_unique$curve4)

curve3<-rbind(curve3, rand_intesect_curve3)
curve4<-rbind(curve4, rand_intesect_curve4)


curve_merge<-rbind(curve3, curve4)


curve_merge_order <- curve_merge[order(curve_merge[, "pt"]), , drop = FALSE]

rownames(curve_merge_order)<-curve_merge_order$cell

merge_t<-curve_merge_order$pt


dev_curve_merge<-subset(dev.pal.2, cells= curve_merge_order$cell)
View(max_cell_curve3_sub_order$rownames.df_curve3.)

library(circlize)
col_fun = colorRamp2(c(-2,0, 6), c("blue", "white", "red"))
col_fun(seq(-3, 3))



heatdata_merge <- as.matrix(dev_curve_merge@assays$RNA@data[rownames(dev_curve_merge@assays$RNA@data) %in% top_markers_dormed_unique, order(merge_t, na.last = NA)])

heatdata_merge_order<- heatdata_merge[max_cell_curve3_sub_order$rownames.df_curve3.,]

heatdata_merge_order_scale = t(scale(t(heatdata_merge_order)))

View(top_markers_dormed_unique_TF)
library("ComplexHeatmap")
library("zoo")

ha = rowAnnotation(foo = anno_mark(at = c(40, 16, 34, 17, 4, 37, 29, 22, 19, 23, 50, 36), labels = top_markers_dormed_unique_TF))


pdf(paste0("2-17heatmap_pseudotime_merge_dormed_v8.pdf"), width=10, height=10)
Heatmap(
        heatdata_merge_order_scale,
        cluster_rows=FALSE, cluster_columns=FALSE, show_row_names = FALSE,
      show_column_names=FALSE, heatmap_height = unit(20, "cm"), right_annotation = ha, 
      column_split =  rep(c("DM Trajectory","VL Trajectory"),c(1048,1189)),
        name = "heatmap")
dev.off()




View(max_cell_curve4_sub_order$rownames.df_curve4.)
View(top_markers_venlat_unique_TF)

heatdata_VL_merge <- as.matrix(dev_curve_merge@assays$RNA@data[rownames(dev_curve_merge@assays$RNA@data) %in% top_markers_venlat_unique, order(merge_t, na.last = NA)])

heatdata_VL_merge_order<- heatdata_VL_merge[max_cell_curve4_sub_order$rownames.df_curve4.,]

heatdata_VL_merge_order_scale = t(scale(t(heatdata_VL_merge_order)))


hv = rowAnnotation(foo = anno_mark(at = c(17, 37, 16,  10, 11, 8, 6), labels = top_markers_venlat_unique_TF))

pdf(paste0("2-17heatmap_pseudotime_merge_Venlat.pdf"), width=10, height=10)
Heatmap(
        heatdata_VL_merge_order_scale,
        cluster_rows=FALSE, cluster_columns=FALSE, show_row_names = FALSE,
        right_annotation = hv, show_column_names=FALSE, heatmap_height = unit(20, "cm"), 
        name = "heatmap", column_split =  rep(c("DM Trajectory","VL Trajectory"),c(1048,1189)))
dev.off()







