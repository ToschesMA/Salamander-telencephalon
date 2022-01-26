# runs on the server

########################################################################
# Setup

.libPaths("/burg/home/mt3353/R/")

library(Seurat, lib.loc="/burg/home/mt3353/R")
library(corrplot, lib.loc="/burg/home/mt3353/R")
library(MatrixGenerics, lib.loc="/burg/home/mt3353/R")
library(ggplot2, lib.loc="/burg/home/mt3353/R")
library(dplyr, lib.loc="/burg/home/mt3353/R")

setwd("/burg/tosches/users/mt3353/comparative/")

# Import Seurat objects
mouse.z0 <- readRDS("zeisel_telencephalic_neurons.RDS")
mouse.y0 <- readRDS("mouse_yao_10x.rds")
pleuro.0 <- readRDS("neurons_fina2_corrected.rds")

metadata <- read.csv("zeisel_metadata.csv", header=F)

# import eggnog files
eggnog_mouse <- read.table("mus_eggnog_pruned.txt" , header=T)
eggnog_salamander <- read.table("pleurodeles_eggnog_pruned.txt" , header=T)

# setup cluster identities for analysis
z_glut <- c(grep("TEGLU",levels(Idents(mouse.z0)), value = TRUE),grep("DGGR", levels(Idents(mouse.z0)), value = TRUE))
y_glut <- c("CA1-ProS","CA2-IG-FC","CA3","Car3","CT SUB","DG","L2 IT ENTl","L2 IT ENTm","L2/3 IT CTX","L2/3 IT ENTl","L2/3 IT PPP","L2/3 IT RHP","L3 IT ENT","L4 RSP-ACA","L4/5 IT CTX","L5 IT CTX","L5 PPP","L5 PT CTX","L5/6 IT TPE-ENT","L5/6 NP CTX"," L6 CT CTX","L6 IT CTX","L6 IT ENTl","L6b CTX","L6b/CT ENT","NP PPP","NP SUB","SUB-ProS")
s_glut <- grep("TEGLU",levels(Idents(pleuro.0)), value = TRUE)

z_gaba <- c(grep("TEINH", levels(Idents(mouse.z0)), value = TRUE),grep("MSN", levels(Idents(mouse.z0)), value = TRUE), grep("OBINH", levels(Idents(mouse.z0)), value = TRUE),"OBDOP","OBDOP2", "TECHO" )
y_gaba <- c("Lamp5","Meis2","Pvalb","Sncg","Sst","Sst Chodl", "Vip")
s_gaba <- grep("TEGABA", levels(Idents(pleuro.0)), value = TRUE)

# variables for saving output files
date = "I7_220117_CCA_"
Sz = "Z"
Sy = "Y"
Ss = "S"

z_idents_keep = levels(Idents(mouse.z0))[levels(Idents(mouse.z0) %in% z_glut)]
y_idents_keep = levels(Idents(mouse.y0))[levels(Idents(mouse.y0) %in% y_glut)]
s_idents_keep = levels(Idents(pleuro.0))[levels(Idents(pleuro.0) %in% s_glut)]
vars_z = c("nFeature_RNA","nCount_RNA","Age","percent.mt","SampleID")
vars_y = c("nFeature_RNA","nCount_RNA","external_donor_name_label")
vars_s = c("nFeature_RNA","nCount_RNA","animal","percent.mt")

# number of dimensions to use for IntegrateData
nDims = 80

# number of variable genes to use to identify integration anchors
ngenes = 2000

########################################################################
### Run integration ###

# Bring the objects in the same gene space

mouse.z0 <- subset(mouse.z0, idents = z_idents_keep)
mouse.y0 <- subset(mouse.y0, idents = y_idents_keep)
pleuro.0 <- subset(pleuro.0, idents = s_idents_keep)

z_data <- GetAssayData(mouse.z0, slot="counts")
z_data <- z_data[rownames(z_data) %in% eggnog_mouse[,1],]

y_data <- GetAssayData(mouse.y0, slot="counts") 
y_data <- y_data[rownames(y_data) %in% eggnog_mouse[,1],]

s_data <- GetAssayData(pleuro.0, slot="counts")
s_data <- s_data[rownames(s_data) %in% eggnog_salamander[,1],]

eggnog_z <- eggnog_mouse[eggnog_mouse[,1] %in% rownames(z_data),]
eggnog_y <- eggnog_mouse[eggnog_mouse[,1] %in% rownames(y_data),]
eggnog_s <- eggnog_salamander[eggnog_salamander[,1] %in% rownames(s_data),]

z_data <- z_data[order(match(rownames(z_data), eggnog_z[,1])),]
y_data <- y_data[order(match(rownames(y_data), eggnog_y[,1])),]
s_data <- s_data[order(match(rownames(s_data), eggnog_s[,1])),]

all.equal(rownames(s_data),as.character(eggnog_s[,1]))
all.equal(rownames(y_data),as.character(eggnog_y[,1]))
all.equal(rownames(z_data),as.character(eggnog_z[,1]))

# replace gene names with eggnog names in expression tables (only one-to-ones)
rownames(z_data) <- eggnog_z$eggnog 
rownames(y_data) <- eggnog_y$eggnog 
rownames(s_data) <- eggnog_s$eggnog

dim(z_data)
dim(y_data)
dim(s_data)

z_data[1:4,1:4]
y_data[1:4,1:4]
s_data[1:4,1:4]

# Create metadata column with species information
mouse.z0[["cluster_label"]] <- Idents(object = mouse.z0)
mouse.y0[["cluster_label"]] <- Idents(object = mouse.y0)
pleuro.0[["cluster_label"]] <- Idents(object = pleuro.0)

meta.seurat <- data.frame(
  cluster_label = c(paste(Sz, as.character(Idents(mouse.z0)), sep="-"), paste(Sy,as.character(Idents(mouse.y0)), sep="-"),paste(Ss,as.character(Idents(pleuro.0)), sep="-")),
  species = c(rep("mouse", length(Idents(mouse.z0))),rep("mouse", length(Idents(mouse.y0))),rep("salamander", length(Idents(pleuro.0)))),
  row.names=c(colnames(z_data),colnames(y_data),colnames(s_data)))


# Running the integration - with CCA

comb.list <- vector("list",3)
comb.list[[1]] <- CreateSeuratObject(counts=z_data, meta.data = meta.seurat[colnames(z_data),])
comb.list[[1]] <- AddMetaData(comb.list[[1]], mouse.z0@meta.data)
Idents(comb.list[[1]]) <- comb.list[[1]]$cluster_label
comb.list[[1]] <- SCTransform(comb.list[[1]], vars.to.regress = vars_z)

comb.list[[2]] <- CreateSeuratObject(counts=y_data, meta.data = meta.seurat[colnames(y_data),])
comb.list[[2]] <- AddMetaData(comb.list[[2]], mouse.y0@meta.data)
Idents(comb.list[[2]]) <- comb.list[[2]]$cluster_label
comb.list[[2]] <- SCTransform(comb.list[[2]], vars.to.regress = vars_y)

comb.list[[3]] <- CreateSeuratObject(counts=s_data, meta.data = meta.seurat[colnames(s_data),])
comb.list[[3]] <- AddMetaData(comb.list[[3]], pleuro.0@meta.data)
Idents(comb.list[[3]]) <- comb.list[[3]]$cluster_label
comb.list[[3]] <- SCTransform(comb.list[[3]], vars.to.regress = vars_s)

options(future.globals.maxSize = 60000 * 1024^2)
comb.features <- SelectIntegrationFeatures(object.list = comb.list, nfeatures = ngenes)
comb.list.prep <- PrepSCTIntegration(object.list = comb.list, anchor.features = comb.features, verbose = FALSE)

comb.anchors <- FindIntegrationAnchors(object.list = comb.list.prep, normalization.method = "SCT", anchor.features = comb.features, verbose = FALSE, reduction="cca")
comb.integrated <- IntegrateData(anchorset = comb.anchors, normalization.method = "SCT", verbose = FALSE, dims=1:nDims)

comb.integrated <- RunPCA(comb.integrated, verbose = FALSE)
comb.integrated <- RunUMAP(comb.integrated, dims = 1:nDims)
comb.integrated <- FindNeighbors(comb.integrated, dims = 1:nDims)
comb.integrated <- FindClusters(comb.integrated)


# Plots
pdf(paste0(date, "_integration_zeisel_yao_salamander","_Dims",nDims,"_features",ngenes,".pdf"))

DimPlot(comb.integrated, group.by="species", raster=F) + ggtitle("Species")
DimPlot(comb.integrated, split.by="species", label = T, raster=F) + ggtitle("Species")
DimPlot(comb.integrated, label = T, raster=F) + NoLegend() + ggtitle("Integrated clusters")
DimPlot(comb.integrated, group.by="cluster_label", label=T, label.size = 1, raster=F) + NoLegend()

dev.off()

# Saving results
saveRDS(comb.integrated, file=paste0(date, "_integration_zeisel_yao_salamander","_Dims",nDims,"_features_",ngenes,".RDS"))
