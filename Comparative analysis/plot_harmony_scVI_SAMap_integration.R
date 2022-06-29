# alternative integration methods

library(Seurat)
 
# load integration object
int <- readRDS("/Users/maria/Documents/Lab_Columbia/Pleurodeles_scRNAseq/Comparative2/Final_3sp_integration.rds")


####################################
########### Harmony ################
####################################

library(harmony)
int <- ScaleData(int)
int <- RunHarmony(int, group.by.vars = "species", assay.use = "SCT", project.dim = F)
# careful because these commands overwrite existing ones
int <- FindNeighbors(int, reduction="harmony", dims=1:80)
# int <- FindClusters(int, res=2.1)
int <- RunUMAP(int, reduction="harmony", dims=1:80)


# save results
harmony.embeddings <- Embeddings(int, reduction="harmony")
harmony.umap <- Embeddings(int, reduction="umap")
saveRDS(list(harmony.embeddings, harmony.umap), file="Final_3sp_harmony.rds")
# saveRDS(list(harmony.embeddings, harmony.umap), file="Final_4sp_harmony.rds")

# save seurat object
saveRDS(int, file="Final_3sp_integration_harmony.rds")

# to re-attach Harmony results to seurat object:
# harmony.data <- readRDS("Final_3sp_harmony.rds")
# int[["harmony"]] <- CreateDimReducObject(embeddings = harmony.data[[2]], key = "harmony_", assay = "RNA")


# plots 

pdf("Final_3sp_harmony_UMAP.pdf")
# species_cols <- c("#EECC66", "#6699CC","#EE99AA","#FFA500")
species_cols <- c("#EECC66", "#EE99AA", "#6699CC")
DimPlot(int, reduction="harmony", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle=T, pt.size=0.5)
DimPlot(int, reduction="harmony", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle = T, pt.size=0.5) + NoLegend()
Idents(int) <- int@meta.data$id_by_cluster_label_merge_sp
mycols.5 <- c("#AA4499", "#004488", "#B2182B", "#2166AC", "#009988", "#762A83", "#9E8F32", 'darkgrey', "#1B7837", rep('darkgrey', 7), "#9E8F32", "#1B7837", "#2166AC", "#004488", "#B2182B", "#009988", "#762A83", "#6773B5", "#EA5F5B", "#458AC9", "#8666A9", "#CCBB44", "#FF6699", "#33CC99", 'darkgrey', "#2DA736", rep('darkgrey', 3))
DimPlot(int, reduction="harmony", cols = adjustcolor(mycols.5, alpha=0.7), shuffle=T, pt.size=0.5)
DimPlot(int, reduction="harmony", cols = adjustcolor(mycols.5, alpha=0.7), shuffle = T, pt.size=0.5) + NoLegend()
DimPlot(int, reduction="harmony", group.by="tosches_annot", label.size = 2, shuffle = T, label=T, pt.size=0.5) + NoLegend()
DimPlot(int, reduction="harmony", group.by="unassigned", label.size = 2, shuffle = T, label=T, pt.size=0.5) + NoLegend()
DimPlot(int, reduction="harmony", group.by="cluster_labels", label.size = 2, shuffle = T, label = T, pt.size=0.5) + NoLegend()
dev.off()

# 4 species, plot once the harmony-based UMAP is attached to the Seurat object
pdf("Final_4sp_Harmony_UMAP.pdf")
species_cols <- c("#EECC66", "#6699CC","#EE99AA","#FFA500")
mycols18 <- c("#00BA38", "#FF61C3", "#F8766D", "#FF699C", "#AE87FF", "#00ADFA", "#F564E3", "#D39200", "#00B9E3","#B79F00", "#00BF74", "#E88526", "#00C19F", "#DB72FB", "#93AA00", "#5EB300", "#619CFF", "#00BFC4")
mycols3 <- c("#F8766D", "#00BA38", "#619CFF")
DimPlot(int, reduction="harmony", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle=T, pt.size=0.3)
DimPlot(int, reduction="harmony", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle = T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", label=T, label.size = 2, shuffle=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", group.by="tosches_annot_final", label.size = 2, shuffle = T, label=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", group.by="unassigned", label.size = 2, shuffle = T, label=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", group.by="superclusters", label.size = 2, shuffle = T, label=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", group.by="cluster_labels", label.size = 2, shuffle = T, label = T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", group.by="subclass_label", label.size = 2, shuffle = T, label = T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="harmony", group.by="high_anno", label.size = 2, shuffle = T, label = F, pt.size=0.3, cols=mycols3, na.value = "grey85") 
DimPlot(int, reduction="harmony", group.by="high_anno", label.size = 2, shuffle = T, label = T, pt.size=0.3, cols=mycols3, na.value = "grey85") + NoLegend()
DimPlot(int, reduction="harmony", group.by="region_label", label.size = 2, shuffle = T, label = T, pt.size=0.3, cols=mycols18, na.value = "grey85") 
DimPlot(int, reduction="harmony", group.by="region_label", label.size = 2, shuffle = T, label = T, pt.size=0.3,cols=mycols18, na.value = "grey85") + NoLegend()
DimPlot(int, reduction="harmony", group.by="region_label", label.size = 2, shuffle = T, label = F, pt.size=0.3,cols=mycols18, na.value = "grey85") + NoLegend()
dev.off()


####################################
############# scVI #################
####################################

# https://docs.scvi-tools.org/en/0.14.1/tutorials/notebooks/scvi_in_R.html
setwd("~/Documents/Lab_Columbia/Pleurodeles_scRNAseq/Comparative2/scVI")

sc3k <- readRDS("scvi-tools_3Sp_2960genes_v1.rds")
scvi.umap <- Embeddings(sc3k, reduction = "umap")

sc4k <- readRDS("4Sp_scvi_v1.rds")
scvi.umap <- Embeddings(sc4k, reduction = "umap")

int[["scVI"]] <- CreateDimReducObject(embeddings = scvi.umap, key="scvi_", assay="RNA")

# plots
pdf("Final_3sp_scVI_UMAP.pdf")
species_cols <- c("#EECC66", "#EE99AA", "#6699CC")
DimPlot(int, reduction="scVI", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle=T, pt.size=0.5)
DimPlot(int, reduction="scVI", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle = T, pt.size=0.5) + NoLegend()
mycols.5 <- c("#AA4499", "#004488", "#B2182B", "#2166AC", "#009988", "#762A83", "#9E8F32", 'darkgrey', "#1B7837", rep('darkgrey', 7), "#9E8F32", "#1B7837", "#2166AC", "#004488", "#B2182B", "#009988", "#762A83", "#6773B5", "#EA5F5B", "#458AC9", "#8666A9", "#CCBB44", "#FF6699", "#33CC99", 'darkgrey', "#2DA736", rep('darkgrey', 3))
DimPlot(int, reduction="scVI", cols = adjustcolor(mycols.5, alpha=0.7), shuffle=T, pt.size=0.5)
DimPlot(int, reduction="scVI", cols = adjustcolor(mycols.5, alpha=0.7), shuffle = T, pt.size=0.5) + NoLegend()
DimPlot(int, reduction="scVI", group.by="tosches_annot", label.size = 2, shuffle = T, label=T, pt.size=0.5) + NoLegend()
DimPlot(int, reduction="scVI", group.by="unassigned", label.size = 2, shuffle = T, label=T, pt.size=0.5) + NoLegend()
DimPlot(int, reduction="scVI", group.by="cluster_labels", label.size = 2, shuffle = T, label = T, pt.size=0.5) + NoLegend()
dev.off()

pdf("Final_4sp_scVI_UMAP.pdf")
species_cols <- c("#EECC66", "#6699CC","#EE99AA","#FFA500")
mycols18 <- c("#00BA38", "#FF61C3", "#F8766D", "#FF699C", "#AE87FF", "#00ADFA", "#F564E3", "#D39200", "#00B9E3","#B79F00", "#00BF74", "#E88526", "#00C19F", "#DB72FB", "#93AA00", "#5EB300", "#619CFF", "#00BFC4")
mycols3 <- c("#F8766D", "#00BA38", "#619CFF")
DimPlot(int, reduction="scVI", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle=T, pt.size=0.3)
DimPlot(int, reduction="scVI", group.by = "species", cols = adjustcolor(species_cols, alpha=0.7), shuffle = T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", label=T, label.size = 2, shuffle=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", group.by="tosches_annot_final", label.size = 2, shuffle = T, label=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", group.by="unassigned", label.size = 2, shuffle = T, label=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", group.by="superclusters", label.size = 2, shuffle = T, label=T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", group.by="cluster_labels", label.size = 2, shuffle = T, label = T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", group.by="subclass_label", label.size = 2, shuffle = T, label = T, pt.size=0.3) + NoLegend()
DimPlot(int, reduction="scVI", group.by="high_anno", label.size = 2, shuffle = T, label = F, pt.size=0.3, cols=mycols3, na.value = "grey85") 
DimPlot(int, reduction="scVI", group.by="high_anno", label.size = 2, shuffle = T, label = T, pt.size=0.3, cols=mycols3, na.value = "grey85") + NoLegend()
DimPlot(int, reduction="scVI", group.by="region_label", label.size = 2, shuffle = T, label = T, pt.size=0.3, cols=mycols18, na.value = "grey85") 
DimPlot(int, reduction="scVI", group.by="region_label", label.size = 2, shuffle = T, label = T, pt.size=0.3,cols=mycols18, na.value = "grey85") + NoLegend()
DimPlot(int, reduction="scVI", group.by="region_label", label.size = 2, shuffle = T, label = F, pt.size=0.3,cols=mycols18, na.value = "grey85") + NoLegend()
dev.off()


####################################
############# SAMap ################
####################################

# load integration object
int <- readRDS("/Users/maria/Documents/Lab_Columbia/Pleurodeles_scRNAseq/Comparative2/Final_3sp_integration.rds")

# SAMap UMAP data
pw.umap <- read.csv("3sp_SAMap_UMAP_pw.csv", row.names = 1, header = T)
pv.umap <- read.csv("3sp_SAMap_UMAP_pv.csv", row.names = 1, header = T)
cp.umap <- read.csv("3sp_SAMap_UMAP_cp.csv", row.names = 1, header = T)
 
samap.umap <- rbind(pw.umap, pv.umap, cp.umap)
samap.umap <- samap.umap[order(match(rownames(samap.umap), Cells(int))),]
samap.umap <- as.matrix(samap.umap)
 
all.equal(Cells(int), rownames(samap.umap))
int[["SAMap"]] <- CreateDimReducObject(embeddings = samap.umap, key = "SAMap_", assay = DefaultAssay(int))

pdf("Final_3sp_SAMap_UMAP.pdf")
species_cols <- c("#EECC66", "#EE99AA", "#6699CC")
DimPlot(int, reduction = "SAMap", group.by = "species", shuffle = T, cols=adjustcolor(species_cols, alpha=0.7), pt.size = 0.5)
DimPlot(int, reduction = "SAMap", group.by = "species", shuffle = T, cols=adjustcolor(species_cols, alpha=0.7), pt.size = 0.5) + NoLegend()
Idents(int) <- int@meta.data$id_by_cluster_label_merge_sp
mycols.5 <- c("#AA4499", "#004488", "#B2182B", "#2166AC", "#009988", "#762A83", "#9E8F32", 'darkgrey', "#1B7837", rep('darkgrey', 7), "#9E8F32", "#1B7837", "#2166AC", "#004488", "#B2182B", "#009988", "#762A83", "#6773B5", "#EA5F5B", "#458AC9", "#8666A9", "#CCBB44", "#FF6699", "#33CC99", 'darkgrey', "#2DA736", rep('darkgrey', 3))
DimPlot(int, reduction="SAMap", cols = adjustcolor(mycols.5, alpha=0.7), shuffle=T, pt.size = 0.5)
DimPlot(int, reduction="SAMap", cols = adjustcolor(mycols.5, alpha=0.7), shuffle = T, pt.size = 0.5) + NoLegend()
DimPlot(int, reduction="SAMap", group.by="tosches_annot", label.size = 2, shuffle = T, label=T, pt.size = 0.5) + NoLegend()
DimPlot(int, reduction="SAMap", group.by="cluster_labels", label.size = 2, shuffle = T, label = T, pt.size = 0.5) + NoLegend()
DimPlot(int, reduction="SAMap", group.by="unassigned", label.size = 2, shuffle = T, label = T, pt.size = 0.5) + NoLegend()
dev.off()
