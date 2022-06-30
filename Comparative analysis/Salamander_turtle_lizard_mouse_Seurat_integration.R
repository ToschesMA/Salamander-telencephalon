.libPaths("/burg/home/ao2721/R/")
library(Seurat)
library(corrplot)
library(MatrixGenerics)
library(ggplot2)
library(dplyr)
library(gplots)
library(Hmisc)
library(stringr)

date = "4Sp_v4_220523_CCA_"

# number of dimensions to use for IntegrateData
nDims = 80

# variable genes to be extracted from each SCTransformed object
ngenes = 5000

# Reading all the individual files that we'll integrate
lizard.0 <- readRDS("/burg/tosches/users/ao2721/integrated/lizard_clean_integration1000.RDS")
turtle.0 <- readRDS("/burg/tosches/users/ao2721/integrated/170528_turtle_neur.RDS")
pleuro.0 <- readRDS("/burg/tosches/users/ao2721/integrated/neurons_final5.rds")
mouse.y0 <- readRDS("/burg/tosches/users/ao2721/integrated/mouse_yao_10x_downsampled3.rds")

# Reading the eggnong files for each species to select one to one orthologues in the downstream analysis
eggnog_turtle <- read.table("/burg/tosches/users/ao2721/integrated/chrysemys_eggnog_pruned.txt" , header=T)
eggnog_lizard <- read.table("/burg/tosches/users/ao2721/integrated/pogona_eggnog_pruned.txt" , header=T)
eggnog_salamander <- read.table("/burg/tosches/users/ao2721/integrated/pleurodeles_eggnog_pruned.txt" , header=T)
eggnog_mouse <- read.table("/burg/tosches/users/ao2721/integrated/mus_eggnog_pruned.txt" , header=T)

noisy.genes <- readRDS("/burg/tosches/users/ao2721/integrated/noisy_list_12.rds") 

# grabbing the cells from each species to integrate and compare
liz_cpal <- c(grep("MCtx", levels(Idents(lizard.0)), value = TRUE),
              grep("DCtx", levels(Idents(lizard.0)), value = TRUE),
              grep("MGE", levels(Idents(lizard.0)), value = TRUE),
              grep("CGE", levels(Idents(lizard.0)), value = TRUE))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
tur_cpal <- c("e07", "e08", "e13", "e14", "e15", "e16", "e24", "e26", "e27",
              "e28", "e29", "e30", "e31", "e32", "e33", "e34", "e35", "e36", 
              "e37", "e38", "i07", "i08", "i09", "i10", "i11", "i12", "i13", 
              "i14", "i15", "i16", "i17", "i18")

ple_cpal <- c("TEGLU1", "TEGLU2", "TEGLU4", "TEGLU5", "TEGLU3", "TEGLU6", 
              "TEGLU7", "TEGLU11", "TEGLU9", "TEGLU10", "TEGLU8", "TEGLU13", "TEGLU20",
              "TEGABA25", "TEGABA26", "TEGABA24", "TEGABA23", "TEGABA27", "TEGABA22", 
              "TEGABA21", "TEGABA20", "TEGABA18", "TEGABA19")

y_glut <- c("CA1-ProS","CA2-IG-FC","CA3","Car3","CT SUB","DG","L2 IT ENTl","L2 IT ENTm",
            "L2/3 IT CTX","L2/3 IT ENTl","L2/3 IT PPP","L2/3 IT RHP","L3 IT ENT","L4 RSP-ACA",
            "L4/5 IT CTX","L5 IT CTX","L5 PPP","L5 PT CTX","L5/6 IT TPE-ENT","L5/6 NP CTX",
            "L6 CT CTX","L6 IT CTX","L6 IT ENTl","L6b CTX","L6b/CT ENT","NP PPP","NP SUB","SUB-ProS")
y_gaba <- c("Lamp5","Pvalb","Sncg","Sst","Sst Chodl", "Vip")
y_idents_keep <- append(y_glut, values = y_gaba)

Sl = "L"
St = "T"
Ss = "S"
Sy = "M"

# variables to regress
vars_l = c("percent.mito", "animal")
vars_t = c("animalident", "percent.mito")
vars_s = c("animal", "percent.mt")
vars_y = c("external_donor_name_label")

# subsetting the objects to include only the comparable cells
lizard.0 <- subset(lizard.0, idents = liz_cpal)
turtle.0 <- subset(turtle.0, idents = tur_cpal)
pleuro.0 <- subset(pleuro.0, idents = ple_cpal)
mouse.y0 <- subset(mouse.y0, idents = y_idents_keep)

# getting raw counts from the original objects and selecting one to one orthologues
l_data <- GetAssayData(lizard.0, slot="counts")
l_data <- l_data[rownames(l_data) %in% eggnog_lizard[,1],]

t_data <- GetAssayData(turtle.0, slot="counts") 
t_data <- t_data[rownames(t_data) %in% eggnog_turtle[,1],]

s_data <- GetAssayData(pleuro.0, slot="counts")
s_data <- s_data[rownames(s_data) %in% eggnog_salamander[,1],]

y_data <- GetAssayData(mouse.y0, slot="counts") 
y_data <- y_data[rownames(y_data) %in% eggnog_mouse[,1],]

eggnog_l <- eggnog_lizard[eggnog_lizard[,1] %in% rownames(l_data),]
eggnog_t <- eggnog_turtle[eggnog_turtle[,1] %in% rownames(t_data),]
eggnog_s <- eggnog_salamander[eggnog_salamander[,1] %in% rownames(s_data),]
eggnog_y <- eggnog_mouse[eggnog_mouse[,1] %in% rownames(y_data),]

l_data <- l_data[order(match(rownames(l_data), eggnog_l[,1])),]
t_data <- t_data[order(match(rownames(t_data), eggnog_t[,1])),]
s_data <- s_data[order(match(rownames(s_data), eggnog_s[,1])),]
y_data <- y_data[order(match(rownames(y_data), eggnog_y[,1])),]

all.equal(rownames(s_data),as.character(eggnog_s[,1]))
all.equal(rownames(t_data),as.character(eggnog_t[,1]))
all.equal(rownames(l_data),as.character(eggnog_l[,1]))
all.equal(rownames(y_data),as.character(eggnog_y[,1]))

rownames(l_data) <- eggnog_l$eggnog 
rownames(t_data) <- eggnog_t$eggnog 
rownames(s_data) <- eggnog_s$eggnog
rownames(y_data) <- eggnog_y$eggnog 

dim(l_data)
dim(t_data)
dim(s_data)
dim(y_data)

common.genes <- intersect(rownames(l_data), rownames(t_data))
common.genes <- intersect(common.genes, rownames(s_data))
common.genes <- intersect(common.genes, rownames(y_data))

l_data[1:4,1:4]
t_data[1:4,1:4]
s_data[1:4,1:4]
y_data[1:4,1:4]

# preparing metadata for integrated object
lizard.0[["cluster_label"]] <- Idents(object = lizard.0)
turtle.0[["cluster_label"]] <- Idents(object = turtle.0)
pleuro.0[["cluster_label"]] <- Idents(object = pleuro.0)
mouse.y0[["cluster_label"]] <- Idents(object = mouse.y0)

meta.seurat <- data.frame(
  ident_track = c(paste(Sl, as.character(Idents(lizard.0)), sep="-"), 
                  paste(St, as.character(Idents(turtle.0)), sep="-"), 
                  paste(Ss, as.character(Idents(pleuro.0)), sep="-"),
                  paste(Sy, as.character(Idents(mouse.y0)), sep = "-")),
  species = c(rep("lizard", length(Idents(lizard.0))),
              rep("turtle", length(Idents(turtle.0))), 
              rep("salamander", length(Idents(pleuro.0))),
              rep("mouse", length(Idents(mouse.y0)))),
  row.names=c(colnames(l_data),
              colnames(t_data),
              colnames(s_data),
              colnames(y_data)))

# we create a list that will contain all the three objects, normalized with SCTransform v2
comb.list <- vector("list",4)
comb.list[[1]] <- CreateSeuratObject(counts=l_data, meta.data = meta.seurat[colnames(l_data),])
comb.list[[1]] <- AddMetaData(comb.list[[1]], lizard.0@meta.data)
Idents(comb.list[[1]]) <- comb.list[[1]]$ident_track
comb.list[[1]] <- SCTransform(comb.list[[1]], vars.to.regress = vars_l, variable.features.n = 5000, vst.flavor = "v2")

comb.list[[2]] <- CreateSeuratObject(counts=t_data, meta.data = meta.seurat[colnames(t_data),])
comb.list[[2]] <- AddMetaData(comb.list[[2]], turtle.0@meta.data)
Idents(comb.list[[2]]) <- comb.list[[2]]$ident_track
comb.list[[2]] <- SCTransform(comb.list[[2]], vars.to.regress = vars_t, variable.features.n = 5000, vst.flavor = "v2")

comb.list[[3]] <- CreateSeuratObject(counts=s_data, meta.data = meta.seurat[colnames(s_data),])
comb.list[[3]] <- AddMetaData(comb.list[[3]], pleuro.0@meta.data)
Idents(comb.list[[3]]) <- comb.list[[3]]$ident_track
comb.list[[3]] <- SCTransform(comb.list[[3]], vars.to.regress = vars_s, variable.features.n = 5000, vst.flavor = "v2")

comb.list[[4]] <- CreateSeuratObject(counts=y_data, meta.data = meta.seurat[colnames(y_data),])
comb.list[[4]] <- AddMetaData(comb.list[[4]], mouse.y0@meta.data)
Idents(comb.list[[4]]) <- comb.list[[4]]$ident_track
comb.list[[4]] <- SCTransform(comb.list[[4]], vars.to.regress = vars_y, variable.features.n = 5000, vst.flavor = "v2")

# taking top 5000 variable genes from each dataset 
lizard.vf <- VariableFeatures(comb.list[[1]])[1:ngenes] 
pleuro.vf <- VariableFeatures(comb.list[[3]])[1:ngenes]
turtle.vf <- VariableFeatures(comb.list[[2]])[1:ngenes]
mouse.vf  <- VariableFeatures(comb.list[[4]])[1:ngenes]

# intersect each list of variable genes among all the species in order to balance the proportion of variable genes for each species
comb.features <- intersect(lizard.vf, pleuro.vf)
comb.features <- intersect(comb.features, turtle.vf)
comb.features <- intersect(comb.features, mouse.vf) 

options(future.globals.maxSize = 60000 * 1024^2)
# removing those variable genes that have a 'salt and pepper' expression
comb.features <- comb.features[comb.features %nin% noisy.genes]

# finding anchors and integrating data
comb.list.prep <- PrepSCTIntegration(object.list = comb.list, anchor.features = comb.features, verbose = FALSE, assay ="SCT")
comb.anchors <- FindIntegrationAnchors(object.list = comb.list.prep, normalization.method = "SCT", anchor.features = comb.features, verbose = FALSE, reduction="cca")
comb.integrated <- IntegrateData(anchorset = comb.anchors, normalization.method = "SCT", verbose = FALSE, dims=1:nDims)

comb.integrated <- RunPCA(comb.integrated, verbose = FALSE, npcs = 200)
comb.integrated <- RunUMAP(comb.integrated, dims = 1:nDims)
comb.integrated <- FindNeighbors(comb.integrated, dims = 1:nDims)
comb.integrated <- FindClusters(comb.integrated)

# adding metadata [['ident_track']] column
meta_df <- comb.integrated@meta.data
meta_df$data_origin <- sub("\\-.*", "", meta_df$ident_track)
comb.integrated <- AddMetaData(comb.integrated, meta_df)

# preparing the dataset of finding marker genes in the SCT space
comb.integrated <- PrepSCTFindMarkers(comb.integrated)

comb.integrated$integrated_clusters <- Idents(comb.integrated)

# confusion matrices
meta.data.table <- comb.integrated@meta.data[,c("species", "cluster_label", "integrated_clusters")]
hClust.order <- ple_cpal

# Integration mouse - salamander
meta.data.table.mouse <- meta.data.table[meta.data.table$species=="mouse",]
meta.data.table.salamander <- meta.data.table[meta.data.table$species=="salamander",]

m_overlap <- table(meta.data.table.mouse$cluster_label, meta.data.table.mouse$integrated_clusters)
s_overlap <- table(meta.data.table.salamander$cluster_label, meta.data.table.salamander$integrated_clusters)

s_ratio <- s_overlap/rowSums(s_overlap)
m_ratio <- m_overlap/rowSums(m_overlap)

overlap.salmou <- s_ratio %*% t(m_ratio)
overlap.salmou <- as.data.frame(overlap.salmou)

overlap.salmou.order <- (match(hClust.order, rownames(overlap.salmou)))
overlap.salmou <- overlap.salmou[overlap.salmou.order, ]

# Integration salamander - reptiles
meta.data.table.lizard <- meta.data.table[meta.data.table$species=="lizard",]
meta.data.table.turtle <- meta.data.table[meta.data.table$species=="turtle",]
meta.data.table.salamander <- meta.data.table[meta.data.table$species=="salamander",]

l_overlap <- table(meta.data.table.lizard$cluster_label, meta.data.table.lizard$integrated_clusters)
t_overlap <- table(meta.data.table.turtle$cluster_label, meta.data.table.turtle$integrated_clusters)
r_overlap <- rbind(l_overlap, t_overlap)
s_overlap <- table(meta.data.table.salamander$cluster_label, meta.data.table.salamander$integrated_clusters)

s_ratio <- s_overlap/rowSums(s_overlap)
r_ratio <- r_overlap/rowSums(r_overlap)

overlap.salrep <- s_ratio %*% t(r_ratio)
overlap.salrep <- as.data.frame(overlap.salrep)

overlap.salrep.order <- (match(hClust.order, rownames(overlap.salrep)))
overlap.salrep <- overlap.salrep[overlap.salrep.order, ]

# Integration salamander - lizard
meta.data.table.lizard <- meta.data.table[meta.data.table$species=="lizard",]
meta.data.table.salamander <- meta.data.table[meta.data.table$species=="salamander",]

l_overlap <- table(meta.data.table.lizard$cluster_label, meta.data.table.lizard$integrated_clusters)
s_overlap <- table(meta.data.table.salamander$cluster_label, meta.data.table.salamander$integrated_clusters)

s_ratio <- s_overlap/rowSums(s_overlap)
l_ratio <- l_overlap/rowSums(l_overlap)

overlap.salliz <- s_ratio %*% t(l_ratio)
overlap.salliz <- as.data.frame(overlap.salliz)

overlap.salliz.order <- (match(hClust.order, rownames(overlap.salliz)))
overlap.salliz <- overlap.salliz[overlap.salliz.order, ]

# Integration salamander - turtle
meta.data.table.turtle <- meta.data.table[meta.data.table$species=="turtle",]
meta.data.table.salamander <- meta.data.table[meta.data.table$species=="salamander",]

t_overlap <- table(meta.data.table.turtle$cluster_label, meta.data.table.turtle$integrated_clusters)
s_overlap <- table(meta.data.table.salamander$cluster_label, meta.data.table.salamander$integrated_clusters)

s_ratio <- s_overlap/rowSums(s_overlap)
t_ratio <- t_overlap/rowSums(t_overlap)

overlap.saltur <- s_ratio %*% t(t_ratio)
overlap.saltur <- as.data.frame(overlap.saltur)

overlap.saltur.order <- (match(hClust.order, rownames(overlap.saltur)))
overlap.saltur <- overlap.saltur[overlap.saltur.order, ]

# renaming clusters starting from 1 instead of 0
DefaultAssay(comb.integrated) <- "integrated"
new.ids <- seq(1, length(levels(Idents(comb.integrated))), by=1)
names(new.ids) <- levels(Idents(comb.integrated))
comb.integrated <- RenameIdents(comb.integrated, new.ids)
species_annot <- comb.integrated@meta.data$species
meta_df <- comb.integrated@meta.data
names(species_annot) <- rownames(comb.integrated@meta.data)
meta_df$speciesTree <- species_annot

comb.integrated <- AddMetaData(comb.integrated, meta_df)

# saving the final integrated object
saveRDS(comb.integrated, file=paste0(date, "Dims", nDims,"_features_", length(comb.features), ".rds"))

library(conos)
library(pagoda2)
library(speciesTree)

# building the Dendrogram with leaves colores by species' mixing
# get the data from integrated assay
expression.matrix <- GetAssayData(comb.integrated, slot="data", assay="integrated")
meta_clusters <- as.integer(Idents(comb.integrated))
names(meta_clusters) <- names(Idents(comb.integrated))
upperlevelinfo = NULL
species_annot <- comb.integrated@meta.data$species
species_annot <- sub("turtle", "reptile", species_annot)
species_annot <- sub("lizard", "reptile", species_annot)
names(species_annot) <- rownames(comb.integrated@meta.data)

#we modified the original function for building the distance matrices using SPearman method
cluster.matrix.expression.distances.S <-
  function (counts, groups = NULL, useVariablegenes = F, variableGenes = NULL,
            dist = "cor", use.single.cell.comparisons = FALSE, use.scaled.data = FALSE,
            min.cluster.size = 1, max.n.cells = Inf, n.cells = 200)
  {
    require(abind)
    require(sccore)
    if (is.null(groups)) {
      stop("no groups specified")
    }
    else {
      groups <- as.factor(groups)
    }
    valid.dists <- c("JS", "cor")
    if (!dist %in% valid.dists)
      stop(paste("only the following distance types are supported:",
                 paste(valid.dists, collapse = ", ")))
    cl <- factor(groups[match(rownames(counts), names(groups))],
                 levels = levels(groups))
    tt <- table(cl)
    empty <- tt < min.cluster.size
    if (any(tt > max.n.cells)) {
      scn <- unlist(tapply(names(cl), cl, function(x) sample(x,
                                                             min(max.n.cells, length(x)))))
      cl[!(names(cl) %in% scn)] <- NA
      tt <- table(cl)
    }
    if (useVariablegenes) {
      if (is.null(variableGenes)) {
        stop("no variable genesets provided")
      }
      else {
        counts <- counts[, colnames(counts) %in% variableGenes]
      }
    }
    if (use.single.cell.comparisons) {
      tcd <- parallel::mclapply(1:n.cells, function(i) {
        scn <- unlist(tapply(names(cl), cl, function(x) sample(x,
                                                               1)))
        tc <- as.matrix(counts[na.omit(as.character(scn)),
        ])
        rownames(tc) <- names(scn)[!is.na(scn)]
        tc <- tc[match(levels(cl), rownames(tc)), ]
        rownames(tc) <- levels(cl)
        if (dist == "JS") {
          if (use.scaled.data) {
            tc <- t(tc)
            tcd <- pagoda2:::jsDist(tc)
            dimnames(tcd) <- list(colnames(tc), colnames(tc))
          }
          else {
            tc <- t(tc/pmax(1, rowSums(tc)))
            tcd <- pagoda2:::jsDist(tc)
            dimnames(tcd) <- list(colnames(tc), colnames(tc))
          }
        }
        else {
          if (use.scaled.data) {
            tcd <- 1 - cor(t(tc), method="s")
          }
          else {
            tc <- log10(t(tc/pmax(1, rowSums(tc))) * 1000 +
                          1)
            tcd <- 1 - cor(tc, method="s")
          }
        }
        tcd[empty, ] <- tcd[, empty] <- NA
        tcd
      }, mc.cores = 10) %>% abind(along = 3) %>% apply(c(1,
                                                         2), median, na.rm = T)
    }
    else {
      tc <- sccore:::colSumByFactor(counts, cl)
      tc <- tc[-1, , drop = F]
      if (dist == "JS") {
        tc <- t(tc/pmax(1, rowSums(tc)))
        tcd <- pagoda2:::jsDist(tc)
        dimnames(tcd) <- list(colnames(tc), colnames(tc))
      }
      else {
        if (use.scaled.data) {
          tc <- t(tc)
        }
        else {
          tc <- log10(t(tc/pmax(1, rowSums(tc))) * 1000 +
                        1)
        }
        tcd <- 1 - cor(tc, method="s")
      }
    }
    return(tcd)
  }

# obtaining distance matrix and performing hierarchical clustering
d <- cluster.matrix.expression.distances.S(t(expression.matrix), groups=meta_clusters, dist="cor",  useVariablegenes=FALSE,  use.scaled.data=TRUE)
dendr <- hclust(as.dist(d), method='ward.D2')
dend <- as.dendrogram(dendr)

# coloring dendrogram leaves according to the proportion of cells from each species in each cluster
dendr <- TransferDend(dend, renameCluster = FALSE, cls.groups = meta_clusters)
cls.groups <- dendr$new.groups
dend <- dendr$dendrogram
leafcontent <- dendr$leafcontent
stability.measurements = NULL
dend <- AddTreeAttribute(dend, species_annot, leafcontent)
dend <- dendSetWidthBysize(dend, scale = 8)
colorpallete <- colorRampPalette(c("blue", "grey", "grey",  "grey", "red"))(101)

upperLevelnodes = NULL
fac <- as.factor(species_annot)
totalCells <- table(fac)

cc2col <- function(cc, rate=15, base=0.001){
  cc <- round((cc[2:4]/totalCells)/sum(cc[2:4]/totalCells), 2)
  cv <- cc
  cv <- dexp(cv, rate)
  cv <- cv/rate * (1-base)
  col <- adjustcolor(rgb(cv[1],cv[2],cv[3], 1), offset = c(0.5, 0.5, 0.5, 0.1))
  return(col)
}

cbm <- function(d,fac) {
  if(is.leaf(d)) {
    #lc <- fac[leafContent[[as.numeric(attr(d,'label'))]]]
    cc <- attr(d, "cc")
    col <- cc2col(cc)
    attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
    return(d);
  } else {
    oa <- attributes(d);
    d <- lapply(d,cbm,fac=fac);
    attributes(d) <- oa;
    cc <- attr(d, "cc")
    col <- cc2col(cc)
    attr(d,"edgePar") <- c(attr(d,"edgePar"),list(col=col))
    return(d);
  }
}

dend <- cbm(dend,species_annot)
tree <- dend

saveRDS(tree, file = paste0(date,"Spearman_tree.rds"))