.libPaths("/burg/home/ao2721/R/")
library(Seurat)
library(corrplot)
library(MatrixGenerics)
library(ggplot2)
library(dplyr)
library(gplots)
library(Hmisc)
library(stringr)

date = "R18_220523_CCA_"

# number of dimensions to use for IntegrateData
nDims = 80

# number of variable genes to use to identify integration anchors
ngenes = 2000

noisy.genes <- readRDS("/burg/tosches/users/ao2721/integrated/noisy_list_11.rds")

# Reading all the individual files that we'll integrate
lizard.0 <- readRDS("/burg/tosches/users/ao2721/integrated/lizard_clean_integration1000.rds")
turtle.0 <- readRDS("/burg/tosches/users/ao2721/integrated/170528_turtle_neur.RDS")
pleuro.0 <- readRDS("/burg/tosches/users/ao2721/integrated/neurons_final5.rds")

# Reading the eggnong files for each species to select one to one orthologues in the downstream analysis
eggnog_turtle <- read.table("/burg/tosches/users/ao2721/integrated/chrysemys_eggnog_pruned.txt" , header=T)
eggnog_lizard <- read.table("/burg/tosches/users/ao2721/integrated/pogona_eggnog_pruned.txt" , header=T)
eggnog_salamander <- read.table("/burg/tosches/users/ao2721/integrated/pleurodeles_eggnog_pruned.txt" , header=T)

# liz no PThE, no non-telencehalic and no OB interneurons
liz_all <- c("Str_PENK_1", "MGE_INs_5", "MCtx_2", "MCtx_4", "MGE_INs_2", "pDCtx", "Str_PENK_2", "Sept_2", 
             "Sept_1", "aDVR_1", "pDVR", "Str_PENK_3", "CGE_INs_2", "OT", "MGE_INs_1", "aDVR_2", "CGE_INs_3", 
             "LCtx_1", "DLA_1", "MCtx_8", "MGE_INs_3", "aDVR_3", "MCtx_7", "aDVR_6", "aDVR_4", "MCtx_1", 
             "aDVR_5", "MCtx_5", "amDVR", "CGE_INs_1", "MCtx_3", "MCtx_6", "ITCs", "Str_TAC1", "MGE_INs_4", 
             "MCtx_10", "Sept_4", "LCtx_2", "Sept_3", "MGE_INs_9", "DLA_2", "MGE_INs_8", 
             "MGE_INs_6", "aDCtx", "Chol_1", "MGE_INs_7", "OT_3", "MCtx_9", "Chol_3", "Chol_2")

tur_all <- c(grep("e",levels(Idents(turtle.0)), value = TRUE),
             grep("i",levels(Idents(turtle.0)), value = TRUE))

# no OB and no MT, no migration from Hypo (TEGLU37, TEGLU38)
s_tele <- c("TEGLU1", "TEGLU2", "TEGLU4", "TEGLU5", "TEGLU3", "TEGLU6", "TEGLU7", "TEGLU14", 
            "TEGLU15", "TEGLU16", "TEGLU12", "TEGLU11", "TEGLU9", "TEGLU10", "TEGLU8", "TEGLU13", 
            "TEGLU20", "TEGLU25", "TEGLU24", "TEGLU23", "TEGLU22", "TEGLU19", "TEGLU17", "TEGLU18", 
            "TEGLU21", "TEGLU31", "TEGLU30", "TEGLU36", "TEGLU35", "TEGLU33", "TEGLU34", "TEGLU39", "TEGLU40", 
            "TEGLU41", "TEGABA52",  "TEGABA53", "TEGABA58", "TEGABA54",  "TEGABA55", "TEGABA57", "TEGABA56", 
            "TEGABA59", "TEGABA61", "TEGABA60", "TEGABA62", "TEGABA63", 
            "TEGABA64", "TEGABA65", "TEGABA34", "TEGABA35", "TEGABA33", "TEGABA31", "TEGABA32", "TEGABA30", 
            "TEGABA29", "TEGLU32", "TEGABA28", "TEGABA37", "TEGABA36", "TEGABA38", "TEGABA51", "TEGABA48", 
            "TEGABA47", "TEGABA46", "TEGABA50", "TEGABA49", "TEGABA39", "TEGABA40", "TEGABA42", "TEGABA43", 
            "TEGABA45", "TEGABA44", "TEGLU41", "TEGABA41", "TEGABA25", "TEGABA26", "TEGABA24", "TEGABA23", 
            "TEGABA27", "TEGABA22", "TEGABA21", "TEGABA20", "TEGABA18", "TEGABA19")
Sl = "L"
St = "T"
Ss = "S"

# variables to regress
vars_l = c("percent.mito", "animal")
vars_t = c("animalident", "percent.mito")
vars_s = c("animal", "percent.mt")

# subsetting the objects to include only the comparable cells
lizard.0 <- subset(lizard.0, idents = liz_all)
turtle.0 <- subset(turtle.0, idents = tur_all)
pleuro.0 <- subset(pleuro.0, idents = s_tele)

# getting raw counts from the original objects and selecting one to one orthologues
l_data <- GetAssayData(lizard.0, slot="counts")
l_data <- l_data[rownames(l_data) %in% eggnog_lizard[,1],]

t_data <- GetAssayData(turtle.0, slot="counts") 
t_data <- t_data[rownames(t_data) %in% eggnog_turtle[,1],]

s_data <- GetAssayData(pleuro.0, slot="counts")
s_data <- s_data[rownames(s_data) %in% eggnog_salamander[,1],]

eggnog_l <- eggnog_lizard[eggnog_lizard[,1] %in% rownames(l_data),]
eggnog_t <- eggnog_turtle[eggnog_turtle[,1] %in% rownames(t_data),]
eggnog_s <- eggnog_salamander[eggnog_salamander[,1] %in% rownames(s_data),]

l_data <- l_data[order(match(rownames(l_data), eggnog_l[,1])),]
t_data <- t_data[order(match(rownames(t_data), eggnog_t[,1])),]
s_data <- s_data[order(match(rownames(s_data), eggnog_s[,1])),]

all.equal(rownames(s_data),as.character(eggnog_s[,1]))
all.equal(rownames(t_data),as.character(eggnog_t[,1]))
all.equal(rownames(l_data),as.character(eggnog_l[,1]))

rownames(l_data) <- eggnog_l$eggnog 
rownames(t_data) <- eggnog_t$eggnog 
rownames(s_data) <- eggnog_s$eggnog

dim(l_data)
dim(t_data)
dim(s_data)

l_data[1:4,1:4]
t_data[1:4,1:4]
s_data[1:4,1:4]

# preparing metadata for integrated object
lizard.0[["cluster_label"]] <- Idents(object = lizard.0)
turtle.0[["cluster_label"]] <- Idents(object = turtle.0)
pleuro.0[["cluster_label"]] <- Idents(object = pleuro.0)

meta.seurat <- data.frame(
  ident_track = c(paste(Sl, as.character(Idents(lizard.0)), sep="-"), 
                  paste(St, as.character(Idents(turtle.0)), sep="-"), 
                  paste(Ss, as.character(Idents(pleuro.0)), sep="-")),
  species = c(rep("lizard", length(Idents(lizard.0))),
              rep("turtle", length(Idents(turtle.0))), 
              rep("salamander", length(Idents(pleuro.0)))),
  row.names=c(colnames(l_data),
              colnames(t_data),
              colnames(s_data)))

# we create a list that will contain all the three objects, normalized with SCTransform v2
comb.list <- vector("list",3)
comb.list[[1]] <- CreateSeuratObject(counts=l_data, meta.data = meta.seurat[colnames(l_data),])
comb.list[[1]] <- AddMetaData(comb.list[[1]], lizard.0@meta.data)
Idents(comb.list[[1]]) <- comb.list[[1]]$ident_track
comb.list[[1]] <- SCTransform(comb.list[[1]], vars.to.regress = vars_l, vst.flavor = "v2")

comb.list[[2]] <- CreateSeuratObject(counts=t_data, meta.data = meta.seurat[colnames(t_data),])
comb.list[[2]] <- AddMetaData(comb.list[[2]], turtle.0@meta.data)
Idents(comb.list[[2]]) <- comb.list[[2]]$ident_track
comb.list[[2]] <- SCTransform(comb.list[[2]], vars.to.regress = vars_t, vst.flavor = "v2")

comb.list[[3]] <- CreateSeuratObject(counts=s_data, meta.data = meta.seurat[colnames(s_data),])
comb.list[[3]] <- AddMetaData(comb.list[[3]], pleuro.0@meta.data)
Idents(comb.list[[3]]) <- comb.list[[3]]$ident_track
comb.list[[3]] <- SCTransform(comb.list[[3]], vars.to.regress = vars_s, vst.flavor = "v2")

# selecting integration features (genes), excluding those in the noisy list, finding integration anchors and integrating data with Seurat pipeline
options(future.globals.maxSize = 60000 * 1024^2)
comb.features <- SelectIntegrationFeatures(object.list = comb.list, nfeatures = ngenes)
comb.features <- comb.features[comb.features %nin% noisy.genes]
comb.list.prep <- PrepSCTIntegration(object.list = comb.list, anchor.features = comb.features, verbose = FALSE)
comb.anchors <- FindIntegrationAnchors(object.list = comb.list.prep, normalization.method = "SCT", anchor.features = comb.features, verbose = FALSE, reduction="cca")
comb.integrated <- IntegrateData(anchorset = comb.anchors, normalization.method = "SCT", verbose = FALSE, dims=1:nDims)

comb.integrated <- RunPCA(comb.integrated, verbose = FALSE, npcs = 200)
comb.integrated <- RunUMAP(comb.integrated, dims = 1:nDims)
comb.integrated <- FindNeighbors(comb.integrated, dims = 1:nDims)
comb.integrated <- FindClusters(comb.integrated)

# adding metadata to track back species identities in each original cluster [['ident_track']]
meta_df <- comb.integrated@meta.data
meta_df$data_origin <- sub("\\-.*", "", meta_df$ident_track)
comb.integrated <- AddMetaData(comb.integrated, meta_df)

# finding marker genes in the SCT space
comb.integrated$integrated_clusters <- Idents(comb.integrated)
comb.integrated <- PrepSCTFindMarkers(comb.integrated)

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
saveRDS(comb.integrated, file=paste0(date,"integration_reptiles","_Dims",nDims,"_features_",ngenes,".RDS"))

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

pdf(paste0(date,"tree.pdf"), height = 10, width = 20)
plot(tree)
dev.off()

# saving the colored dendrogram
saveRDS(tree, file = paste0(date,"Spearman_tree.rds"))
