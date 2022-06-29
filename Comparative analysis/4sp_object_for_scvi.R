library(Seurat)
library(corrplot)
library(MatrixGenerics)
library(ggplot2)
library(dplyr)
library(gplots)
library(Hmisc)
library(stringr)

lizard.0 <- readRDS("lizard_clean_integration1000.RDS")
turtle.0 <- readRDS("170528_turtle_neur.RDS")
pleuro.0 <- readRDS("neurons_final5.rds")
mouse.y0 <- readRDS("mouse_yao_10x_downsampled3.rds")

eggnog_turtle <- read.table("chrysemys_eggnog_pruned.txt" , header=T)
eggnog_lizard <- read.table("pogona_eggnog_pruned.txt" , header=T)
eggnog_salamander <- read.table("pleurodeles_eggnog_pruned.txt" , header=T)
eggnog_mouse <- read.table("mus_eggnog_pruned.txt" , header=T)

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

vars_l = c("percent.mito", "animal")
vars_t = c("animalident", "percent.mito")
vars_s = c("animal", "percent.mt")
vars_y = c("external_donor_name_label")

lizard.0 <- subset(lizard.0, idents = liz_cpal)
turtle.0 <- subset(turtle.0, idents = tur_cpal)
pleuro.0 <- subset(pleuro.0, idents = ple_cpal)
mouse.y0 <- subset(mouse.y0, idents = y_idents_keep)

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

lizard.1 <- CreateSeuratObject(counts=l_data, meta.data = meta.seurat[colnames(l_data),])
lizard.1 <- AddMetaData(lizard.1, lizard.0@meta.data)
Idents(lizard.1) <- lizard.1$ident_track

turtle.1 <- CreateSeuratObject(counts=t_data, meta.data = meta.seurat[colnames(t_data),])
turtle.1 <- AddMetaData(turtle.1, turtle.0@meta.data)
Idents(turtle.1) <- turtle.1$ident_track

pleuro.1 <- CreateSeuratObject(counts=s_data, meta.data = meta.seurat[colnames(s_data),])
pleuro.1 <- AddMetaData(pleuro.1, pleuro.0@meta.data)
Idents(pleuro.1) <- pleuro.1$ident_track

mouse.y1 <- CreateSeuratObject(counts=y_data, meta.data = meta.seurat[colnames(y_data),])
mouse.y1 <- AddMetaData(mouse.y1, mouse.y0@meta.data)
Idents(mouse.y1) <- mouse.y1$ident_track

scvi.int <- merge(lizard.1, y = list(turtle.1, pleuro.1, mouse.y1), project = "scvi-tools")





