## Preparing next object for scvi integration: REPTILE  

eggnog_turtle <- read.table("chrysemys_eggnog_pruned.txt" , header=T)
eggnog_lizard <- read.table("pogona_eggnog_pruned.txt" , header=T)
eggnog_salamander <- read.table("pleurodeles_eggnog_pruned.txt" , header=T)

liz.all <- c("Str_PENK_1", "MGE_INs_5", "MCtx_2", "MCtx_4", "MGE_INs_2", "pDCtx", "Str_PENK_2", "Sept_2", 
             "Sept_1", "aDVR_1", "pDVR", "Str_PENK_3", "CGE_INs_2", "OT", "MGE_INs_1", "aDVR_2", "CGE_INs_3", 
             "LCtx_1", "DLA_1", "MCtx_8", "MGE_INs_3", "aDVR_3", "MCtx_7", "aDVR_6", "aDVR_4", "MCtx_1", 
             "aDVR_5", "MCtx_5", "amDVR", "CGE_INs_1", "MCtx_3", "MCtx_6", "ITCs", "Str_TAC1", "MGE_INs_4", 
             "MCtx_10", "Sept_4", "LCtx_2", "Sept_3", "MGE_INs_9", "DLA_2", "MGE_INs_8", 
             "MGE_INs_6", "aDCtx", "Chol_1", "MGE_INs_7", "OT_3", "MCtx_9", "Chol_3", "Chol_2")

tur_all <- c(grep("e",levels(Idents(turtle.0)), value = TRUE),
             grep("i",levels(Idents(turtle.0)), value = TRUE))

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

vars_l = c("percent.mito", "animal")
vars_t = c("animalident", "percent.mito")
vars_s = c("animal", "percent.mt")

lizard.0 <- subset(lizard.0, idents = liz.all)
turtle.0 <- subset(turtle.0, idents = tur_all)
pleuro.0 <- subset(pleuro.0, idents = s_tele)

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

lizard.1 <- CreateSeuratObject(counts=l_data, meta.data = meta.seurat[colnames(l_data),])
lizard.1 <- AddMetaData(lizard.1, lizard.0@meta.data)
Idents(lizard.1) <- lizard.1$ident_track

turtle.1 <- CreateSeuratObject(counts=t_data, meta.data = meta.seurat[colnames(t_data),])
turtle.1 <- AddMetaData(turtle.1, turtle.0@meta.data)
Idents(turtle.1) <- turtle.1$ident_track

pleuro.1 <- CreateSeuratObject(counts=s_data, meta.data = meta.seurat[colnames(s_data),])
pleuro.1 <- AddMetaData(pleuro.1, pleuro.0@meta.data)
Idents(pleuro.1) <- pleuro.1$ident_track

scvi <- merge(lizard.1, y = list(turtle.1, pleuro.1), project = "scvi-tools-3sp")

Idents(lizard.0) <- lizard.0$animal
Idents(pleuro.0) <- pleuro.0$animal
Idents(turtle.0) <- turtle.0$animalident

animal.metadata <- c(Idents(lizard.0), Idents(pleuro.0), Idents(turtle.0))

scvi[["animal.ident"]] <- animal.metadata
scvi[["animal"]] <- NULL
scvi[["animalident"]] <- NULL


Idents(lizard.0) <- lizard.0$percent.mito
Idents(pleuro.0) <- pleuro.0$percent.mt
Idents(turtle.0) <- turtle.0$percent.mito

mito.metadata <- c(Idents(lizard.0), Idents(pleuro.0), Idents(turtle.0))

scvi[["mt.percent"]] <- mito.metadata
scvi[["percent.mt"]] <- NULL
scvi[["percent.mito"]] <- NULL


saveRDS(scvi, file = "scvi_tools_3sp_merged.rds")

lizard.vf <- VariableFeatures(comb.list[[1]])[1:1000] 
pleuro.vf <- VariableFeatures(comb.list[[3]])[1:1000]
turtle.vf <- VariableFeatures(comb.list[[2]])[1:1000]

comb.features <- intersect(lizard.vf, pleuro.vf)
comb.features <- intersect(comb.features, turtle.vf)





