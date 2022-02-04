library(gplots)

# integrated is the integrated Seurat object. 
# save the clusters as a metadata column
integrated$integrated_clusters <- Idents(integrated)

# cluster_labels and ClusterName are the metadata columns where the original cluster identities of each dataset are stored; replace with whatever relevant metadata column
meta.data.table <- integrated@meta.data[,c("species","cluster_labels", "ClusterName","integrated_clusters")]

meta.data.table.mouse <- meta.data.table[meta.data.table$species=="mouse",]
meta.data.table.salamander <- meta.data.table[meta.data.table$species=="salamander",]

m_overlap <- table(meta.data.table.mouse$ClusterName, meta.data.table.mouse$integrated_clusters)
s_overlap <- table(meta.data.table.salamander$cluster_labels, meta.data.table.salamander$integrated_clusters)

s_ratio <- s_overlap/rowSums(s_overlap)
m_ratio <- m_overlap/rowSums(m_overlap)

overlap <- s_ratio %*% t(m_ratio)
 
# plot
heatmap.2(as.matrix(t(overlap)),trace="none", col=colorRampPalette(c("white","grey","black")), scale="none", Colv=F, Rowv=F, dendrogram="none")
