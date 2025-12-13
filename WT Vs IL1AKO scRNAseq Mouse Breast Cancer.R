library(Signac)
library(Seurat)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(patchwork)
set.seed(1234)
library(SeuratDisk)
library(patchwork)
library(harmony)
library(rliger)
library(reshape2)
library(RColorBrewer)
library(dplyr)
library(limma)
library(NMF)
library(monocle3)
library(CellChat)
library(ggplot2)
library(ggplot2)
library(ggpubr)
library(grid)
library(ComplexHeatmap)
library(openxlsx)
library(ggalluvial)
library(bsicons)

#Load Raw data
wt1 <- Read10X("WT1/")
wt2 <- Read10X("WT2/")
ko1 <- Read10X("KO1/")
ko2 <- Read10X("KO2/")

#Removing prefix mm10___ from features
wt1@Dimnames[1] <- lapply(wt1@Dimnames[1], function(x) gsub("mm10___", "", x))
wt2@Dimnames[1] <- lapply(wt2@Dimnames[1], function(x) gsub("mm10___", "", x))
ko1@Dimnames[1] <- lapply(ko1@Dimnames[1], function(x) gsub("mm10___", "", x))
ko2@Dimnames[1] <- lapply(ko2@Dimnames[1], function(x) gsub("mm10___", "", x))


#Create Seurat objects
WT1 <- CreateSeuratObject(wt1, project = "WT")
WT2 <- CreateSeuratObject(wt2, project = "WT")
KO1 <- CreateSeuratObject(ko1, project = "KO")
KO2 <- CreateSeuratObject(ko2, project = "KO")

#Remove unnecessary objects
rm(wt1, wt2, ko1, ko2)

# Calculate percent.mt and perform filtering for each sample
WT1[["percent.mt"]] <- PercentageFeatureSet(WT1, pattern = "^MT-")
WT2[["percent.mt"]] <- PercentageFeatureSet(WT2, pattern = "^MT-")
KO1[["percent.mt"]] <- PercentageFeatureSet(KO1, pattern = "^MT-")
KO2[["percent.mt"]] <- PercentageFeatureSet(KO2, pattern = "^MT-")

WT1 <- subset(WT1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
WT2 <- subset(WT2, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
KO1 <- subset(KO1, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)
KO2 <- subset(KO2, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 15)

# Merge the Seurat objects and rename cell identities
WT <- merge(WT1, WT2)
KO <- merge(KO1, KO2)

# Update cell identities
Idents(WT) <- "WT"
Idents(KO) <- "KO"

# Combine the merged Seurat objects
pan.integrated <- merge(WT, KO)

# Update cell identities in pan.integrated
pan.integrated_merged <- subset(pan.integrated, idents = c("WT", "KO"))

pan.integrated_merged <- NormalizeData(pan.integrated_merged, verbose = F, normalization.method = "LogNormalize")
pan.integrated_merged <- FindVariableFeatures(pan.integrated_merged, selection.method = "vst", nfeatures = 2000, verbose = F)

DefaultAssay(pan.integrated_merged) <- "RNA"
pan.integrated_merged <- ScaleData(pan.integrated_merged, verbose = FALSE)
pan.integrated_merged <- RunPCA(pan.integrated_merged, npcs = 30, verbose = FALSE)
pan.integrated_merged <- RunUMAP(pan.integrated_merged, reduction = "pca", dims = 1:30)
pan.integrated_merged <- FindNeighbors(pan.integrated_merged, reduction = "pca", dims = 1:30)
pan.integrated_merged <- FindClusters(pan.integrated_merged, resolution = 0.5)

pan.integrated_merged <- JoinLayers(pan.integrated_merged)
pan.integrated_merged <- RenameIdents(pan.integrated_merged,`2` = "Mono/mac" ,`10` = "DC",`9` = "pDC",`12` = "Neutro", `13` = "Gamma-delta cells", `7` = "Other T cells",`3` = "CD8 cells",`6` = "Fibro", `4` = "CD4 cells",`1` = "NKT",`0` = "NK Cells", `14` = "B cells",`5` = "T Reg", `11` = "Endo",`8` = "Tumor")

#WT Vs IL1aKO whole tumor UMAP is made
DimPlot(pan.integrated_merged, reduction = "umap", label = T,split.by = "orig.ident" )

#Reclustering Mono/Mac
DefaultAssay(pan.integrated_merged) <- "RNA"
list <- c("Mono/mac")
Myeloid <- subset(pan.integrated_merged, idents = list)

DefaultAssay(Myeloid) <- "RNA"
Myeloid <- ScaleData(Myeloid, verbose = FALSE)
Myeloid <- RunPCA(Myeloid, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
Myeloid <- RunUMAP(Myeloid, reduction = "pca", dims = 1:10)
Myeloid <- FindNeighbors(Myeloid, reduction = "pca", dims = 1:10)
Myeloid <- FindClusters(Myeloid, resolution = 0.3)
#WT Vs IL1aKO reclustered Mono/Mac UMAP is made
DimPlot(Myeloid, reduction = "umap",label = T, split.by = "orig.ident")

#Genes used to identify unique clusters
VlnPlot(Myeloid,features = c("Ly6c2","Itgax","Adgre1","Retnla","Mrc1","Cx3cr1","C1qa","Ccr2","H2-Aa"),stack = T,flip = T)+theme(axis.title.y = element_text(size = 14, color = "black", face = "bold"),axis.text.y = element_text(size = 12,colour = "black",face = "bold"),axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1,size = 12,colour = "black",face = "bold"),axis.line = element_line(size = 0.8, color = "black"))
Myeloid <- RenameIdents(Myeloid, Myeloid <- RenameIdents(Myeloid, `2` = "Mono",`1` = "Mac1", `3` = "Mac2",`4` = "Mac3",`6`= "Mac4",`0` = "Mac5",`5` = "Mac6"))

#Moncole3 traj analysis
DefaultAssay(Myeloid) <- "RNA"

x <- SplitObject(Myeloid, split.by = "orig.ident")


b<- as.cell_data_set(x$WT)
c<- as.cell_data_set(x$KO)
b <- cluster_cells(cds = b, reduction_method = "UMAP")
c <- cluster_cells(cds = c, reduction_method = "UMAP")

b <- learn_graph(b, use_partition = TRUE,learn_graph_control = list(ncenter= 200))
c <- learn_graph(c, use_partition = TRUE,learn_graph_control = list(ncenter= 200))

WT.cds <- order_cells(b)
KO.cds <- order_cells(c)


plot_cells(
  cds =  WT.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points = F,
  label_roots = F,
  label_leaves = F,
)+ theme(axis.title.y = element_text(size = 14, color = "black", face = "bold"),axis.text.y = element_text(size = 12,colour = "black",face = "bold"),axis.text.x = element_text(vjust = 0.5, hjust = 1,size = 12,colour = "black",face = "bold"),axis.line = element_line(size = 0.8, color = "black"))


plot_cells(
  cds = KO.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
)+ theme(axis.title.y = element_text(size = 14, color = "black", face = "bold"),axis.text.y = element_text(size = 12,colour = "black",face = "bold"),axis.text.x = element_text(vjust = 0.5, hjust = 1,size = 12,colour = "black",face = "bold"),axis.line = element_line(size = 0.8, color = "black"))


#WT Vs IL1aKO genes were 

markers <- c("Itgam","Ccr2","Ly6c2","Itgax","Adgre1","Retnla","Mrc1","Gatm","C1qa","Mafb","Apoe","H2-Aa","H2-K1","Cx3cr1","Cadm1","Trem2","Spp1","Folr2","Vcam1","Nos2","Arg1","Chil3","Hilpda","Bnip3","Ndrg1","Csf1r","Csf2ra","Csf2rb","Csf3r")
DotPlot(Myeloid, features = markers, split.by = "orig.ident",   cols = c("Blue", "Red")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,size = 12))

#******CellChat******#


b <- SplitObject(pan.integrated_merged, split.by = "orig.ident")

# WT and KO objects assigned
WT_obj <- b$WT
KO_obj <- b$KO


# Create CellChat objects for WT and KO
# for WT
meta_WT <- data.frame(
  labels = Idents(WT_obj),
  samples = WT_obj$orig.ident,
  row.names = names(Idents(WT_obj))
)
WT <- createCellChat(object = WT_obj, meta = meta_WT, group.by = "labels")

# for KO
meta_KO <- data.frame(
  labels = Idents(KO_obj),
  samples = KO_obj$orig.ident,
  row.names = names(Idents(KO_obj))
)
KO <- createCellChat(object = KO_obj, meta = meta_KO, group.by = "labels")

# Save CellChat objects
saveRDS(WT, file = "WT_chat.rds")
saveRDS(KO, file = "KO_chat.rds")


# Update CellChat objects
WT <- updateCellChat(WT)
KO <- updateCellChat(KO)


# Merge CellChat objects
object.list <- list(WT = WT, KO = KO)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
save(cellchat, file = "cellchat_merged.RData")

execution.time <- Sys.time() - ptm
cat("Execution time (seconds):", as.numeric(execution.time, units = "secs"), "\n")


# Compare signaling patterns

i <- 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)

ht1 <- netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing",
                                         signaling = pathway.union, title = names(object.list)[i],
                                         width = 5, height = 6)
ht2 <- netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing",
                                         signaling = pathway.union, title = names(object.list)[i+1],
                                         width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# Bubble plots
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = 1:15,
                        comparison = c(1,2), max.dataset = 1,
                        title.name = "MonoMac to All", angle.x = 45, remove.isolate = TRUE)
gg2 <- netVisual_bubble(cellchat, sources.use = 1:15, targets.use = 1,
                        comparison = c(1,2), max.dataset = 1,
                        title.name = "All to MonoMac", angle.x = 45, remove.isolate = TRUE)
gg1 + gg2

# Identify unique and shared pathways
unique_WT <- setdiff(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)
unique_KO <- setdiff(object.list[[2]]@netP$pathways, object.list[[1]]@netP$pathways)
shared_pathways <- intersect(object.list[[1]]@netP$pathways, object.list[[2]]@netP$pathways)

cat("Unique pathways in WT:\n"); print(unique_WT)
cat("Unique pathways in KO:\n"); print(unique_KO)
cat("Shared pathways:\n"); print(shared_pathways)

# Chord diagrams

pdf("CellChat_ChordPlots_WT_KO.pdf", width = 8, height = 6)

# Loop through unique WT pathways
if (length(unique_WT) > 0) {
  for (pathway in unique_WT) {
    netVisual_chord_cell(WT, targets.use = "Mono/mac", signaling = pathway,
                         title.name = paste0("Unique Pathway in WT: ", pathway))
  }
}
dev.off()

# Loop through unique KO pathways
if (length(unique_KO) > 0) {
  for (pathway in unique_KO) {
    netVisual_chord_cell(KO, targets.use = "Mono/mac", signaling = pathway,
                         title.name = paste0("Unique Pathway in KO: ", pathway))
  }
}
dev.off()


# Example: To inspect TGFb signaling in WT

TGFb_genes <- cellchat@LR$WT$LRsig[cellchat@LR$WT$LRsig$pathway_name == "TGFb", ]
print(TGFb_genes)
print(TGFb_genes$receptor.symbol)

for (gene in genes) {
  
  # Generate violin plot for the current gene
  plot <- VlnPlot(
    a, features = gene, split.by = "orig.ident", split.plot = TRUE, 
    flip = TRUE, combine = TRUE, cols = c("Blue", "Red")
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12))
  
  # Find max y-value and adjust significance star positioning dynamically
  y_max <- max(FetchData(a, vars = gene), na.rm = TRUE) * 1.2  # Increase 20% buffer to avoid cropping
  
  # Add statistical significance with adjusted y.position
  plot <- plot + stat_compare_means(
    method = "t.test",
    label = "p.signif",
    hide.ns = TRUE,
    size = 8
  )
  
  # Save the plot as a PDF with a larger height to accommodate stars
  pdf(paste0(gene, "_Ttestviolin_plot.pdf"), width = 11, height = 9)  
  print(plot)
  dev.off()
}

#WT Vs IL1aKO cellchat


