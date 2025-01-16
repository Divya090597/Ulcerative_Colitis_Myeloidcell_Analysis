



# Load the seurat objects -------------------------------------------------


obj1 = readRDS("myeloid_PMC6662628_seurat.rds")
obj2 = readRDS("myeloid_GSE250498_subsetted.rds")


# Normalize 
obj1 <- NormalizeData(obj1)

obj2 <- NormalizeData(obj2)

# Find variable features for obj1 and obj2
obj1 <- FindVariableFeatures(obj1, selection.method = "vst", nfeatures = 2000)
obj2 <- FindVariableFeatures(obj2, selection.method = "vst", nfeatures = 2000)

# Select integration features
features <- SelectIntegrationFeatures(object.list = list(obj1, obj2))

# Perform PCA to prepare for integration
obj1 <- ScaleData(obj1, features = features)
obj1 <- RunPCA(obj1, features = features)

obj2 <- ScaleData(obj2, features = features)
obj2 <- RunPCA(obj2, features = features)

# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = list(obj1, obj2), anchor.features = features)

# Perform integration 
integrated_data <- IntegrateData(anchorset = anchors, dims = 1:30)

# Perform a integration analysis
DefaultAssay(integrated_data) <- "integrated"

integrated_data <- ScaleData(integrated_data, verbose = FALSE)
integrated_data <- RunPCA(integrated_data, npcs = 30, verbose = FALSE)
integrated_data <- FindNeighbors(integrated_data, reduction = "pca", dims = 1:30)
integrated_data <- FindClusters(integrated_data, resolution = 0.5)
integrated_data <- RunUMAP(integrated_data, reduction = "pca", dims = 1:30)

# Dimplot

DimPlot(integrated_data, reduction = "umap",group.by = "fine_annotations_MK_V4",label = TRUE)

saveRDS(integrated_data, file = "Integrated_myeoild.rds")

# Check how many cells are missing annotations
sum(is.na(integrated_data@meta.data$fine_annotations_MK_V4))

# Create a mapping of cluster IDs to annotations
annotation <- c(
  "0" = "08A-MNP macrophage",
  "1" = "08C-MNP mDC",
  "2" = "08B-MNP monocyte",
  "3" = "08B-MNP monocyte",
  "4" = "08A-MNP macrophage",
  "5" = "08C-MNP mDC",
  "6" = "08C-MNP mDC",
  "7" = "10F-Epithelial DCS",
  "8" = "08C-MNP mDC",
  "9" = "08D-MNP-pDC",
  "10" = "08A-MNP macrophage",
  "11" = "08C-MNP mDC"
  )


annotation <- c(
  "0" = "Macrophages",
  "1" = "Macrophages",
  "2" = "Macrophages",
  "3" = "Inflammatory Monocytes",
  "4" = "Macrophages",
  "5" = "MT-hi",
  "6" = "DC1",
  "7" = "?",
  "8" = "Cycling Monocytes",
  "9" = "DC2",
  "10" = "Macrophages",
  "11" = "DC2"
)



# Fill missing annotations based on cluster IDs
integrated_data@meta.data$fine_annotations_MK_V4[is.na(integrated_data@meta.data$fine_annotations_MK_V4)] <- 
  annotation[as.character(integrated_data@meta.data$seurat_clusters[is.na(integrated_data@meta.data$fine_annotations_MK_V4)])]

# Verify updated annotations
table(integrated_data@meta.data$fine_annotations_MK_V4)



# SingleR -----------------------------------------------------------------


library(SingleR)

# Use SingleR for automated annotation

reference_HPCA <- celldex::HumanPrimaryCellAtlasData()
singleR_result_HPCA <- SingleR(test = integrated_data@assays$integrated$data, ref = reference_HPCA, labels = reference_HPCA$label.main)
integrated_data$SingleR.labels.HPCA <- singleR_result_HPCA$labels

# Plot SingleR annotations

Idents(integrated_data) = "SingleR.labels.HPCA"
DimPlot(integrated_data, reduction = "umap") + ggtitle("SingleR annotation HPCA")

DimPlot(seu_obj_processed,group.by = c("CellOntology_name","Immune_UCell", "Lymphoid_UCell", "PanBcell_UCell" ,"Bcell_UCell", "APC_UCell") ,label = T, label.size = 3) 

genes = c("ITGAM", "ITGAX", "LYZ" ,"IL1B","IL6" , "HECTD3", "XCR1", "CLEC9A", "IRF8", "BATF3", "IRF4", "CD1C", "FCER1A", "IRF7", "GZMB", "TCF4", "CD14", "FCGR3A", "S100A8", "S100A9", "IL3RA" , "CCR2", "ZBTB46", "CD16")
DotPlot(integrated_data, features = genes, group.by = "fine_annotations_MK_V4")+
  coord_fixed()+
  RotatedAxis()

DefaultAssay(integrated_data) = "integrated"

# Find all marker genes for all clusters
all_markers <- FindAllMarkers(
  obj3,
  only.pos = TRUE,              
  min.pct = 0.25,               
  logfc.threshold = 0.25        
)

# View the top markers
head(all_markers)

library(data.table)
fwrite(all_markers, file = "~/data/HECTD3/Integrated_myeoild_upregulated_genes.csv")

# Time to annotate by myself

# Create a contingency table
annotation_overlap <- table(
  integrated_data@meta.data$Cluster,
  integrated_data@meta.data$fine_annotations_MK_V4
)

# View the table
print(annotation_overlap)

library(pheatmap)

# Plot a heatmap of the overlap
pheatmap(
  annotation_overlap,
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = TRUE, 
  color = colorRampPalette(c("white", "blue"))(100)
)

# Create a unified annotation column
integrated_data@meta.data$final_annotation <- ifelse(
  !is.na(integrated_data@meta.data$Cluster),
  integrated_data@meta.data$Cluster,
  integrated_data@meta.data$fine_annotations_MK_V4
)

# Check the distribution of the final annotations
table(integrated_data@meta.data$final_annotation)

DimPlot(integrated_data, group.by = "final_annotation", label = TRUE)


