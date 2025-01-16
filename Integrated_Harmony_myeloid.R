library(harmony)

# Merge the Seurat objects
merged_data <- merge(obj1, y = obj2)

# Normalize and find variable features
merged_data <- NormalizeData(merged_data)
merged_data <- FindVariableFeatures(merged_data)
merged_data <- ScaleData(merged_data)
merged_data <- RunPCA(merged_data)

# Add a batch column based on `orig.ident`
merged_data$batch <- merged_data$orig.ident

# Verify the batch column
table(merged_data$batch)

# Ensure the batch column is a factor
merged_data$batch <- as.factor(merged_data$batch)

# Verify the levels of the batch factor
levels(merged_data$batch)



# Run Harmony for batch correction
merged_data <- RunHarmony(
  merged_data,
  group.by.vars = "batch"  # Specify the metadata column indicating batch
)
# Verify Harmony embeddings
head(Embeddings(merged_data, reduction = "harmony"))

# Run UMAP and clustering on Harmony embeddings
merged_data <- RunUMAP(merged_data, reduction = "harmony", dims = 1:30)
merged_data <- FindNeighbors(merged_data, reduction = "harmony", dims = 1:30)
merged_data <- FindClusters(merged_data, resolution = 0.5)

# Visualize the integrated data
DimPlot(merged_data, reduction = "umap", label = TRUE)

genes = c("ITGAM", "ITGAX", "LYZ" ,"IL1B","IL6" , "HECTD3", "XCR1", "CLEC9A", "IRF8", "BATF3", "IRF4", "CD1C", "FCER1A", "IRF7", "GZMB", "TCF4", "CD14", "FCGR3A", "S100A8", "S100A9", "IL3RA" , "CCR2", "ZBTB46","CD1A","HLA-DRB5")
DotPlot(merged_data, features = genes, group.by = "RNA_snn_res.0.5")+
  coord_fixed()+
  RotatedAxis()

# Check how many cells are missing annotations
sum(is.na(merged_data@meta.data$fine_annotations_MK_V4))

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

saveRDS(merged_data, file = "~/data/HECTD3/Integraged_Harmony.rds" )

