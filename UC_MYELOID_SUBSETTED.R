

# Load Data ---------------------------------------------------------------

seu_obj = readRDS("~/data/HECTD3/UC_biopsy.rds")

genes = c("ITGAM", "ITGAX", "LYZ" ,"IL1B","IL6" , "HECTD3", "XCR1", "CLEC9A", "IRF8", "BATF3", "IRF4", "CD1C", "FCER1A", "IRF7", "GZMB", "TCF4", "CD14", "FCGR3A", "S100A8", "S100A9", "IL3RA" , "CCR2", "ZBTB46", "CD16")
DotPlot(myeloid_resubset, features = genes, group.by = "fine_annotations_MK_V4")+
  coord_fixed()

FeaturePlot(seu_obj, features = c("ITGAM", "ITGAX"))

# Extract expression data for ITGAM and ITGAX
expression_data <- FetchData(seu_obj, vars = c("ITGAM", "ITGAX"))

# summary
summary(expression_data)

# double-positive cells
double_positive_cells <- WhichCells(seu_obj, expression = ITGAM > 0 & ITGAX > 0)

length(double_positive_cells)

DimPlot(seu_obj, cells.highlight = double_positive_cells, cols.highlight = "red")

# Subset Seurat object
double_positive_subset <- subset(seu_obj, cells = double_positive_cells)

# Plot UMAP of the subset
DimPlot(double_positive_subset, reduction = "X_umap")

# Normalize the data
double_positive_subset <- NormalizeData(double_positive_subset)
double_positive_subset <- FindVariableFeatures(double_positive_subset, selection.method = "vst", nfeatures = 2000)
double_positive_subset <- ScaleData(double_positive_subset, features = rownames(double_positive_subset))
double_positive_subset <- RunPCA(double_positive_subset, features = VariableFeatures(object = double_positive_subset))
double_positive_subset <- FindNeighbors(double_positive_subset, dims = 1:10)  
double_positive_subset <- FindClusters(double_positive_subset, resolution = 0.5)  
double_positive_subset <- RunUMAP(double_positive_subset, dims = 1:10)
DimPlot(double_positive_subset, reduction = "umap", group.by = "fine_annotations_MK_V4", label = T)



all_features <- rownames(double_positive_subset)
c("DPP4", "FCGR1A") %in% all_features

VlnPlot(double_positive_subset, features = c("DPP4", "FCGR1A"), pt.size = 0.1)

# Fetch expression data
expr_data <- FetchData(myeloid_subset, vars = c("DPP4", "FCGR1A"))

# Summary
summary(expr_data$DPP4)
summary(expr_data$FCGR1A)

# Co-expression scatter plot
FeatureScatter(myeloid_subset, feature1 = "DPP4", feature2 = "FCGR1A") +
  ggtitle("Co-Expression of CD26 (DPP4) and CD64 (FCGR1A)")

# thresholds for high expression
cd26_threshold <- 0  
cd64_threshold <- 0 

# Categorize cells
myeloid_subset$CD26_CD64_Group <- ifelse(
  FetchData(myeloid_subset, vars = "DPP4")[, 1] > cd64_threshold &
    FetchData(myeloid_subset, vars = "FCGR1A")[, 1] > cd26_threshold, "CD26_High_CD64_High",
  ifelse(FetchData(myeloid_subset, vars = "DPP4")[, 1] > cd64_threshold &
           FetchData(myeloid_subset, vars = "FCGR1A")[, 1] <= cd26_threshold, "CD26_High_CD64_Low",
         ifelse(FetchData(myeloid_subset, vars = "DPP4")[, 1] <= cd64_threshold &
                  FetchData(myeloid_subset, vars = "FCGR1A")[, 1] > cd26_threshold, "CD26_Low_CD64_High",
                "CD26_Low_CD64_Low"))
)



# Table of distribution 
table(myeloid_subset$CD26_CD64_Group)

FeatureScatter(myeloid_subset, feature1 = "DPP4", feature2 = "FCGR1A") +
  ggtitle("Co-Expression of CD26 (DPP4) and CD64 (FCGR1A)")

# Proportion of cells with non-zero expression
cd26_positive <- sum(FetchData(double_positive_subset, vars = "DPP4")[, 1] > 0) / nrow(double_positive_subset)
cd64_positive <- sum(FetchData(double_positive_subset, vars = "FCGR1A")[, 1] > 0) / nrow(double_positive_subset)

cat("Proportion of cells expressing CD26 (DPP4):", cd26_positive, "\n")
cat("Proportion of cells expressing CD64 (FCGR1A):", cd64_positive, "\n")

# Set group identities
Idents(double_positive_subset) <- "CD26_CD64_Group"

# Find markers for CD26_High_CD64_Low
markers_high_low <- FindMarkers(double_positive_subset, ident.1 = "CD26_High_CD64_Low", ident.2 = NULL)

# View top markers
head(markers_high_low)

# Define myeloid clusters
myeloid_subclusters <- c("08B-MNP monocyte", "08A-MNP macrophage", "08C-MNP mDC", "08D-MNP-pDC", "10F-Epithelial DCS")

# Subset the Seurat object
myeloid_resubset <- subset(seu_obj, subset = fine_annotations_MK_V4 %in% myeloid_subclusters)

# Normalize the data
myeloid_resubset <- NormalizeData(myeloid_resubset)
myeloid_resubset <- FindVariableFeatures(myeloid_resubset, selection.method = "vst", nfeatures = 2000)
myeloid_resubset <- ScaleData(myeloid_resubset, features = rownames(myeloid_resubset))
myeloid_resubset <- RunPCA(myeloid_resubset, features = VariableFeatures(object = myeloid_resubset))
myeloid_resubset <- FindNeighbors(myeloid_resubset, dims = 1:10)  
myeloid_resubset <- FindClusters(myeloid_resubset, resolution = 0.5)  
myeloid_resubset <- RunUMAP(myeloid_resubset, dims = 1:10)
DimPlot(myeloid_resubset, reduction = "umap", group.by = "originalexp_snn_res.0.5", label = T)

DimPlot(myeloid_resubset, reduction = "umap", group.by = "fine_annotations_MK_V4", label = T)

# Plot cDC1 markers
FeaturePlot(myeloid_resubset, features = c("CLEC9A", "CD1E", "CCR2", "IL1B"), reduction = "umap", cols = c("lightgray", "blue"))
# Plot cDC2 markers
FeaturePlot(myeloid_resubset, features = c("CD68", "CD14", "FCGR3A"), reduction = "umap", cols = c("lightgray", "blue"))
# Plot pDC markers
FeaturePlot(myeloid_resubset, features = c("CD123", "IRF7", "GZMB", "TCF4"), reduction = "umap", cols = c("lightgray", "blue"))


