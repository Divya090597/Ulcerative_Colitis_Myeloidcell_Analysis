library(ComplexHeatmap)
library(dplyr)
library(readr)
library(here)
library(DESeq2)
library(ggplot2)


# Set Idents 
Idents(myeloid_resubset) <- "fine_annotations_MK_V4"

# Check the Idents
levels(Idents(myeloid_resubset))

# Differential expression between two clusters
dge_results <- FindMarkers(myeloid_resubset, ident.1 = "08B-MNP monocyte")

# View the results
head(dge_results)

# Define thresholds for significance
logfc_threshold <- 0.25  # Log2 fold-change threshold
pval_threshold <- 0.05   # P-value threshold

# Add a column for significance
dge_results$Significance <- with(dge_results, ifelse(p_val_adj < pval_threshold & abs(avg_log2FC) > logfc_threshold, 
                                                     ifelse(avg_log2FC > 0, "Upregulated", "Downregulated"), 
                                                     "Not Significant"))



# Create the volcano plot
volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot: 08A-MNP macrophage vs 08B-MNP monocyte",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

# Display the plot
print(volcano_plot)

# Sort by adjusted p-value and get top genes for each category
top_upregulated <- dge_results %>%
  filter(Significance == "Upregulated") %>%
  arrange(p_val_adj) %>%
  head(5)  # Top 5 upregulated genes

top_downregulated <- dge_results %>%
  filter(Significance == "Downregulated") %>%
  arrange(p_val_adj) %>%
  head(5)  # Top 5 downregulated genes

# Combine top genes into a single list
top_genes <- c(rownames(top_upregulated), rownames(top_downregulated))

# Add gene labels to the plot
library(ggrepel)


# Create a subset of the data containing only the top genes
label_data <- dge_results[rownames(dge_results) %in% top_genes, ]

# Ensure that label_data has rownames as a column for use in ggplot
label_data$gene <- rownames(label_data)

# Add labels for the top genes
volcano_plot <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_text_repel(data = label_data,
                  aes(label = gene),
                  size = 3, max.overlaps = 10) +
  labs(title = "Volcano Plot for MNP-monocytes",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

# Display the plot
print(volcano_plot)
