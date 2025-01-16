library(ComplexHeatmap)
library(dplyr)
library(readr)
library(here)
library(DESeq2)
library(ggplot2)


# Set Idents 
Idents(integrated_data) <- "integrated_snn_res.0.5"

# Check the Idents
levels(Idents(integrated_data))

# Differential expression 
gene_results <- FindMarkers(integrated_data, ident.1 = "1")

# View the results
head(gene_results)

# Define thresholds for significance
logfc_threshold <- 0.25  # Log2 fold-change threshold
pval_threshold <- 0.05   # P-value threshold

# Add a column for significance
gene_results$Significance <- with(gene_results, ifelse(p_val_adj < pval_threshold & abs(avg_log2FC) > logfc_threshold, 
                                                     ifelse(avg_log2FC > 0, "Upregulated", "Downregulated"), 
                                                     "Not Significant"))



# Create the volcano plot
volcano_plot <- ggplot(gene_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  labs(title = "Volcano Plot: cluster8",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

# Display the plot
print(volcano_plot)

# Sort by adjusted p-value and get top genes for each category
top_upregulated <- gene_results %>%
  filter(Significance == "Upregulated") %>%
  arrange(p_val_adj) %>%
  head(15)  # Top 5 upregulated genes

top_downregulated <- gene_results %>%
  filter(Significance == "Downregulated") %>%
  arrange(p_val_adj) %>%
  head(15)  # Top 5 downregulated genes

# Combine top genes into a single list
top_genes <- c(rownames(top_upregulated), rownames(top_downregulated))

# Add gene labels to the plot
library(ggrepel)


# Create a subset of the data containing only the top genes
label_data <- gene_results[rownames(gene_results) %in% top_genes, ]

# Ensure that label_data has rownames as a column for use in ggplot
label_data$gene <- rownames(label_data)

# Add labels for the top genes
volcano_plot <- ggplot(gene_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  geom_text_repel(data = label_data,
                  aes(label = gene),
                  size = 3, max.overlaps = 10) +
  labs(title = "Volcano Plot for cluster8",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

# Display the plot
print(volcano_plot)
