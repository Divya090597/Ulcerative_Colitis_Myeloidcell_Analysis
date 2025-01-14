
library(Seurat)
library(affy)
library(R.utils)
library(hgu133plus2.db)

# Install the annotation package
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)


# Load Data ---------------------------------------------------------------

seu_obj = readRDS("~/data/HECTD3/UC_biopsy.rds")

untar("~/data/HECTD3/GSE38713_RAW.tar", exdir = "GSE38713_RAW")

setwd("GSE38713_RAW")
data <- ReadAffy()

summary(data)

normalized_data <- rma(data)

# View the normalized expression matrix
expr_matrix <- exprs(normalized_data)
head(expr_matrix)
dim(expr_matrix)

# Extract probe IDs from the expression matrix
probe_ids <- rownames(expr_matrix)

# Map probe IDs to gene symbols
gene_symbols <- mapIds(hgu133plus2.db, 
                       keys = probe_ids, 
                       column = "SYMBOL", 
                       keytype = "PROBEID", 
                       multiVals = "first")

# Add gene symbols as a new column
expr_matrix_with_symbols <- data.frame(GeneSymbol = gene_symbols, expr_matrix[, 1])

# View the first few rows
head(expr_matrix_with_symbols)

# Define genes of interest
genes_of_interest <- c("IL6", "IL1B","DPP4", "FCGR1A", "ITGAM", "ITGAX", "HECTD3") 

# Extract expression values for these genes
expr_values <- expr_matrix_with_symbols[expr_matrix_with_symbols$GeneSymbol %in% genes_of_interest, ]

# Print the expression values
print(expr_values)

# Bar plot for selected genes
barplot(as.numeric(expr_values$expr_matrix...1.),
        names.arg = expr_values$GeneSymbol,
        main = "Expression Levels",
        ylab = "Expression Level",
        col = "red",
        las = 2)  # Rotate axis labels

# Rename the column for clarity
colnames(expr_values)[2] <- "Expression"

# Verify the data frame
print(expr_values)


# Aggregate expression values by GeneSymbol
expr_aggregated <- aggregate(Expression ~ GeneSymbol, data = expr_values, FUN = mean)

# Bar plot with aggregated values
barplot(expr_aggregated$Expression,
        names.arg = expr_aggregated$GeneSymbol,
        main = "Expression Levels(Aggregated)",
        ylab = "Expression Level",
        col = "red",
        las = 2)



