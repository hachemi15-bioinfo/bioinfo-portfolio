# ==============================================================================
# Title: Differential Gene Expression Analysis with DESeq2
# Author: Bioinfo Portfolio
# Description: Identification of DEGs from RNA-Seq dataset using DESeq2
# ==============================================================================

# 1. Load Libraries
library(DESeq2)
library(tidyverse)
library(pheatmap)
library(EnhancedVolcano)

# 2. Load Dataset (Airway dataset)
if (!requireNamespace("airway", quietly = TRUE)) {
    BiocManager::install("airway")
}
library(airway)

data("airway")
se <- airway

# Extract count matrix and metadata
count_matrix <- assay(se)
sample_metadata <- as.data.frame(colData(se))

# 3. Build DESeq2 Object
dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = ~ dex
)

# Filter low-count genes (row sum >= 10)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 4. Differential Expression Analysis
dds <- DESeq(dds)
res <- results(dds)

# 5. Format and Export Results
res_df <- as.data.frame(res) %>% 
    rownames_to_column(var = "gene_id") %>% 
    filter(!is.na(padj))

# Save results to the results directory
write.csv(res_df, "../results/differential_expression_results.csv", row.names = FALSE)

print("Analysis complete. Results exported to results/differential_expression_results.csv")