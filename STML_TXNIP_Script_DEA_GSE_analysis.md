# Load Required Libraries
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
library(ggnewscale)
library(biomaRt)

# Set output directory
output_dir <- "./"

# Load Input Data
data <- fread("TXNIP_raw_counts.csv", header = TRUE)

# Extract 'Gene-ID' column containing Ensembl IDs for annotation
ensembl_ids <- unique(data$`Gene-ID`)

# Set Up Ensembl Mart for Human Genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Convert Ensembl IDs to Gene Symbols and Entrez IDs using biomaRt
conversion_results <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "entrezgene_id"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# Ensure there are no duplicate Ensembl IDs in conversion results
conversion_results <- conversion_results[!duplicated(conversion_results$ensembl_gene_id), ]

# Merge conversion results back with the original data
data <- merge(data, conversion_results, by.x = "Gene-ID", by.y = "ensembl_gene_id", all.x = TRUE)

# Remove rows with missing gene symbols
data <- data[!is.na(data$external_gene_name), ]

# Prepare data for aggregation
data_for_aggregation <- data[, c(
  "external_gene_name", 
  names(data)[!(names(data) %in% c("Gene-ID", "entrezgene_id", "Gene", "Gene-Description", "external_gene_name"))]
)]

# Ensure all selected columns are numeric
if (!all(sapply(data_for_aggregation[, -1], is.numeric))) {
  stop("Not all selected columns are numeric. Check your data.")
}

# Aggregate by Mean Expression for Duplicate Gene Symbols
data_aggregated <- aggregate(. ~ external_gene_name, 
                             data = data_for_aggregation, 
                             FUN = mean)

# Set row names to gene symbols
rownames(data_aggregated) <- data_aggregated$external_gene_name
data_aggregated$external_gene_name <- NULL

# Pre-process Data
## Remove Low Counts Across All Samples
min_counts <- 10
data_filtered <- data_aggregated[rowSums(data_aggregated) >= min_counts, ]

## Remove Genes with Low Variation
variances <- apply(data_filtered, 1, var)
var_threshold <- quantile(variances, 0.10)
data_filtered <- data_filtered[variances > var_threshold, ]

## Remove Genes with Too Many Zeros
data_filtered <- data_filtered[rowSums(data_filtered == 0) <= (0.5 * ncol(data_filtered)), ]

# Create Sample Metadata
sample_info <- data.frame(
  Batch = factor(rep(c("Rep1", "Rep2", "Rep3"), each = 4)),  # Replicate batches
  Group = factor(rep(c("WT-DMSO", "WT-Thaps", "KO-DMSO", "KO-Thaps"), times = 3))  # Replicate groups
)
rownames(sample_info) <- colnames(data_filtered)

# Print Sample Info
print(sample_info)

# Normalize Data
## Convert raw counts to log-CPM
dge <- DGEList(counts = data_filtered)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)

# Save normalized data
write.csv(logCPM, "normalized_logCPM_data.csv", row.names = TRUE)

# Create design matrix for limma
design <- model.matrix(~ 0 + Group + Batch, data = sample_info)
colnames(design) <- make.names(gsub("Group", "", colnames(design)))

# Fit linear model
fit <- lmFit(logCPM, design)

# Define contrasts
contrast_matrix <- makeContrasts(
  TXNIP_WT_Thaps_vs_DMSO = WT.Thaps - WT.DMSO,
  TXNIP_KO_Thaps_vs_DMSO = KO.Thaps - KO.DMSO,
  TXNIP_KO_vs_WT_DMSO = KO.DMSO - WT.DMSO,
  TXNIP_KO_vs_WT_Thaps = KO.Thaps - WT.Thaps,
  levels = design
)

# Fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract and save DEA results
## WT Thaps vs DMSO
results_WT_Thaps_vs_DMSO <- topTable(fit2, coef = "TXNIP_WT_Thaps_vs_DMSO", adjust.method = "BH", number = Inf)
results_WT_Thaps_vs_DMSO$Gene <- rownames(results_WT_Thaps_vs_DMSO)
write.csv(results_WT_Thaps_vs_DMSO, file.path(output_dir, "DEA_results_WT_Thaps_vs_DMSO.csv"), row.names = TRUE)

## KO Thaps vs DMSO
results_KO_Thaps_vs_DMSO <- topTable(fit2, coef = "TXNIP_KO_Thaps_vs_DMSO", adjust.method = "BH", number = Inf)
results_KO_Thaps_vs_DMSO$Gene <- rownames(results_KO_Thaps_vs_DMSO)
write.csv(results_KO_Thaps_vs_DMSO, file.path(output_dir, "DEA_results_KO_Thaps_vs_DMSO.csv"), row.names = TRUE)

## KO vs WT DMSO
results_KO_vs_WT_DMSO <- topTable(fit2, coef = "TXNIP_KO_vs_WT_DMSO", adjust.method = "BH", number = Inf)
results_KO_vs_WT_DMSO$Gene <- rownames(results_KO_vs_WT_DMSO)
write.csv(results_KO_vs_WT_DMSO, file.path(output_dir, "DEA_results_KO_vs_WT_DMSO.csv"), row.names = TRUE)

## KO vs WT Thaps
results_KO_vs_WT_Thaps <- topTable(fit2, coef = "TXNIP_KO_vs_WT_Thaps", adjust.method = "BH", number = Inf)
results_KO_vs_WT_Thaps$Gene <- rownames(results_KO_vs_WT_Thaps)
write.csv(results_KO_vs_WT_Thaps, file.path(output_dir, "DEA_results_KO_vs_WT_Thaps.csv"), row.names = TRUE)

# Display top 50 genes for each contrast
print("Top 50 Genes for WT Thaps vs DMSO:")
print(head(results_WT_Thaps_vs_DMSO[order(results_WT_Thaps_vs_DMSO$P.Value), ], 50))

print("Top 50 Genes for KO Thaps vs DMSO:")
print(head(results_KO_Thaps_vs_DMSO[order(results_KO_Thaps_vs_DMSO$P.Value), ], 50))

print("Top 50 Genes for KO vs WT DMSO:")
print(head(results_KO_vs_WT_DMSO[order(results_KO_vs_WT_DMSO$P.Value), ], 50))

print("Top 50 Genes for KO vs WT Thaps:")
print(head(results_KO_vs_WT_Thaps[order(results_KO_vs_WT_Thaps$P.Value), ], 50))

# Enrichment Analysis GSE

#
#

# First Step (Down and Up together)

#

# Load necessary libraries

library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)

# Set outputs directory

output_dir <- "~/Desktop/carlos/leonardo/REFRESH"

# WT GO Plot

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_GSEA_GO.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE GO - WT Thaps vs DMSO",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_GSEA_GO.png"), plot = gsea_plot, width = 8.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)

# KO GO Plot

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_GSEA_GO.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE GO - KO Thaps vs DMSO",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_GSEA_GO.png"), plot = gsea_plot, width = 8.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)
#
#
#
#
#

# WT KEGG Plot

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_GSEA_KEGG.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE KEGG - WT Thaps vs DMSO",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_GSEA_KEGG.png"), plot = gsea_plot, width = 6.7, height = 3.8, dpi = 600, units = "in")
print(gsea_plot)

# KO KEGG Plot

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_GSEA_KEGG.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE KEGG - KO Thaps vs DMSO",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_GSEA_KEGG.png"), plot = gsea_plot, width = 6.5, height = 3.9, dpi = 600, units = "in")
print(gsea_plot)
#
#
#

# WT Reactome Plot

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_GSEA_Reactome.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE Reactome - WT Thaps vs DMSO",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_GSEA_Reactome.png"), plot = gsea_plot, width = 7.7, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)

# KO Reactome Plot

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_GSEA_Reactome.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE Reactome - KO Thaps vs DMSO",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_GSEA_Reactome.png"), plot = gsea_plot, width = 7, height = 3.7, dpi = 600, units = "in")
print(gsea_plot)
#
#
#

#  2nd Step - Enrichment Analysis separated in two datasets (Down and UP)

#
#
#

# Separate upregulated and downregulated genes using logFC thresholds

upregulated_WT <- results_WT_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_WT <- results_WT_Thaps_vs_DMSO %>% filter(logFC < -1)
#
upregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC < -1)
#
#
#

# Rank genes for GSEA

ranked_up_WT <- rank_genes(upregulated_WT)
ranked_down_WT <- rank_genes(downregulated_WT)

ranked_up_KO <- rank_genes(upregulated_KO)
ranked_down_KO <- rank_genes(downregulated_KO)

# Convert to Entrez IDs

ranked_up_WT_entrez <- convert_to_entrez(ranked_up_WT)
ranked_down_WT_entrez <- convert_to_entrez(ranked_down_WT)

ranked_up_KO_entrez <- convert_to_entrez(ranked_up_KO)
ranked_down_KO_entrez <- convert_to_entrez(ranked_down_KO)
#
#
#

# Perform GSEA for WT Upregulated Genes

perform_gsea(ranked_up_WT_entrez, "WT_Thaps_vs_DMSO_Up")

# Perform GSEA for WT Downregulated Genes

perform_gsea(ranked_down_WT_entrez, "WT_Thaps_vs_DMSO_Down")

# Perform GSEA for KO Upregulated Genes

perform_gsea(ranked_up_KO_entrez, "KO_Thaps_vs_DMSO_Up")

# Perform GSEA for KO Downregulated Genes

perform_gsea(ranked_down_KO_entrez, "KO_Thaps_vs_DMSO_Down")
#
#
#

# WT Upregulated GO Plot

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_Up_GSEA_GO.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE GO - WT Upregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_Up_GSEA_GO.png"), plot = gsea_plot, width = 6.5, height = 3.6, dpi = 600, units = "in")
print(gsea_plot)

# Repeat the same for Downregulated GO genes

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_Down_GSEA_GO.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE GO - WT Downregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_Down_GSEA_GO.png"), plot = gsea_plot, width = 6.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)


# #

# Separate upregulated and downregulated genes for KO Thaps vs DMSO

upregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC < -1)

# Rank upregulated genes for GSEA

ranked_up_KO <- rank_genes(upregulated_KO)
ranked_up_KO_entrez <- convert_to_entrez(ranked_up_KO)

# Perform GSEA for upregulated genes

perform_gsea(ranked_up_KO_entrez, "KO_Thaps_vs_DMSO_Up")

# Rank downregulated genes for GSEA

ranked_down_KO <- rank_genes(downregulated_KO)
ranked_down_KO_entrez <- convert_to_entrez(ranked_down_KO)

# Perform GSEA for downregulated genes

perform_gsea(ranked_down_KO_entrez, "KO_Thaps_vs_DMSO_Down")

# Plot GSEA results for upregulated genes (GO)

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_Up_GSEA_GO.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE GO - KO Upregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_Up_GSEA_GO.png"), plot = gsea_plot, width = 8, height = 4, dpi = 600, units = "in")
print(gsea_plot)

# Plot GSEA results for downregulated genes (GO)

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_Down_GSEA_GO.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE GO - KO Downregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_Down_GSEA_GO.png"), plot = gsea_plot, width = 7.5, height = 4, dpi = 600, units = "in")
print(gsea_plot)


# #

# Separate upregulated and downregulated genes for WT Thaps vs DMSO

upregulated_WT <- results_WT_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_WT <- results_WT_Thaps_vs_DMSO %>% filter(logFC < -1)

# Rank upregulated genes for GSEA

ranked_up_WT <- rank_genes(upregulated_WT)
ranked_up_WT_entrez <- convert_to_entrez(ranked_up_WT)

# Perform GSEA for upregulated genes

perform_gsea(ranked_up_WT_entrez, "WT_Thaps_vs_DMSO_Up")

# Rank downregulated genes for GSEA

ranked_down_WT <- rank_genes(downregulated_WT)
ranked_down_WT_entrez <- convert_to_entrez(ranked_down_WT)

# Perform GSEA for downregulated genes

perform_gsea(ranked_down_WT_entrez, "WT_Thaps_vs_DMSO_Down")

# Plot GSEA results for upregulated genes (KEGG)

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_Up_GSEA_KEGG.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE KEGG - WT Upregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_Up_GSEA_KEGG.png"), plot = gsea_plot, width = 6.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)

# Plot GSEA results for downregulated genes (KEGG)

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_Down_GSEA_KEGG.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE KEGG - WT Downregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_Down_GSEA_KEGG.png"), plot = gsea_plot, width = 8.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)
#
#
#

# Separate upregulated and downregulated genes for KO Thaps vs DMSO

upregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC < -1)

# Rank upregulated genes for GSEA

ranked_up_KO <- rank_genes(upregulated_KO)
ranked_up_KO_entrez <- convert_to_entrez(ranked_up_KO)

# Perform GSEA for upregulated genes

perform_gsea(ranked_up_KO_entrez, "KO_Thaps_vs_DMSO_Up")

# Rank downregulated genes for GSEA

ranked_down_KO <- rank_genes(downregulated_KO)
ranked_down_KO_entrez <- convert_to_entrez(ranked_down_KO)

# Perform GSEA for downregulated genes

perform_gsea(ranked_down_KO_entrez, "KO_Thaps_vs_DMSO_Down")

# Plot GSEA results for upregulated genes (KEGG)

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_Up_GSEA_KEGG.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE KEGG - KO Upregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_Up_GSEA_KEGG.png"), plot = gsea_plot, width = 6.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)

# Plot GSEA results for downregulated genes (KEGG)

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_Down_GSEA_KEGG.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE KEGG - KO Downregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_Down_GSEA_KEGG.png"), plot = gsea_plot, width = 8.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)


# #

# Separate upregulated and downregulated genes for KO Thaps vs DMSO

upregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_KO <- results_KO_Thaps_vs_DMSO %>% filter(logFC < -1)

# Rank upregulated genes for GSEA

ranked_up_KO <- rank_genes(upregulated_KO)
ranked_up_KO_entrez <- convert_to_entrez(ranked_up_KO)

# Perform GSEA for upregulated genes

perform_gsea(ranked_up_KO_entrez, "KO_Thaps_vs_DMSO_Up")

# Rank downregulated genes for GSEA

ranked_down_KO <- rank_genes(downregulated_KO)
ranked_down_KO_entrez <- convert_to_entrez(ranked_down_KO)

# Perform GSEA for downregulated genes

perform_gsea(ranked_down_KO_entrez, "KO_Thaps_vs_DMSO_Down")

# Plot GSEA results for upregulated genes (Reactome)

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_Up_GSEA_Reactome.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE Reactome - KO Upregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_Up_GSEA_Reactome.png"), plot = gsea_plot, width = 6.4, height = 4, dpi = 600, units = "in")
print(gsea_plot)

# Plot GSEA results for downregulated genes (Reactome)

gsea_data <- read.csv(file.path(output_dir, "KO_Thaps_vs_DMSO_Down_GSEA_Reactome.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE Reactome - KO Downregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "KO_Thaps_vs_DMSO_Down_GSEA_Reactome.png"), plot = gsea_plot, width = 8.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)
#
#
#

# Separate upregulated and downregulated genes for WT Thaps vs DMSO

upregulated_WT <- results_WT_Thaps_vs_DMSO %>% filter(logFC > 1)
downregulated_WT <- results_WT_Thaps_vs_DMSO %>% filter(logFC < -1)

# Rank upregulated genes for GSEA

ranked_up_WT <- rank_genes(upregulated_WT)
ranked_up_WT_entrez <- convert_to_entrez(ranked_up_WT)

# Perform GSEA for upregulated genes

perform_gsea(ranked_up_WT_entrez, "WT_Thaps_vs_DMSO_Up")

# Rank downregulated genes for GSEA

ranked_down_WT <- rank_genes(downregulated_WT)
ranked_down_WT_entrez <- convert_to_entrez(ranked_down_WT)

# Perform GSEA for downregulated genes

perform_gsea(ranked_down_WT_entrez, "WT_Thaps_vs_DMSO_Down")

# Plot GSEA results for upregulated genes (Reactome)

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_Up_GSEA_Reactome.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE Reactome - WT Upregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_Up_GSEA_Reactome.png"), plot = gsea_plot, width = 6.4, height = 4, dpi = 600, units = "in")
print(gsea_plot)

# Plot GSEA results for downregulated genes (Reactome)

gsea_data <- read.csv(file.path(output_dir, "WT_Thaps_vs_DMSO_Down_GSEA_Reactome.csv"), stringsAsFactors = FALSE)
gsea_data <- gsea_data %>% filter(!is.na(NES) & !is.na(p.adjust))
gsea_data$EnrichmentFold <- gsea_data$setSize * abs(gsea_data$enrichmentScore)  # Calculate Enrichment Fold
top10_gsea <- gsea_data %>% arrange(p.adjust) %>% slice_head(n = 10)
top10_gsea$Description <- str_to_title(top10_gsea$Description)
gsea_plot <- ggplot(top10_gsea, aes(x = EnrichmentFold, y = reorder(Description, EnrichmentFold))) +
  geom_point(aes(size = setSize, color = -log10(p.adjust)), shape = 16) +
  scale_color_gradient(low = "blue", high = "red", name = "-log10(adj-pvalue)") +
  scale_size_continuous(name = "Gene Count", range = c(2, 4)) +
  theme_minimal() +
  labs(
    title = "GSE Reactome - WT Downregulated (Thaps vs DMSO)",
    x = "Enrichment Fold",
    y = "Pathways"
  ) +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    axis.text.y = element_text(size = 8, color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  guides(
    color = guide_colorbar(title = "-log10(adj-pvalue)", order = 1),
    size = guide_legend(title = "Gene Count", order = 2)
  )
ggsave(filename = file.path(output_dir, "WT_Thaps_vs_DMSO_Down_GSEA_Reactome.png"), plot = gsea_plot, width = 8.5, height = 3.5, dpi = 600, units = "in")
print(gsea_plot)
#
#
#





