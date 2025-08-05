# ---- Setup chunk: Configure knitr options ----
# (In R scripts, this is typically unnecessary unless using knitr for reports)
knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

# ---- Load Required Libraries ----
# Install required Bioconductor and CRAN packages if not already installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("GO.db")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("biomaRt")
#BiocManager::install("ReactomePA")
#BiocManager::install("enrichplot")
#BiocManager::install("clusterProfiler")

# Load libraries
library(GO.db)
library(org.Hs.eg.db)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(tidyverse)

# ---- Load Datasets ----

# Load list of human protein-coding genes
protein_coding_genes <- read_delim(
  "C:/Users/HP-ssd/Desktop/Short term project/protein coding genes/gene_with_protein_product.txt", 
  delim = "\t", escape_double = FALSE, trim_ws = TRUE
)
protein_coding_genes_list <- protein_coding_genes$symbol

# Load FUSIL annotations (gene essentiality categories)
fusil_m_gene <- read_delim("C:/Users/HP-ssd/Desktop/Short term project2/fusil.csv")

# Load human gene paralogues
human_gene_paralogues <- read.csv("C:/Users/HP-ssd/Desktop/Short term project2/paralogues/human_gene_paralogues.csv")

# Clean and rename gene symbol column
human_gene_paralogues <- human_gene_paralogues %>%
  select(-1, -2, -4) %>%
  rename(gene_symbol = external_gene_name)

# ---- Prepare FUSIL-Paralogue Matrix ----

# Keep relevant columns from FUSIL file
Fusil_genes <- fusil_m_gene %>% select(-1, -2)

# Join paralogues with FUSIL annotations based on paralogue names
Fusil_genes_paralogues <- human_gene_paralogues %>%
  left_join(Fusil_genes, by = c("hsapiens_paralog_associated_gene_name" = "gene_symbol")) %>%
  filter(hsapiens_paralog_associated_gene_name %in% protein_coding_genes_list) %>%
  rename("fusil_paralogue" = "fusil")

# Join FUSIL annotations for the main gene too
Fusil_genes_paralogues <- Fusil_genes_paralogues %>%
  left_join(Fusil_genes, by = c("gene_symbol" = "gene_symbol"))

# Reorder columns and remove rows with missing values
Fusil_genes_paralogues <- Fusil_genes_paralogues %>%
  relocate(fusil, .after = "gene_symbol") %>%
  relocate(fusil_paralogue, .after = "hsapiens_paralog_associated_gene_name") %>%
  na.omit()

# Annotate whether gene and paralogue have matching FUSIL bins and their similarity bin
fusil_match <- Fusil_genes_paralogues %>%
  mutate(FUSIL_match = ifelse(fusil == fusil_paralogue, "Match", "Mismatch")) %>%
  mutate(SIMILARITY_bin = case_when( 
    hsapiens_paralog_perc_id >= 80 ~ "High >80% ",
    hsapiens_paralog_perc_id >= 60 ~ "Medium-High 60-80%",
    hsapiens_paralog_perc_id >= 40 ~ "Medium 40-60%",
    hsapiens_paralog_perc_id >= 20 ~ "Medium-Low 20-50%",
    TRUE ~ "Low <20%"))

# ---- GSE (Gene Set Enrichment) Analysis per FUSIL Category ----

# Setup connection to Ensembl BioMart
hs_mart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

# Define gene universe: all FUSIL-annotated genes
genes <- unique(fusil_match$gene_symbol)
gene_entrez_id <- getBM(
  attributes = c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id'),
  filters = 'hgnc_symbol',
  values = genes,
  mart = hs_mart
)

# Filter valid Entrez IDs for universe
reference_set_entrez <- unique(gene_entrez_id$entrezgene_id)
reference_set_entrez <- reference_set_entrez[!is.na(reference_set_entrez)]
reference_set_entrez <- as.character(reference_set_entrez)

# Prepare to store enrichment results for each FUSIL category
categories <- unique(fusil_match$fusil)
enrichment_results_list <- list()

# Loop through each FUSIL bin
for (cat in categories) {
  message(paste("Processing category:", cat))
  
  # Subset data for current category
  matching_subset <- fusil_match %>%
    filter(fusil == cat)
  
  gene_match <- unique(matching_subset$hsapiens_paralog_associated_gene_name)
  
  # Map paralogue symbols to Entrez IDs
  gene_mapping <- getBM(
    attributes = c('hgnc_symbol', 'ensembl_gene_id', 'entrezgene_id'),
    filters = 'hgnc_symbol',
    values = gene_match,
    mart = hs_mart
  )
  
  # Filter valid Entrez IDs
  gene_set_entrez <- unique(gene_mapping$entrezgene_id)
  gene_set_entrez <- gene_set_entrez[!is.na(gene_set_entrez)]
  gene_set_entrez <- gene_set_entrez[gene_set_entrez %in% reference_set_entrez]
  
  message(paste("Genes in test set:", length(gene_set_entrez)))
  
  # Skip if too few genes
  if (length(gene_set_entrez) < 5) {
    message(paste("Too few genes for enrichment in category:", cat, "- skipping."))
    next
  }
  
  # Perform GO enrichment using clusterProfiler
  enrichment_result <- enrichGO(
    gene = gene_set_entrez,
    universe = reference_set_entrez,
    OrgDb = org.Hs.eg.db,
    ont = "BP",  # Biological Process
    pAdjustMethod = "BH",
    readable = TRUE
  )
  
  # Skip empty enrichment
  if (is.null(enrichment_result) || nrow(as.data.frame(enrichment_result)) == 0) {
    message(paste("No enrichment found for category:", cat, "- skipping plot."))
    next
  }
  
  enrichment_results_list[[cat]] <- enrichment_result
  
  # Visualize enrichment
  p <- dotplot(enrichment_result, showCategory = 8) + 
    ggtitle(paste("GO BP Enrichment:", cat))
  
  message(paste("Analysis done for category:", cat))
  print(p)
}

# ---- Compare Enrichment Across All Categories ----

# Group paralogues by FUSIL bin
paralogue_lists <- fusil_match %>%
  group_by(fusil) %>%
  summarise(paralogues = list(unique(hsapiens_paralog_associated_gene_name))) %>%
  deframe()

# Convert gene symbols to Entrez IDs for each category
entrez_sets <- lapply(paralogue_lists, function(gene_symbols) {
  gene_mapping <- getBM(
    attributes = c('hgnc_symbol', 'entrezgene_id'),
    filters = 'hgnc_symbol',
    values = gene_symbols,
    mart = hs_mart
  )
  entrez_ids <- unique(gene_mapping$entrezgene_id)
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  entrez_ids[entrez_ids %in% reference_set_entrez]
})


# Remove categories with <5 genes
entrez_sets <- entrez_sets[sapply(entrez_sets, length) >= 5]

# Perform cluster-based GO enrichment comparison
compare_result <- compareCluster(
  geneCluster = entrez_sets,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  universe = reference_set_entrez,
  pAdjustMethod = "BH",
  readable = TRUE
)

# Visualize comparative enrichment across categories
dotplot(compare_result, showCategory = 4) + 
  ggtitle("GO BP Comparison Across Categories")



# Perform cluster-based MF enrichment comparison
compare_result <- compareCluster(
  geneCluster = entrez_sets,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  ont = "MF",
  universe = reference_set_entrez,
  pAdjustMethod = "BH",
  readable = TRUE
)

# Visualize comparative enrichment across categories
dotplot(compare_result, showCategory = 5) + 
  ggtitle("MF BP Comparison Across Categories")





