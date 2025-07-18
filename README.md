# READ ME

<a href="https://creativecommons.org">Untitled</a> © 1999 by <a href="https://creativecommons.org">Jane Doe</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">


# GSE Analysis by FUSIL Category and Paralogues

This repository contains an R Markdown (`.Rmd`) analysis pipeline for performing Gene Set Enrichment (GSE) using gene paralogues grouped by **FUSIL essentiality categories**. The analysis integrates multiple bioinformatics packages from Bioconductor and uses GO enrichment to compare biological processes across essentiality groups.

---

## Files

- `GSE_analysis_with_comments.Rmd` – Main annotated R Markdown script performing the entire analysis.

---

## Key Features

- ✅ Loads and processes protein-coding gene lists, paralogues, and FUSIL annotations.
- ✅ Integrates gene essentiality (FUSIL) and paralog identity via BioMart.
- ✅ Performs GO Biological Process (BP) enrichment using `clusterProfiler`.
- ✅ Compares enrichment across all FUSIL categories via `compareCluster`.
- ✅ Fully annotated and reproducible R Markdown for transparency.

---

## Required Packages

- `biomaRt`
- `clusterProfiler`
- `org.Hs.eg.db`
- `ReactomePA`
- `enrichplot`
- `GO.db`
- `tidyverse`

Install packages (if needed):

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("biomaRt", "clusterProfiler", "org.Hs.eg.db", "ReactomePA", "enrichplot", "GO.db"))
install.packages("tidyverse")


