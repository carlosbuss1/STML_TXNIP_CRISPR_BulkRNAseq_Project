# STML TXNIP CRISPR Bulk-RNAseq Project

## Overview

This project investigates the role of Thioredoxin-Interacting Protein (TXNIP) in human pluripotent stem cells using **CRISPR/Cas12 genome editing**. 
The goal is to evaluate the impact of TXNIP deficiency on glucose metabolism and the differentiation of hepatocyte-like cells and insulin-producing islet-like aggregates.

**Reference:** Traini, Negueruela, (...) Buss et al., 2025.
https://stemcellres.biomedcentral.com/articles/10.1186/s13287-025-04314-5

------

## Key Objectives

1. **CRISPR/Cas12 Editing**: Target TXNIP in human pluripotent stem cells.

2. **Bulk RNA Sequencing**: Analyze transcriptomic changes following TXNIP knockdown.

3. Differential Expression Analysis (DEA)

   :

   - Compare conditions (e.g., WT vs KO under DMSO and Thapsigargin treatments).

4. Pathway Enrichment Analysis

   :

   - Identify biological pathways impacted by TXNIP knockdown using Gene Set Enrichment Analysis (GSEA).

------

## Analysis Workflow

1. **Data Preprocessing**:
   - Load and annotate raw count data (`TXNIP_raw_counts.csv`).
   - Normalize data using log-CPM transformation.
   - Filter genes based on expression and variability thresholds.
2. **Differential Expression Analysis (DEA)**:
   - Perform DEA using `limma`.
   - Define contrasts for various conditions:
     - WT Thaps vs DMSO
     - KO Thaps vs DMSO
     - KO vs WT (DMSO and Thaps)
3. **Gene Set Enrichment Analysis (GSEA)**:
   - Perform pathway analysis (GO, KEGG, and Reactome).
   - Separate analyses for upregulated and downregulated genes.
4. **Visualization**:
   - Generate dot plots for GSEA results.
   - Export DEA and GSEA outputs.

------

## Repository Structure

```
graphqlCopiar código.
├── Dockerfile                    # Docker configuration
├── combined_txnip_dea.R          # Main R script for DEA and GSEA
├── TXNIP_raw_counts.csv          # Raw count data (input)
├── DEA_results_WT_Thaps_vs_DMSO.csv  # Example DEA output
├── README.md                     # Project overview and instructions
```

------

## Docker Workflow

### Prerequisites

- Install Docker.

### Running the Project

1. **Build the Docker Image**:

   ```
   bash
   
   
   Copiar código
   docker build -t txnip_project .
   ```

2. **Run the Analysis**:

   ```
   bash
   
   
   Copiar código
   docker run -v $(pwd):/app -w /app txnip_project
   ```

3. **Outputs**:

   - DEA results (`DEA_results_*.csv`) and GSEA visualizations (`*.png`) will be generated in the current directory.

------

## Citation

If you use this repository or any part of the analysis in your research, please cite:
**Traini, L., Negueruela, J., et al. (2025)**. CRISPR/Cas12 genome editing of TXNIP in human pluripotent stem cells to generate hepatocyte-like cells and insulin-producing islet-like aggregates. *Manuscript in preparation*.

------

## Contact

For inquiries, please contact:
**Carlos Buss** - Bioinformatician at STML laboratory - ULB Erasme Campus, Brussels, Belgium.
GitHub: [@carlosbuss1](https://github.com/carlosbuss1)
