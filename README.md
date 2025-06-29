
# Lignature

<!-- badges: start -->
<!-- badges: end -->

Lignature is a full-featured R package for analyzing cell-cell communication from bulk & single-cell RNA-seq data

<hr>

**Lignature: A Comprehensive Database of Ligand Signatures to predict cell-cell communication** <br />
Ligand-receptor interactions mediate intercellular communication, inducing transcriptional changes that regulate physiological and pathological processes. Ligand-induced transcriptomic signatures can be used to predict active ligands; however, the absence of a comprehensive set of ligand-response signatures has limited their practical application in predicting ligand-receptor interactions. To bridge this gap, we developed Lignature, a curated [database](https://github.com/yingxinac/LignatureData/) encompassing intracellular transcriptomic signatures for human ligands. Lignature compiles signatures from published transcriptomic datasets and established resources, generating both gene- and pathway-based signatures for each ligand. <br />
Using the Lignature database as a reference, we developed a companion computational tool Lignature that infers the ligands responsible for transcriptomic changes within receiver cells, and the corresponding cell-cell interaction networks between multiple cell types or clusters.

<hr>

<div  align="center">
<img src="Figures/fig1v3.png" width = "850" height = "500" alt="Lignature" align=center />
</div>

<hr>


### Installation

**Install devtools**
``` r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```

**Install dependencies**
``` r
install.packages(c("tidyverse", "Seurat", "RColorBrewer", "pheatmap", "lsa", "data.table", "circlize", "randomcoloR", "hrbrthemes", "ggrepel"))
bio_pkgs = c("fgsea", "graphite")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(bio_pkgs)
library(devtools)
devtools::install_github("cdesterke/geneconverter")
```

**Install Lignature**
``` r
devtools::install_github("yingxinac/Lignature")
```

<hr>

### Database
Lignature database is available at
[LignatureData](https://github.com/yingxinac/LignatureData/)

<hr>

### Using Lignature
Example:
[Perform Lignature analysis](vignettes/Lignature.md)
Example dataset can be downloaded at
[google drive](https://drive.google.com/drive/folders/15G2RAnpb5wBJqNuY5oEF-5SAYg_xtXo-?usp=drive_link)

