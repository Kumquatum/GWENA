# GWENA

## Overview
**GWENA** (Gene Whole co-Expression Network Analysis) is an R package to perform gene co-expression network analysis in a single pipeline. This pipeline includes functional enrichment of modules of co-expressed genes, phenotypcal association, topological analysis and comparisons of networks between conditions.

> ![Figure 1. Analysis pipeline of GWENA, from expression data to characterization of the modules and comparison of conditions.](vignettes/figure_pipeline_schema.png)  

Using transcriptomics data from either RNA-seq  or microarray, the package follows the steps displayed in Figure 1:

1. **Input**: data is provided as a data.frame or a matrix of expression intensities (pre-normalized).
2. **Gene filtering**: data is filtered according to the transcriptomic technology used.
3. **Network building**: a matrix of similarity score is computed between each gene with Spearman correlation, then transformed into an adjacency matrix, and finally into a topological overlap matrix.
4. **Modules detection**: groups of genes with closest similarity scores are detected as modules.
5. **Biological integration**: gene set enrichment analysis and phenotypic association (if phenotypes are provided) are performed on modules.
6. **Graph and topological analysis**: identification of hub genes is made available, as well as modules visualization.
7. **Networks comparison**: if multiple conditions are available (time points, treatments, phenotype, etc.), analysis of modules preservation/non-preservation between conditions can be performed.


## Installation
```R
# Prerequisites
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly=TRUE))
  install.packages("devtools")
```

Installation can either be from:

1. ~~the official version of the last Bioconductor release (recommended).~~ Not yet available
```R
# BiocManager::install("GWENA") # Not yet available
```
2. the last stable version from the Bioc Devel branch.
```R
# 2. From Bioconductor devel
BiocManager::install("GWENA", version = "devel")
```
3. the day-to-day development version from the [Github repository](https://github.com/Kumquatum/GWENA).
```R
# 3. From Github repository
BiocManager::install("Kumquatum/GWENA")
# OR
devtools::install_github("Kumquatum/GWENA")
```

## Package tutorial
A vignette is available at [vignette/GWENA_guide.Rmd](vignette/GWENA_guide.Rmd). To see the html version, use `vignette("GWENA_guide")`.
Note: if you want to install from Github, please specify `build_vignettes = TRUE` in `install_github`.


## Licence and use
GWENA is a free software; you can redistribute it and/or modify it under the terms of the [GNU General Public License 3 (GPL 3)](LICENSE.md) as published by the Free Software Foundation. 
It is distributed in the hope that it will be useful, but without any warranty and without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details.


## Contact
* If you have any question, feel free to send a mail to the author and maintainer Gwenaëlle Lemoine at lemoine.gwenaelle[@t)gmail{d0t]com.
* To report a bug, please use the [Github issues system](https://github.com/Kumquatum/GWENA/issues).
