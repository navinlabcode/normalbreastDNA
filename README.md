# normalbreastDNA
This repository contains the scripts used in the manuscript: *Normal Breast Tissues Harbor Rare Populations of Aneuploid Epithelial Cells*.


#### scDNA data analysis 
- Demultiplexed fastq files were first pre-processed by the pipeline: [_CNV_pipeline_](https://github.com/navinlabcode/CNV_pipeline).
- Then R scripts in _scDNA_analysis_ were used for downstream analysis and visualization. 

#### Statistical analysis 
- R scripts in _stat_analysis_ were used for statistical testing and visualization.

#### Spatial transcriptomic analysis
- 

### _Dependencies_
------------
Session info:
```
R version 4.1.2 (2021-11-01)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblasp-r0.3.3.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] ape_5.7                     ggtree_3.2.1                ggpubr_0.4.0                magrittr_2.0.3             
 [5] BiocParallel_1.28.3         glue_1.6.2                  forcats_1.0.0               stringr_1.5.0              
 [9] dplyr_1.1.0                 purrr_1.0.1                 readr_2.1.4                 tidyr_1.3.0                
[13] tibble_3.1.8                ggplot2_3.4.1               tidyverse_1.3.1             copykit_0.1.2              
[17] DNAcopy_1.68.0              Rsubread_2.8.2              SingleCellExperiment_1.16.0 SummarizedExperiment_1.24.0
[21] Biobase_2.54.0              GenomicRanges_1.46.1        GenomeInfoDb_1.30.1         IRanges_2.28.0             
[25] S4Vectors_0.32.4            BiocGenerics_0.40.0         MatrixGenerics_1.6.0        matrixStats_0.63.0    
```

### _Data source_
------------
The original sequencing data from this study has been deposited to the Sequence Read Archive (SRA): [PRJNA1003661](https://dataview.ncbi.nlm.nih.gov/object/PRJNA1003661).

### _Contact_
------------
For any additional information, please [email](mailto:nnavin@mdanderson.org) corresponding author.

