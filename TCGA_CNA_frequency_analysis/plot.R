setwd("./published_data/TCGA/IBC_v2/")
source("SWITCHplusforFrequencyplots.R")

plot.freq.SW(data = "BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt",
             grps = "sample_cat.tsv",
             species = "chr_length_hg19.txt",
             threshold = 0.3,
             returnSummary = T)
