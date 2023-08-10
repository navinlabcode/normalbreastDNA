library(copykit)
library(SummarizedExperiment)
library(magrittr)
library(BiocParallel)

register(MulticoreParam(workers = 50, progressbar = F),default = T)

bcmhbca_41l<-runVarbin('./correct_bam_input/BCMHBCA41L_bam/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca_47r<-runVarbin('./BCMHBCA_47R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca01l<-runVarbin('./BCMHBCA01L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca04r<-runVarbin('./BCMHBCA04R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca09l<-runVarbin('./BCMHBCA09L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca10r<-runVarbin('./BCMHBCA10R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca15l<-runVarbin('./BCMHBCA15L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca21l<-runVarbin('.BCMHBCA21L/PM2429/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca22l<-runVarbin('./BCMHBCA22L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca24r<-runVarbin('./BCMHBCA24R/PM2429/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca35l<-runVarbin('./BCMHBCA35L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca51r<-runVarbin('./correct_bam_input/BCMHBCA51R_bam/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca60l<-runVarbin('./correct_bam_input/BCMHBCA60L_bam/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca63r<-runVarbin('./BCMHBCA63R/renamed_output/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca66l<-runVarbin('./correct_bam_input/BCMHBCA66L_bam/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca68r<-runVarbin('./BCMHBCA68R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca71r<-runVarbin('./BCMHBCA71R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca72r<-runVarbin('/./BCMHBCA72R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca75l_1<-runVarbin('/./BCMHBCA75L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca75l_2<-runVarbin('./correct_bam_input/BCMHBCA75L_2_bam/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca82r<-runVarbin('./BCMHBCA82R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca83l<-runVarbin('./BCMHBCA83L/PM2408/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca85l<-runVarbin('./BCMHBCA85L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca89l<-runVarbin('./BCMHBCA89L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca03l<-runVarbin('./BCMHBCA/BCMHBCA03L/output_PM2453/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca06r<-runVarbin('./BCMHBCA06R/epi_selection/output/BCMHBCA06R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca08r<-runVarbin('./BCMHBCA08R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca20r<-runVarbin('./BCMHBCA20R/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca27l<-runVarbin('./BCMHBCA27L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca36l<-runVarbin('./BCMHBCA36L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca42l<-runVarbin('./BCMHBCA42L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=F,min_bincount=1)
bcmhbca44l<-runVarbin('./BCMHBCA44L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)
bcmhbca67l<-runVarbin('./BCMHBCA67L/output/sort/',genome='hg19',remove_Y=T,is_paired_end=T,min_bincount=1)

cbs_list <- list(bcmhbca_41l,
                 bcmhbca_47r,
                 bcmhbca15l,
                 bcmhbca10r, 
                 bcmhbca66l, 
                 bcmhbca35l, 
                 bcmhbca82r, 
                 bcmhbca71r, 
                 bcmhbca51r, 
                 bcmhbca60l, 
                 bcmhbca22l, 
                 bcmhbca01l, 
                 bcmhbca63r, 
                 bcmhbca83l, 
                 bcmhbca04r, 
                 bcmhbca09l, 
                 bcmhbca75l_1,
                 bcmhbca75l_2,
                 bcmhbca72r,
                 bcmhbca24r,
                 bcmhbca21l,
                 bcmhbca68r,
                 bcmhbca85l,
                 bcmhbca89l,
                 bcmhbca03l,
                 bcmhbca06r,
                 bcmhbca08r,
                 bcmhbca20r,
                 bcmhbca27l,
                 bcmhbca36l,
                 bcmhbca42l,
                 bcmhbca44l,
                 bcmhbca67l)

merged_ckobj_cbs <- do.call(cbind, cbs_list)
saveRDS(merged_ckobj_cbs, file = "./cbs_0501_v1/32pts_merged_ckobj_cbs.RDS")

## Rscript Step0_gather_samples_cbs.R 2> log
