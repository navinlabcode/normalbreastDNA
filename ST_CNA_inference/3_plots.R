
#load package
library(dplyr)
library(gsDesign)
library(readr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(copykit, lib.loc = "/opt/R/4.1.2/lib/R/library")#0.1.0
library(pwr)
library(Seurat)
source('./Copykit_copykat_functions.R')

svdir = '/volumes/USR1/yiyun/Project/HBCA/copykat_ST/'


#select patient
para_df = read_tsv('ST_sample_keep.tsv')

#set color
col = list(CNA_prediction = c('Aneuploid'='#E41A1B','Diploid' = '#4EAF4A'),
           pred_celltype = c('LumSec'='#E41A1C', 'LumHR' = '#43997A', 'Basal' = '#F781BF',
                             'Fibroblasts' = '#B7B7B7', 'Vascular' = '#BFCCB5',
                             'Pericytes' = '#E5BEEC', 'Lymphatic' = '#E4D0D0',
                             'T.cells' = '#B0DAFF','B.cells'='#FFD422','Myeloid' = '#98D8AA'))

#plot on ST H&E image: celltype, aneuploid prediction
if(T){
  plist = list()
  for(sample.name in para_df$sample.name){
    # sample.name = 'BCMHBCA52R1'
    
    copykat_dir = paste0(svdir,'/improve_copykat/',sample.name,'/');
    dir.create(copykat_dir);setwd(copykat_dir)
    
    seu_ST = read_rds(paste0(copykat_dir,sample.name,'_seu_ST.rds'))
    
    table(seu_ST@meta.data$CNA_prediction)
    is.na(seu_ST@meta.data$CNA_prediction) %>% sum()
    
    
    seu_ST$pred_celltype = seu_ST$pred.id
    col = list(CNA_prediction = c('Aneuploid'='#E41A1B','Diploid' = '#4EAF4A'),
               pred_celltype = c('LumSec'='#E41A1C', 'LumHR' = '#43997A', 'Basal' = '#F781BF',
                                 'Fibroblasts' = '#B7B7B7', 'Vascular' = '#BFCCB5',
                                 'Pericytes' = '#E5BEEC', 'Lymphatic' = '#E4D0D0',
                                 'T.cells' = '#B0DAFF','B.cells'='#FFD422','Myeloid' = '#98D8AA'))
    
    # if(sample.name == 'BCMHBCA04L'){RT = 1.1/1.3}
    # if(sample.name == 'BCMHBCA52R'){RT = 1.1/1.08}
    # if(sample.name == 'BCMHBCA60L1_061223'){RT = 1.1/0.98}
    RT = 1.1/1.05
    
    plist[['CNA_prediction']] = SpatialDimPlot(seu_ST, group.by = "CNA_prediction", pt.size.factor = 1,
                                               stroke = NA, combine = F)[[1]]+
      scale_fill_manual(breaks = names(col$CNA_prediction), values = col$CNA_prediction)+
      theme(aspect.ratio=RT);
    plist[['pred_celltype']] = SpatialDimPlot(seu_ST, group.by = "pred_celltype", pt.size.factor = 1,
                                              stroke = NA, combine = F)[[1]]+
      scale_fill_manual(breaks = names(col$pred_celltype), values = col$pred_celltype)+
      theme(aspect.ratio=RT);
    seu_ST$pred_celltype_sub = seu_ST$pred_celltype
    seu_ST$pred_celltype_sub[!grepl('Aneuploid',seu_ST$CNA_prediction)]=NA
    plist[['pred_celltype_sub']] = SpatialDimPlot(seu_ST, group.by = "pred_celltype_sub", pt.size.factor = 1,
                                                  stroke = NA, combine = F,alpha = c(NA))[[1]]+
      scale_fill_manual(breaks = names(col$pred_celltype), values = col$pred_celltype)+
      theme(aspect.ratio=RT);
    plot = lapply(plist, function(x) x+theme(legend.position = 'None'))
    lg = lapply(plist, function(x) get_legend(x))
    
    pdf(paste0(sample.name, '_copykat_map_HE.pdf'),width = 8,height = 8)
    plot_grid(plot_grid(plotlist = plot,ncol=1), plot_grid(plotlist = lg,ncol=1)) %>% print()
    dev.off()
  }
}

#plot copy number heatmap: celltype, aneuploid prediction
if(T){
  plist = list()
  for(sample.name in 'BCMHBCA65L2_090723'){
    
    # sample.name = 'BCMHBCA52R1'
    copykat_dir = paste0(svdir,'/improve_copykat/',sample.name,'/');
    dir.create(copykat_dir);setwd(copykat_dir)
    
    cpk_ST = read_rds(paste0(sample.name,'_cpk_ST.rds'))
    colData(cpk_ST)$pred_celltype = colData(cpk_ST)$pred.id
    
    label =c('CNA_prediction','pred_celltype')
    p1 = plotHeatmap_copykat(cpk_ST, label = label,label_colors = col[label],
                             row_split = 'CNA_prediction',order_cells = NULL, highscale = 0.3, lowscale = -0.3)
    pdf(paste0(sample.name,'_heatmap_improved.pdf'),width = 8,height = 10)
    print(p1)
    dev.off()
  }
}