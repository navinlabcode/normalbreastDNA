# Rscript ./copykat_ST/2.Runmin_improve_ST.R
## Notes ##
# source('/volumes/USR1/yiyun/Script/Rum_lib.R')

ordir ='/volumes/USR1/yiyun/Project/HBCA/'
setwd(ordir)

# source('/volumes/USR2/rumwei/tools/copykat_knn_dn.R')
`%!in%` = Negate(`%in%`)
source('./plot_0701/upload_code/ST_CNA_infer_functions.R')
library(Seurat)
require(Hmisc)
require(zoo)
require(RColorBrewer)
require(mclust)
require(grid)
require(copynumber)
library(dplyr)
library(readr)
require(copykat)
require(magrittr)
require(Matrix)
require(ggpubr)




options(future.globals.maxSize = 20000 * 1024^2)
plan("multiprocess", workers = 100)

#data dir
dir_vec <- c(paste0(grep('_061223',list.files('/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/',full.names = T),value = T),'/outs/'),
             '/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/BCMHBCA52R1/BCMHBCA52R1/outs/',
             '/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/BCMHBCA04_spara1-2_ref382020/outs/')
smp_vec <-c(grep('_061223',list.files('/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/'),value = T),'BCMHBCA52R','BCMHBCA04L')

#run copykat  
for (i in 1:length(dir_vec)) {
  print(i)
  obj_i <- Load10X_Spatial(dir_vec[i])
  count_i <- as.matrix(obj_i@assays$Spatial@counts)
  print('running CopyKat...')
  for(w in c(25)){
    for(k in c(0.05)){
      wkdir =paste0('./copykat_ST/win_',w,'_ks_',k);dir.create(wkdir)
      setwd(wkdir)
      tryCatch({
        copykat::copykat(rawmat=count_i, id.type='S', ngene.chr=5, LOW.DR=0.05, UP.DR=0.1, win.size=w, norm.cell.names = "", KS.cut = k, sam.name=smp_vec[i], n.cores=50)}, 
        error=function(e) {cat("ERROR :",conditionMessage(e), "\n")})
    }
  }
}

#define spots with CNAs
if(T){
  # Preparation -------------------------------------------------------------
  st_out_dir <- c(paste0(grep('_061223',list.files('/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/',full.names = T),value = T),'/outs/'),
                  '/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/BCMHBCA52R1/BCMHBCA52R1/outs/',
                  '/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/BCMHBCA04_spara1-2_ref382020/outs/')
  all_st_sample <- c(grep('_061223',list.files('/volumes/seq/projects/10X_ST/Fresh/BCMHBCA/spaceranger-1.2.0/'),value = T),'BCMHBCA52R','BCMHBCA04L')
  sample = c("BCMHBCA52R","BCMHBCA04L","BCMHBCA60L1_061223")
  pick_out_dir <-st_out_dir[all_st_sample %in% sample]
  pick_sample <-all_st_sample[all_st_sample %in% sample]
  
  if(T){
    hbca_sc_int_sub <- readRDS('./copykat_ST/ST_ref/HBCA_SC_int_sub_3hr_020922.rds')
    p_emp <- ggplot() + theme_void()
    
    # DEG ---------------------------------------------------------------------
    hbca_sc_int_deg <- FindAllMarkers(hbca_sc_int_sub)
    write_rds(hbca_sc_int_deg,'./copykat_ST/ST_ref/hbca_sc_int_deg.rds')
    
    hbca_sc_int_deg = read_rds('./copykat_ST/ST_ref/hbca_sc_int_deg.rds')
    hbca_sc_int_deg$pct_log2FC <- log(hbca_sc_int_deg$pct.1/hbca_sc_int_deg$pct.2, base = 2)
    hbca_sc_int_deg_top <- hbca_sc_int_deg %>% dplyr::filter(pct.1>.2, pct_log2FC>1, avg_log2FC>.5) %>% group_by(gene) %>% top_n(1, avg_log2FC) %>% group_by(cluster) %>% top_n(20, avg_log2FC)
    hbca_sc_int_deg_top_genelist <- as.character(hbca_sc_int_deg_top$cluster) %>% set_names(hbca_sc_int_deg_top$gene) %>% to_list()
    
    # Loop over samples -------------------------------------------------------
    for (i in 3) {
      print(i)
      sample_temp <- pick_sample[i]
      st_dir_temp <- pick_out_dir[i]
      wkdir = paste0(ordir,'/copykat_ST/improve_copykat/',sample_temp,'/');dir.create(wkdir,recursive = T)
      setwd(wkdir)
      
      # Load CopyKat results  ---------------------------------------------------------
      cpkt_res <- cpkt_load(mat_dir = paste0(ordir,'/copykat_ST/win_25_ks_005/', sample_temp, '_copykat_CNA_results.txt'), 
                            prd_dir = paste0(ordir,'/copykat_ST/win_25_ks_005/', sample_temp, '_copykat_prediction.txt'))
      cpkt_res$cpkt_mat[1:5,1:5]
      # MultiPCF used later in plotting ----------------------------------------------------------------
      cpkt_multipcf <- multipcf_wrp(mat=cpkt_res$cpkt_mat) %>% as.data.frame %>% dplyr::mutate(cell.names=rownames(.))
      cpkt_multipcf[1:5,1:5]
      hg19 = read_tsv(paste0(ordir,'/copykat_ST/hg19rg.tsv'))
      index = table(hg19 %>% mutate(chr_arm = paste0(chr,arm)) %>% pull(chr_arm))
      names(index) = gsub('chrX','chr23',names(index))
      cpkt_multipcf_mat <- do.call(cbind,
                                   lapply(colnames(cpkt_multipcf)[1:41],
                                          function(x){
                                            print(x)
                                            matrix(rep(cpkt_multipcf[,x], index[x]),
                                                   ncol = index[x],
                                                   byrow = F)
                                          }))
      cpkt_multipcf_mat[1:5,1:5]
      colnames(cpkt_multipcf_mat) = colnames(cpkt_res$cpkt_mat)
      rownames(cpkt_multipcf_mat) = rownames(cpkt_res$cpkt_mat)
      cpkt_res$cpkt_pcf_mat = cpkt_multipcf_mat
      
      
      # CNV score ---------------------------------------------------------------
      cnvscore_df <- score_cut(mat=cpkt_res$cpkt_mat, G=3) %>% dplyr::mutate(cell.names=rownames(.))
      
      # UMAP and Cluster -----------------------------------------------------------------
      umap_df <- umap_clust(mat=cpkt_res$cpkt_mat, metric = 'correlation', n_threads=15, n_neighbors=15, min_dist=0, seed=42, grp_cut=10, k=15, eps=.5, minPts=10) %>% 
        dplyr::mutate(cell.names=rownames(.)) %>% as.data.frame %>% set_rownames(.$cell.names)
      
      # Merge -------------------------------------------------------------------
      umap_df <- purrr::reduce(.f = inner_join, list(umap_df, cnvscore_df, cpkt_res$cpkt_pred), by='cell.names') %>% as.data.frame %>% set_rownames(.$cell.names)
      umap_df %>% dplyr::group_by(cnv_grp) %>% dplyr::summarise(cnvscore = mean(cnvscore)) %>% arrange(-cnvscore) 
      # Tumor identification ----------------------------------------------------
      ## This step requires sanity check ##
      if(sample_temp == 'BCMHBCA04L'){
        cnv_grp_pick = umap_df %>% dplyr::group_by(cnv_grp) %>% dplyr::summarise(cnvscore = mean(cnvscore)) %>% arrange(-cnvscore) %>% head(2) %>% pull(cnv_grp);print(cnv_grp_pick)
        umap_df <- tumor_assign(umap_df = umap_df, cnv_grp_ = cnv_grp_pick[2], cnvscore_cut_=sort(unique(umap_df$cnvscore_cut), decreasing = T)[1], 
                                cpkt.pred_=c('dipd'), grp_otl_=c(1), ncell_=10, prop_=70) %>% as.data.frame %>% set_rownames(.$cell.names)
      }else{
        cnv_grp_pick = umap_df %>% dplyr::group_by(cnv_grp) %>% dplyr::summarise(cnvscore = mean(cnvscore)) %>% arrange(-cnvscore) %>% head(2) %>% pull(cnv_grp);print(cnv_grp_pick)
        umap_df <- tumor_assign(umap_df = umap_df, cnv_grp_ = cnv_grp_pick[1], cnvscore_cut_=sort(unique(umap_df$cnvscore_cut), decreasing = T)[1], 
                                cpkt.pred_=c('dipd'), grp_otl_=c(1), ncell_=10, prop_=70) %>% as.data.frame %>% set_rownames(.$cell.names)
      }
      
      # Plot-1 ----------------------------------------------------------------
      p_cnv <- cnv_plot(umap_df, mat=cpkt_res$cpkt_pcf_mat, out_path=paste0(sample_temp,'test0.pdf'))
      
      # metadata generate -------------------------------------------------------
      metadata <- purrr::reduce(.f = inner_join, list(umap_df %>% dplyr::select(cell.names, cnv_X1=X1, cnv_X2=X2, cnv_grp, cnvscore, cnvscore_cut, cpkt.pred, tumor), 
                                                      cpkt_multipcf), by='cell.names') %>% as.data.frame %>% set_rownames(.$cell.names)
      
      # ST Analyis ----------------------------------------------------
      st_res <- st_analysis(st_dir=st_dir_temp, clust_res=.5, nFeature_Spatial_=.05, 
                            sc_srt=hbca_sc_int_sub, pred_='celltype', metadata = metadata)
      srt_temp <- st_res$srt_out
      
      # Plot-2 ----------------------------------------------------------------
      p_st <- st_plot(st_res=st_res, pt.size.factor=1.2, cnvevent_plot=c('chr1q'), out_path=paste0(sample_temp,'test1.pdf'))
      
      # Plot-all ----------------------------------------------------------------
      p_all <- plot_grid(p_cnv, p_st, nrow=2, rel_heights = c(.4, 1))
      
      pdf(paste0(sample_temp, '_cnv_st_summary.pdf'), width=15, height=30)
      p_all
      dev.off()
      
      # Outputdata --------------------------------------------------------------
      print('Saving data...')
      srt_temp@assays$Spatial@scale.data <- matrix(NA, 1, 1)
      saveRDS(srt_temp, paste0(sample_temp, '_Seurat_ST.rds'))
    }
  }
  
}
