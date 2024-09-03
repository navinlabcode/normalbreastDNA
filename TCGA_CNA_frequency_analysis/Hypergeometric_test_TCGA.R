# Rscript /volumes/USR1/yiyun/Project/HBCA/Revision/R3.C9_correlation_of_TCGA_frequency_plot.R
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(GenomicRanges)
library(magrittr)
library(BiocParallel)
library(ggplot2)


data = "./TCGA files/BRCA.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt"
grps = "./TCGA/IBC_v2/sample_cat.tsv"
cytoBand = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz"
threshold = 0.3

svdir = '/volumes/USR1/yiyun/Project/HBCA'


#read TCGA freq and adding arm info---------
if(F){
  seg_TCGA = read.table(data,sep="\t",header=T)
  seg_TCGA$Segment_Mean %>% hist(.,breaks = 100)

  
  #construct chr arm level ranges
  hg19_rg <- fread(cytoBand, 
                   col.names = c("chr","start","end","name","gieStain"))%>% data.frame() %>% 
    mutate(arm = substring(name, 1, 1),
           length = end-start) %>% 
    mutate(chr_arm = paste0(chr,arm))

  hg19_rg_byarm_list = split(hg19_rg,hg19_rg$chr_arm) 
  hg19_rg_byarm_df = do.call(rbind,
                             lapply(hg19_rg_byarm_list,function(x){
                               hg19_rg_byarm_df_i = x[1,] %>% select(chr,start,end,arm,chr_arm)
                               hg19_rg_byarm_df_i$start =  min(x$start)
                               hg19_rg_byarm_df_i$end =  max(x$end)
                               return(hg19_rg_byarm_df_i)
                             })
  )
  gr_varbin_full <-makeGRangesFromDataFrame(hg19_rg_byarm_df,keep.extra.columns = T)
  GenomeInfoDb::seqlevelsStyle(gr_varbin_full) <- "Ensembl"
  gr_varbin_full <- GenomeInfoDb::renameSeqlevels(gr_varbin_full, 
                                                  c(X = 23, Y = 24))
  
  #add chr arm by patient 
  seg_TCGA_list = split(seg_TCGA,seg_TCGA$Sample)
  do.call(c,lapply(seg_TCGA_list, nrow)) %>% hist(.,breaks = 100)
  do.call(c,lapply(seg_TCGA_list, nrow)) %>% median
  
  seg_TCGA_df = do.call(rbind,
                        bplapply(names(seg_TCGA_list),
                                 function(x){
                                   rg_TCGA = seg_TCGA_list[[x]] %>% 
                                     mutate(chr = Chromosome,start=Start,
                                            end = End,chr.start.end = paste0(Chromosome,'_',Start,'_',End),
                                            segratio = Segment_Mean) %>% 
                                     select(chr,start,end,chr.start.end,segratio)
                                   
                                   rg_TCGA_form = makeGRangesFromDataFrame(rg_TCGA,keep.extra.columns = T)
                                   overlaps <- findOverlaps(query = rg_TCGA_form, subject = gr_varbin_full)
                                   overlaps %>% data.frame()
                                   olp_df = data.frame(chr.start.end = rg_TCGA_form$chr.start.end[queryHits(overlaps)], 
                                                       chr_arm = gr_varbin_full$chr_arm[subjectHits(overlaps)]) 
                                   olp_df = left_join(olp_df,rg_TCGA)
                                   olp_df$Sample = x
                                   return(olp_df)
                                 },
                                 BPPARAM = bpparam()))

  #set cutoff for gain/loss
  cat_TCGA_df_chrarm = seg_TCGA_df %>% dplyr::group_by(Sample,chr_arm) %>% dplyr::summarise(segratio = median(segratio))
  hist(cat_TCGA_df_chrarm$segratio,breaks = 100)
  cat_TCGA_df_chrarm$cat = if_else(cat_TCGA_df_chrarm$segratio>(1*threshold), 1, if_else(cat_TCGA_df_chrarm$segratio<(-1*threshold), (-1), 0))
  cat_TCGA_df_chrarm$cat %>% table()
  
  
  #add ER group
  ERgrp = read_tsv(grps) %>% set_colnames(c('Sample','ER'))
  cat_TCGA_df_chrarm =  left_join(cat_TCGA_df_chrarm,ERgrp)
  lvl = lapply(paste0('chr',c(1:22,"X")),function(x) paste0(x,c('p','q'))) %>% unlist()
  cat_TCGA_df_chrarm$chr_arm = factor(cat_TCGA_df_chrarm$chr_arm,levels = lvl)
  cat_TCGA_df_chrarm = cat_TCGA_df_chrarm %>% filter(!is.na(ER))%>% arrange(Sample,chr_arm)
  cat_TCGA_df_chrarm %>% filter(chr_arm=='chr1q',ER == 'ERpos') %>% group_by(cat) %>% count() %>% ungroup() %>% mutate(pct = n/sum(n))
  cat_TCGA_df_chrarm %>% filter(chr_arm=='chr16q',ER == 'ERpos') %>% group_by(cat) %>% count() %>% ungroup() %>% mutate(pct = n/sum(n))
  cat_TCGA_df_chrarm %>% filter(chr_arm=='chr1q',ER == 'ERneg') %>% group_by(cat) %>% count() %>% ungroup() %>% mutate(pct = n/sum(n))
  cat_TCGA_df_chrarm %>% filter(chr_arm=='chr5q',ER == 'ERneg') %>% group_by(cat) %>% count() %>% ungroup() %>% mutate(pct = n/sum(n))
  
}

#read Normal breast data and summarize it into same format as TCGA is---------
if(T){
  meta_chr = read_rds(paste0(svdir,'/event_sum/filtered_cells_w_event.rds')) %>% filter(Sample_ID!='BCMHBCA06R')
  head(meta_chr)
  meta_chr$cell %>% unique() %>% length()
  
  lvl = lapply(paste0('chr',c(1:22,"X")),function(x) paste0(x,c('p','q'))) %>% unlist()
  meta_chr$chr_arm = factor(meta_chr$chrom,levels = lvl)
  meta_chr$count = 1
  meta_chr_w =meta_chr  %>% 
    select(chrom_event_new,count,cell,Sample_ID)%>% 
    pivot_wider(names_from = chrom_event_new,values_from = count,values_fill = 0)
  meta_chr_w %>% filter()
  
  cat_Normal_df_chrarm = do.call(rbind,
                                 bplapply(lvl,
                                          function(chr_arm_){
                                            event = grep(chr_arm_,colnames(meta_chr_w),value = T)
                                            if(length(event)>0){
                                              print(event)
                                              cat_ = rep(0, nrow(meta_chr_w))
                                              gain_event = grep("Gain",event,value =T)
                                              if(length(gain_event)>0){
                                                cat_[meta_chr_w[,gain_event]==1] = 1
                                              }
                                              loss_event = grep("Loss",event,value =T)
                                              if(length(loss_event)>0){
                                                cat_[meta_chr_w[,loss_event]==1] = -1
                                              }
                                              cat_Normal_df_chrarm = data.frame(Cell =meta_chr_w$cell, chr_arm = chr_arm_, cat = cat_, Project = 'Normal')
                                            }
                                          },BPPARAM = bpparam()))
  
}

#permutation analysis(null distribution generated based on sampling "gain/loss/neutral" over the whole genome)---------
if(T){
 
  #TCGA
  if(T){
    pval_df_all=data.frame()
    plist = list()
    for(ER_ in c('ERpos','ERneg')){
      ### permutation ###
      #null: the frequency of gain/loss is the just the noise 
      #construct null distribution: randomly sample the gain/loss/neutral from the 46 chromosome arms * number of samples(ERpos = 763, ERneg = 217) in the cohort
      #H1: reject null, meaning the signal detected from this specific chromosome arm is not a noise.
      set.seed(40)
      nperm = 1000
      cat_TCGA_df_chrarm_i = cat_TCGA_df_chrarm  %>% data.frame() %>% filter(ER == ER_)
      size_ = cat_TCGA_df_chrarm %>% data.frame() %>% 
        filter(ER == ER_) %>% pull(Sample) %>% unique() %>% length()
      freqdf_perm = do.call(rbind,bplapply(1:nperm,
                                           function(i){
                                             perm_i = sample(x = cat_TCGA_df_chrarm_i$cat,size = size_,replace = T)
                                             freqdf_i = table(perm_i) %>% data.frame() %>% set_colnames(c('cat','n')) %>% mutate(pct = n/sum(n))
                                           },BPPARAM = bpparam()))
      plist[[ER_]] <- ggplot(freqdf_perm, aes(x=pct,color=cat)) + 
        geom_histogram(fill="white", bins=500)+
        scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),breaks = c(0,1,-1))+
        facet_wrap(~cat,scales = 'free_x')+
        ylab(ER_)+
        theme_light()
      
      for(chr_arm_ in levels(cat_TCGA_df_chrarm$chr_arm)){
        print(chr_arm_)
        obsv_i = cat_TCGA_df_chrarm %>% data.frame() %>% 
          filter(ER == ER_, chr_arm == chr_arm_) 
        if(nrow(obsv_i)>0){
          obsv_freqdf = table(obsv_i$cat) %>% data.frame() %>% set_colnames(c('cat','n')) %>% mutate(pct = n/sum(n))
          for(cat_ in c(1,-1)){
            ref = obsv_freqdf %>% filter(cat==cat_) %>% pull(pct)
            if(length(ref)>0){
              null = freqdf_perm %>% filter(cat==cat_) %>% pull(pct)
              mean_null = mean(null)
              sd_null = sd(null)
              zi = (ref-mean_null)/sd_null
              pi = pnorm(q=zi,lower.tail = F)
              
              #write the result
              pval_df_i = data.frame(ER = ER_,chr_arm = chr_arm_, cat = cat_, z= zi, pval = pi, pavlr = round(pi,digits = 4))
              pval_df_all = rbind(pval_df_all,pval_df_i)
            }
          }
        }
      }
    }
  }
  pval_df_TCGA = pval_df_all
  #Normal
  if(T){
    ### permutation ###
    #null: the increase frequence of gain/loss is the just the noise 
    #construct null distribution: 
    # if I randomly sample the gain/loss/neutral from the 46 chromosome arms * number of samples(ERpos = 763, ERneg = 217) in the cohort,
    # the gain and loss should be in a very low frequency
    #H1: reject null, meaning the signal detected from this specific chromosome arm is not a noise
    
    set.seed(40)
    nperm = 1000
    size_ = cat_Normal_df_chrarm %>% data.frame() %>% pull(Cell) %>% unique() %>% length()
    freqdf_perm = do.call(rbind,bplapply(1:nperm,
                                         function(i){
                                           perm_i = sample(x = cat_Normal_df_chrarm$cat,size = size_,replace = T)
                                           freqdf_i = table(perm_i) %>% data.frame() %>% set_colnames(c('cat','n')) %>% mutate(pct = n/sum(n))
                                         },BPPARAM = bpparam()))
    
    plist[['Normal']] <- ggplot(freqdf_perm, aes(x=pct,color=cat)) + 
      geom_histogram(fill="white", bins=500)+
      scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),breaks = c(0,1,-1))+
      facet_wrap(~cat,scales = 'free_x')+
      ylab('Normal')+
      theme_light()
    
    
    lvl = lapply(paste0('chr',c(1:22,"X")),function(x) paste0(x,c('p','q'))) %>% unlist()
    pval_df_all=data.frame()
    for(chr_arm_ in lvl){
      print(chr_arm_)
      obsv_i = cat_Normal_df_chrarm %>% data.frame() %>% filter(chr_arm == chr_arm_) 
      if(nrow(obsv_i)>0){
        obsv_freqdf = table(obsv_i$cat) %>% data.frame() %>% set_colnames(c('cat','n')) %>% mutate(pct = n/sum(n))
        for(cat_ in c(1,-1)){
          ref = obsv_freqdf %>% filter(cat==cat_) %>% pull(pct)
          if(length(ref)>0){
            null = freqdf_perm %>% filter(cat==cat_) %>% pull(pct)
            mean_null = mean(null)
            sd_null = sd(null)
            zi = (ref-mean_null)/sd_null
            pi = pnorm(q=zi,lower.tail = F)
            
            #write the result
            pval_df_i = data.frame(Project = 'Normal',chr_arm = chr_arm_, cat = cat_, z= zi, pval = pi, pavlr = round(pi,digits = 4))
            pval_df_all = rbind(pval_df_all,pval_df_i)
          }
        }
      }
    }
  }
  pval_df_normal = pval_df_all
  
  #check null distribution
  cowplot::plot_grid(plotlist = plist,ncol = 1) %>% print()
}

#enrichment analysis
if(T){
  #select significant events for TCGA
  TCGA_sum = cat_TCGA_df_chrarm %>% dplyr::group_by(ER,chr_arm,cat) %>% 
    dplyr::count() %>% dplyr::group_by(ER,chr_arm)%>% mutate(pct = n/sum(n))%>% 
    mutate(event = if_else(cat==1,'Gain',if_else(cat == -1, 'Loss','Neu'))) %>% 
    mutate(chr_arm_event = paste0(chr_arm,'_',event))
  pval_df_TCGA = left_join(pval_df_TCGA,TCGA_sum)
  
  p_threshold = 0.0005
  pct_threshold =0.15
  pval_df_TCGA %>% filter(ER == 'ERpos',pval<p_threshold,pct>pct_threshold) %>% arrange(-pct,pval)
  pval_df_TCGA %>% filter(ER == 'ERneg',pval<p_threshold,pct>pct_threshold) %>% arrange(-pct,pval)
  ERpos_event = pval_df_TCGA %>% filter(ER == 'ERpos',pval<p_threshold,pct>pct_threshold) %>% arrange(-pct,pval) %>% pull(chr_arm_event);ERpos_event#ERpos sig event
  ERneg_event = pval_df_TCGA %>% filter(ER == 'ERneg',pval<p_threshold,pct>pct_threshold) %>% arrange(-pct,pval) %>% pull(chr_arm_event);ERneg_event#ERneg sig event
  
  #select significant events for normal
  Noraml_sum = cat_Normal_df_chrarm %>%dplyr:: group_by(chr_arm,cat) %>% dplyr::count() %>% dplyr::group_by(chr_arm)%>% mutate(pct = n/sum(n)) %>% 
    mutate(event = if_else(cat==1,'Gain',if_else(cat == -1, 'Loss','Neu'))) %>% 
    mutate(chr_arm_event = paste0(chr_arm,'_',event))
  pval_df_normal = left_join(pval_df_normal,Noraml_sum)
  pval_df_normal %>% filter(pval<p_threshold) %>% arrange(-pct,pval) 
  Normal_event = pval_df_normal %>% filter(pval<p_threshold) %>% arrange(-pct,pval) %>% pull(chr_arm_event);Normal_event
  Sporadic_event = c('chr20p_Gain','chr20q_Gain','chr10p_Loss','chr3q_Loss')
  #extract overlap
  intersect(Normal_event,ERpos_event)
  intersect(Normal_event,ERneg_event)
  intersect(Sporadic_event,ERpos_event) %>% print()
  intersect(Sporadic_event,ERneg_event) %>% print()
  
  
  #hypergeometric testp
  library(fgsea)
  #is a list
  pathways_ = list(ERpos_event = ERpos_event,
                   ERneg_event = ERneg_event)
  stats_ = pval_df_normal %>% filter(pval<p_threshold) %>% arrange(z) %>% pull(z,chr_arm_event)
  gene_ = pval_df_normal %>% filter(pval<p_threshold) %>% arrange(z) %>% pull(chr_arm_event)
  lvl = lapply(paste0('chr',c(1:22,"X")),function(x) paste0(x,c('p','q'))) %>% unlist()
  universe_ = do.call(c,lapply(c('Gain','Neu','Loss'),function(x) paste0(lvl,'_',x)))
  
  set.seed(101)
  fgseaRes1 <- fora(pathways = pathways_, 
                    genes = gene_[!grepl('chrX',gene_)],
                    universe = universe_)
  fgseaRes1
  #expanded
  # pathway       pval       padj overlap size                       overlapGenes
  # 1: ERpos_event 0.02457497 0.04914994       3   11 chr1q_Gain,chr16q_Loss,chr22q_Loss
  # 2: ERneg_event 0.05890118 0.05890118       3   15 chr1q_Gain,chr10q_Loss,chr16q_Loss
  
  #is a list
  pathways_ = list(ERpos_event = ERpos_event,
                   ERneg_event = ERneg_event)
  stats_ = pval_df_normal %>% filter(chr_arm_event %in% Sporadic_event) %>% arrange(z) %>% pull(z,chr_arm_event)
  gene_ = pval_df_normal %>% filter(chr_arm_event %in% Sporadic_event) %>% arrange(z) %>% pull(chr_arm_event)
  lvl = lapply(paste0('chr',c(1:22,"X")),function(x) paste0(x,c('p','q'))) %>% unlist()
  universe_ = do.call(c,lapply(c('Gain','Neu','Loss'),function(x) paste0(lvl,'_',x)))
  
  set.seed(101)
  fgseaRes2 <- fora(pathways = pathways_, 
                    genes = gene_[!grepl('chrX',gene_)],
                    universe = universe_)
  fgseaRes2
  #sporadic
  # pathway      pval      padj overlap size overlapGenes
  # 1: ERpos_event 0.2854483 0.5708967       1   11  chr20q_Gain
  # 2: ERneg_event 1.0000000 1.0000000       0   15    
}


#same analysis was used to test the frequency in CNA events from LumSec, LumHR and basal cells, 
#just need to redo this step: read Normal breast data and summarize it into same format as TCGA is--------- 