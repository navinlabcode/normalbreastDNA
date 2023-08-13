library(readr)
library(dplyr)
library(tidyr)
library(copykit)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(BiocParallel)
################# getEventMat function #######################
getEventMat <- function(
    consensus,          # consensus CN matrix of which will be converted to event matrix
    bin_adj = 2,    # number of bins allowed to be adjusted to consider as the same breakpoint
    ploidy_trunc = 8,  # maximum integer value, all integer value larger than this will be set to this
    rmY= F
){
  
  ## trunc integer matrix
  seg_df = t(consensus)
  seg_df[seg_df>=ploidy_trunc] = ploidy_trunc
  
  range = read_rds('/volumes/USR1/yiyun/Script/CNA/hg19_range.rds')
  range = range %>%
    dplyr::as_tibble() %>%
    dplyr::select(seqnames, start, end, arm)
  dim(range)
  if(rmY){
    range = range %>% filter(seqnames != 'chrY')
  }
  intmat <- range %>% 
    cbind(seg_df) %>%
    tibble::remove_rownames()
  intmat[1:5,1:5]
  tail(intmat)
  
  ## merge segments
  res_int <- as.data.frame(intmat[1,])
  for(i in 2:nrow(intmat)){
    if(identical(as.character(intmat[i,-c(1:4)]), as.character(intmat[i-1,-c(1:4)]))){
      next
    }else{
      res_int <- rbind(res_int, intmat[i,])
    }
  }
  res_int[,10:15] %>% head()
  res_int[,10:15] %>% tail()
  res_int$bin <- as.numeric(rownames(res_int))
  
  ## finding common breakpoints
  res_int_cbp <- as.data.frame(res_int[1,])
  for(i in 2:(nrow(res_int)-1)){
    if(res_int$bin[i+1]-res_int$bin[i]<=bin_adj){
      next
    }else{
      res_int_cbp <- rbind(res_int_cbp, res_int[i,])
    }
  }
  res_int_cbp <- rbind(res_int_cbp, res_int[nrow(res_int),])
  
  res_int_cbp$n.bins=c(res_int_cbp$bin[-1], ncol(consensus)+1) - res_int_cbp$bin
  tail(res_int_cbp)
  res_int_cbp$end.pos = intmat$end[as.numeric(res_int_cbp$bin)+res_int_cbp$n.bins-1]
  res_int_cbp$end.chr = intmat$seqnames[as.numeric(res_int_cbp$bin)+res_int_cbp$n.bins-1]
  res_df <- res_int_cbp %>%
    dplyr::rename(start.chr=seqnames, 
                  start.pos=start) %>%
    dplyr::select(starts_with("start"), end.chr, end.pos, bin, n.bins, everything(), -end) %>%
    tibble::remove_rownames()
  
  return(res_df)
}

#subset copykit object with the picked co-occurred events
if(T){
  filtered_aneu_obj = read_rds('./rds/filtered_aneu_obj.rds')
  filtered_aneu_obj = filtered_aneu_obj[,colData(filtered_aneu_obj)$Sample_ID!='BCMHBCA06R']
  lvls = colData(filtered_aneu_obj)$Sample_ID %>% table() %>% sort(decreasing = T) %>% names()
  colData(filtered_aneu_obj)$Sample_ID = factor(colData(filtered_aneu_obj)$Sample_ID,lvls)
  colData(filtered_aneu_obj) %>% colnames
  
  #cluster by chr events
  meta_df_w = read_tsv('./event_sum/cell_events_summary.tsv')
  meta_df_w$index = str_extract(meta_df_w$combined_events,pattern = '^chr\\d{1,2}|^chrX')
  meta_df_w$index = factor(meta_df_w$index, levels = paste0('chr',c(1:22,'X'))[paste0('chr',c(1:22,'X')) %in% unique(meta_df_w$index) ])
  meta_df_w %>% arrange(index)
  
  #add cluster
  head(colData(filtered_aneu_obj)$sample)
  (colData(filtered_aneu_obj)$sample %in% meta_df_w$cell) %>% sum()
  meta_df_w$sample = meta_df_w$cell
  meta = colData(filtered_aneu_obj) %>% data.frame()
  meta = left_join(meta,meta_df_w %>% select(sample, combined_events,index))
  
  colData(filtered_aneu_obj)$combined_events = meta$combined_events
  colData(filtered_aneu_obj)$index = meta$index

  
  #plot phy for picked combined events
  result_df = read_rds('./plot_0701/Fig3/picked_occurr_event')
  result_df %>% arrange(col_chr)# check events on the list and also in figure 3d
  # chr1q_Gain == 1|chr16q_Loss == 1|chr7q_Loss == 1|chr5p_Gain ==1|chr10q_Loss ==1|chr14q_Gain ==1|chr17p_Loss==1|chr12p_Loss ==1|chr22q_Loss ==1||chr23p_Loss==1|chr23p_Gain==1|chr21q_Loss==1
  event_chr1 = paste0('chr1q_Gain.*.',c('chr10q_Loss','chr5p_Gain','chr16q_Loss.*','chr3p_Loss'))
  event_chr7 = paste0('chr7q_Loss.*.',c('chr16q_Loss','chr22q_Loss'))
  event_chr8 = paste0('chr8p_Loss.*.',c('chr16q_Loss'))
  event_chr14 = paste0('chr14q_Gain.*.',c('chr22q_Loss','chr18p_Loss'))
  event_chr16 =paste0('chr16q_Loss.*.',c('chr17p_Loss'))
  event_chr17 = paste0('chr17p_Loss.*.',c('chr19p_Gain.*','chr22q_Loss'))
  event = c('^chr1q_Gain$','^chr8p_Loss$','^chr14q_Gain$','^chr16q_Loss$','^chr17p_Loss$',
            event_chr1,event_chr7,event_chr8,event_chr14,event_chr16,event_chr17)
  meta= colData(filtered_aneu_obj) %>% data.frame() %>% 
    filter(grepl(paste(event,collapse = '|'),combined_events)) 
  meta$combined_events %>% table()
  colData(filtered_aneu_obj)$combined_events %>% table()
  filtered_aneu_obj_chr1 <- filtered_aneu_obj[,rownames(meta)]#most interacted
  colData(filtered_aneu_obj_chr1)$combined_events %>% table()
  
  ct = colData(filtered_aneu_obj_chr1)$combined_events %>% table() %>% sort();ct
  level = c(ct[ct>=5]%>% names(),'other');level
  
  colData(filtered_aneu_obj_chr1)$combined_events = factor(colData(filtered_aneu_obj_chr1)$combined_events, levels = level)
  colData(filtered_aneu_obj_chr1)$combined_events[is.na(colData(filtered_aneu_obj_chr1)$combined_events)]='other'
}

#calculate integer, run medicc2 in concensus
if(T){
  colData(filtered_aneu_obj_chr1)$combined_events %>% table()
  (colData(filtered_aneu_obj_chr1)$combined_events == 'other') %>% sum()
  filtered_aneu_obj_chr1_sub <- filtered_aneu_obj_chr1[,colData(filtered_aneu_obj_chr1)$combined_events != 'other']
  colData(filtered_aneu_obj_chr1_sub)$combined_events %>% unique()
  colData(filtered_aneu_obj_chr1_sub)$combined_events %>% table()
  
  colData(filtered_aneu_obj_chr1_sub)$ploidy=2
  seg <- SummarizedExperiment::assay(filtered_aneu_obj_chr1_sub, "segment_ratios")
  int_values <- round(as.matrix(seg) %*% diag(colData(filtered_aneu_obj_chr1_sub)$ploidy)) %>% 
    as.data.frame() %>% set_colnames(colnames(seg)) %>% set_rownames(rownames(seg))
  int_values[1:5,1:5]
  long_list <- split(as.data.frame(t(int_values)), colData(filtered_aneu_obj_chr1_sub)$combined_events)
  consensus_list <- BiocParallel::bplapply(long_list, function(x) {
    apply(x, 2, median)}, BPPARAM = bpparam())
  cs_df <- as.data.frame(t(do.call(rbind, consensus_list)))
  names(cs_df) <- names(consensus_list)
  cs_df[1:5,1:5]
  dim(cs_df)
  cluster_consensus = cs_df %>% select(-other) %>% t() %>% data.frame();cluster_consensus[,1:5]
  
  eventmat <- getEventMat(cluster_consensus, bin_adj = 2, ploidy_trunc = 4,rmY = T)
  popseg_long <- as.data.frame(apply(as.data.frame(t(eventmat %>% dplyr::select(rownames(cluster_consensus)))), 1, 
                                     function(m) {rep.int(m, eventmat$n.bins)}))
  popseg_long %>% head()
  
  medicc_input <- cbind(SummarizedExperiment::rowRanges(filtered_aneu_obj_chr1_sub) %>% dplyr::as_tibble() %>% 
                          dplyr::select(seqnames, start, end), popseg_long) %>%
    mutate(diploid=2) %>%          ## here we manually add a diploid cell CN profile to be used as root in the medicc tree
    dplyr::rename(chrom=seqnames) %>% gather('sample_id', 'CN', -chrom, -start, -end) %>% dplyr::select(sample_id,chrom, everything())
  medicc_input %>% head()
  write_tsv(medicc_input, file = ("./plot_0701/Fig3/medicc2/medicc2_input.tsv"))
  
  
  #----run medicc2 in terminal----------
  # cd ./plot_0701/Fig3/medicc2/
  # conda activate medicc_env
  # medicc2 -a CN --total-copy-numbers -j 40 -vv medicc2_input.tsv output
}

#read medicc2 distance and run fastme.ols tree
if(T){
  dis = read_tsv('./plot_0701/Fig3/medicc2/output/medicc2_input_pairwise_distances.tsv') %>% data.frame()
  dis_mat = dis %>% set_rownames(dis$sample_id) %>% select(-sample_id) %>% as.matrix()
  tree = ape::fastme.ols(dis_mat)
  # tree = ape::nj(dis_mat)
  # tree = ape::fastme.bal(dis_mat)
  plot(tree)
  tree <- ape::root.phylo(tree, outgroup = 'diploid', resolve.root = TRUE)
  plot(tree)
  
  library(ggtree)
  group_data = colData(filtered_aneu_obj) %>% data.frame() %>% select(combined_events,index) %>% dplyr::distinct()
  list_samples <- split(group_data$combined_events, group_data$index)
  # tree <- ggtree::groupOTU(medic, list_samples)
  tree <- ggtree::groupOTU(tree, list_samples)
  plot(tree)
  
  #add size frequency and branch color
  freq_df = colData(filtered_aneu_obj) %>% data.frame() %>% 
    dplyr::group_by(combined_events) %>% dplyr::count()  %>% 
    data.frame() %>% mutate(prop = n/sum(n)*100)
  freq_df=freq_df %>% filter(combined_events %in% tree$tip.label)
  freq_df$color = if_else(grepl('chr14q',freq_df$combined_events),'#FF9006',
                          if_else(grepl('^chr17p',freq_df$combined_events),'#2A9CFF',
                                  if_else(grepl('^chr8p_Loss$',freq_df$combined_events),'#F781BF',
                                          if_else(grepl('^chr1q_Gain',freq_df$combined_events),'#E41A1C', '#5E985E')
                                  )
                          )
  )
  
  freq_df[nrow(freq_df)+1,] = c('diploid',30,1.5,'grey')
  freq_df$n10 = as.numeric(freq_df$n)/10
  p <- ggtree::ggtree(tree, ladderize = FALSE) %<+% 
    freq_df +  geom_tippoint(aes(size=n10, color=color,fill = color)) + 
    scale_color_manual(values = freq_df$color, breaks = freq_df$color, limits = force) +
    scale_size_continuous(name = 'Events count',
                          breaks = c(1,3,6,10),
                          labels = c(10,30,60,100))+
    theme(legend.position = "right") + 
    geom_treescale(x = 4)+
    geom_tiplab(aes(color = color), size = 3, hjust = -0.1)+ 
    theme_tree2();p
  pdf('./plot_0701/Fig3/phylogeny_me_ols_on_medicc2_root_diploid.pdf',width = 4.5,height = 6)
  p %>% print()
  dev.off()
}