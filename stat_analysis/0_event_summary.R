# Rscript /volumes/USR1/yiyun/Project/HBCA/event_sum/Check_sporatic_events_new_expand_vs_sporadic.R
if(T){
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(igraph)
  library(copykit)#V0.1.0
  library(readr)
  library(base)
  library(reader)
  library(stringr)
  library(magrittr)
  library(ggrepel)
  library(BiocParallel)
  library(ggimage)
  library(cowplot)
  library(stringr)
  library(ComplexHeatmap)
  # source("/volumes/USR1/yiyun/Script/IndexART.R")
  # source("/volumes/USR1/yiyun/Script/CNApackage.R")
  source('/volumes/USR1/yiyun/Project/HBCA/code_deposit/event_summary/chr_arm_event_check_function.R')
}



svdir = '/volumes/USR1/yiyun/Project/HBCA/'
setwd(svdir)
#make summary of patient and event from Kris's cell name, add it into filtered_aneu_obj colData
if(T){
  filtered_aneu_obj = read_rds(paste0(svdir,'/rds/filtered_aneu_obj_junke.rds'))
  meta_chr = read_rds(paste0(svdir,'/event_sum/filtered_cells_w_event.rds'))
  head(meta_chr)
  meta_chr$count = 1
  meta_chr_w =meta_chr %>% select(chrom_event_new,count,cell,Sample_ID)%>% 
    pivot_wider(names_from = chrom_event_new,values_from = count,values_fill = 0)
  combined_events_index = meta_chr_w %>% select(-cell,-Sample_ID) %>% apply(.,1,sum)>1
  
  meta_chr_w$combined_events = meta_chr_w %>% apply(.,1,
                                                    function(x){
                                                      if(length(colnames(meta_chr_w)[x==1])>1){
                                                        combined_events = paste(colnames(meta_chr_w)[x==1],collapse = '.')
                                                      }else{
                                                        combined_events = colnames(meta_chr_w)[x==1]
                                                      }
                                                      return(combined_events)
                                                    })
  meta_chr_w$combined_events %>% table() %>% sort(.,decreasing = T)
  if(grepl(pattern = 'chr23',meta_chr_w$combined_events) %>% sum()>0){
    meta_chr_w$combined_events = gsub('chr23','chrX',meta_chr_w$combined_events)
  }
  meta_chr_w$Sample_ID %>% table()
  unique(meta_chr_w$Sample_ID) %>% length()
  
  # first 1-2 columns are cellname and samplename, 
  # 3-46 columns are individual events detected per cells, 
  # column 47 are combined events eg. 1q_amp+16q_del per cells
  unique(meta_chr_w$combined_events)
  unique(grep('chrX',meta_chr_w$combined_events,value = T))
  meta_chr_w$if_chrX = if_else(grepl('chrX',meta_chr_w$combined_events), 'yes','no')
  # meta_chr_w$if_chrX = if_else((meta_chr_w$combined_events %in% c('chrXp_AMP_chrXq_AMP','chrXq_AMP','chrXp_DEL_chrXq_DEL')),
  #                              'yes','no')
  meta_chr_w$if_chrX %>% table()
  
  meta_chr_w%>% write_rds(.,paste0(svdir,'/event_sum/cell_events_summary.rds'))
  meta_chr_w%>% write_tsv(.,paste0(svdir,'/event_sum/cell_events_summary.tsv'))
  meta_chr_w$cell %>% duplicated() %>% table()
  table(colData(filtered_obj)$sample %in% meta_chr_w$cell) #filtered_obj (all cells object)
  
  # prepare all cells metadata
  meta = filtered_obj %>% colData() %>% data.frame() #filtered_obj (all cells object)
  meta_pt = read_tsv(paste0(svdir,'/metadata/table1_metadata_01092024.tsv'))
  meta = left_join(meta,meta_pt)
  is.na(meta$paper_id) %>% table()
  meta %>% write_tsv(.,paste0(svdir,'/event_sum/all_cells_copykit_metadata.tsv'))
  dim(meta)
  
  # add patient info to metadata
  df = colData(filtered_aneu_obj) %>% data.frame()
  
  dim(df)
  df$Sample_ID %>% unique()
  meta_pt = read_tsv(paste0(svdir,'/metadata/table1_metadata_01092024.tsv'))
  meta_pt$DNAdata_used_assay %>% table()
  
  df = left_join(df,meta_pt)
  dim(df)
  df$paper_id %>% is.na() %>% table()
  colData(filtered_aneu_obj) = cbind(colData(filtered_aneu_obj),
                                     df %>% select(-colnames(colData(filtered_aneu_obj) %>% data.frame())))
  
  
  
  #add events to metadata
  meta_chr_w = read_rds(paste0(svdir,'/event_sum/cell_events_summary.rds'))
  meta_chr_w$sample=meta_chr_w$cell
  colnames(meta_chr_w)
  
  df = colData(filtered_aneu_obj) %>% data.frame()
  df = left_join(df,meta_chr_w %>% select(-cell,-Sample_ID))
  dim(df)
  colData(filtered_aneu_obj) = cbind(colData(filtered_aneu_obj),
                                     df %>% select(-colnames(colData(filtered_aneu_obj) %>% data.frame())))
  write_rds(filtered_aneu_obj,paste0(svdir,'/rds/filtered_aneu_obj.rds'))
  
  meta = colData(filtered_aneu_obj) %>% data.frame()
  meta %>% write_rds(paste0(svdir,'/event_sum/CNA_cells_copykit_metadata.rds'))
  
  
  data = meta_chr_w %>% dplyr::group_by(Sample_ID,combined_events) %>% dplyr::count()
  # data$combined_events[data$combined_events %!in% keep] = 'others'
  # data$combined_events = factor(data$combined_events, levels = c(keep,'others'))
  data$patient = str_extract(string = data$Sample_ID,pattern = 'BCMHBCA\\d{1,2}')
  dim(meta_pt)
  meta_pt = meta_pt %>% filter(!duplicated(patient))
  data = left_join(data,meta_pt) %>% arrange(-Age)
  head(data)
  data$paper_id = factor(data$paper_id,levels = c(unique(data$paper_id)))
  write_rds(data,paste0(svdir,'/event_sum/cell_events_summary_long.rds'))
}

#plot aneuploid heatmap with chr_arm counts
if(T){
  filtered_aneu_obj = read_rds(paste0(svdir,'/rds/filtered_aneu_obj.rds'))
  dim(filtered_aneu_obj)
  col = list(
    paper_id = colormaker(colData(filtered_aneu_obj)$paper_id,'Paired'),
    assay = colormaker(colData(filtered_aneu_obj)$assay,'Set2'),
    chrom = colormaker(colData(filtered_aneu_obj)$chrom,'Set1'),
    is_aneuploid = c('TRUE'='#19376D','FALSE'='gray'),
    paper_id = colormaker(colData(filtered_aneu_obj)$paper_id,'Paired'),
    if_sporadic = c('expanded'='#F07857','sporadic'='#CBD6E2')
  )
  
  ## quantify sporatic events by patient
  cutoff = 3
  pt_event_df_all = data.frame()
  cell_event_df_all = data.frame()
  wkdir =paste0(svdir,'/event_sum/cutoff_',cutoff,'/');dir.create(wkdir);setwd(wkdir)
  # samplename= 'BCMHBCA89L'
  
  for(samplename in unique(colData(filtered_aneu_obj)$Sample_ID)){
    obj = filtered_aneu_obj[,colData(filtered_aneu_obj)$Sample_ID == samplename]
    obj= runPhylo(obj,method = 'nj',metric = 'manhattan',n_threads = 100)
    h = if_else(ncol(obj)*0.1 >3, ncol(obj)*0.1,3)
    
    #events sort by sporadic or expanded? define sporadic, cutoff>=3 is expanded?
    meta = colData(obj) %>% data.frame()
    
    chr_event_sum = meta[,grepl('chr\\d{1,2}|chrXp|chrXq',colnames(meta))] %>% colSums()
    names(chr_event_sum) = gsub('chr23','chrX',names(chr_event_sum))
    colnames(meta) = gsub('chr23','chrX',colnames(meta))
    
    pt_event_df = data.frame(event = names(chr_event_sum), n=unname(chr_event_sum), if_sporadic = if_else(chr_event_sum>=cutoff,'expanded','sporadic'),Sample_ID = samplename)
    pt_event_df_all = rbind(pt_event_df_all,pt_event_df)
    exp_event = pt_event_df %>% filter(if_sporadic == 'expanded') %>% pull(event)
    # Create a function to apply the filter for each events
    filter_func <- function(var) {
      var_sym <- sym(var)
      meta %>% filter(!!var_sym > 0)
    }
    if(sum(pt_event_df$if_sporadic == 'expanded')>0){
      exp_cell <- map_dfr(exp_event, filter_func) %>% pull(sample)
      meta$if_sporadic = if_else(meta$sample %in% exp_cell, 'expanded','sporadic')
    }else{
      meta$if_sporadic  = 'sporadic'
    }
    
    cell_event_df_all = rbind(cell_event_df_all,meta)
    colData(obj)$if_sporadic = meta$if_sporadic 
    
    # what events were counted in this patient?
    p = plotHeatmap_chrarm(scCNA = obj,n_threads = 100,label = 'if_sporadic',order_cells = 'phylogeny',label_colors = col['if_sporadic'],row_split = 'if_sporadic')
    pdf(paste0('heatmap_',samplename,'_CNAs.pdf'),height = h,width = 15)
    draw(p)
    dev.off()
  }
  dim(pt_event_df_all)
  
  svdir = '/volumes/USR1/yiyun/Project/HBCA/'
  setwd(svdir)
  meta_pt = readxl::read_xlsx(paste0(svdir,'/metadata/table1_metadata_01092024.xlsx'))
  meta_pt$Sample_ID  = meta_pt$sample
  pt_event_df_all = left_join(pt_event_df_all,meta_pt %>% select(-sample))
  
  wkdir = paste0(svdir,'/event_sum/cutoff_',cutoff,'/');dir.create(wkdir);setwd(wkdir)
  write_rds(pt_event_df_all,'events_sporatic_summary_pt.rds')
  dim(cell_event_df_all)
  write_rds(cell_event_df_all,'events_sporatic_summary_cell.rds')
}

#quantify by chrX, auto | expanded, nonexpaned
if(T){
  #read event matrix
  cutoff = 3
  svdir = '/volumes/USR1/yiyun/Project/HBCA/'
  setwd(svdir)
  
  wkdir = paste0(svdir,'/event_sum/cutoff_',cutoff,'/')
  cell_event_df_all = read_rds(paste0(wkdir,'events_sporatic_summary_cell.rds'))
  allevent = grep('^chr',colnames(cell_event_df_all),value = T)
  meta_event = cell_event_df_all %>% select(Sample_ID,all_of(allevent),if_chrX,if_sporadic)
  event_mat = meta_event %>% select(all_of(allevent))
  event_mat$combined_events = event_mat %>% apply(.,1,
                                                  function(x){
                                                    if(length(colnames(event_mat)[x==1])>1){
                                                      combined_events = paste(colnames(event_mat)[x==1],collapse = '.')
                                                    }else{
                                                      combined_events = colnames(event_mat)[x==1]
                                                    }
                                                    return(combined_events)
                                                  })
  event_mat$combined_events %>% table()
  meta_event$combined_events = event_mat$combined_events
  
  meta_event$if_sporadic
  meta_event$chrX_sporadic = meta_event$if_sporadic 
  meta_event$chrX_sporadic[meta_event$if_chrX == 'yes'] = 'chrX'
  table(meta_event$chrX_sporadic)
  
  meta_allcells = read_tsv( paste0(svdir,'/event_sum/all_cells_copykit_metadata.tsv'))
  allCells_df = table(meta_allcells$Sample_ID) %>% data.frame() %>% set_colnames(c('Sample_ID','allCells'))
  df_all = data.frame(Sample_ID = names(table(meta_event$Sample_ID)))
  df_all = left_join(df_all,allCells_df)
  for(i in unique(meta_event$if_sporadic)){
    meta_event_i = meta_event %>% filter(if_sporadic == i)
    
    #number of aneuploid cells
    nCells_df = meta_event_i %>% dplyr::group_by(Sample_ID) %>% dplyr::count() %>% set_colnames(c('Sample_ID','nCells'))
    nCells_df = left_join(nCells_df,allCells_df) %>% mutate(prop = nCells/allCells)
    
    #chromosomal arm event
    meta_df_w = meta_event_i %>% select(Sample_ID,all_of(allevent)) %>% group_by(Sample_ID) %>% 
      summarize_if(is.numeric, sum) %>% data.frame()
    meta_df_w = meta_df_w %>% set_rownames(meta_df_w$Sample_ID) %>% select(-Sample_ID)
    meta_df_w$n_events = apply(meta_df_w, 1, function(x){sum(x>0)})
    
    #subclone event quantify by patient
    meta_df_w_comb = table(meta_event_i$Sample_ID,meta_event_i$combined_events) %>% 
      as.matrix.data.frame() %>% 
      set_colnames(names(table(meta_event_i$combined_events))) %>% 
      set_rownames(names(table(meta_event_i$Sample_ID))) %>% 
      data.frame()
    meta_df_w_comb$n_subcl_events = apply(meta_df_w_comb, 1, function(x){sum(x>0)})
    
    #summary 
    df_i = data.frame(Sample_ID = rownames(meta_df_w_comb),
                      nCells = nCells_df$nCells,
                      prop = nCells_df$prop,
                      n_events = meta_df_w$n_events,
                      n_subcl_events =  meta_df_w_comb$n_subcl_events) %>% set_colnames(c('Sample_ID',paste0(c('nCells_','prop_','n_events_','n_subcl_events_'),i)))
    df_all = left_join(df_all,df_i)
  }
  
  for(i in unique(meta_event$if_chrX)){
    meta_event_i = meta_event %>% filter(if_chrX == i)
    
    if(i == 'yes'){
      chr_anno = 'chromoX'
    }else{
      chr_anno = 'chromoAuto'
    }
    #number of aneuploid cells
    nCells_df = meta_event_i %>% dplyr::group_by(Sample_ID) %>% dplyr::count() %>% set_colnames(c('Sample_ID','nCells'))
    nCells_df = left_join(nCells_df,allCells_df) %>% mutate(prop = nCells/allCells)
    
    #chromosomal arm event
    meta_df_w = meta_event_i %>% select(Sample_ID,all_of(allevent)) %>% group_by(Sample_ID) %>% 
      summarize_if(is.numeric, sum) %>% data.frame()
    meta_df_w = meta_df_w %>% set_rownames(meta_df_w$Sample_ID) %>% select(-Sample_ID)
    meta_df_w$n_events = apply(meta_df_w, 1, function(x){sum(x>0)})
    
    #subclone event quantify by patient
    meta_df_w_comb = table(meta_event_i$Sample_ID,meta_event_i$combined_events) %>% 
      as.matrix.data.frame() %>% 
      set_colnames(names(table(meta_event_i$combined_events))) %>% 
      set_rownames(names(table(meta_event_i$Sample_ID))) %>% 
      data.frame()
    meta_df_w_comb$n_subcl_events = apply(meta_df_w_comb, 1, function(x){sum(x>0)})
    
    #summary 
    df_i = data.frame(Sample_ID = rownames(meta_df_w_comb),
                      nCells = nCells_df$nCells,
                      prop = nCells_df$prop,
                      n_events = meta_df_w$n_events,
                      n_subcl_events =  meta_df_w_comb$n_subcl_events) %>% set_colnames(c('Sample_ID',paste0(c('nCells_','prop_','n_events_','n_subcl_events_'),chr_anno)))
    df_all = left_join(df_all,df_i)
  }
  df_all %>% head()
  df_all[is.na(df_all)]=0
  
  #add metadata
  meta_pt = read_tsv(paste0(svdir,'/metadata/table1_metadata_01092024.tsv'))#patient,Age,BMI,hyperplasia_bin,metaplasia_bin,Meno_short,Ethnicity,parity_cat,
  df_all = left_join(df_all,meta_pt)
  df_all$nAcells = df_all$nCells_chromoAuto+df_all$nCells_chromoX
  df_all$propAcells = df_all$nAcells/df_all$allCells
  df_all$age_10 = df_all$Age*0.1
  
  write_rds(df_all, paste0(svdir,'/event_sum/cutoff_',cutoff,'/events_prop_summary.rds'))
}


