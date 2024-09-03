library(readr)
library(dplyr)
library(tidyr)
library(copykit)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cowplot)

svdir = '/volumes/USR1/yiyun/Project/HBCA'
# generate matrix for association analysis
if(T){
  #events count, type of combined events (subclone)
  meta_df_w <-read_tsv(paste0(svdir,'/event_sum/cell_events_summary.tsv')) %>% filter(Sample_ID !='BCMHBCA06R')
  df_cells = meta_df_w$Sample_ID %>% table() %>% data.frame() %>% set_colnames(c('Sample_Id','nCells'))
  meta_event = table(meta_df_w$Sample_ID,meta_df_w$combined_events) %>% 
    as.matrix.data.frame() %>% 
    set_colnames(names(table(meta_df_w$combined_events))) %>% 
    set_rownames(names(table(meta_df_w$Sample_ID))) %>% 
    data.frame()
  meta_event$n_subcl_events = apply(meta_event, 1, function(x){sum(x>0)})
  dim(meta_event)
  #get number of events(chr arm level)
  df_events = meta_df_w %>% dplyr::group_by(Sample_ID) %>% dplyr::summarise_if(is.numeric,sum)
  meta_event$n_events = apply(df_events,1,function(x){sum(x>0)})
  meta_event$Sample_ID =rownames(meta_event)
  #number of cells, proportion of cells
  meta_allcells = read_tsv(paste0(svdir,'/event_sum/all_cells_copykit_metadata.tsv'))
  meta_aneu = table(meta_allcells$is_aneuploid,meta_allcells$Sample_ID) %>% 
    data.frame() %>% 
    set_colnames(c('is_aneuploid','Sample_ID','nCells')) %>% 
    dplyr::group_by(Sample_ID) %>% mutate(prop = nCells/sum(nCells)) %>% filter(is_aneuploid == T)
  meta_event = left_join(meta_event,meta_aneu)
  head(meta_event)
  meta_event$patient = str_extract(meta_event$Sample_ID,'BCMHBCA\\d{1,3}')
  meta_event %>% filter(Sample_ID == 'BCMHBCA104R')
  
  #patient info
  meta_pt = read_tsv(paste0(svdir,'/metadata/table1_metadata_01092024.tsv'))
  meta_pt %>% filter(Sample_ID == 'BCMHBCA104R')
  dim(meta_pt)
  head(meta_pt)
  meta_event = left_join(meta_event,meta_pt)
  head(meta_event)
  meta_event %>% arrange(-n_subcl_events) 
  
  meta_event %>% select(paper_id,prop,Age) %>% arrange(Age)
  median(meta_event$Age)
  meta_event %>% select(paper_id,prop,Age) %>% arrange(Age) %>% filter(Age==43) %>% pull(prop) %>% mean()
  write_tsv(meta_event,'./event_meta_stat_perpt.tsv')
}
