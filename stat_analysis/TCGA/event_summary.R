library(readr)
library(dplyr)
library(tidyr)
library(copykit)
library(stringr)

setwd('/volumes/USR1/yiyun/Project/HBCA/')


#make summary of patient and event from
load('./cbs_0501_v1/20230608_step2_plot.Rdata')
if(T){
  meta_chr = read_rds('./event_sum/filtered_cells_w_event.rds')
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
  
  meta_chr_w%>% write_rds(.,'./event_sum/cell_events_summary.rds')
  meta_chr_w%>% write_tsv(.,'./event_sum/cell_events_summary.tsv')
  meta_chr_w$cell %>% duplicated() %>% table()
  table(colData(filtered_obj)$sample %in% meta_chr_w$cell)
  
  # prepare all cells metadata
  meta = filtered_obj %>% colData() %>% data.frame()
  meta_pt = read_tsv('./metadata/table1_metadata_05302023.tsv')
  meta = left_join(meta,meta_pt)
  is.na(meta$paper_id) %>% table()
  meta %>% write_tsv(.,'./event_sum/all_cells_copykit_metadata.tsv') 
  
  # add patient info to metadata
  df = colData(filtered_aneu_obj) %>% data.frame()
  
  dim(df)
  df$Sample_ID %>% unique()
  meta_pt = read_tsv('./metadata/table1_metadata_05302023.tsv')
  meta_pt$assay %>% table()
  
  df = left_join(df,meta_pt)
  dim(df)
  df$paper_id %>% is.na() %>% table()
  colData(filtered_aneu_obj) = cbind(colData(filtered_aneu_obj),
                                     df %>% select(-colnames(colData(filtered_aneu_obj) %>% data.frame())))
  
  meta_chr_w = read_rds('./event_sum/cell_events_summary.rds')
  meta_chr_w$sample=meta_chr_w$cell
  colnames(meta_chr_w)
  
  data = meta_chr_w %>% dplyr::group_by(Sample_ID,combined_events) %>% dplyr::count()
  # data$combined_events[data$combined_events %!in% keep] = 'others'
  # data$combined_events = factor(data$combined_events, levels = c(keep,'others'))
  data$patient = str_extract(string = data$Sample_ID,pattern = 'BCMHBCA\\d{1,2}')
  dim(meta_pt)
  meta_pt = meta_pt %>% filter(!duplicated(patient))
  data = left_join(data,meta_pt) %>% arrange(-Age)
  head(data)
  data$paper_id = factor(data$paper_id,levels = c(unique(data$paper_id)))
  write_rds(data,'./event_sum/cell_events_summary_long.rds')
  
}