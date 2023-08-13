library(readr)
library(dplyr)
library(tidyr)
library(copykit)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cowplot)
# setwd('./')

col = read_rds('./rds/col.rds')
color = read_rds('./rds/meta_color.rds')

# all events co-occurrence p value by cell
if(T){
  meta_df_w = read_tsv('./event_sum/cell_events_summary.tsv')
  dim(meta_df_w)
  meta_df_w = meta_df_w %>% filter(Sample_ID != 'BCMHBCA06R')
  colnames(meta_df_w)
  df_count = meta_df_w %>% select(-cell,-Sample_ID,-combined_events,-if_chrX) %>% data.frame()
  
  #co-occurrence, fish.exact
  co_occur_mat <- matrix(nrow = ncol(df_count), ncol = ncol(df_count)) %>% 
    set_colnames(colnames(df_count)) %>% set_rownames(colnames(df_count))
  prop_pt <- matrix(nrow = ncol(df_count), ncol = ncol(df_count)) %>% 
    set_colnames(colnames(df_count)) %>% set_rownames(colnames(df_count))
  pmat <- matrix(nrow = ncol(df_count), ncol = ncol(df_count)) %>% 
    set_colnames(colnames(df_count)) %>% set_rownames(colnames(df_count))
  omat <- matrix(nrow = ncol(df_count), ncol = ncol(df_count)) %>% 
    set_colnames(colnames(df_count)) %>% set_rownames(colnames(df_count))
  df_count[1:5,1:5]
  # Double loop to count co-occurrences for each pair of events
  for(i in 1:ncol(df_count)) {
    for(j in 1:ncol(df_count)) {
      co_occur_mat[i, j] <- sum(df_count[,i] & df_count[,j])/sum(df_count[,i])
      prop_pt[i, j] <- (meta_df_w$Sample_ID[df_count[,i] & df_count[,j]] %>% unique() %>% length())/(meta_df_w$Sample_ID %>% unique %>% length())
      test = fisher.test(df_count[,i], df_count[,j] )
      print(paste0(colnames(df_count)[i],'_',colnames(df_count)[j]))
      print(table(df_count[,i], df_count[,j] ))
      pmat[i, j]<- test$p.value
      print(test$p.value)
      omat[i, j]<- test$estimate
    }
  }
  padjmat = pmat
  p_adjust<-p.adjust(pmat,method="BH")
  for (i in 1:ncol(padjmat)){
    for (j in 1:nrow(padjmat)){
      padjmat[j,i]=p_adjust[(i-1)*nrow(padjmat)+j]
    }
  }
  p_sigmat = apply(padjmat,2,
                   function(x){
                     x[x<=0.05] ='*'
                     x[x>0.05]=''
                     return(x)})
  p_sigmat[co_occur_mat==0] = ''
  diag(p_sigmat)=''
  #basic co-occur heatmap
  library(circlize)
  col_fun = colorRamp2(c(0,0.5), c("#FFFFFF",'#DC4869FF'))
  library(ComplexHeatmap)
  Heatmap(co_occur_mat,
          col = col_fun,
          cluster_rows = F,
          cluster_columns = F,
          show_row_names = T,
          show_column_names = T,
          # row_split = rowsp,
          # column_split = colsp,
          # top_annotation = top_ha,
          # left_annotation = row_ha,
          show_row_dend = F,
          show_column_dend = F,
          row_gap = unit(2, "mm"), column_gap = unit(2, "mm"), border = T,
          name = 'co-occurrence')
  
  
  #dot heatmap
  prop_pt_trunc = apply(prop_pt,2,
                        function(x){
                          x[x>0.4]=0.4
                          return(x)})
  lvl = lapply(paste0('chr',c(1:22,"X")),function(x) paste0(x,c('p','q'))) %>% unlist()
  lvl_event = lapply(lvl,function(x) paste0(x,c('_Gain','_Loss'))) %>% unlist()
  events = lvl_event[lvl_event %in% rownames(co_occur_mat)]
  
  co_occur_mat = co_occur_mat[events,events]
  prop_pt = prop_pt[events,events]
  prop_pt_trunc = prop_pt_trunc[events,events]
  padjmat = padjmat[events,events]
  p_sigmat = p_sigmat[events,events]
  
  col_fun = colorRamp2(c(0,0.5), c("#FFFFFF",'#DC4869FF'))
  layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = 'grey', fill = NA))
    grid.circle(x=x,y=y,r= pindex(prop_pt_trunc, i, j)/0.4 * unit(2, "mm"),
                gp = gpar(fill = col_fun(pindex(co_occur_mat, i, j)), col = NA))
  }
  cell_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = 'grey', fill = NA))
    grid.circle(x=x,y=y,r= pindex(prop_pt_trunc, i, j)/0.4 * unit(2, "mm"),
                gp = gpar(fill = col_fun(pindex(co_occur_mat, i, j)), col = NA))
    grid.text(p_sigmat[i, j], x = x, y = y)
  }
  
  
  lgd_list = list(
    Legend( labels = paste0(c(0,0.1,0.2,0.3,0.4)*100, '%'), title = "patient",
            graphics = list(
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0 * unit(2, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2, "mm"),
                                               gp = gpar(fill = "black")),
              function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                               gp = gpar(fill = "black")))
    ))
  
  set.seed(123) 
  hc = hclust(dist(co_occur_mat,method = 'euclidean'),method = 'ward.D2')
  or = cutree(hc,8) %>% data.frame() %>% set_colnames('cluster')
  
  hp<- Heatmap(co_occur_mat,
               heatmap_legend_param=list(title="co-occurrence"),
               column_title = "events co-occurrence across cells", 
               col=col_fun,
               rect_gp = gpar(type = "none"),
               # layer_fun = layer_fun,
               cell_fun = cell_fun,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               cluster_rows = F,
               cluster_columns = F,
               # row_split = or$cluster,
               # column_split = or$cluster,
               border = "black")
  draw(hp, annotation_legend_list = lgd_list)
  pdf('./plot_0701/Fig3/all_events_co-occurence.pdf',width = 12,height = 10)
  draw(hp, annotation_legend_list = lgd_list)
  dev.off()
  
  save(co_occur_mat,prop_pt,padjmat,p_sigmat, 
       file = './plot_0701/Fig3/Events_co-occurence.Rdata')
}

#pick events shown in the main figure
if(T){
  # select events that is occurred in more than 3 patient
  load('./plot_0701/Fig3/Events_co-occurence.Rdata')
  prop_pt[1:5,1:5]
  co_occur_mat[1:5,1:5]
  co_occur_mat[upper.tri(co_occur_mat)]=NA
  co_occur_mat= as.matrix(co_occur_mat);dim(co_occur_mat)
  prop_pt = as.matrix(prop_pt);dim(prop_pt)
  prop_pt[upper.tri( prop_pt)]=NA
  colnames(co_occur_mat)
  colnames(prop_pt)
  
  co_occur_mat[prop_pt < (3/31)]=NA
  # co_occur_mat[padjmat > 0.05] = NA
  co_occur_mat[co_occur_mat==1]=NA
  
  matrix_data <- co_occur_mat
  dim(matrix_data)
  
  # Flatten the matrix and get the indices of ordered occurrence value
  flatten_matrix <- as.vector(matrix_data)
  indices <- order(flatten_matrix, decreasing = TRUE)
  num_rows <- nrow(matrix_data)
  values <- flatten_matrix[indices]
  row_indices <- (indices - 1) %% num_rows + 1
  column_indices <- (indices - 1) %/% num_rows + 1
  row_names <-rownames(matrix_data)[row_indices]
  column_names <-colnames(matrix_data)[column_indices] 
  
  
  # Combine the results into a data frame
  result_df <- data.frame(co_occur = values,
                          Row_Name = row_names,
                          Column_Name = column_names)
  
  # Print the result
  print(result_df)
  result_df <- result_df %>% mutate(row_chr = str_extract(Row_Name,'chr\\d{1,2}|chrX'),
                                    col_chr = str_extract(Column_Name,'chr\\d{1,2}|chrX'))
  result_df = result_df %>% filter(row_chr!=col_chr,!is.na(co_occur)) %>% head(15)
  print(result_df)
  write_rds(result_df,'./plot_0701/Fig3/picked_occurr_event')
  result_df$Row_Name %>% unique()
  result_df$Column_Name %>% unique()
}

# picked eventsco-occurrence heatmap 
if(T){
  load('./plot_0701/Fig3/Events_co-occurence.Rdata')
  result_df = read_rds('./plot_0701/Fig3/picked_occurr_event')
  
  meta_df_w = read_tsv('./event_sum/cell_events_summary.tsv')
  dim(meta_df_w)
  meta_df_w = meta_df_w %>% filter(Sample_ID != 'BCMHBCA06R')
  colnames(meta_df_w)
  df_count = meta_df_w %>% select(-cell,-Sample_ID,-combined_events,-if_chrX) %>% data.frame()
  
  fil = c(unique(c(result_df$Row_Name,result_df$Column_Name)))
  # double check: 
  # grep('chr1q_Gain|chr16q_Loss|chr12p_Loss|chr10q_Loss|chr5p_Gain|chr7q_Loss|chr8p_Loss|chr22q_Loss|chr11q_Loss|chr14q_Gain|chr15q_Gain|chr17p_Loss|chr21q_Loss|chrXp_Loss|chrXp_Gain|chr3p_Loss',colnames(df_count),value = T)
  
  df_count = df_count[fil,fil]
  
  #co-occurrence, fish.exact
  co_occur_mat_fil <- co_occur_mat[fil,fil]
  prop_pt_fil <- prop_pt[fil,fil]
  p_sigmat_fil <- p_sigmat[fil,fil]
  
  #dot heatmap
  col_fun = colorRamp2(c(0,0.5), c("#FFFFFF",'#DC4869FF'))
  prop_pt_trunc_fil = apply(prop_pt_fil,2,
                            function(x){
                              x[x>0.4]=0.4
                              return(x)})
  
  #subset
  from = c('chr1q_Gain','chr7q_Loss','chr8p_Loss','chr14q_Gain','chr17p_Loss')
  to = colnames(co_occur_mat_fil)[colnames(co_occur_mat_fil) %!in% from]
  co_occur_mat = co_occur_mat_fil[to,from] %>% t()
  prop_pt_trunc = prop_pt_trunc_fil[to,from]%>% t()
  p_sigmat = p_sigmat_fil[to,from]%>% t()
  
  #set heatmap parameter
  col_fun = colorRamp2(c(0,0.5), c("#FFFFFF",'#DC4869FF'))
  layer_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = 'grey', fill = NA))
    grid.circle(x=x,y=y,r= pindex(prop_pt_trunc, i, j)/0.4 * unit(2, "mm"),
                gp = gpar(fill = col_fun(pindex(co_occur_mat, i, j)), col = NA))
  }
  cell_fun = function(j, i, x, y, w, h, fill){
    grid.rect(x = x, y = y, width = w, height = h, 
              gp = gpar(col = 'grey', fill = NA))
    grid.circle(x=x,y=y,r= pindex(prop_pt_trunc, i, j)/0.4 * unit(2, "mm"),
                gp = gpar(fill = col_fun(pindex(co_occur_mat, i, j)), col = NA))
    grid.text(p_sigmat[i, j], x = x, y = y)
  }
  
  
  lgd_list = list(
    Legend(labels = paste0(c(0.01,0.1,0.2,0.3,0.4)*100, '%'), title = "patient",
           graphics = list(
             function(x, y, w, h) grid.circle(x = x, y = y, r = 0.025 * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x = x, y = y, r = 0.25 * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x = x, y = y, r = 0.5 * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x = x, y = y, r = 0.75 * unit(2, "mm"),
                                              gp = gpar(fill = "black")),
             function(x, y, w, h) grid.circle(x = x, y = y, r = 1 * unit(2, "mm"),
                                              gp = gpar(fill = "black")))
    ))
  
  set.seed(123) 
  # hc = hclust(dist(co_occur_mat,method = 'euclidean'),method = 'ward.D2')
  # or = cutree(hc,8) %>% data.frame() %>% set_colnames('cluster')
  
  hp<- Heatmap(co_occur_mat,
               heatmap_legend_param=list(title="co-occurrence"),
               column_title = "events co-occurrence across cells", 
               col=col_fun,
               rect_gp = gpar(type = "none"),
               # layer_fun = layer_fun,
               cell_fun = cell_fun,
               row_names_gp = gpar(fontsize = 5),
               column_names_gp = gpar(fontsize = 5),
               cluster_rows = F,
               cluster_columns = F,
               # row_split = or$cluster,
               # column_split = or$cluster,
               border = "black")
  draw(hp, annotation_legend_list = lgd_list)
  pdf('./plot_0701/Fig3/C_picked_Events_co-occurence.pdf',width = 5,height = 2.6)
  draw(hp, annotation_legend_list = lgd_list)
  dev.off()
}

