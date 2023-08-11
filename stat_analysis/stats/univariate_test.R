library(readr)
library(dplyr)
library(tidyr)
library(copykit)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cowplot)

setwd('/volumes/USR1/yiyun/Project/HBCA/')

#age,BMI as numeric variable association with number of events, number of cells, proportion of cells
if(T){
  #events count, type of combined events (subclone)
  meta_df_w <-read_tsv('./event_sum/cell_events_summary.tsv')
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
  meta_allcells = read_tsv('./event_sum/all_cells_copykit_metadata.tsv')
  meta_aneu = table(meta_allcells$is_aneuploid,meta_allcells$Sample_ID) %>% 
    data.frame() %>% 
    set_colnames(c('is_aneuploid','Sample_ID','nCells')) %>% 
    dplyr::group_by(Sample_ID) %>% mutate(prop = nCells/sum(nCells)) %>% filter(is_aneuploid == T)
  meta_event = left_join(meta_event,meta_aneu)
  head(meta_event)
  meta_event$patient = str_extract(meta_event$Sample_ID,'BCMHBCA\\d{1,2}')
  
  #patient info
  meta_pt = read_tsv('./metadata/table1_metadata_05302023.tsv')
  dim(meta_pt)
  head(meta_pt)
  meta_event = left_join(meta_event,meta_pt)
  head(meta_event)
  write_tsv(meta_event,'./plot_0701/Fig2/event_meta_stat_perpt.tsv')
  write_rds(meta_event,'./plot_0701/Fig2/event_meta_stat_perpt.rds')
  
  #Age,BMI correlation with events and cell proportion 
  p1 = ggscatter(meta_event%>% filter(paper_id!='P33'), x = 'Age', y = "n_subcl_events",
                 color = 'black', palette = "jco",
                 add = "reg.line",                                 
                 conf.int = TRUE,
                 add.params = list(color = "black",
                                   fill = "lightgray"))+ 
    ylab('number of events (subclone)')+
    stat_cor(method = "spearman");p1
  p2 = ggscatter(meta_event%>% filter(paper_id!='P33'), x = 'Age', y = "prop",
                 color = 'black', palette = "jco",
                 add = "reg.line",                                 
                 conf.int = TRUE,
                 add.params = list(color = "black",
                                   fill = "lightgray"))+ 
    ylab('proportion of aneuploid cells')+
    stat_cor(method = "spearman");p2
  p_pick = plot_grid(p1,p2,ncol = 2)
  
  pdf('./plot_0701/Fig2/B_num_events_association_to_age_pick.pdf',width = 6, height = 3)
  print(p_pick)
  dev.off()
  
  
  p3 = ggscatter(meta_event%>% filter(paper_id!='P33'), x = 'BMI', y = "n_subcl_events",
                 color = 'black', palette = "jco",
                 add = "reg.line",                                 
                 conf.int = TRUE,
                 add.params = list(color = "black",
                                   fill = "lightgray"))+ 
    ylab('number of events (subclone)')+
    stat_cor(method = "spearman");p3
  p4 = ggscatter(meta_event%>% filter(paper_id!='P33'), x = 'BMI', y = "prop",
                 color = 'black', palette = "jco",
                 add = "reg.line",                                 
                 conf.int = TRUE,
                 add.params = list(color = "black",
                                   fill = "lightgray"))+ 
    ylab('proportion of aneuploid cells')+
    stat_cor(method = "spearman");p4
  p_pick_s = plot_grid(p3,p4,ncol = 2)
  
  pdf('./plot_0701/FigS2/A_num_events_association_to_BMI_pick.pdf',width = 6.5, height = 3)
  print(p_pick_s)
  dev.off()
}

#age as a categorical variable association with number of events, number of cells, proportion of cells
if(T){
  meta_event = read_rds('./plot_0701/Fig2/event_meta_stat_perpt.rds')
  meta_event$age_cat = factor(meta_event$age_cat,levels = c('<40','40-50','>50'))
  
  continuous_var = meta_event %>% select(n_subcl_events, prop) %>% colnames() #n_events
  cat = 'age_cat'
  plist = list()
  library(rstatix)
  library(ggthemes)
  temp = meta_event %>% filter(paper_id != 'P33') %>% select(all_of(cat),all_of(continuous_var),patient) %>% 
    filter(!is.na(!!sym(cat))) %>% pivot_longer(all_of(continuous_var),names_to = 'Continous',values_to='Observance')%>%
    filter(!is.na(Observance)) 
  stat.test <- temp  %>% 
    group_by(Continous) %>%
    wilcox_test(as.formula(paste0('Observance ~', cat))) 
  y = temp %>% dplyr::group_by(Continous) %>% dplyr::summarise(y.position = max(Observance))
  stat.test = left_join(stat.test,y)
  
  for(con in continuous_var){
    temp_sub = temp %>% filter(Continous == con,!is.na(!!sym(cat)))
    p_temp = ggboxplot(temp_sub, x = all_of(cat), y = 'Observance',
                       color = 'black', palette = "jco",
                       add = "jitter")+
      xlab(cat)+
      ylab(con)+
      theme_clean()+
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'None',
            panel.border = element_blank())+
      stat_pvalue_manual(stat.test %>% filter(Continous == con), label = "p",size = 4,
                         step.increase = 0.15)+
      facet_wrap(~Continous,nrow = 1)
    plist[[paste0(con,'_',cat)]]=p_temp
  }
  length(plist)

  pdf('./plot_0701/Fig2/C_num_events_association_to_age_cat.pdf',width = 6.5, height = 4)
  print(plot_grid(plotlist = plist[1:2], nrow = 1))
  dev.off()
}

#categorical variables association with number of events, number of cells, proportion of cells
if(T){
  meta_event=read_tsv('./plot_0701/Fig2/event_meta_stat_perpt.tsv')
  meta_event %>% head()
  colnames(meta_event)
  meta_event = meta_event %>% select(paper_id,n_subcl_events:parity_cat)
  
  #set levels for categorical variables
  meta_event$hyperplasia_bin = factor(meta_event$hyperplasia_bin,levels = c('Not_observed','HP'))
  meta_event$metaplasia_bin = factor(meta_event$metaplasia_bin,levels = c('Not_observed','MP'))
  meta_event$Meno_short = factor(meta_event$Meno_short,levels = c('Pre','On','Post'))#'Meno','Unknown' removed
  meta_event$Ethnicity = factor(meta_event$Ethnicity,levels = c('African-American','Caucasian','Hispanic'))
  meta_event$parity_cat = factor(meta_event$parity_cat,levels =c('P0','P>0'))#,'Unknown'removed
  
  #for numeric variable and metadata
  colnames(meta_event)
  category_var = meta_event %>% select(hyperplasia_bin:parity_cat) %>% colnames()
  continuous_var = meta_event %>% select(n_subcl_events, prop) %>% colnames() #n_events
  
  plist = list()
  library(rstatix)
  library(ggthemes)
  for(cat in category_var){
    temp = meta_event %>% filter(paper_id != 'P33') %>% select(all_of(continuous_var),all_of(cat),patient) %>% 
      filter(!is.na(!!sym(cat))) %>% pivot_longer(all_of(continuous_var),names_to = 'Continous',values_to='Observance')%>%
      filter(!is.na(Observance)) 
    stat.test <- temp  %>% 
      group_by(Continous) %>%
      wilcox_test(as.formula(paste0('Observance ~', cat))) %>% 
      adjust_pvalue(method = "bonferroni") %>%
      add_significance("p.adj") %>% 
      mutate(print = paste0('p=',round(p,3),', (', p.adj.signif,', adjusted)'))
    y = temp %>% dplyr::group_by(Continous) %>% dplyr::summarise(y.position = max(Observance))
    stat.test = left_join(stat.test,y)
    
    for(con in continuous_var){
      temp_sub = temp %>% filter(Continous == con,!is.na(!!sym(cat)))
      p_temp = ggboxplot(temp_sub, x = all_of(cat), y = 'Observance',
                         color = 'black', palette = "jco",
                         add = "jitter")+
        xlab(cat)+
        ylab(con)+
        theme_clean()+
        stat_pvalue_manual(stat.test %>% filter(Continous == con), label = "p",size = 4,
                           step.increase = 0.15)+
        facet_wrap(~Continous,nrow = 1)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = 'None',
              panel.border = element_blank(),
              strip.background = element_rect(fill = NA, colour = NA),
              strip.text = element_text(size = 12, color = "black", face = "bold"))
      plist[[paste0(con,'_',cat)]]=p_temp
    }
  }
  length(plist)
  
  pdf('./plot_0701/FigS2/B1_num_events_association_to_allcategorical.pdf',width = 10, height = 3.5)
  print(plot_grid(plist$n_subcl_events_hyperplasia_bin,plist$prop_hyperplasia_bin, 
                  plist$n_subcl_events_metaplasia_bin,plist$prop_metaplasia_bin, 
                  nrow = 1))
  dev.off()
  pdf('./plot_0701/FigS2/B2_num_events_association_to_allcategorical.pdf',width = 10, height = 3)
  print(plot_grid(plist$n_subcl_events_parity_cat,plist$prop_parity_cat, plist$n_subcl_events_Meno_short,plist$prop_Meno_short,
                  nrow = 1))
  dev.off()
  pdf('./plot_0701/FigS2/B3_num_events_association_to_allcategorical.pdf',width = 5, height = 3.7)
  print(plot_grid(plist$n_subcl_events_Ethnicity,plist$prop_Ethnicity,
                  nrow = 1))
  dev.off()
}