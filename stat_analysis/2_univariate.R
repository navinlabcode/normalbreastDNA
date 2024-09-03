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
# Spearman's rank association: association of age,BMI and number of events, proportion of cells
if(T){
  meta_event=read_tsv('./event_meta_stat_perpt.tsv') %>% filter(paper_id != 'P50')
  
  #Age,BMI correlation with events and cell proportion 
  continuous_var = meta_event %>% select(Age,BMI) %>% colnames()
  plist = list()
  for(con in continuous_var){

    p1 = ggscatter(meta_event, x = con, y = "prop",
                   color = 'black', palette = "jco",
                   add = "reg.line",                                 
                   conf.int = TRUE,
                   add.params = list(color = "black",
                                     fill = "lightgray"))+ 
      ylab('proportion of aneuploid cells')+
      stat_cor(method = "spearman");
    
    
    p2 = ggscatter(meta_event, x = con, y = "n_subcl_events",
                   color = 'black', palette = "jco",
                   add = "reg.line",                                 
                   conf.int = TRUE,
                   add.params = list(color = "black",
                                     fill = "lightgray"))+ 
      ylab('number of events (subclone)')+
      stat_cor(method = "spearman");
    
    plist[[con]] =  plot_grid(p1,p2,ncol = 2)
  }
  
  #print association plot
  plist$Age
  plist$BMI
}

# Spearman's rank association: Autosomal vs chrX; Expanded vs Sporadic to Age association
if(T){
  df_all = read_rds(paste0(svdir,'/event_sum/cutoff_3/events_prop_summary.rds')) %>% filter(paper_id!='P50')
  df_all %>% filter(Sample_ID == 'BCMHBCA104R')
  df_all %>% head()
  
  plist = list()
  for(i in c('chromoX','chromoAuto','sporadic','expanded')){
    continuous_var_x = df_all %>% select(Age) %>% colnames()
    continuous_var_y = paste0(c('prop_','n_subcl_events_'),i)
    for(x_ in continuous_var_x){
      p1 = ggscatter(df_all, x = x_, y = continuous_var_y[1],
                     color = 'black', palette = "jco",
                     add = "reg.line",                                 
                     conf.int = TRUE,
                     add.params = list(color = "black",
                                       fill = "lightgray"))+ 
        ylab(paste0('number of aneuploid cells ',i))+
        stat_cor(method = "spearman");p1
      
      p2 = ggscatter(df_all, x = x_, y = continuous_var_y[2],
                     color = 'black', palette = "jco",
                     add = "reg.line",                                 
                     conf.int = TRUE,
                     add.params = list(color = "black",
                                       fill = "lightgray"))+ 
        ylab(paste0('proportion of aneuploid cells ',i))+
        stat_cor(method = "spearman");p2
      
      
      plist[[paste0(i,'_',x_)]] = plot_grid(p1,p2,ncol = 2)  
    }
  }
  
  plot_grid(plist$chromoX_Age,plist$chromoAuto_Age,ncol=1)
  plot_grid(plist$expanded_Age,plist$sporadic_Age,ncol=1)
}

# wilcoxon rank sum test: association of other categorical clinical metadata and number of events, proportion of cells
if(T){
  meta_event=read_tsv('./event_meta_stat_perpt.tsv') %>% filter(paper_id != 'P50')
  meta_event %>% head()
  meta_event$Ethnicity %>% table()
  colnames(meta_event)
  meta_event = meta_event %>% select(paper_id,n_subcl_events:parity_cat) %>% data.frame()
  meta_event %>% head()
  
  
  #set levels for categorical variables;
  #category with only one patient will be removed
  meta_event$age_cat = factor(meta_event$age_cat,levels = c('<40','40-50','>50'))
  meta_event$bmi_cat = factor(meta_event$bmi_cat,levels = c('18.5-25','25-30','>30')) #remove '<18.5'
  meta_event$hyperplasia_bin = factor(meta_event$hyperplasia_bin,levels = c('Not_observed','HP'))
  meta_event$metaplasia_bin = factor(meta_event$metaplasia_bin,levels = c('Not_observed','MP'))
  meta_event$Meno_short = factor(meta_event$Meno_short,levels = c('Pre','Post'))#remove 'Unknown','On'
  meta_event$Ethnicity = factor(meta_event$Ethnicity,levels = c('African-American','Caucasian','Hispanic')) #remove 'Asian','AA-Hispanic'
  meta_event$parity_cat = factor(meta_event$parity_cat,levels = c('P0','P>0'))# remove 'Unknown'
  meta_event %>% select(age_cat,bmi_cat,hyperplasia_bin, metaplasia_bin, Meno_short, Ethnicity, parity_cat) %>% apply(.,2, table)
  
  #for numeric variable and metadata
  colnames(meta_event)
  category_var = meta_event %>% select(hyperplasia_bin:parity_cat) %>% colnames()
  continuous_var = meta_event %>% select(n_subcl_events, prop) %>% colnames() 
  
  
  plist = list()
  library(rstatix)
  library(ggthemes)
  for(con in continuous_var){
    temp = meta_event %>% select(all_of(category_var),all_of(con),patient) %>% 
      pivot_longer(all_of(category_var),names_to = 'category',values_to='cat_value')%>%
      filter(!is.na(cat_value)) %>% set_colnames(c('con_value','pattient','category','cat_value'))
    
    stat.test <- temp  %>% 
      group_by(category) %>%
      wilcox_test(as.formula('con_value~ cat_value')) %>% 
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj") %>% 
      mutate(print = paste0('p=',round(p,3),', (', p.adj.signif,', adjusted)'), p.adj =round(p.adj,3))
    y = temp %>% dplyr::group_by(category) %>% dplyr::summarise(y.position = max(con_value))
    stat.test = left_join(stat.test,y)
    
    for(cat in category_var){
      temp_sub = temp %>% filter(category == cat,!is.na(cat_value))
      p_temp = ggboxplot(temp_sub, x = 'cat_value', y = 'con_value',
                         color = 'black', palette = "jco",
                         add = "jitter")+
        xlab(cat)+
        ylab(con)+
        theme_clean()+
        stat_pvalue_manual(stat.test %>% filter(category == cat), label = "p.adj",size = 4,
                           step.increase = 0.15)+
        # facet_wrap(~category,nrow = 1)+
        ggtitle(cat)+
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = 'None',
              panel.border = element_blank(),
              strip.background = element_rect(fill = NA, colour = NA),
              strip.text = element_text(size = 12, color = "black", face = "bold"))
      plist[[paste0(con,'_',cat)]]=p_temp
    }
  }
  names(plist)
  plot_grid(plotlist = plist[grepl('prop_',names(plist))],nrow = 1)
  plot_grid(plotlist = plist[grepl('n_subcl_',names(plist))],nrow = 1)
}   
