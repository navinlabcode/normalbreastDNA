library(readr)
library(dplyr)
library(tidyr)
library(copykit)
library(stringr)
library(magrittr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(rstatix)
library(ggthemes)


### multivariate regression of n_events, prop and all other meta---------------------------

if(T){
  meta_event=read_tsv('./plot_0701/Fig2/event_meta_stat_perpt.tsv')
  meta_event %>% head()
  dim(meta_event)
  meta_event = meta_event %>% select(paper_id,n_subcl_events:BMI)
  
  dim(meta_event)
  
  #set levels for categorical variables
  meta_event$age_cat = factor(meta_event$age_cat,levels = c('<40','40-50','>50'))
  meta_event$bmi_cat = factor(meta_event$bmi_cat,levels = c('18.5-25','25-30','>30'))
  meta_event$hyperplasia_bin = factor(meta_event$hyperplasia_bin,levels = c('Not_observed','HP'))
  meta_event$metaplasia_bin = factor(meta_event$metaplasia_bin,levels = c('Not_observed','MP'))
  meta_event$Meno_short = factor(meta_event$Meno_short,levels = c('Pre','On','Post'))
  meta_event$Ethnicity = factor(meta_event$Ethnicity,levels = c('African-American','Caucasian','Hispanic'))
  meta_event$parity_cat = factor(meta_event$parity_cat,levels =c('P0','P>0'))
  
  colnames(meta_event)
  continuous_var = meta_event %>% select(n_subcl_events, prop) %>% colnames() #n_events


 
  #n_events, prop - multivariate, 
  # remove metaplasia and hyperplasia since the only 2 cases HP,MP.

  # con = 'n_events'
  meta_event = meta_event %>% filter(paper_id != 'P33')
  meta_event$age_10 = meta_event$Age*0.1
  p_all =list()
  for(con in continuous_var){
    dependent = con
    explanatory =colnames(meta_event)[c(13:17,19:20)]
    form = as.formula(paste0(dependent,' ~',paste(explanatory,collapse = '+')))
    print(form)
    if(con %in% c('n_subcl_events','n_events')){
      model <- glm(formula = form,
                   family = 'poisson',
                   data = meta_event,na.action = 'na.omit')
      summary(model) 
    }else{
      model <- glm(formula = form,
                   family = 'gaussian',
                   data = meta_event,na.action = 'na.omit')
      summary(model) 
    }
    
    tidy_glm_Multi = broom::tidy(model,conf.int = TRUE, conf.level = 0.95) %>% 
      mutate(sig = case_when( p.value < .001 ~ "***", p.value < .01 ~ "**", p.value < .05 ~ "*", TRUE ~ "" )) %>% 
      filter(term!='(Intercept)') %>% 
      mutate(type = 'Multivariate',variable = term)
    tidy_glm_Multi %>% print()
    tidy_glm_ggforest = tidy_glm_Multi %>%
      mutate(name = term, beta = estimate,se = std.error)
    
    #create plot
    library(ggforestplot)
    p_all[[con]] = tidy_glm_ggforest%>% 
      ggplot(aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
      geom_pointrange() +
      geom_text(aes(label = sig,size = 4,fontface = 'bold'), vjust = -0.5) + # Display stars next to the points
      geom_hline(yintercept = 0, linetype = "dashed") +
      coord_flip() +
      theme_bw() +
      theme(legend.position = 'None')+
      xlab("") + # Remove x-axis label
      ylab("log(OR)") + # Add y-axis label
      ggtitle(con) # Add title
    
    library(forestmodel)
    p_all[[con]] = forestmodel::forest_model(model,
                                             exponentiate = F,
                                             format_options = forest_model_format_options( colour="black", 
                                                                                           color=NULL, shape=20, 
                                                                                           text_size=5, point_size=3, 
                                                                                           banded=TRUE ))
    
  }
  
  pdf('./plot_0701/FigS2/Forestplot_n_events_prop_association_to_all_multivariables.pdf',width = 12, height = 9)
  plot_grid(plotlist = p_all,ncol=1) %>% print()
  dev.off()
  
}