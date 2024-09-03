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

# multivariate regression of n_events, prop to all collected clinical metadata 
if(T){
  meta_event=read_tsv('./event_meta_stat_perpt.tsv')%>% filter(paper_id != 'P50')
  meta_event %>% head()
  dim(meta_event)
  meta_event = meta_event %>% select(paper_id,n_subcl_events:BMI)
  meta_event %>% head()
  dim(meta_event)
  
  
  #set levels for categorical variables;
  #category with only one patient will be removed
  meta_event$hyperplasia_bin = factor(meta_event$hyperplasia_bin,levels = c('Not_observed','HP'))
  meta_event$metaplasia_bin = factor(meta_event$metaplasia_bin,levels = c('Not_observed','MP'))
  meta_event$Meno_short = factor(meta_event$Meno_short,levels = c('Pre','Post'))#remove 'Unknown','On'
  meta_event$Ethnicity = factor(meta_event$Ethnicity,levels = c('Caucasian','African-American','Hispanic')) #remove 'Asian','AA-Hispanic'
  meta_event$parity_cat = factor(meta_event$parity_cat,levels = c('P0','P>0'))# remove 'Unknown'
  meta_event %>% select(age_cat,bmi_cat,hyperplasia_bin, metaplasia_bin, Meno_short, Ethnicity, parity_cat) %>% apply(.,2, table)
  meta_event$age_10 = meta_event$Age*0.1
  
  colnames(meta_event)
  continuous_var = meta_event %>% select(n_subcl_events, prop) %>% colnames() #n_events
  
  #n_events, prop - multivariate, 
  plist =list()
  # con = 'n_events'
  for(con in continuous_var){
    dependent = con
    explanatory =c("hyperplasia_bin","metaplasia_bin","Meno_short","Ethnicity","parity_cat","BMI","age_10")
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
    
    
    library(forestmodel)
    plist[[con]] = forestmodel::forest_model(model,
                                              exponentiate = F,
                                              format_options = forest_model_format_options( colour="black", 
                                                                                            color=NULL, shape=20, 
                                                                                            text_size=5, point_size=3, 
                                                                                            banded=TRUE ))
  }
  
  plist$n_subcl_events
  plist$prop
}
