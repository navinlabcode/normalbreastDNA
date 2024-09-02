
#load package
library(dplyr)
library(gsDesign)
library(readr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(copykit)#V0.1.0
library(pwr)
library(Seurat)
source('./Copykit_copykat_functions.R')

svdir = '/volumes/USR1/yiyun/Project/HBCA/copykat_ST/'


#make copykit object and merge objects for all ST patients 
if(T){
  setwd(paste0(svdir,'/improve_copykat/'))
  ST_df = read.csv(paste0(svdir,'all_ST_dir.csv'),header = T)
  
  for(sample.name in ST_df$sample){
    tryCatch({
      #set object saving directory
      # sample.name = 'BCMHBCA60L1_061223'
      cpk_svdir = paste0(svdir,'/improve_copykat/',sample.name)
      
      #read copykat result
      copykatdata <- read.table(paste0(svdir,'/win_25_ks_0.05/all_ST_copykat_result/',sample.name,'_copykat_CNA_results.txt'),header = T)
      copykatdata[1:5,1:5]
      dim(copykatdata)
      message(sample.name)
      
      #make copykit object
      cpkobj <- Copykit_copykat(CNA_result = copykatdata, genome_version = 'hg19',resolution = '200k',method = 'multipcf')
      colData(cpkobj) %>% head()
      colData(cpkobj)$sample.name = sample.name
      message(ncol(cpkobj))
      
      #merge samples
      colnames(cpkobj) = paste0(sample.name,colnames(cpkobj))
      #BCMHBCA01L_15L_061223 is the first sample name
      if(sample.name == 'BCMHBCA01L_15L_061223'){
        cpkobj_int = cpkobj
      }else{
        cpkobj_int = cbind(cpkobj_int,cpkobj)
      }
      message(ncol(cpkobj_int))
    }, 
    error=function(e) {cat("ERROR :",conditionMessage(e), "\n")})
  }
  setwd(paste0(svdir,'/improve_copykat/int/'))
  write_rds(cpkobj_int, 'cpkobj_int_allpatients.rds')
}


#check distribution of segmentation for chromosomes by sample 
if(T){
  setwd(paste0(svdir,'/improve_copykat/int/'))
  cpkobj_int = read_rds('cpkobj_int_allpatients.rds')
  
  meta = colData(cpkobj_int) %>% data.frame()
  chr_arm = c('chr1q','chr10q','chr16q','chr7q','chr22q','chrXp','chrXq') #most frequent event
  df_state = data.frame(chr_arm = c('chr1q','chr10q','chr16q','chr7q','chr22q','chrXp','chrXq'),
                        state = c('Gain','Loss','Loss','Loss','Loss','Loss','Loss'))
  hist_list = list()
  para_list = list()
  # arm_i = 'chr1q'
  for(arm_i in chr_arm){
    message(arm_i)
    ghist = NULL
    hist_data_all = meta %>% dplyr::select(all_of(arm_i),sample.name) %>%
      set_colnames(c('logRatio','sample.name')) %>% filter(!is.na(logRatio))
    hist_data_all %>% head() %>% print()
    
    hist_data_list = split(hist_data_all,hist_data_all$sample.name)
    # hist_data = hist_data_list$BCMHBCA52R1 #check data here
    
    hist_list_armi = lapply(hist_data_list,function(hist_data){
      #fit model 
      fit_2 <- Mclust(hist_data$logRatio, G = 2)
      message(paste0('2 gaussion model bic is ',fit_2$bic))
      fit = fit_2
      params =fit$parameters
      if(length(params$variance$sigmasq)==1){
        params$variance$sigmasq=c(params$variance$sigmasq,params$variance$sigmasq)
      }
      #parameters of mix distribution
      para_df = data.frame(miu1 = params$mean[1],sd1 = sqrt(params$variance$sigmasq[1]),
                           miu2 = params$mean[2],sd2 = sqrt(params$variance$sigmasq[2]),
                           pro1 = params$pro[1], pro2 = params$pro[2],
                           diff = abs(params$mean[1]-params$mean[2]),
                           chr_arm = arm_i)
      para_df = left_join(para_df,df_state)
      #the miu of distribution whco is more far away from 0, then that distribution is the H1 
      if(abs(para_df$miu1)>abs(para_df$miu2)){
        #gain
        z0_0.5 = qnorm(0.5, mean = para_df$miu2, sd = para_df$sd2)#H0_z0.5
        z1_0.5 = qnorm(0.5, mean = para_df$miu1, sd = para_df$sd1)#H1_z0.5
        z1_0.6 = qnorm(0.6, mean = para_df$miu1, sd = para_df$sd1)#H1_z0.6
        sd1 = para_df$sd1
        pro1 = para_df$pro1
        H0_95CI = qnorm(c(0.05), mean = para_df$miu2, sd = para_df$sd2)
        H1_95CI = qnorm(c(0.05), mean = para_df$miu1, sd = para_df$sd1)
        H0_pro = para_df$pro2
        H1_pro = para_df$pro1
      }else{
        #loss
        z0_0.5 = qnorm(0.5, mean = para_df$miu1, sd = para_df$sd1)#H0_z0.5
        z1_0.5 = qnorm(0.5, mean = para_df$miu2, sd = para_df$sd2)#H1_z0.5
        z1_0.6 = qnorm(0.6, mean = para_df$miu2, sd = para_df$sd2)#H1_z0.6
        sd1 = para_df$sd2
        pro1 = para_df$pro2
        H0_95CI = qnorm(c(0.95), mean = para_df$miu1, sd = para_df$sd1)
        H1_95CI = qnorm(c(0.95), mean = para_df$miu2, sd = para_df$sd2)
        H0_pro = para_df$pro1
        H1_pro = para_df$pro2
      }
      para_df$z1_0.6 = z1_0.6
      
      #plot
      ghist <- ggplot(hist_data, aes(x = logRatio)) +
        geom_histogram(aes(y = ..density..), binwidth = 0.002, color = "black", fill = "white")+
        xlim(-0.2,0.2)+
        ylim(0,25)+
        ylab(paste0('cell density (',arm_i,')'))+
        theme_classic()
      
      
      #reject conditions: 
      #1.H095CI<H1_miu: check1>0.001
      #2.H0_pro>H1_pro 0.6,0.4: check2>0.2
      #3.H195CI>0.1
      #4.state & miu
      check1 = abs(abs(z1_0.5)-abs(H0_95CI))
      check2 = H0_pro-H1_pro
      if(z1_0.5>0){
        check3 = abs(H1_95CI)>abs(0.1)
      }else{
        check3 = abs(H1_95CI)>abs(-0.1)
      }
      if(para_df$state =='Gain'){
        check4 = para_df$pro2<para_df$pro1 
      }else{
        check4 = para_df$pro1<para_df$pro2
      }
      
      if(check1>0.001 & check2>0.2 & check3 & check4){
        para_df$G2_dist = T
        if(z1_0.5>0){
          check = try(z1_intsct<- uniroot(function(x) {
            dnorm(x, para_df$miu1, para_df$sd1)*para_df$pro1 - dnorm(x, para_df$miu2, para_df$sd2)*para_df$pro2
          }, interval = c(0,0.2))$root)# intersect of these two distributions
          if(grepl('Error',check)){
            z1_intsct = 0.3
          }else{
            z1_intsct<- uniroot(function(x) {
              dnorm(x, para_df$miu1, para_df$sd1)*para_df$pro1 - dnorm(x, para_df$miu2, para_df$sd2)*para_df$pro2
            }, interval = c(0,0.2))$root
          }
          para_df$intsct = z1_intsct
          
          ghist_line <- ghist +
            geom_vline(xintercept = params$mean[1], linetype = "dashed", color = "blue") +
            geom_text(x=params$mean[1]+0.05,y=20,label = paste0('u1 = ',round(params$mean[1],4)), color = "blue")+
            geom_vline(xintercept = params$mean[2], linetype = "dashed", color = "red") +
            geom_text(x=params$mean[2]+0.05,y=15,label = paste0('u2 = ',round(params$mean[2],4)), color = "red")+
            geom_vline(xintercept = z1_intsct, linetype = "dashed", color = "#186F65") +
            geom_text(x=z1_intsct+0.05,y=24,label = paste0('intersect = ',round(z1_intsct,4)), color = "#186F65")+
            stat_function(fun = function(x) {
              dnorm(x, mean = params$mean[1], sd = sqrt(params$variance$sigmasq[1])) * params$pro[1]
            }, color = 'blue')+
            stat_function(fun = function(x) {
              dnorm(x, mean = params$mean[2], sd = sqrt(params$variance$sigmasq[2])) * params$pro[2]
            }, color = 'red')+
            stat_function(fun = ~dnorm(.x, mean = z1_0.5, sd = sd1)*pro1*(.x >= z1_0.5), 
                          fill='red', geom = 'polygon', alpha = 0.5)
          
        }
        
        if(z1_0.5<0){
          check = try(z1_intsct<- uniroot(function(x) {
            dnorm(x, para_df$miu1, para_df$sd1)*para_df$pro1 - dnorm(x, para_df$miu2, para_df$sd2)*para_df$pro2
          }, interval = c(-0.2,0))$root)# intersect of these two distributions
          if(grepl('Error',check)){
            z1_intsct = -0.3
          }else{
            z1_intsct<- uniroot(function(x) {
              dnorm(x, para_df$miu1, para_df$sd1)*para_df$pro1 - dnorm(x, para_df$miu2, para_df$sd2)*para_df$pro2
            }, interval = c(-0.2,0))$root
          }
          para_df$intsct = z1_intsct
          
          ghist_line <- ghist +
            geom_vline(xintercept = params$mean[1], linetype = "dashed", color = "blue") +
            geom_text(x=params$mean[1]+0.05,y=20,label = paste0('u1 = ',round(params$mean[1],4)), color = "blue")+
            geom_vline(xintercept = params$mean[2], linetype = "dashed", color = "red") +
            geom_text(x=params$mean[2]+0.05,y=15,label = paste0('u2 = ',round(params$mean[2],4)), color = "red")+
            geom_vline(xintercept = z1_intsct, linetype = "dashed", color = "#186F65") +
            geom_text(x=z1_intsct+0.05,y=24,label = paste0('intersect = ',round(z1_intsct,4)), color = "#186F65")+
            stat_function(fun = function(x) {
              dnorm(x, mean = params$mean[1], sd = sqrt(params$variance$sigmasq[1])) * params$pro[1]
            }, color = 'blue')+
            stat_function(fun = function(x) {
              dnorm(x, mean = params$mean[2], sd = sqrt(params$variance$sigmasq[2])) * params$pro[2]
            }, color = 'red')+
            stat_function(fun = ~dnorm(.x, mean = z1_0.5, sd = sd1)*pro1*(.x <= z1_0.5), 
                          fill='blue', geom = 'polygon', alpha = 0.5)
        }
      }else{
        para_df$G2_dist = F
        para_df$intsct = NA
        
        ghist_line <- ghist +
          geom_vline(xintercept = params$mean[1], linetype = "dashed", color = "blue") +
          geom_text(x=params$mean[1]+0.05,y=20,label = paste0('u1 = ',round(params$mean[1],4)), color = "blue")+
          geom_vline(xintercept = params$mean[2], linetype = "dashed", color = "red") +
          geom_text(x=params$mean[2]+0.05,y=15,label = paste0('u2 = ',round(params$mean[2],4)), color = "red")+
          geom_vline(xintercept = H0_95CI[1], color = "black") +
          geom_vline(xintercept = H0_95CI[2], color = "black") +
          geom_text(x=-0.15,y=23,label = paste0('reject'), color = "black")+
          stat_function(fun = function(x) {
            dnorm(x, mean = params$mean[1], sd = sqrt(params$variance$sigmasq[1])) * params$pro[1]
          }, color = 'blue')+
          stat_function(fun = function(x) {
            dnorm(x, mean = params$mean[2], sd = sqrt(params$variance$sigmasq[2])) * params$pro[2]
          }, color = 'red')
      }
      res = list(plot = ghist_line, para = para_df)
      return(res)
    })
    hist_list[[arm_i]] = plot_grid(plotlist = lapply(hist_list_armi,function(x) x$plot),ncol = 7,
                                   labels = names(hist_list_armi),label_size = 8)
    para_list[[arm_i]] = do.call(rbind,lapply(hist_list_armi,function(x) x$para))
  }
  write_rds(para_list, './mixgaussian_para_select_chr_all_pt.rds')
  pdf('multiPCF_logRatio_distribution_select_chr_armlevel_byALLpatient.pdf',width = 5*7,height = 4*7)
  lapply(hist_list, function(x) print(x))
  dev.off()
}


#make copykit and seurat object for selected patients who passed the distribution analysis
if(T){
  setwd(paste0(svdir,'/improve_copykat/int/'))
  
  #select patient
  para_list = read_rds('./mixgaussian_para_select_chr_all_pt.rds')
  para_list = lapply(para_list,
                     function(x){
                       x$sample.name = rownames(x)
                       return(x)
                     })
  para_df = do.call(rbind,para_list) %>% filter(G2_dist==T)
  
  
  #read copykit object
  cpkobj_int = read_rds('./cpkobj_int_allpatients.rds')
  
  
  #make seurat object
  options(future.globals.maxSize = 20000 * 1024^2)
  plan("multiprocess", workers = 100)
  
  ST_df = read.csv(paste0(svdir,'/all_ST_dir.csv'),header = T)
  ST_df_sub = ST_df 
  
  # sample.name = 'BCMHBCA52R1'
  hbca_sc_int_sub <- readRDS('HBCA_SC_int_sub_3hr_020922.rds') #subset HBCA data for cell type prediction
  for(sample.name in unique(para_df$sample.name)){
    dir.create(paste0(svdir,sample.name))
    setwd(paste0(svdir,sample.name))
    st_dir_temp = ST_df_sub %>% filter(sample == sample.name) %>% pull(outdir) %>% paste0(.,'/')
    p_emp <- ggplot() + theme_void()
    
    cpkobj = cpkobj_int[,colData(cpkobj_int)$sample.name %in% sample.name]
    write_rds(cpkobj, paste0(sample.name,'_cpk_ST.rds'))
    
    metadata = colData(cpkobj) %>% data.frame() %>% select(chr1q,chr10q)
    st_res <- st_analysis(st_dir=st_dir_temp,samplename =sample.name, clust_res=.5, nFeature_Spatial_=.05, 
                          sc_srt=hbca_sc_int_sub, pred_='celltype', metadata = metadata)
    seuobj <- st_res$srt_out
    write_rds(seuobj,paste0(sample.name,'_seu_ST.rds'))
  }
}