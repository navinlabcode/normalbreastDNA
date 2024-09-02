
#load packages
library(Seurat)
require(copykat)
require(magrittr)
require(Matrix)

svdir = '/volumes/USR1/yiyun/Project/HBCA/copykat_ST/'

#prepare the Spaceranger output directory
ST_df = read.csv(paste0(svdir,'all_ST_dir.csv'),header = T)
ST_df_sub = ST_df %>% dplyr::filter(if_chr1q == T)
dir_vec = ST_df_sub$outdir
smp_vec <-ST_df_sub$sample
packageVersion("copykat") #1.1.0


for (i in 1:length(dir_vec)) {
  w = 25
  k = 0.05
  wkdir =paste0(svdir,'/win_',w,'_ks_',k);dir.create(wkdir)
  setwd(wkdir)
  tryCatch(
    {print(i)
      obj_i <- Load10X_Spatial(dir_vec[i])
      count_i <- as.matrix(obj_i@assays$Spatial@counts)
      print('running CopyKat...')
      copykat::copykat(rawmat=count_i, id.type='S', ngene.chr=5, LOW.DR=0.05, UP.DR=0.1, win.size=w, norm.cell.names = "", KS.cut = k, sam.name=smp_vec[i], n.cores=50)
    }, error=function(e) {cat("ERROR :",conditionMessage(e), "\n")})
}
