library(copykit)
library(SummarizedExperiment)
library(magrittr)
library(BiocParallel)
library(tidyverse)
library(glue)

register(MulticoreParam(workers = 50, progressbar = F),default = T)
brkp_cutoff <- 15
res=0.3
merged_ckobj_cbs <- readRDS("./cbs_0501_v1/32pts_merged_ckobj_cbs.RDS")
outputdir <- "./cbs_0501_v1/"
col_fun=circlize::colorRamp2(breaks = c(-1,0,1), 
                             c("blue", "white", "red"))

colData(merged_ckobj_cbs)$Sample_ID <- stringr::str_extract(string = colData(merged_ckobj_cbs)$sample, pattern = "BCMHBCA[_]?[0-9]+[RL]")
colData(merged_ckobj_cbs)$matchname <- str_split_fixed(tolower(colData(merged_ckobj_cbs)$sample), pattern = fixed("."), n = 10)[,1]
merged_ckobj_cbs <- merged_ckobj_cbs[, !colData(merged_ckobj_cbs)$Sample_ID %>% is.na()]

merged_ckobj_cbs <- runMetrics(merged_ckobj_cbs)
merged_ckobj_cbs <- findAneuploidCells(merged_ckobj_cbs, remove_XY = F, simul = F)
pdf(glue("{outputdir}/allaneuploid_ht.pdf", height = 10, width = 8))
print(plotHeatmap(merged_ckobj_cbs[,colData(merged_ckobj_cbs)$is_aneuploid==T], col = col_fun, order_cells = "hclust", n_threads = 100))
dev.off()
pdf(glue("{outputdir}/alldiploid_ht.pdf", height = 10, width = 8))
print(plotHeatmap(merged_ckobj_cbs[,colData(merged_ckobj_cbs)$is_aneuploid==F], col = col_fun, order_cells = "hclust", n_threads = 100))
dev.off()

merged_ckobj_cbs_aneu <- merged_ckobj_cbs[, colData(merged_ckobj_cbs)$is_aneuploid==T]


cons_mat <- segment_ratios(merged_ckobj_cbs_aneu[,colData(merged_ckobj_cbs_aneu)$breakpoint_count<brkp_cutoff])
cons_mat_cat <- apply(cons_mat, c(1,2), function(x){
  ifelse(x<1-res, "DEL", ifelse(x>1+res, "AMP", "NEU"))
})

cons_mat_count <- apply(cons_mat_cat, 2, function(x){
  r <- rle(x)
  a <- sapply(which(r$values!="NEU"), function(i){
    if(i==1){
      c(1, r$lengths[i], r$values[i])
    }else{
      c(sum(r$lengths[1:i-1]), sum(r$lengths[1:i]), r$values[i])
    }
  }) %>% t() %>% dplyr::as_tibble()
  colnames(a) <-c("start", "end", "event")
  a

})


mat <- cons_mat_count[[1]]
mat$cell=names(cons_mat_count)[1]
sum_cna <- mat
for (i in 2:length(cons_mat_count)) {
  mat <- cons_mat_count[[i]]
  if(nrow(mat)==0 | ncol(mat)==0){
    next
  }
  mat$cell=names(cons_mat_count)[i]
  sum_cna <- rbind(sum_cna, mat)
}

erged_ckobj_cbs_aneu_filt <- findOutliers(merged_ckobj_cbs_aneu[, colData(merged_ckobj_cbs_aneu)$sample %in% unique(sum_cna$cell)], resolution = 0.7, k = 2)
pdf(glue("{outputdir}/findcellwCNA_bp15_res03_knn07_ht.pdf"), height = 10, width = 8)
print(plotHeatmap(merged_ckobj_cbs_aneu_filt, label = "outlier", row_split = "outlier", order_cells = "hclust", n_threads = 100, col = col_fun))
dev.off()
merged_ckobj_cbs_aneu_filt <- merged_ckobj_cbs_aneu_filt[,colData(merged_ckobj_cbs_aneu_filt)$outlier==FALSE]
merged_ckobj_cbs_aneu_filt <- runUmap(merged_ckobj_cbs_aneu_filt)
merged_ckobj_cbs_aneu_filt <- findClusters(merged_ckobj_cbs_aneu_filt, k_subclones = 60)
plotHeatmap(merged_ckobj_cbs_aneu_filt, label = "subclones", row_split = "subclones",order_cells = "hclust", n_threads = 100, col = col_fun)
merged_ckobj_cbs_aneu_filt <- merged_ckobj_cbs_aneu_filt[, !colData(merged_ckobj_cbs_aneu_filt)$subclones%in%c("c18","c20")]


ref <- hg19_rg %>%
  dplyr::as_tibble() %>%
  dplyr::select(chr, start, end, arm) %>%
  dplyr::mutate(seqnames2=as.character(chr)) %>%
  dplyr::mutate(chr = if_else(seqnames2=="chrX", "chr23", if_else(seqnames2=="chrY", "chr24", seqnames2))) %>%
  dplyr::select(chr, start, end, arm) %>%
  dplyr::filter(chr!="chr24")
ref <- cbind(rownames(ref), ref)
colnames(ref)[1] <- "bin"

f <- left_join(sum_cna %>% dplyr::filter(cell%in%colData(merged_ckobj_cbs_aneu_filt[,colData(merged_ckobj_cbs_aneu_filt)$outlier==FALSE])$sample), ref[,c(1,2,3,5)], by=c("start"="bin")) %>%
  dplyr::rename(start.chr=chr, start.arm=arm) %>%
  left_join(ref[,c(1,2,4,5)],  by=c("end"="bin")) %>%
  dplyr::rename(end.chr=chr, end.arm=arm)
chrm_level_cna_cell_count <- f %>%
  dplyr::filter(cell%in%colData(merged_ckobj_cbs_aneu_filt[,colData(merged_ckobj_cbs_aneu_filt)$outlier==FALSE])$sample) %>%
  dplyr::group_by(cell) %>%
  unite(chrcomb, start.chr, end.chr,sep="+") %>%
  dplyr::summarise(n=length(unique(chrcomb)))



## filter chromosome by chromosome
chr = "chr1"
# for example, chr1
selected_cells <- sort(c(unique(f$cell[f$end.chr==chr]), unique(f$cell[f$start.chr==chr])))
colData(merged_ckobj_cbs_aneu_filt)$dup <- ifelse(colData(merged_ckobj_cbs_aneu_filt)$sample %in% selected_cells[duplicated(selected_cells)], TRUE, FALSE)
plotHeatmap(merged_ckobj_cbs_aneu_filt[, (colData(merged_ckobj_cbs_aneu_filt)$sample %in% selected_cells)] , label = c("Sample_ID","dup"), row_split = "dup", order="hclust", n_threads=100)
tmp <- runUmap(merged_ckobj_cbs_aneu_filt[, (colData(merged_ckobj_cbs_aneu_filt)$sample %in% selected_cells[duplicated(selected_cells)])], n_neighbors = 10)
tmp <- findClusters(tmp, k_subclones=5)
plotHeatmap(tmp,label = c("subclones"), row_split = "subclones", order="hclust", n_threads=100)
excluded_cells <- colData(tmp)$sample[colData(tmp)$subclones%in%c("c15","c35","c36","c37")]
del <- colData(tmp)$sample[colData(tmp)$subclones%in%c("c34")]
tmp <- runUmap(tmp[, colData(tmp)$subclones%in%c("c0")], n_neighbors = 2)
tmp <- findClusters(tmp, k_subclones=2)
plotHeatmap(tmp,label = c("Sample_ID","subclones"), row_split = "subclones", order="hclust", n_threads=100, col = col_fun )
excluded_cells <- c(colData(tmp)$sample[colData(tmp)$subclones%in%c("c2","c1","c7","c15")], excluded_cells)
del <- c(colData(tmp)$sample[colData(tmp)$subclones%in%c("c5","c6")], c(del,"BCMHBCA10R_r26c58_yes_C1_Good_S13835_R1_001.sort.markdup"))

selected_cells <- colData(merged_ckobj_cbs_aneu_filt)$sample[colData(merged_ckobj_cbs_aneu_filt)$sample %in% selected_cells[duplicated(selected_cells)]]
selected_cells <- selected_cells[!selected_cells%in%excluded_cells]

pdf(glue("{outputdir}/32pts_bp15_res03_knn07_{chr}_ht.pdf"), width = 11, height = 3)
plotHeatmap(merged_ckobj_cbs_aneu_filt[, (colData(merged_ckobj_cbs_aneu_filt)$sample %in% c(del,unique(selected_cells)))], label = c("Sample_ID"), order="hclust", col=col_fun, n_threads=100)
dev.off()

final_samples <- data.frame(chrom="chr1q", event = "AMP", cell = selected_cells)
final_samples <- rbind(final_samples,data.frame(chrom="chr1q", event = "DEL", cell = del))


save(list = c("final_samples",
              ".Random.seed",
              "chr",
              "brkp_cutoff",
              "col_fun",
              "res",
              "f",
              "merged_ckobj_cbs",
              "outputdir"),
     file = glue("{outputdir}/find_cell_wCNA_simp.Rdata"))
              