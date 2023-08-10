library(copykit)
library(SummarizedExperiment)
library(magrittr)
library(BiocParallel)
library(tidyverse)
library(glue)

register(MulticoreParam(workers = 100, progressbar = F),default = T)
# brkp_cutoff <- 15
# res=0.3
load("./cbs_0501_v1/find_cell_wCNA_simp.Rdata")
outputdir <- "./cbs_0501_v1/"

filtered_obj <- cbind(merged_ckobj_cbs[,colData(merged_ckobj_cbs)$is_aneuploid==F], merged_ckobj_cbs[,colData(merged_ckobj_cbs)$sample%in%unique(final_samples$cell)])

for(id in unique(colData(filtered_obj)$Sample_ID)){
  
  pdf(glue("{outputdir}/filternoisy/{id}_filternoisy_allcells_ht.pdf", height = 10, width = 8))
  print(plotHeatmap(filtered_obj[,colData(filtered_obj)$Sample_ID==id], col = col_fun, order_cells = "hclust", n_threads = 100, label = "is_aneuploid", row_split = "is_aneuploid"))
  dev.off()
  
}

colData(merged_ckobj_cbs)$filter_tag <- 
  colData(merged_ckobj_cbs) %>% as_tibble() %>%
  mutate(filter_tag = if_else(sample%in%unique(final_samples$cell), "CNAevent", if_else(is_aneuploid, "noise","diploid"))) %>%
  dplyr::pull(filter_tag)
for(id in unique(colData(merged_ckobj_cbs)$Sample_ID)){
  
  pdf(glue("{outputdir}/filter_aneu_tag/{id}_filt_aneu_allcells_ht.pdf", height = 10, width = 8))
  print(plotHeatmap(merged_ckobj_cbs[,colData(merged_ckobj_cbs)$Sample_ID==id], col = col_fun, order_cells = "hclust", n_threads = 100, label = "filter_tag", row_split = "filter_tag"))
  dev.off()
  
}

pdf(glue("{outputdir}/32pts_bp15_res03_knn07_freqplot.pdf"), height = 3, width = 10)
plotFreq(merged_ckobj_cbs[,colData(merged_ckobj_cbs)$sample%in%unique(final_samples$cell)], high_threshold = 1.25, low_threshold = 0.75) + ylim(-0.5,0.5)
dev.off()


filtered_aneu_obj <- merged_ckobj_cbs[,colData(merged_ckobj_cbs)$sample%in%unique(final_samples$cell)]
filtered_aneu_obj_wo_X <- merged_ckobj_cbs[,colData(merged_ckobj_cbs)$sample%in%unique(final_samples$cell[!final_samples$chrom%in%c("chr23p","chr23q")])]


colData(filtered_aneu_obj)$ID_fac <- factor(x = colData(filtered_aneu_obj)$Sample_ID,
                                            levels = colData(filtered_aneu_obj) %>% as_tibble %>%
                                                      dplyr::group_by(Sample_ID) %>%
                                                      dplyr::summarise(n=n()) %>%
                                                      arrange(-n) %>%
                                                      dplyr::pull(Sample_ID))

### Plot by patient heatmap
pdf(glue("{outputdir}/32pts_bp15_res03_knn07_bypt_ht.pdf", height = 10, width = 8))
print(plotHeatmap(filtered_aneu_obj, label = "ID_fac", row_split = "ID_fac", col = col_fun, order_cells = "hclust", n_threads = 100))
dev.off()


final_samples$sample <- stringr::str_extract(string = final_samples$cell, pattern = "BCMHBCA[_]?[0-9]+[RL]")
final_samples %<>% filter(cell!="BCMHBCA66L_r30c41_yes_D1_Good_S7386") ## remove 1 mannually confirmed diploid cell
final_samples %>% dplyr::group_by(chrom, event) %>% 
  dplyr::summarise(nSample = length(unique(sample)),
            nCell = length(unique(cell))) %>%
  arrange(-nCell) %>%
  write_tsv(glue("{outputdir}/32pts_bp15_res03_knn07_event_sumstat.tsv"))

chrs <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chr23")
a <- final_samples
final_samples$chr <- factor(stringr::str_extract(string = final_samples$chrom, pattern = "chr[0-9]+"), levels = chrs)
a <- final_samples %>% arrange(chr, cell) %>%
  dplyr::select(chr,cell) %>%
  distinct(chr,cell)
num <- sapply(a$cell, function(x){which(colData(filtered_aneu_obj)$sample==x)})
cons_event_filt_final <- filtered_aneu_obj[,as.numeric(num)]
colData(cons_event_filt_final)$chrom <- a$chr

colnames(cons_event_filt_final) <- colData(cons_event_filt_final) %>% as_tibble() %>%
  unite(col = name, sample, chrom, sep = "+") %$%
  name

### plot by chromosome level heatmap
pdf(glue("{outputdir}/32pts_bp15_res03_knn07_bychr_ht.pdf"), width = 15, height = 20)
print(plotHeatmap(cons_event_filt_final, order_cells = "hclust", n_threads = 100, row_split = "chrom", label = c("ID_fac","chrom")))
dev.off()


## all events statistics
sum_cna_tbl <- colData(filtered_obj) %>% as_tibble() %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::summarise(nCellCNA = sum(is_aneuploid==T),
            nCellDiploid = sum(is_aneuploid==F)) %>%
  mutate(PropCellwCNA=nCellCNA/(nCellCNA+nCellDiploid))


## autosomes events statistics
filtered_obj_wo_X <- cbind(merged_ckobj_cbs[,colData(merged_ckobj_cbs)$is_aneuploid==F], merged_ckobj_cbs[,colData(merged_ckobj_cbs)$sample%in%unique(final_samples$cell[!final_samples$chrom%in%c("chr23p","chr23q")])]) 
Xonly_cellid <- final_samples %>% 
  dplyr::distinct(cell, event, chr, .keep_all = T) %>%
  dplyr::group_by(cell) %>%
  dplyr::filter(n() ==1) %>%
  dplyr::ungroup() %>%
  dplyr::filter(chrom %in%c("chr23p","chr23q")) %>%
  dplyr::pull(cell) %>%
  unique()
tmp <- merged_ckobj_cbs[,merged_ckobj_cbs$sample%in%Xonly_cellid]
tmp$is_aneuploid = F
filtered_obj_wo_X <- cbind(filtered_obj_wo_X, tmp)
sum_cna_tbl_wo_X <- colData(filtered_obj_wo_X) %>%
  as_tibble() %>%
  dplyr::group_by(Sample_ID) %>%
  dplyr::summarise(nCellCNA_woX = sum(is_aneuploid==T),
            nCellDiploid_woX = sum(is_aneuploid==F)) %>%
  mutate(PropCellwCNA_woX=nCellCNA_woX/(nCellCNA_woX+nCellDiploid_woX))

## cell frequency bar plots
library(viridis)
library(hrbrthemes)
ggplot(sum_cna_tbl %>% dplyr::select(-PropCellwCNA) %>% gather(key = "type",value = "count", -Sample_ID), 
       aes(fill=type, y=count, x=Sample_ID)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Pastel2") +
  geom_text(aes(label=count),color="gray30",size=3,position=position_stack(vjust=0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        legend.position = "bottom")
ggsave(filename = glue("{outputdir}/32pts_bp15_res03_knn07_count_barplot.pdf"), width = 8, height = 5)

ggplot(sum_cna_tbl_wo_X %>% dplyr::select(-PropCellwCNA_woX) %>% gather(key = "type",value = "count", -Sample_ID), 
       aes(fill=type, y=count, x=Sample_ID)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(palette = "Pastel2") +
  geom_text(aes(label=count),color="gray30",size=3,position=position_stack(vjust=0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
        legend.position = "bottom")
ggsave(filename = glue("{outputdir}/32pts_bp15_res03_knn07_count_woX_barplot.pdf"), width = 8, height = 5)

final_samples$Sample_ID <- stringr::str_extract(string = final_samples$cell, pattern = "BCMHBCA[_]?[0-9]+[RL]")

sum_cna_tbl <- final_samples %>% dplyr::group_by(Sample_ID) %>%
  dplyr::summarise(nCNAEvent = length(unique(chrom))) %>%
  full_join(sum_cna_tbl)
sum_cna_tbl <- final_samples %>% filter(!chrom%in%c("chr23p","chr23q")) %>% dplyr::group_by(Sample_ID) %>%
  dplyr::summarise(nCNAEvent_woX = length(unique(chrom))) %>%
  full_join(sum_cna_tbl_wo_X) %>%
  full_join(sum_cna_tbl)

sumstat_bychr <- final_samples %>% group_by(chrom, event) %>% 
  dplyr::summarise(nSample = length(unique(sample)),
            nCell = length(unique(cell))) %>%
  arrange(-nCell)

save(list = ls(all.names = TRUE),
     file = glue("{outputdir}/step2_plot.Rdata"))

write_tsv(final_samples, glue("{outputdir}/filtered_cells_w_event.tsv"))
