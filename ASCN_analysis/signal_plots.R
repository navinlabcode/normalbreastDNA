library(signals)
library(glue)
library(tidyverse)
library(copykit)
library(BiocParallel)
library(ggtree)
library(ape)
register(MulticoreParam(workers = 100, progressbar = F), default = T)

load("20230608_find_cell_wCNA_corrected_simp.Rdata") ## load DNA data
source("plotHeatmap_state.R")
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
chrlevels <- c("1",  "2",  "3",  "4",  "5",  "6",  "7",  "8",  "9",  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y" )
# cytoband <- read_tsv("/volumes/USR1/junke/Projects/HBCA/ASCN/scripts/cytoBand.txt.gz", col_names = F) %>% dplyr::filter(X1!="chrY")
# cytoband$chr <- factor(gsub("chr","",cytoband$X1), levels = chrlevels)
# cytoband <- cytoband %>% arrange(chr, X2, X3)
# cytoband.gr <- makeGRangesFromDataFrame(cytoband, seqnames.field = "chr", start.field = "X2", end.field = "X3", keep.extra.columns = T)
# cytoband.chr.gr <- makeGRangesFromDataFrame(cytoband, seqnames.field = "X1", start.field = "X2", end.field = "X3", keep.extra.columns = T)

cytoband <- hg19chrom_coordinates[-c(1:24,71,72),]
cytoband$seqname <- sapply(cytoband$chr, function(x){paste0("chr",x)})
cytoband.gr <- makeGRangesFromDataFrame(cytoband%>%dplyr::filter(chr!="Y"), seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = T)
cytoband.chr.gr <- makeGRangesFromDataFrame(cytoband%>%dplyr::filter(chr!="Y"), seqnames.field = "seqname", start.field = "start", end.field = "end", keep.extra.columns = T)

firsttime=F
samp <- "BCMHBCA63R"
cutoff <- "min3"

if(cutoff=="min3"){
  load(glue("{samp}/signals.Rdata"))
}else if(cutoff=="min10"){
  signals_res <- readRDS(glue("{samp}/signals/signals.Rdata"))
  hscn <- signals_res$hscn
  
  if(firsttime){
    pdf(glue("{samp}/signals/signals_hscn_state_A.pdf"), width = 10, height=10)
    print(signals::plotHeatmap(signals_res$hscn, plotcol = "A", normalize_ploidy = T, plottree = F, plotfrequency = F, show_library_label = F, maxCNcol = 4))
    dev.off()
    
    pdf(glue("{samp}/signals/signals_hscn_state_B.pdf"), width = 10, height=10)
    print(signals::plotHeatmap(signals_res$hscn, plotcol = "B", normalize_ploidy = T, plottree = F, plotfrequency = F, show_library_label = F, maxCNcol = 4))
    dev.off()
  }
}

hscn <- filtercn(hscn)
cell_keeped <- (str_split_fixed(unique(hscn$data$cell_id),"_S", n = Inf))[,1]

obj <- merged_ckobj_cbs[, (str_split_fixed(as.character(merged_ckobj_cbs$sample), pattern = "_S", n = Inf))[,1] %in% cell_keeped]
obj <- calcInteger(obj, method = "fixed", ploidy_value = 2)

state_p <- hscn$data %>% 
  as_tibble() %>% 
  mutate(state_p = case_when(state_phase == "A-Gained" ~ 3,
                             state_phase == "Balanced" ~ 2,
                             state_phase == "B-Gained" ~ 4,
                             state_phase == "B-Hom" ~ 0,
                             state_phase == "A-Hom" ~ 1,
                             state_phase == NA ~ 2)) %>%
  select(chr, start, end, cell_id, state_p) %>%
  spread(key="cell_id", value = "state_p")

state_p$chr <- factor(state_p$chr, levels = chrlevels)
state_p <- state_p %>% arrange(chr, start, end)

## remove Y
state_p <- state_p[1:(which(state_p$chr=="Y")[1]-1),]

# hc <- fastcluster::hclust(dist(t(assay(obj, "integer")), method = "euclidean"), method = "ward.D2")
hc <- fastcluster::hclust(dist(t(state_p[,-c(1:3)]), method = "manhattan"))
obj_samp <- (str_split_fixed(as.character(obj$sample), pattern = "_S", n = Inf))[,1]

if(firsttime){
  h1 <- copykit::plotHeatmap(obj[,match((str_split_fixed(hc$labels[hc$order],"_S", n = Inf))[,1], obj_samp)], assay = "integer", n_threads = 100)
  
  cl <- data.frame(cell_id=hc$labels[hc$order], clone_id = seq(1, length(hc$order), 1))
  h2 <- signals::plotHeatmap(hscn$data,
                             plottree = FALSE,
                             plotcol = "state_phase",
                             reorderclusters = T, 
                             clusters = cl,
                             show_legend = F, 
                             show_clone_label = F, 
                             show_library_label = F)
  
  pdf(glue("{samp}/tot_cn_state_phased_{cutoff}.pdf"), width = 25, height = 12)
  print(ComplexHeatmap::draw(h1 + h2,
                             ht_gap = unit(0.6, "cm"),
                             column_title = paste0("Total number of cells: ", length(unique(hscn$data$cell_id))),
                             column_title_gp = grid::gpar(fontsize = 20),
                             heatmap_legend_side = "bottom",
                             annotation_legend_side = "bottom",
                             show_heatmap_legend = TRUE))
  dev.off()
}


hc <- fastcluster::hclust(dist(t(state_p[,-c(1:3)]), method = "manhattan"))
plot(x = hc, labels =  row.names(hc), cex = 0.5)
cl_members <- cutree(tree = hc, h = 4000)
rect.hclust(tree = hc, h=4000,cluster = cl_members)
table(cl_members)
keep_samples <- names(cl_members)[cl_members%in%which(table(cl_members)>10)]

ranges_sig <- makeGRangesFromDataFrame(as_tibble(state_p[,c(1:3)]), seqnames.field = "chr")
ranges_sig <- sortSeqlevels(ranges_sig)
ranges_sig <- sort.GenomicRanges(ranges_sig)



#############################################
### make copykit object from signals call ###
#############################################

ckobj_sig <- CopyKit(
  assays = list(integer = state_p[,-c(1:3)],
                segment_ratios = state_p[,-c(1:3)]/2),
  rowRanges = ranges_sig
)
ckobj_sig <- logNorm(ckobj_sig, assay = "segment_ratios")
metadata(ckobj_sig)$genome <- "hg19"
colData(ckobj_sig)$sample <- names(state_p[,-c(1:3)])
colData(ckobj_sig)$ploidy <- 2
# plotHeatmap(ckobj_sig, assay = "integer", row_split = "outlier", label = "outlier",order_cells = "hclust")

ckobj_sig_filt <- ckobj_sig[, ckobj_sig$sample%in%keep_samples]
# ckobj_sig_filt <- findOutliers(ckobj_sig_filt, assay = "logr", k = 5, res = 0.1)


hc_filt <- fastcluster::hclust(dist(t(assay(ckobj_sig_filt,"integer")), method = "manhattan"))
# plotHeatmap(ckobj_sig_filt[, hc_filt$order], assay = "integer")

obj_filt <- obj[,match((str_split_fixed(hc_filt$labels[hc_filt$order],"_S", n = Inf))[,1], obj_samp)]
h1 <- copykit::plotHeatmap(obj[,match((str_split_fixed(hc_filt$labels[hc_filt$order],"_S", n = Inf))[,1], obj_samp)], assay = "integer", n_threads = 100)

state_chr <- hscn$data %>% 
  as_tibble() %>% 
  mutate(state_p = case_when(state_phase == "A-Gained" ~ "AAB",
                             state_phase == "Balanced" ~ "AB",
                             state_phase == "B-Gained" ~ "ABB",
                             state_phase == "B-Hom" ~ "B",
                             state_phase == "A-Hom" ~ "A",
                             state_phase == NA ~ "AB")) %>%
  select(chr, start, end, cell_id, state_p) %>%
  spread(key="cell_id", value = "state_p")

state_chr$chr <- factor(state_chr$chr, levels = chrlevels)
state_chr <- state_chr %>% arrange(chr, start, end)
state_chr <- state_chr[1:(which(state_chr$chr=="Y")[1]-1), keep_samples]

h2 <- plotHeatmap_state(df=state_chr[, hc_filt$order], chr_ranges = as.data.frame(ranges_sig), raster_quality = 5)
pdf(glue("{samp}/tot_cn_state_phased_hcfilt_4000_{cutoff}_3.pdf"), width = 15, height = 6)
print(ComplexHeatmap::draw(h1 + h2,
                           ht_gap = unit(0.6, "cm"),
                           column_title = paste0("Sample: ", samp),
                           column_title_gp = grid::gpar(fontsize = 20),
                           heatmap_legend_side = "bottom",
                           annotation_legend_side = "bottom",
                           show_heatmap_legend = TRUE))
dev.off()

## map 2 matrices to the cytoband level
source("getcytolevelmat.R")
ck_cyto_mat <- getcytolevelmat(origin.gr = rowRanges(obj_filt), cytoband.gr = cytoband.chr.gr, origin.mat = assay(obj_filt,"integer"))
sig_cyto_mat <-  getcytolevelmat(origin.gr = ranges_sig, cytoband.gr = cytoband.gr, origin.mat = assay(ckobj_sig_filt[,hc_filt$order],"integer"))

cytoband.ind <- (ck_cyto_mat$V1>10)&(sig_cyto_mat$V1>10)
cbind(
  as.data.frame(cytoband.gr), 
  ck_cyto_mat[,-1]
  )[cytoband.ind,] %>% 
  gather(key="cellname", 
         value = "tot_cn", 
         -names(as.data.frame(cytoband.gr))) %>%
  cbind(
    cbind(
      as.data.frame(cytoband.gr), 
      sig_cyto_mat[,-1]
    )[cytoband.ind,] %>%
      gather(key="cellname_sig", 
             value = "state_p", 
             -names(as.data.frame(cytoband.gr))) %>%
      dplyr::select(cellname_sig, state_p)
  ) -> ck_sig_cyto_mat_long

# ck_sig_cyto_mat_long$state_p[is.na(ck_sig_cyto_mat_long$state_p)] <- 2

## merge cell ck and signals CN state
apply(ck_sig_cyto_mat_long[, c("tot_cn", "state_p")], 1, function(x){
  if(x[1]==2){
    if(x[2]==0){
      aB=2;aA=0
    }else if(x[2]==1){
      aA=2;aB=0
    }else{
      aA=1;aB=1
    }
  }else if(x[1]>2){
    if(x[2]==3){
      aA=x[1]-1;aB=1
    }else if(x[2]==4){
      aB=x[1]-1;aA=1
    }else{
      aA=NA;aB=NA
    }
  }else if(x[1]==1){
    if(x[2]==0){
      aB=1;aA=0
    }else if(x[2]==1){
      aA=1;aB=0
    }else{
      aA=NA;aB=NA
    }
  }else if(x[1]==0){
    aA=0;aB=0
  }
  return(c(aA,aB))
}) -> ascn_cyto

ck_sig_cyto_mat_long$aA <- as_tibble(t(ascn_cyto))$V1
ck_sig_cyto_mat_long$aB <- as_tibble(t(ascn_cyto))$V2
ck_sig_cyto_mat <- ck_sig_cyto_mat_long %>% unite(ascn, aA, aB, sep = "/") %>% dplyr::select(seqnames,start,end,cellname, ascn) %>%
  spread(key=cellname, value=ascn)
# ck_sig_cyto_mat_aA <- ck_sig_cyto_mat_long  %>% dplyr::select(seqnames,start,end,cellname, aA) %>%
#   spread(key=cellname, value=aA)
# ck_sig_cyto_mat_aB <- ck_sig_cyto_mat_long  %>% dplyr::select(seqnames,start,end,cellname, aB) %>%
#   spread(key=cellname, value=aB)
df <- as.data.frame(t(ck_sig_cyto_mat[,-c(1:3)]))

## agregate cell subclonal profiles
min_cell_clone <- 1 ## 44L set to 2
agre_mat <- cbind(aggregate(cbind(df,numdup=1), df, length)$numdup,aggregate(cbind(df,numdup=1), df, length)[,c(1:sum(cytoband.ind))])
agre_mat <- agre_mat[agre_mat[,1]>min_cell_clone,]
names(agre_mat) <- c("numcell", cytoband$chrarm[cytoband.ind])

agre_mat$sample_id <- apply(agre_mat[,-1], 1, function(x){
  a <- x[x!="1/1"]
  names(a) <- names(x)[x!="1/1"]
  sapply(names(a), function(chrarm){
    y=a[chrarm]
    sub=unlist(strsplit(y,"/")); 
    names(sub) <- c("a","b")
    sapply(names(sub), function(an){
      acn <- sub[an]
      if(acn>1){
        o = paste0(chrarm,"AMP",an)
        if(acn>2){
          o = paste0(chrarm,"AMP+",an)
        }
      }else if(acn<1){
        o = paste0(chrarm,"DEL",an)
      }else{
        o=NA
      }
      o
    }) -> alist
    alist[!is.na(alist)]
  }) -> ablist
  paste0(unlist(ablist), collapse = "_")
})



agre_mat_long <- agre_mat %>% gather(key = "chrarm", value = "state", -numcell, -sample_id) %>% 
  left_join(cytoband%>%dplyr::select(chr, start, end, chrarm))

agre_mat_long$cn_a <- sapply(agre_mat_long$state, function(x){
  as.integer(unlist(strsplit(x,"/"))[1])
})
agre_mat_long$cn_b <- sapply(agre_mat_long$state, function(x){
  as.integer(unlist(strsplit(x,"/"))[2])
})

system(glue("mkdir -p {samp}/medicc"))
write_tsv(agre_mat_long,  glue("{samp}/medicc/ascn_armlevel_long.tsv"))

ascn_armlevel_long <- agre_mat_long %>% dplyr::rename(chrom=chr) %>% dplyr::select(chrom, start, end, sample_id,cn_a, cn_b)
write_tsv(ascn_armlevel_long, glue("{samp}/medicc/input.tsv"))
glue("cd {samp}/medicc; conda activate snakemake; medicc2 -j 60 -vv input.tsv output;")


tree <- read.tree(glue("{samp}/medicc/output/input_final_tree.new"))
freq_df <- agre_mat %>% dplyr::select(sample_id, numcell)
freq_df$group <- c(1,2,3,4,2)
freq_df$group <- as.character(freq_df$group)
summary_freq <- summary(freq_df$numcell)

ggtree::ggtree(tree, ladderize = FALSE, size = .3) %<+% 
  freq_df +  geom_tippoint(aes(size=numcell, color=group, fill = group)) + 
  scale_color_brewer(palette = "Set1") +
  scale_size_continuous(name = 'Cell count',
                        breaks = c(3,10,30),
                        labels = c(3,10,30))+
  theme(legend.position = "right") + 
  geom_treescale(x = 4)+
  geom_tiplab(aes(color = group), size = 3, hjust = -0.1)+ 
  geom_text(aes(x=branch, label=round(branch.length, 0), vjust=-.5), size = 3) +
  theme_tree2() -> treeplt

cowplot::save_plot(plot = treeplt, filename = glue("{samp}/medicc/output/treeplt.pdf"), base_width = 4, base_height = 4)

