library(dplyr)
#readcopykat
library(BiocParallel)
Copykit_copykat <- function(CNA_result,
         genome_version = 'hg19',
         method = 'copykat',
         resolution = '200k'){
  
  # read copykat result
  CNA_result[1:5,1:5]
  varbin_counts_df = CNA_result[,4:ncol(CNA_result)]
  varbin_counts_df[1:5,1:5]
  library(copykit, lib.loc = "/opt/R/4.1.2/lib/R/library")
  
  # make rowRanges
  if(genome_version == 'hg19'){
    gr_varbin_full <-makeGRangesFromDataFrame(hg19_rg,keep.extra.columns = T)
  }else{
    gr_varbin_full <-makeGRangesFromDataFrame(hg38_rg,keep.extra.columns = T)
  }
  GenomeInfoDb::seqlevelsStyle(gr_varbin_full) <- "Ensembl"
  gr_varbin_full <- GenomeInfoDb::renameSeqlevels(gr_varbin_full,  c(X = 23, Y = 24))
  rg <- CNA_result %>% dplyr::select(c(chrom, chrompos, abspos)) %>% as.data.frame()
  IRanges::start(gr_varbin_full) %>% length()
  length(rg$chrompos)
  key_ref <- paste0(GenomicRanges::seqnames(gr_varbin_full), 
                    "_", IRanges::start(gr_varbin_full))
  idx = grepl('24_',x=key_ref)
  g <- gr_varbin_full[!idx, ]
  g$abspos <- rg$abspos
  g <- GenomeInfoDb::renameSeqlevels(g, c(`23` = "X", `24` = "Y"))
  GenomeInfoDb::seqlevelsStyle(g) <- "UCSC"
  
  # generate seg_ratios
  ref = hg19_rg[1:nrow(varbin_counts_df),]
  ref_chrarm <- ref %>%
    dplyr::mutate(chrarm = paste0(gsub("chr", "", chr), arm))
  levels_chrarm <- gtools::mixedsort(unique(ref_chrarm$chrarm))
  ref_chrarm <- ref_chrarm %>%
    dplyr::mutate(chrarm = as.factor(chrarm)) %>%
    dplyr::mutate(chrarm = forcats::fct_relevel(chrarm, levels_chrarm))
  if (method == "CBS") {
    seg_list <-
      BiocParallel::bplapply(
        varbin_counts_df,
        FUN = function(x) {
          CNA_object <-
            DNAcopy::CNA(x,
                         ref_chrarm$chrarm,
                         ref$start,
                         data.type = "logratio",
                         sampleid = names(x)
            )
          
          withr::with_seed(seed = 17,
                           segment_smoothed_CNA_object <-
                             .quiet(
                               DNAcopy::segment(
                                 CNA_object,
                                 alpha = 1e-5,
                                 min.width = 5,
                                 undo.splits = 'prune'
                               )
                             )
          )
          short_cbs <- segment_smoothed_CNA_object[[2]]
          log_seg_mean_LOWESS <-
            rep(short_cbs$seg.mean, short_cbs$num.mark)
        },
        BPPARAM = bpparam()
    )
    seg_df <- dplyr::bind_cols(seg_list) %>%
      as.data.frame() %>%
      round(2)
  }
  if (method == 'multipcf'){
    dat=CNA_result
    mpcf <- copynumber::multipcf(dat %>% select(-abspos),
                                 arms = vapply(
                                   regmatches(ref_chrarm$chrarm,
                                              regexec("[pq]",
                                                      ref_chrarm$chrarm)),
                                   FUN =  "[",
                                   1,
                                   FUN.VALUE = character(1)
                                 ),
                                 normalize = T,
                                 return.est = F
    )
    
    seg_df <- apply(mpcf[, 6:ncol(mpcf)], 2, function(x) {
      rep.int(x, mpcf$n.probes)
    })
    # seg_df <- round(as.data.frame(seg_df), 2)
    seg_df[1:5,1:5]
  }
  if(method == 'copykat'){
    seg_df = varbin_counts_df
  }
  
  
  # read copykit
  cna_obj <- CopyKit(
    assays = list(bincounts = varbin_counts_df,
                  segment_ratios = seg_df),
    rowRanges = g
  )
  S4Vectors::metadata(cna_obj)$genome <- genome_version
  S4Vectors::metadata(cna_obj)$resolution <- resolution
  colData(cna_obj)$cell.name = colnames(seg_df)
  if (method == "multipcf") {
    mpcf_df = mpcf[, 6:ncol(mpcf)] %>% t() %>% data.frame() %>% set_colnames(paste0('chr',levels(ref_chrarm$chrarm)))
    mpcf_df$cell.name = rownames(mpcf_df)
    meta = colData(cna_obj) %>% data.frame()
    meta = left_join(meta,mpcf_df)
    colData(cna_obj) = cbind(colData(cna_obj),meta %>% select(-colnames(colData(cna_obj))))
  }
  
  SummarizedExperiment::colData(cna_obj)$sample = colnames(cna_obj)
  return(cna_obj)
}

st_analysis <- function(st_dir='',samplename, nFeature_Spatial_=.05, clust_res=0.3, sc_srt, pred_='celltype', metadata) {
  ## Test input ##
  # st_dir=st_dir_temp
  # sc_srt=hbca_sc_int_sub
  ###
  cat('Seurat loading...\n')
  srt_temp <- Load10X_Spatial(st_dir)
  img_hires <- png::readPNG(paste0(st_dir, 'spatial/tissue_hires_image.png'))
  srt_temp$id <- names(srt_temp$orig.ident)
  srt_temp$mt_percent <- PercentageFeatureSet(srt_temp, pattern = "^MT-")
  srt_temp$rp_percent <- PercentageFeatureSet(srt_temp, pattern = "^RPL|^RPS")
  srt_temp <- subset(srt_temp, subset=nFeature_Spatial>quantile(srt_temp$nFeature_Spatial, probs = nFeature_Spatial_))
  srt_temp$sample <- samplename
  srt_temp <- RenameCells(srt_temp, new.names=make.names(Cells(srt_temp)))
  
  cat('Clustering and DE analysis...\n')
  srt_temp <- NormalizeData(srt_temp) %>% ScaleData(features=rownames(.))
  srt_temp <- srt_temp %>% FindVariableFeatures %>% RunPCA %>% RunUMAP(dims=1:20)
  srt_temp <- FindNeighbors(srt_temp) %>% FindClusters(resolution=clust_res)
  srt_temp_deg <- FindAllMarkers(srt_temp, logfc.threshold=.5, min.pct=.3, only.pos=T, densify=T)
  srt_temp_deg_top <- srt_temp_deg %>% group_by(cluster) %>% top_n(10, avg_log2FC)
  
  
  cat('Label transfer...\n')
  srt_temp_anchors <- FindTransferAnchors(reference=sc_srt, query=srt_temp, dims=1:20, reduction='cca')
  srt_temp_pred <- TransferData(anchorset=srt_temp_anchors, refdata=c(sc_srt[[pred_]]), dims = 1:20, weight.reduction = 'cca') %>% 
    set_colnames(., gsub('prediction.score', 'pred.score', colnames(.))) %>% set_colnames(., gsub('predicted.id', 'pred.id', colnames(.)))
  srt_temp@meta.data[, grep('predicted|prediction.score', colnames(srt_temp@meta.data))] <- NULL
  srt_temp <- AddMetaData(srt_temp, srt_temp_pred)
  
  
  cat('Add metadata...\n')
  cat('# Cells overlapped:', length(intersect(Cells(srt_temp), rownames(metadata))), '\n')
  srt_temp <- AddMetaData(srt_temp, metadata = metadata)
  
  return(list(srt_out=srt_temp, 
              img_hires=img_hires,
              deg=srt_temp_deg))
}



multipcf_wrp <- function(mat, gamma=40) {
  cat('Run multipcf...\n')
  mat %>% t %>% data.frame(chrom=gsub('chr|\\_.*', '', rownames(.)) %>% as.numeric,
                           chrompos=gsub('.*\\_(.*)', '\\1', rownames(.)) %>% as.numeric, ., check.names = F) %>%
    copynumber::multipcf(., gamma = gamma, normalize = F, return.est = T) %>%
    extract2('segments') %>% set_rownames(., paste0('chr', .$chrom, .$arm) %>% make.unique()) %>% .[, -c(1:5)] %>% t
}

require(mclust)
score_cut <- function(mat, G=3) {
  mat_out <- data.frame(cnvscore=apply(mat, 1, function(x) sqrt(mean(x^2)))) %>% set_rownames(rownames(mat))
  mat_out$cnvscore_cut <- paste0('cnvscore_', mclust::Mclust(mat_out$cnvscore, G=G)$classification)
  mat_out$cnvscore_cut <- factor(mat_out$cnvscore_cut, levels = paste0('cnvscore_',1:length(unique(mat_out$cnvscore_cut))))
  return(mat_out)
}

grp_dbsc_otl <- function (data, group, eps, minPts) {
  data_out <- data.frame(group, dbsc_group=NA) %>% set_rownames(rownames(data))
  group <- as.character(group)
  grp_unq <- unique(group)
  for (i in 1:length(grp_unq)) {
    grp_i <- grp_unq[i]
    data_i <- data[group==grp_i, ]
    dbsc_i <- dbscan::dbscan(data_i, eps=eps, minPts=minPts)
    data_out$dbsc_group[group==grp_i] <- paste0('g',dbsc_i$cluster)
  }
  data_out
}

umap_clust <- function(mat, n_threads=10, metric=NULL, n_neighbors = 30, min_dist=0, seed=42, grp_cut=10, k=15, eps=.5, minPts=10) {
  cat('UMAP...\n')
  
  if(is.null(metric)){
    umap_res <- uwot::umap(mat,n_threads=n_threads, n_neighbors=n_neighbors, min_dist=min_dist, ret_extra='fgraph')
  }else{
    umap_res <- uwot::umap(mat, metric=metric, n_threads=n_threads, n_neighbors=n_neighbors, min_dist=min_dist, ret_extra='fgraph')
  }
  
  head(umap_res)
  umap_out <- umap_res$embedding %>% data.frame %>% set_rownames(rownames(mat))
  head(umap_out)
  
  cat('Clustering...\n')
  umap_graph <- umap_res$fgraph %>% set_rownames(rownames(mat)) %>% set_colnames(rownames(mat)) %>% 
    igraph::graph_from_adjacency_matrix(., mode = 'undirected', weighted = T)
  set.seed(seed)
  umap_graph_clust <- igraph::cluster_louvain(umap_graph)
  umap_graph_clust$membership %>% table %>% print()
  umap_out$cnv_grp <- paste0('C', umap_graph_clust$membership)
  umap_out
  
  cat('KNN reassign...\n')
  ns_g <- names(table(umap_out$cnv_grp))[table(umap_out$cnv_grp)<grp_cut]
  train_data <- umap_out %>% dplyr::filter(cnv_grp %!in% ns_g)
  test_data <- umap_out %>% dplyr::filter(cnv_grp %in% ns_g)
  if (nrow(test_data)>0 & nrow(train_data)>0) {
    knn_model <- FNN::knn(train=train_data[, c('X1', 'X2')], test=test_data[, c('X1', 'X2')], cl=train_data$cnv_grp, k=k)
    umap_out$cnv_grp[umap_out$cnv_grp %in% ns_g] <- as.character(knn_model)
  }
  umap_out$cnv_grp <- umap_out$cnv_grp %>% as.factor %>% as.numeric %>% paste0('C', .)
  
  cat('Outlier detection...\n')
  umap_out$grp_otl <- grp_dbsc_otl(data=umap_out[, c('X1', 'X2')], group=umap_out$cnv_grp, eps = eps, minPts=minPts)$dbsc_group
  
  umap_out
}


tumor_assign <- function(meta,
                         cnv_grp_=NULL,
                         cnvscore_cut_=sort(unique(meta$cnvscore_cut), decreasing = T)[1],
                         cpkt.pred_=c('aneuploid'), ncell_ = 3,
                         grp_otl_=NULL, corr_=0, cnvscore_=0.06) {

  cat('copykat predicted tumor cell number is ', meta %>% filter(cpkt.pred == cpkt.pred_) %>% nrow, '\n')
  if (is.null(cnv_grp_)) {
    if(is.null(grp_otl_)){
      tum_grp = meta %>%
        dplyr::group_by(cnv_grp) %>%
        dplyr::summarise(corr = mean(corr),cnvscore = mean(cnvscore)) %>%
        dplyr::filter(corr<corr_, cnvscore>cnvscore_)%>% extract2('cnv_grp')
      cat('1st filter: tumor group narrowed down to:', tum_grp, '\n')
      cat('2nd filter: cnvscore_cut set to:', cnvscore_cut_, '\n',
          'group aneuploid cell number set to be more than 3', '\n')
      meta$tumor <- 'TN'
      tum_idx <- which(meta$cnvscore_cut %in% cnvscore_cut_ &
                         meta$cnv_grp %in% tum_grp &
                         meta$copykat.pred %in% cpkt.pred_)
      meta$tumor[tum_idx] <- paste0('T', gsub('C', '', meta$cnv_grp[tum_idx]))
      tum_grp_up = meta %>%
        dplyr::group_by(cnv_grp,tumor) %>% dplyr::count() %>%
        dplyr::filter(tumor!='TN',n>=ncell_)%>% extract2('cnv_grp')
      cat('update: tumor group narrowed down to:', tum_grp_up, '\n')
      meta$tumor <- 'TN'
      tum_idx_up <- which(meta$cnvscore_cut %in% cnvscore_cut_ &
                         meta$cnv_grp %in% tum_grp_up &
                         meta$copykat.pred %in% cpkt.pred_)
      meta$tumor[tum_idx_up] <- paste0('T', gsub('C', '', meta$cnv_grp[tum_idx_up]))
      cat('now filtered predicted tumor cell number is:', length (tum_idx_up), '\n')
    }else{
      tum_grp = meta %>% filter(grp_otl %in% grp_otl_) %>%
        dplyr::group_by(cnv_grp) %>%
        dplyr::summarise(corr = mean(corr),cnvscore = mean(cnvscore)) %>%
        dplyr::filter(corr<corr_, cnvscore>cnvscore_)%>% extract2('cnv_grp')
      cat('1st filter: tumor group narrowed down to:', tum_grp, '\n')

      cat('2nd filter: cnvscore_cut set to:', cnvscore_cut_, '\n',
          'group aneuploid cell number set to be more than 3', '\n')
      meta$tumor <- 'TN'
      tum_idx <- which(meta$grp_otl %in% grp_otl_ &
                         meta$cnvscore_cut %in% cnvscore_cut_ &
                         meta$cnv_grp %in% tum_grp &
                         meta$copykat.pred %in% cpkt.pred_)
      meta$tumor[tum_idx] <- paste0('T', gsub('C', '', meta$cnv_grp[tum_idx]))
      tum_grp_up = meta %>% filter(grp_otl %in% grp_otl_) %>%
        dplyr::group_by(cnv_grp,tumor) %>% dplyr::count() %>%
        dplyr::filter(tumor!='TN',n>=ncell_)%>% extract2('cnv_grp')
      cat('update: tumor group narrowed down to:', tum_grp_up, '\n')
      meta$tumor <- 'TN'
      tum_idx_up <- which(meta$cnvscore_cut %in% cnvscore_cut_ &
                            meta$cnv_grp %in% tum_grp_up &
                            meta$copykat.pred %in% cpkt.pred_)
      meta$tumor[tum_idx_up] <- paste0('T', gsub('C', '', meta$cnv_grp[tum_idx_up]))
      cat('now filtered predicted tumor cell number is:', length (tum_idx_up), '\n')
    }
  } else {
    tum_grp <- cnv_grp_
    cat('Find tumor group:', tum_grp, '\n')
    meta$tumor <- 'TN'
    tum_idx <- which(meta$grp_otl %in% grp_otl_ & meta$cnvscore_cut %in% cnvscore_cut_ & meta$cnv_grp %in% c(tum_grp))
    meta$tumor[tum_idx] <- paste0('T', gsub('C', '', meta$cnv_grp[tum_idx]))
  }
  return(meta)
}

tumor_assign <- function(umap_df, cnv_grp_=NULL, cnvscore_cut_=sort(unique(umap_df$cnvscore_cut), decreasing = T)[1], cpkt.pred_=c('anpd'), grp_otl_=c(1), ncell_=10, prop_=70) {
  if (is.null(cnv_grp_)) {
    tum_grp <- umap_df %>% dplyr::group_by(cnv_grp, grp_otl, cnvscore_cut, cpkt.pred) %>% dplyr::summarize(ncell=n()) %>% group_by(cnv_grp) %>% mutate(prop=ncell*100/sum(ncell)) %>% ungroup() %>% 
      dplyr::filter(grp_otl %in% grp_otl_, cnvscore_cut %in% cnvscore_cut_, cpkt.pred %in% cpkt.pred_, ncell>ncell_, prop>prop_) %>% extract2('cnv_grp')
    cat('Find tumor group:', tum_grp, '\n')
    umap_df$tumor <- 'TN'
    tum_idx <- which(umap_df$grp_otl %in% grp_otl_ & umap_df$cnvscore_cut %in% cnvscore_cut_ & umap_df$cnv_grp %in% tum_grp & umap_df$cpkt.pred %in% cpkt.pred_)
    umap_df$tumor[tum_idx] <- paste0('T', gsub('C', '', umap_df$cnv_grp[tum_idx]))
  } else {
    tum_grp <- cnv_grp_
    cat('Find tumor group:', tum_grp, '\n')
    umap_df$tumor <- 'TN'
    tum_idx <- which(umap_df$grp_otl %in% grp_otl_ & umap_df$cnvscore_cut %in% cnvscore_cut_ & umap_df$cnv_grp %in% c(tum_grp))
    umap_df$tumor[tum_idx] <- paste0('T', gsub('C', '', umap_df$cnv_grp[tum_idx]))
  }
  umap_df
}


to_list <- function (in_vec) {
  res <- list()
  grps <- sort(unique(in_vec))
  for (i in 1:length(grps)) {
    grp_i <- grps[i]
    res[[i]] <- names(which(in_vec==grp_i))
  }
  # names(res) <- paste0('clust_', 1:length(grps))
  names(res) <- make.names(grps)
  res
}

library(RColorBrewer)
colormaker <- function(x, pattern = "Set2"){
  n_color = length(unique(x)); 
  if(n_color>8){
    minor_col<-setNames(colorRampPalette(brewer.pal(8, pattern))(n_color),sort(unique(x),decreasing = F)) 
  }else{minor_col<-setNames(brewer.pal(n_color, pattern),sort(unique(x),decreasing = F))}
  return(minor_col[1:n_color])
}

library(stringr)
plotHeatmap_copykat = function (scCNA, assay = "segment_ratios", pt.name = '',
                                order_cells = c("consensus_tree", "hclust", "phylogeny"), 
                                label = NULL, label_colors = NULL, lowscale = -0.5,highscale = 0.5,
                                consensus = FALSE, rounding_error = FALSE, row_split = NULL, n_threads = 48) 
{
  order_cells <- match.arg(order_cells)
  if (is.null(label) & !is.null(label_colors)) {
    stop("Please provide a label argument if colors are being specified for it.")
  }
  if (!is.null(label_colors) & !is.list(label_colors)) {
    stop("label_colors argument must be a named list.")
  }
  if (!is.null(label_colors)) {
    if (length(label_colors) != length(label)) {
      stop("Label and Label colors arguments must have the same length.")
    }
  }
  if (!is.null(label) & !is.character(label)) {
    stop("Label must be a character vector.")
  }
  if (is.null(SummarizedExperiment::colData(scCNA)$subclones) && 
      order_cells != "phylogeny") {
    message("Ordering by consensus requires cluster information.\nSwitching to hclust")
    order_cells <- "hclust"
  }
  if (rounding_error == TRUE && assay != "integer") {
    stop("Rounding error argument must be used with assay 'integer'.")
  }
  if (consensus == TRUE) {
    if (attr(consensus(scCNA), "consensus_assay") == "integer") {
      assay <- "integer"
    }
  }
  
  color_heat = circlize::colorRamp2(breaks = c(lowscale,0,highscale), c("dodgerblue3", "white", "firebrick3"))
  if (assay == "integer") {
    mean_ploidy <- mean(SummarizedExperiment::colData(scCNA)$ploidy)
    ploidy_trunc <- 2 * round(mean_ploidy)
    color_heat <- structure(pals::ocean.balance(length(0:ploidy_trunc)), 
                            names = 0:ploidy_trunc)
    if (round(mean_ploidy) == 2) {
      color_heat <- structure(c("#3787BA", "#95B8C5", "#F0ECEB", 
                                "#D7A290", "#BF583B", "#8D1128", "#3C0912"), 
                              names = c("0", "1", "2", "3", "4", "5", "6"))
    }
  }
  
  color_heat = circlize::colorRamp2(breaks = c(lowscale,0,highscale), c("dodgerblue3", "white", "firebrick3"))
  
  seg_data <- t(SummarizedExperiment::assay(scCNA, assay))
  seg_data[1:5,1:5]
  chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
  if (any(chr_ranges$seqnames == "24") || any(chr_ranges$seqnames == 
                                              "Y") || any(chr_ranges$seqnames == "chrY")) {
    chr_binary <- rep(c(2, 1), length(chr_lengths)/2)
  }
  else {
    chr_binary <- c(rep(c(2, 1), (length(chr_lengths)/2)), 
                    2)
  }
  chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  chr_rl_c <- c(1, cumsum(chr_lengths))
  chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], 
                       b = chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))
  chrom.names <- c(1:22, "X", "Y")
  v <- vector(length = sum(chr_lengths), mode = "character")
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""
  chr_bar <- ComplexHeatmap::HeatmapAnnotation(chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)], 
                                                                                    gp = grid::gpar(fontsize = 14)), df = as.character(chr[1:nrow(chr), 
                                                                                    ]), show_legend = FALSE, show_annotation_name = FALSE, 
                                               which = "column", col = list(df = c(`1` = "grey88", `2` = "black")))
  if (consensus == FALSE) {
    if (order_cells == "phylogeny") {
      if (ape::Ntip(phylo(scCNA)) == 0) {
        stop("No phylogeny detected in scCNA object. Use runPhylo")
      }
      tree <- phylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      seg_data_ordered <- seg_data[tree_tips_order, ]
    }
    if (order_cells == "hclust") {
      if (length(copykit::distMat(scCNA)) == 0) {
        message("No distance matrix detected in the scCNA object.")
        scCNA <- runDistMat(scCNA, metric = "euclidean",n_threads = n_threads)
      }
      if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
        stop("Number of samples in the distance matrix different\n          from number of samples in the scCNA object.\n        Perhaps you filtered your dataset? use copykit::runDistMat() to update it.")
      }
      hc <- fastcluster::hclust(distMat(scCNA), method = "ward.D2")
      seg_data_ordered <- seg_data[hc$order, ]
    }
    if (order_cells == "consensus_tree") {
      if (nrow(consensus(scCNA)) == 0) {
        scCNA <- calcConsensus(scCNA)
      }
      scCNA <- runConsensusPhylo(scCNA)
      consensus_by <- attr(consensus(scCNA), "consensus_by")
      meta <- as.data.frame(colData(scCNA)) %>% dplyr::select(sample, 
                                                              !!consensus_by)
      meta_info <- as.character(dplyr::pull(meta, !!consensus_by))
      tree <- consensusPhylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      meta_o <- meta[order(match(meta_info, tree_tips_order)), 
      ]
      seg_data_ordered <- seg_data[meta_o$sample, ]
    }else{
      seg_data_ordered <- seg_data
    }
  }
  else {
    scCNA <- runConsensusPhylo(scCNA)
    tree <- consensusPhylo(scCNA)
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
      rev()
    seg_data <- as.matrix(t(consensus(scCNA)))
    seg_data_ordered <- seg_data[tree_tips_order, ]
  }
  if (assay == "integer") {
    seg_data_int <- seg_data_ordered
    seg_data_int[seg_data_int > ploidy_trunc] <- ploidy_trunc
    if (rounding_error == TRUE) {
      int_nr <- as.matrix(SummarizedExperiment::assay(scCNA, 
                                                      "segment_ratios")) %*% diag(SummarizedExperiment::colData(scCNA)$ploidy)
      names(int_nr) <- names(SummarizedExperiment::assay(scCNA, 
                                                         "segment_ratios"))
      err <- abs(int_nr - assay(scCNA, "integer"))
      err <- err[rownames(seg_data_int)]
    }
  }
  message("Plotting Heatmap.")
  complex_args <- list(use_raster = TRUE, column_title = "genomic coordinates", 
                       column_title_gp = grid::gpar(fontsize = 18), column_title_side = "bottom", 
                       row_title = paste0(pt.name, '             ',nrow(seg_data_ordered), " samples"), 
                       row_title_gp = grid::gpar(fontsize = 18), top_annotation = chr_bar, 
                       cluster_rows = FALSE, border = TRUE, cluster_columns = FALSE, 
                       show_column_names = FALSE, show_row_names = FALSE, show_heatmap_legend = TRUE)
  if (is.null(label)) {
    if (assay != "integer") {
      suppressMessages(do.call(ComplexHeatmap::Heatmap, 
                               c(list(matrix = seg_data_ordered, 
                                      heatmap_legend_param = list(title = "log2 (segratio)"),
                                      col = color_heat), 
                                 complex_args)))
    }
    else {
      if (rounding_error == FALSE) {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                heatmap_legend_param = list(title = "copy number"), 
                                                col = color_heat), complex_args))
      }
      else {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                heatmap_legend_param = list(title = "rounding error"), 
                                                col = viridis::viridis(200)), complex_args))
      }
    }
  }
  else {
    metadata <- SummarizedExperiment::colData(scCNA) %>% 
      as.data.frame()
    if (consensus == FALSE) {
      metadata <- metadata[rownames(seg_data_ordered), 
      ]
    }
    metadata_anno_df <- metadata %>% dplyr::select(dplyr::all_of(label))
    if (consensus == TRUE) {
      cons_attr <- attr(consensus(scCNA), "consensus_by")
      if (length(label) > 1) {
        stop("Label must be of length 1 for consensus heatmap annotation.")
      }
      if (cons_attr != label) {
        stop("Consensus heatmap can only be annotated with the same metadata element\n             used for generating the consensus matrix.")
      }
      metadata_anno_df <- metadata_anno_df[label] %>% dplyr::distinct()
      rownames(metadata_anno_df) <- metadata_anno_df %>% 
        dplyr::pull(!!cons_attr)
      metadata_anno_df <- metadata_anno_df[rownames(seg_data_ordered), 
                                           , drop = FALSE]
    }
    if (is.null(label_colors)) {
      h <- 15
      l <- 65
      cont_options <- c("D", "A", "E", "C", "B")
      cont_i <- 1
      label_colors <- list()
      label_colors <- c(list(superclones = superclones_pal(), 
                             subclones = subclones_pal(), filtered = c(removed = "#DA614D", 
                                                                       kept = "#5F917A"), is_normal = c(`TRUE` = "#396DB3", 
                                                                                                        `FALSE` = "#11181D")), label_colors)
      default_labels <- c("superclones", "subclones", "filtered", 
                          "is_normal")
      for (i in 1:length(label)) {
        if (any(str_detect(label[i], default_labels))) {
          label_colors[i] <- label_colors[default_labels[stringr::str_detect(label[i], 
                                                                             default_labels)]]
          names(label_colors)[i] <- label[i]
        }
        else if (is.numeric(dplyr::pull(metadata_anno_df, 
                                        label[i]))) {
          n = 300
          min_v = min(dplyr::pull(metadata_anno_df, label[i]))
          max_v = max(dplyr::pull(metadata_anno_df, label[i]))
          label_colors[i] <- list(circlize::colorRamp2(seq(min_v, 
                                                           max_v, length = n), viridis::viridis(n, option = cont_options[cont_i])))
          names(label_colors)[i] <- label[i]
          cont_i <- cont_i + 1
        }
        else {
          elements <- metadata_anno_df %>% dplyr::pull(label[i]) %>% 
            unique() %>% as.character() %>% sort()
          n <- length(elements)
          hex <- (scales::hue_pal(h = c(0, 360) + h, 
                                  l = 65))(n)
          col <- structure(hex, names = elements)
          label_colors[i] <- list(col)
          names(label_colors)[i] <- label[i]
          l <- l - 10
          h <- h + 15
        }
      }
    }
    label_colors[sapply(label_colors, is.null)] <- NULL
    cluster_anno <- ComplexHeatmap::rowAnnotation(df = metadata_anno_df, 
                                                  col = label_colors, show_annotation_name = T)
    if (!is.null(row_split)) {
      if (length(row_split) > 1) {
        stop("row_split length must be 1")
      }
      else {
        if (assay != "integer") {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_ordered, row_split = dplyr::pull(metadata_anno_df, 
                                                                                                     row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "log2 (segratio)"),col = color_heat), 
                                             complex_args))
        }
        else {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  row_split = dplyr::pull(metadata_anno_df, 
                                                                          row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
      }
    }
    else {
      if (assay != "integer") {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_ordered, left_annotation = cluster_anno, heatmap_legend_param = list(title = "log2 (segratio)"),
                                                col = color_heat), 
                                           complex_args))
      }
      else {
        if (rounding_error == FALSE) {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
        else {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "rounding error"), 
                                                  col = viridis::viridis(200)), complex_args))
        }
      }
    }
  }
}


plotHeatmap_copykatST = function (scCNA, assay = "segment_ratios", pt.name = '',
                                order_cells = c("no_order","consensus_tree", 
                                                           "hclust", "phylogeny"), label = NULL, label_colors = NULL, 
          group = NULL, consensus = FALSE, rounding_error = FALSE, 
          genes = NULL, col = NULL, row_split = NULL, use_raster = TRUE, 
          raster_quality = 5, n_threads = 1) 
{
  order_cells <- match.arg(order_cells)
  group_value <- NULL
  if (is.null(label) & !is.null(label_colors)) {
    stop("Please provide a label argument.")
  }
  if (!is.null(label_colors) & !is.list(label_colors)) {
    stop("label_colors argument must be a named list.")
  }
  if (!is.null(label_colors)) {
    if (length(label_colors) != length(label)) {
      stop("Label and Label colors arguments must have the same length.")
    }
  }
  if (!is.null(label) & !is.character(label)) {
    stop("Label must be a character vector.")
  }
  
  if(order_cells == 'no_order'){
    order_cells <- 'no_order'
  }else{
    if (is.null(SummarizedExperiment::colData(scCNA)$subclones) && 
        order_cells != "phylogeny") {
      message("Ordering by consensus requires cluster information.")
      message("Switching to hclust.")
      order_cells <- "hclust"
    }
  }
  
  if (rounding_error == TRUE && assay != "integer") {
    stop("Rounding error argument must be used with assay 'integer'.")
  }
  if (consensus == TRUE) {
    consensus_assay <- attr(consensus(scCNA), "consensus_assay")
    if (assay == "integer" & consensus_assay != "integer") {
      stop("Consensus must be calculated from the integer assay.")
    }
    if (consensus_assay == "integer") {
      assay <- "integer"
    }
  }
  if (assay == "integer") {
    mean_ploidy <- mean(SummarizedExperiment::colData(scCNA)$ploidy)
    ploidy_trunc <- 2 * round(mean_ploidy)
    color_heat <- structure(ocean.balance(length(0:ploidy_trunc)), 
                            names = 0:ploidy_trunc)
    if (round(mean_ploidy) == 2) {
      color_heat <- structure(c("#3787BA", "#95B8C5", "#F0ECEB", 
                                "#D7A290", "#BF583B", "#8D1128", "#3C0912"), 
                              names = c("0", "1", "2", "3", "4", "5", "6"))
    }
  }
  seg_data <- t(SummarizedExperiment::assay(scCNA, assay))
  if (any(duplicated(row.names(seg_data)))) {
    row.names(seg_data) <- make.names(row.names(seg_data), 
                                      unique = TRUE)
  }
  chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
  if (any(chr_ranges$seqnames == "24") || any(chr_ranges$seqnames == 
                                              "Y") || any(chr_ranges$seqnames == "chrY")) {
    chr_binary <- rep(c(2, 1), length(chr_lengths)/2)
  }
  else {
    chr_binary <- c(rep(c(2, 1), (length(chr_lengths)/2)), 
                    2)
  }
  chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  chr_rl_c <- c(1, cumsum(chr_lengths))
  chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], 
                       b = chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))
  chrom.names <- c(1:22, "X", "Y")
  v <- vector(length = sum(chr_lengths), mode = "character")
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""
  chr_bar <- ComplexHeatmap::HeatmapAnnotation(chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)], 
                                                                                    gp = grid::gpar(fontsize = 14)), df = as.character(chr[1:nrow(chr), 
                                                                                    ]), show_legend = FALSE, show_annotation_name = FALSE, 
                                               which = "column", col = list(df = c(`1` = "grey88", `2` = "black")))
  if (consensus == FALSE) {
    if (order_cells == "phylogeny") {
      if (ape::Ntip(phylo(scCNA)) == 0) {
        stop("No phylogeny detected in scCNA object. Use runPhylo")
      }
      if (ape::Ntip(consensusPhylo(scCNA)) == 0) {
        stop("No consensus phylogeny detected in scCNA object.")
      }
      tree <- phylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      seg_data_ordered <- seg_data[tree_tips_order, ]
    }
    if (order_cells == "hclust") {
      if (length(copykit::distMat(scCNA)) == 0) {
        message("No distance matrix detected in the scCNA object.")
        scCNA <- runDistMat(scCNA, metric = "euclidean", 
                            n_threads = n_threads)
      }
      if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
        stop("Number of samples in the distance matrix different from number\n                 of samples in the scCNA object. Perhaps you filtered your\n                 dataset?.")
      }
      hc <- fastcluster::hclust(distMat(scCNA), method = "ward.D2")
      seg_data_ordered <- seg_data[hc$order, ]
    }
    if (order_cells == "consensus_tree") {
      if (nrow(consensus(scCNA)) == 0) {
        scCNA <- calcConsensus(scCNA)
      }
      consensus_by <- attr(consensus(scCNA), "consensus_by")
      meta <- as.data.frame(colData(scCNA)) %>% dplyr::select(sample, 
                                                              !!consensus_by)
      meta_info <- as.character(dplyr::pull(meta, !!consensus_by))
      if (ape::Ntip(consensusPhylo(scCNA)) == 0) {
        stop("No consensus phylogeny in the CopyKit object.")
      }
      else {
        tree <- consensusPhylo(scCNA)
      }
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      meta_o <- meta[order(match(meta_info, tree_tips_order)), 
      ]
      seg_data_ordered <- seg_data[row.names(meta_o), ]
    }
    if(order_cells == 'no_order'){
      seg_data_ordered <- seg_data
    }
  }
  else {
    if (ape::Ntip(consensusPhylo(scCNA)) == 0) {
      stop("Build a consensus with calcConsensus and use runConsensusPhylo\n        to store the order of the consensus matrix.")
    }
    else {
      tree <- consensusPhylo(scCNA)
    }
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
      rev()
    seg_data <- as.matrix(t(consensus(scCNA)))
    seg_data_ordered <- seg_data[tree_tips_order, ]
  }
  if (assay == "integer") {
    seg_data_int <- seg_data_ordered
    seg_data_int[seg_data_int > ploidy_trunc] <- ploidy_trunc
    if (rounding_error == TRUE) {
      int_nr <- as.matrix(SummarizedExperiment::assay(scCNA, 
                                                      "segment_ratios")) %*% diag(SummarizedExperiment::colData(scCNA)$ploidy)
      names(int_nr) <- names(SummarizedExperiment::assay(scCNA, 
                                                         "segment_ratios"))
      err <- abs(int_nr - assay(scCNA, "integer"))
      err <- err[rownames(seg_data_int)]
    }
  }
  if (!is.null(genes)) {
    df <- find_scaffold_genes(scCNA, genes)
    mk <- ComplexHeatmap::columnAnnotation(foo = anno_mark(at = df$pos, 
                                                           labels = df$gene, side = "bottom", labels_gp = grid::gpar(fontsize = 12)), 
                                           show_annotation_name = FALSE, show_legend = FALSE)
  }
  else {
    mk <- NULL
  }
  if (consensus == TRUE) {
    metadata <- as.data.frame(colData(scCNA))
    cons_attr <- attr(consensus(scCNA), "consensus_by")
    if (!is.null(group)) {
      if (!is.null(group) && !(group %in% colnames(metadata))) {
        stop("Group ", group, " is not a column of colData().")
      }
      metadata_gr <- metadata %>% droplevels() %>% dplyr::select(!!cons_attr, 
                                                                 !!group) %>% tidyr::gather(key = "group_key", 
                                                                                            value = "group_value", -!!cons_attr)
      names(metadata_gr)[1] <- "cons_attr"
      metadata_counts <- metadata_gr %>% dplyr::group_by(cons_attr) %>% 
        dplyr::count(group_value) %>% dplyr::mutate(n = n/sum(n)) %>% 
        tidyr::pivot_wider(names_from = group_value, 
                           values_from = n, id_cols = cons_attr, values_fill = 0) %>% 
        as.data.frame()
      rownames(metadata_counts) <- metadata_counts[, 1]
      metadata_counts <- metadata_counts[, -1]
      elements_groups <- sort(unique(as.character(metadata_gr$group_value)))
      metadata_counts <- metadata_counts[tree_tips_order, 
                                         elements_groups]
      n_groups <- length(elements_groups)
      hex <- (scales::hue_pal())(n_groups)
      col_group <- structure(hex, names = elements_groups)
      ha_barplot <- rowAnnotation(foo = anno_barplot(metadata_counts, 
                                                     gp = grid::gpar(fill = col_group)), show_annotation_name = FALSE)
    }
    else {
      ha_barplot <- NULL
    }
  }
  else {
    ha_barplot <- NULL
  }
  message("Plotting Heatmap.")
  complex_args <- list(use_raster = use_raster, raster_quality = raster_quality, 
                       col = col, bottom_annotation = mk, right_annotation = ha_barplot, 
                       column_title = "genomic coordinates", column_title_gp = grid::gpar(fontsize = 18), 
                       column_title_side = "bottom", 
                       row_title = paste0(pt.name, '             ',nrow(seg_data_ordered), " spots"), 
                       row_title_gp = grid::gpar(fontsize = 18), 
                       top_annotation = chr_bar, cluster_rows = FALSE, border = TRUE, 
                       cluster_columns = FALSE, show_column_names = FALSE, show_row_names = FALSE, 
                       show_heatmap_legend = TRUE)
  if (ncol(seg_data_ordered) > 16000) {
    complex_args$raster_by_magick <- FALSE
  }
  if (is.null(label)) {
    if (assay != "integer") {
      suppressMessages(do.call(ComplexHeatmap::Heatmap, 
                               c(list(matrix = log2(seg_data_ordered + 0.001), 
                                      heatmap_legend_param = list(title = "log2 (segratio)")), 
                                 complex_args)))
    }
    else {
      complex_args <- complex_args[which(names(complex_args) != 
                                           "col")]
      if (rounding_error == FALSE) {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                heatmap_legend_param = list(title = "copy number"), 
                                                col = color_heat), complex_args))
      }
      else {
        complex_args <- complex_args[which(names(complex_args) != 
                                             "col")]
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                heatmap_legend_param = list(title = "rounding error"), 
                                                col = viridis::viridis(200)), complex_args))
      }
    }
  }
  else {
    metadata <- SummarizedExperiment::colData(scCNA) %>% 
      as.data.frame()
    if (consensus == FALSE) {
      metadata <- metadata[rownames(seg_data_ordered), 
      ]
    }
    metadata_anno_df <- metadata %>% dplyr::select(dplyr::all_of(label))
    if (consensus == TRUE) {
      cons_attr <- attr(consensus(scCNA), "consensus_by")
      if (length(label) > 1) {
        stop("Label must be of length 1 for consensus heatmap.")
      }
      if (cons_attr != label) {
        stop("Consensus heatmap can only be annotated with the same\n                     metadata element used for generating the consensus matrix.")
      }
      metadata_anno_df <- metadata_anno_df[label] %>% dplyr::distinct()
      rownames(metadata_anno_df) <- metadata_anno_df %>% 
        dplyr::pull(!!cons_attr)
      metadata_anno_df <- metadata_anno_df[rownames(seg_data_ordered), 
                                           , drop = FALSE]
    }
    if (is.null(label_colors)) {
      h <- 15
      l <- 65
      cont_options <- c("D", "A", "E", "C", "B")
      cont_i <- 1
      label_colors <- list()
      label_colors <- c(list(superclones = superclones_pal(), 
                             subclones = subclones_pal(), outlier = c(`TRUE` = "#DA614D", 
                                                                      `FALSE` = "#5F917A"), is_aneuploid = c(`TRUE` = "#396DB3", 
                                                                                                             `FALSE` = "#11181D")), label_colors)
      default_labels <- c("superclones", "subclones", "outlier", 
                          "is_aneuploid")
      for (i in 1:length(label)) {
        if (any(grepl(paste(default_labels, collapse = "|"), 
                      label[i]))) {
          label_colors[i] <- label_colors[default_labels[vapply(default_labels, 
                                                                function(x) grepl(x, label[i]), logical(1))]]
          names(label_colors)[i] <- label[i]
        }
        else if (is.numeric(dplyr::pull(metadata_anno_df, 
                                        label[i]))) {
          n <- 300
          min_v <- min(dplyr::pull(metadata_anno_df, 
                                   label[i]))
          max_v <- max(dplyr::pull(metadata_anno_df, 
                                   label[i]))
          label_colors[i] <- list(circlize::colorRamp2(seq(min_v, 
                                                           max_v, length = n), viridis::viridis(n, option = cont_options[cont_i])))
          names(label_colors)[i] <- label[i]
          cont_i <- cont_i + 1
        }
        else {
          elements <- metadata_anno_df %>% dplyr::pull(label[i]) %>% 
            unique() %>% as.character() %>% sort()
          n <- length(elements)
          hex <- (scales::hue_pal(h = c(0, 360) + h, 
                                  l = 65))(n)
          col_lab <- structure(hex, names = elements)
          label_colors[i] <- list(col_lab)
          names(label_colors)[i] <- label[i]
          l <- l - 10
          h <- h + 15
        }
      }
    }
    label_colors[sapply(label_colors, is.null)] <- NULL
    cluster_anno <- ComplexHeatmap::rowAnnotation(df = metadata_anno_df, 
                                                  col = label_colors, show_annotation_name = FALSE)
    if (!is.null(row_split)) {
      if (length(row_split) > 1) {
        stop("row_split length must be 1")
      }
      else {
        if (assay != "integer") {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = log2(seg_data_ordered + 
                                                                  0.001), row_split = dplyr::pull(metadata_anno_df, 
                                                                                                  row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "log2 (segratio)")), 
                                             complex_args))
        }
        else {
          complex_args <- complex_args[which(names(complex_args) != 
                                               "col")]
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  row_split = dplyr::pull(metadata_anno_df, 
                                                                          row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
      }
    }
    else {
      if (assay != "integer") {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = log2(seg_data_ordered + 
                                                                0.001), left_annotation = cluster_anno, heatmap_legend_param = list(title = "log2 (segratio)")), 
                                           complex_args))
      }
      else {
        if (rounding_error == FALSE) {
          complex_args <- complex_args[which(names(complex_args) != 
                                               "col")]
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
        else {
          complex_args <- complex_args[which(names(complex_args) != 
                                               "col")]
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "rounding error"), 
                                                  col = viridis::viridis(200)), complex_args))
        }
      }
    }
  }
}


plotHeatmap2 <- function (scCNA, assay = "segment_ratios", n_threads = 40, 
                          order_cells = c("consensus_tree", "hclust", "phylogeny", 'none'), label = NULL, label_colors = NULL, 
                          consensus = FALSE, rounding_error = FALSE, row_split = NULL) 
{
  order_cells <- match.arg(order_cells)
  if (is.null(label) & !is.null(label_colors)) {
    stop("Please provide a label argument if colors are being specified for it.")
  }
  if (!is.null(label_colors) & !is.list(label_colors)) {
    stop("label_colors argument must be a named list.")
  }
  if (!is.null(label_colors)) {
    if (length(label_colors) != length(label)) {
      stop("Label and Label colors arguments must have the same length.")
    }
  }
  if (!is.null(label) & !is.character(label)) {
    stop("Label must be a character vector.")
  }
  if (is.null(SummarizedExperiment::colData(scCNA)$subclones) && 
      order_cells != "phylogeny") {
    message("Ordering by consensus requires cluster information.\nSwitching to hclust")
    order_cells <- "hclust"
  }
  if (rounding_error == TRUE && assay != "integer") {
    stop("Rounding error argument must be used with assay 'integer'.")
  }
  if (consensus == TRUE) {
    if (attr(copykit::consensus(scCNA), "consensus_assay") == "integer") {
      assay <- "integer"
    }
  }
  if (assay == "integer") {
    mean_ploidy <- mean(SummarizedExperiment::colData(scCNA)$ploidy)
    ploidy_trunc <- 2 * round(mean_ploidy)
    color_heat <- structure(pals::ocean.balance(length(0:ploidy_trunc)), 
                            names = 0:ploidy_trunc)
    if (round(mean_ploidy) == 2) {
      color_heat <- structure(c("#3787BA", "#95B8C5", "#F0ECEB", 
                                "#D7A290", "#BF583B", "#8D1128", "#3C0912"), 
                              names = c("0", "1", "2", "3", "4", "5", "6"))
    }
  }else{
    color_heat = circlize::colorRamp2(breaks = c(-1,0,1), c("dodgerblue3", "white", "firebrick3"))
  }
  
  seg_data <- t(SummarizedExperiment::assay(scCNA, assay))
  chr_ranges <- as.data.frame(SummarizedExperiment::rowRanges(scCNA))
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
  if (any(chr_ranges$seqnames == "24") || any(chr_ranges$seqnames == 
                                              "Y") || any(chr_ranges$seqnames == "chrY")) {
    chr_binary <- rep(c(2, 1), length(chr_lengths)/2)
  }else {
    chr_binary <- c(rep(c(2, 1), (length(chr_lengths)/2)), 
                    2)
  }
  chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  chr_rl_c <- c(1, cumsum(chr_lengths))
  chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1], 
                       b = chr_rl_c[2:length(chr_rl_c)])
  chr_l_means <- round(rowMeans(chr_df))
  chrom.names <- c(1:22, "X", "Y")
  v <- vector(length = sum(chr_lengths), mode = "character")
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""
  chr_bar <- ComplexHeatmap::HeatmapAnnotation(chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data)], 
                                                                                    gp = grid::gpar(fontsize = 14)), df = as.character(chr[1:nrow(chr), 
                                                                                    ]), show_legend = FALSE, show_annotation_name = FALSE, 
                                               which = "column", col = list(df = c(`1` = "grey88", `2` = "black")))
  
  
  if (consensus == FALSE) {
    if (order_cells == "none") {
      seg_data_ordered <- seg_data
    }
    if (order_cells == "phylogeny") {
      if (ape::Ntip(phylo(scCNA)) == 0) {
        stop("No phylogeny detected in scCNA object. Use runPhylo")
      }
      tree <- phylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      seg_data_ordered <- seg_data[tree_tips_order, ]
    }
    if (order_cells == "hclust") {
      if (length(copykit::distMat(scCNA)) == 0) {
        message("No distance matrix detected in the scCNA object.")
        scCNA <- runDistMat(scCNA, metric = "euclidean", 
                            n_threads = n_threads)
      }
      if (nrow(as.matrix(copykit::distMat(scCNA))) != ncol(scCNA)) {
        stop("Number of samples in the distance matrix different\n          from number of samples in the scCNA object.\n        Perhaps you filtered your dataset? use copykit::runDistMat() to update it.")
      }
      hc <- fastcluster::hclust(distMat(scCNA), method = "ward.D2")
      seg_data_ordered <- seg_data[hc$order, ]
    }
    if (order_cells == "consensus_tree") {
      if (nrow(copykit::consensus(scCNA)) == 0) {
        scCNA <- calcConsensus(scCNA)
      }
      scCNA <- runConsensusPhylo(scCNA)
      consensus_by <- attr(copykit::consensus(scCNA), "consensus_by")
      meta <- as.data.frame(colData(scCNA)) %>% dplyr::select(sample, 
                                                              !!consensus_by)
      meta_info <- as.character(dplyr::pull(meta, !!consensus_by))
      tree <- consensusPhylo(scCNA)
      is_tip <- tree$edge[, 2] <= length(tree$tip.label)
      ordered_tips_index <- tree$edge[is_tip, 2]
      tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
        rev()
      meta_o <- meta[order(match(meta_info, tree_tips_order)),]
      n=1
      df = table(meta_o$subclones) %>% data.frame()
      vall =c()
      for(i in unique(meta_o$subclones)){
        length = df[df$Var1 == i,'Freq']
        v = rep(paste0('t',n),length)
        vall = c(vall,v)
        n=n+1
      }
      meta_o$treeorder = vall
      
      seg_data_ordered <- seg_data[meta_o$sample, ]
    }
  }else {
    scCNA <- runConsensusPhylo(scCNA)
    tree <- consensusPhylo(scCNA)
    is_tip <- tree$edge[, 2] <= length(tree$tip.label)
    ordered_tips_index <- tree$edge[is_tip, 2]
    tree_tips_order <- tree$tip.label[ordered_tips_index] %>% 
      rev()
    seg_data <- as.matrix(t(copykit::consensus(scCNA)))
    seg_data_ordered <- seg_data[tree_tips_order, ]
  }
  if (assay == "integer") {
    seg_data_int <- seg_data_ordered
    seg_data_int[seg_data_int > ploidy_trunc] <- ploidy_trunc
    if (rounding_error == TRUE) {
      int_nr <- as.matrix(SummarizedExperiment::assay(scCNA, 
                                                      "segment_ratios")) %*% diag(SummarizedExperiment::colData(scCNA)$ploidy)
      names(int_nr) <- names(SummarizedExperiment::assay(scCNA, 
                                                         "segment_ratios"))
      err <- abs(int_nr - assay(scCNA, "integer"))
      err <- err[rownames(seg_data_int)]
    }
  }
  message("Plotting Heatmap.")
  pt = unique(colData(scCNA)$patient)
  if(length(pt)==1){
    complex_args <- list(use_raster = T, column_title = "genomic coordinates", raster_quality = 5, 
                         column_title_gp = grid::gpar(fontsize = 18), column_title_side = "bottom", 
                         row_title = paste0(nrow(seg_data_ordered), " cells in ", pt),
                         # row_title_rot = 0,
                         row_title_gp = grid::gpar(fontsize = 18), top_annotation = chr_bar, 
                         cluster_rows = FALSE, border = TRUE, cluster_columns = FALSE, 
                         col= color_heat,
                         show_column_names = FALSE, show_row_names = FALSE, show_heatmap_legend = TRUE)
  }else{
    complex_args <- list(use_raster = T, column_title = "genomic coordinates", raster_quality = 5, 
                         column_title_gp = grid::gpar(fontsize = 18), column_title_side = "bottom", 
                         row_title = paste0(nrow(seg_data_ordered), " cells"),
                         # row_title_rot = 0,
                         row_title_gp = grid::gpar(fontsize = 18), top_annotation = chr_bar, 
                         cluster_rows = FALSE, border = TRUE, cluster_columns = FALSE, 
                         col= color_heat,
                         show_column_names = FALSE, show_row_names = FALSE, show_heatmap_legend = TRUE)
  }
  
  if (is.null(label)) {
    if (assay != "integer") {
      suppressMessages(do.call(ComplexHeatmap::Heatmap, 
                               c(list(matrix = log2(seg_data_ordered + 0.001), 
                                      heatmap_legend_param = list(title = "log2 (segratio)")), 
                                 complex_args)))
    }
    else {
      if (rounding_error == FALSE) {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                heatmap_legend_param = list(title = "copy number")), complex_args))
      }
      else {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                heatmap_legend_param = list(title = "rounding error"), 
                                                col = viridis::viridis(200)), complex_args))
      }
    }
  }
  else {
    metadata <- SummarizedExperiment::colData(scCNA) %>% 
      as.data.frame()
    if (consensus == FALSE) {
      metadata <- metadata[rownames(seg_data_ordered), 
      ]
    }
    metadata_anno_df <- metadata %>% dplyr::select(dplyr::all_of(label))
    if (consensus == TRUE) {
      cons_attr <- attr(copykit::consensus(scCNA), "consensus_by")
      metadata_anno_df <- metadata_anno_df[label] %>% dplyr::distinct()
      rownames(metadata_anno_df) <- metadata_anno_df %>% 
        dplyr::pull(!!cons_attr)
      metadata_anno_df <- metadata_anno_df[rownames(seg_data_ordered), 
                                           , drop = FALSE]
    }
    if (is.null(label_colors)) {
      j <- 15
      l <- 35
      cont_options <- c("D", "A", "E", "C", "B")
      cont_i <- 1
      label_colors <- list()
      label_colors <- c(list(superclones = superclones_pal(), 
                             subclones = subclones_pal(), filtered = c(removed = "#DA614D", 
                                                                       kept = "#5F917A"), is_normal = c(`TRUE` = "#396DB3", 
                                                                                                        `FALSE` = "#11181D")), label_colors)
      default_labels <- c("superclones", "subclones", "filtered", 
                          "is_normal")
      for (i in 1:length(label)) {
        if (any(str_detect(label[i], default_labels))) {
          label_colors[i] <- label_colors[default_labels[stringr::str_detect(label[i], 
                                                                             default_labels)]]
          names(label_colors)[i] <- label[i]
        }
        else if (is.numeric(dplyr::pull(metadata_anno_df, 
                                        label[i]))) {
          n = 300
          min_v = min(dplyr::pull(metadata_anno_df, label[i]))
          max_v = max(dplyr::pull(metadata_anno_df, label[i]))
          label_colors[i] <- list(circlize::colorRamp2(seq(min_v, 
                                                           max_v, length = n), viridis::viridis(n, option = cont_options[cont_i])))
          names(label_colors)[i] <- label[i]
          cont_i <- cont_i + 1
        }
        else {
          elements <- metadata_anno_df %>% dplyr::pull(label[i]) %>% 
            unique() %>% as.character() %>% sort()
          n <- length(elements)
          hues <- seq(j, 375, length = n + 1)
          hex <- hcl(h = hues, l = l, c = 100)[1:n]
          col <- structure(hex, names = elements)
          label_colors[i] <- list(col)
          names(label_colors)[i] <- label[i]
          j <- j + 15
          l <- l + 5
        }
      }
    }
    label_colors[sapply(label_colors, is.null)] <- NULL
    cluster_anno <- ComplexHeatmap::rowAnnotation(df = metadata_anno_df, 
                                                  col = label_colors, show_annotation_name = FALSE)
    if (!is.null(row_split)) {
      if (length(row_split) > 1) {
        stop("row_split length must be 1")
      }
      else {
        if (assay != "integer") {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = log2(seg_data_ordered + 0.001), 
                                                  row_split = dplyr::pull(metadata_anno_df, row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "log2 (segratio)")), 
                                             complex_args))
        }
        else {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  row_split = dplyr::pull(metadata_anno_df, row_split), left_annotation = cluster_anno, 
                                                  heatmap_legend_param = list(title = "copy number")), 
                                             complex_args))
        }
      }
    }
    else {
      if (assay != "integer") {
        do.call(ComplexHeatmap::Heatmap, c(list(matrix = log2(seg_data_ordered + 
                                                                0.001), left_annotation = cluster_anno, heatmap_legend_param = list(title = "log2 (segratio)")), 
                                           complex_args))
      }
      else {
        if (rounding_error == FALSE) {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = seg_data_int, 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "copy number"), 
                                                  col = color_heat), complex_args))
        }
        else {
          do.call(ComplexHeatmap::Heatmap, c(list(matrix = t(err), 
                                                  left_annotation = cluster_anno, heatmap_legend_param = list(title = "rounding error"), 
                                                  col = viridis::viridis(200)), complex_args))
        }
      }
    } 
  }
}
