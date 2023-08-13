# infer ST CNAs with copykat result
# Functions ---------------------------------------------------------------
cpkt_load <- function(mat_dir='', prd_dir='') {

  cat('Loading CopyKat results...\n')
  cpkt_raw <- read.table(mat_dir, header=TRUE, check.names=FALSE)
  cpkt_mat <- t(cpkt_raw[, -c(1:3)]) %>% set_colnames(paste0('chr', cpkt_raw$chrom, '_', cpkt_raw$chrompos))
  cpkt_pred <- read.table(prd_dir, header=TRUE) %>% set_colnames(c('cell.names', 'cpkt.pred'))
  cpkt_pred <- cpkt_pred %>% dplyr::filter(cpkt.pred!='not.defined')
  cpkt_pred$cpkt.pred <- gsub('aneuploid', 'anpd', cpkt_pred$cpkt.pred)
  cpkt_pred$cpkt.pred <- gsub('diploid', 'dipd', cpkt_pred$cpkt.pred)
  cpkt_pred$cell.names <- make.names(cpkt_pred$cell.names)
  rownames(cpkt_pred) <- cpkt_pred$cell.names
  return(list(cpkt_mat=cpkt_mat, cpkt_pred=cpkt_pred))
}
multipcf_wrp <- function(mat, gamma=40) {
  cat('Run multipcf...\n')
  mat %>% t %>% data.frame(chrom=gsub('chr|\\_.*', '', rownames(.)) %>% as.numeric,
                           chrompos=gsub('.*\\_(.*)', '\\1', rownames(.)) %>% as.numeric, ., check.names = F) %>%
    copynumber::multipcf(., gamma = gamma, normalize = F, return.est = T) %>%
    extract2('segments') %>% set_rownames(., paste0('chr', .$chrom, .$arm) %>% make.unique()) %>% .[, -c(1:5)] %>% t
}
score_cut <- function(mat, G=3) {
  mat_out <- data.frame(cnvscore=apply(mat, 1, function(x) sqrt(mean(x^2)))) %>% set_rownames(rownames(mat))
  mat_out$cnvscore_cut <- paste0('cnvscore_', mclust::Mclust(mat_out$cnvscore, G=G)$classification)
  mat_out
}
umap_clust <- function(mat, metric='correlation', n_threads=10, n_neighbors = 15, min_dist=0, seed=42, grp_cut=10, k=15, eps=.5, minPts=10) {
  cat('UMAP...\n')
  umap_res <- uwot::umap(mat, metric=metric, n_threads=n_threads, n_neighbors=n_neighbors, min_dist=min_dist, ret_extra='fgraph')
  umap_out <- umap_res$embedding %>% data.frame %>% set_rownames(rownames(mat))
  
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
tm_frm <- rremove('xylab') + rremove('xy.text') + rremove('ticks') + rremove('grid') + theme(plot.title = element_text(hjust = 0.5))
tt_cen <- theme(plot.title = element_text(hjust = 0.5))
small_legend <- function(myPlot, pointSize = 2, textSize = 6, spaceLegend = 0.1) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
cnv_plot <- function(umap_df, mat, out_path=NULL) {
  ## Test input ##
  # mat=cpkt_res$cpkt_mat
  ####
  cat('CNV Plotting...\n')
  cnv_p0 <- ggplot(umap_df, aes(X1, X2)) + geom_point(aes(color=cnv_grp), size=.5)  +
    scale_color_brewer(palette = 'Set1') + scale_fill_brewer(palette = 'Set1') + theme_bw() + NoLegend() + tm_frm
  cnv_p1 <- ggplot(umap_df, aes(X1, X2)) + geom_point(aes(color=cpkt.pred), size=.5) +
    scale_color_brewer(palette = 'Set1') + scale_fill_brewer(palette = 'Set1') + theme_bw() + NoLegend() + tm_frm
  cnv_p2 <- ggplot(umap_df, aes(X1, X2)) + geom_point(aes(color=cnvscore_cut), size=.5) +
    scale_fill_viridis_d() + scale_color_viridis_d() + theme_bw() + NoLegend() + tm_frm
  cnv_p3 <- (ggplot(umap_df, aes(X1, X2)) + geom_point(aes(color=cnvscore), size=.5) + scale_color_viridis() + theme_bw() + tm_frm) %>% small_legend()
  cnv_p4 <- ggplot(umap_df, aes(X1, X2)) + geom_point(aes(color=tumor), size=.5)  +
    scale_color_brewer(palette = 'Set1') + scale_fill_brewer(palette = 'Set1') + theme_bw() + NoLegend() + tm_frm
  
  ## Heatmap ##
  rowannot <- umap_df %>% dplyr::select(-c(X1, X2)) %>% 
    group_by(cnv_grp) %>% sample_n(if(n()<200) n() else 200) %>% 
    arrange(desc(cnvscore)) %>% ungroup %>% data.frame %>% set_rownames(.$cell.names)
  row_idx <- rowannot$cell.names
  colannot <- data.frame(chrom=gsub('chr|_.*', '', colnames(mat)) %>% as.numeric %>% as.factor) %>% set_rownames(colnames(mat))
  
  ann_col_1 <- c(anpd='#E41A1C', dipd='#377EB8')
  ann_col_3 <- viridis(length(unique(rowannot$cnvscore_cut))) %>% set_names(unique(rowannot$cnvscore_cut) %>% sort)
  ann_col_4 <- ggpubr::get_palette('Set1', length(unique(rowannot$cnv_grp))) %>% set_names(unique(rowannot$cnv_grp) %>% sort)
  ann_col_5 <- ggpubr::get_palette('Set1', length(unique(rowannot$tumor))) %>% set_names(unique(rowannot$tumor) %>% sort)
  ann_col_all <- list(cpkt_pred=c(ann_col_1), cnvscore_cut=c(ann_col_3), cnv_grp=c(ann_col_4), tumor=c(ann_col_5))
  
  mat_temp <- as.matrix(mat[rownames(rowannot), ])
  ht <- ComplexHeatmap::pheatmap(mat_temp, cluster_cols = F, cluster_rows = F,
                                 show_rownames = F, show_colnames = F, 
                                 annotation_row = rowannot[, c(1, 4:7)], 
                                 row_split=factor(rowannot$tumor, levels=sort(unique(rowannot$tumor))),
                                 annotation_col = colannot, column_split=colannot$chrom, 
                                 col=colorRampPalette(c("navy", "white", "firebrick3"))(60), 
                                 annotation_colors = ann_col_all, 
                                 row_title_gp = gpar(fontsize = 6), row_title_rot = 0,
                                 breaks=col_break(mat_temp, len = 60, cut_ = c(-.1, .1)), 
                                 use_raster=TRUE, raster_quality=2)
  
  cnv_hm=grid::grid.grabExpr(ComplexHeatmap::draw(ht))
  p_cnv <- plot_grid(cnv_hm, cowplot::plot_grid(cnv_p0, cnv_p1, cnv_p2, cnv_p3, cnv_p4, ncol = 2), rel_widths = c(.6, .4))
  if (!is.null(out_path)) {
    pdf(file = out_path, width=18, height=8)
    print(p_cnv)
    dev.off()
  }
  p_cnv
}

st_analysis <- function(st_dir='', nFeature_Spatial_=.05, clust_res=0.3, sc_srt, pred_='celltype', metadata) {
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
  srt_temp$sample <- sample_temp
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
st_plot <- function(st_res, pt.size.factor=1.2, cnvevent_plot=c('chr1q', 'chr10q'), out_path=NULL) {
  srt_temp <- st_res$srt_out
  img_hires <- ggplot() + background_image(st_res$img_hires) + 
    coord_fixed(ratio=1, 
                xlim=c(0, dim(st_res$img_hires)[2]), 
                ylim=c(0, dim(st_res$img_hires)[1]), expand = FALSE, clip = 'on')
  
  sptdimplot <- SpatialDimPlot(srt_temp, group.by='seurat_clusters', label=T, pt.size.factor = pt.size.factor) + scale_fill_brewer(palette = 'Paired')
  dimplot <- DimPlot(srt_temp, reduction = 'umap', label=T) + scale_color_brewer(palette = 'Paired')
  srt_deg_top <- st_res$deg %>% group_by(cluster) %>% top_n(7, avg_log2FC)
  rna_htmap_0 <- DoHeatmap(subset(srt_temp, downsample=50), features=srt_deg_top$gene, size=4, raster=T) + theme(text = element_text(size=8)) + NoLegend()
  
  lbtf_dimplot <- DimPlot(srt_temp, group.by='pred.id', reduction='umap', label=T) + scale_color_brewer(palette = 'Paired')
  lbtf_sptdimplot <- SpatialDimPlot(srt_temp, group.by='pred.id', pt.size.factor = pt.size.factor) + scale_fill_brewer(palette = 'Paired')
  
  spftplot_cnvscore <- SpatialFeaturePlot(srt_temp, 'cnvscore', pt.size.factor=pt.size.factor, stroke = NA) + scale_fill_viridis() + theme_void()
  spdimplot_cnvscore_cut <- SpatialDimPlot(srt_temp, 'cnvscore_cut', pt.size.factor=pt.size.factor, stroke = NA, alpha = .8) + scale_fill_viridis_d() + theme_void() + NoLegend() + ggtitle('cnvscore_cut') + tt_cen
  spdimplot_cpkt_pred <- SpatialDimPlot(srt_temp, 'cpkt.pred', pt.size.factor=pt.size.factor, stroke = NA) + scale_fill_brewer(palette = 'Set1') + theme_void() + NoLegend() + ggtitle('cpkt.pred') + tt_cen
  spdimplot_cnv_group <- SpatialDimPlot(srt_temp, 'cnv_grp', pt.size.factor=pt.size.factor, stroke = NA, label.size = 3) +
    scale_fill_manual(values = get_palette('Set1', length(unique(na.omit(srt_temp$cnv_grp))))) + theme_void()
  spdimplot_tumor <- SpatialDimPlot(srt_temp, 'tumor', pt.size.factor=pt.size.factor, stroke = NA, label.size = 3) +
    scale_fill_manual(values = get_palette('Set1', length(unique(na.omit(srt_temp$tumor))))) + theme_void() + NoLegend() + ggtitle('aneuploidy') + tt_cen
  
  cnvevent_plist <- list()
  for (i in 1:length(cnvevent_plot)) {
    event_i <- cnvevent_plot[i]
    p_i <- SpatialFeaturePlot(srt_temp, features = event_i, pt.size.factor=pt.size.factor, stroke = NA) + 
      scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                           midpoint = 0, limits=c(-.1, .1), oob=scales::squish) + NoLegend() + ggtitle(event_i) + tt_cen
    cnvevent_plist[[i]] <- p_i
  }
  cnvevent_p <- plot_grid(plotlist = cnvevent_plist)
  
  p1 <- plot_grid(plot_grid(sptdimplot, lbtf_sptdimplot, nrow=2), rna_htmap_0, nrow=1)
  p2 <- plot_grid(plot_grid(spdimplot_cpkt_pred, spdimplot_cnvscore_cut, spdimplot_cnv_group, nrow=1, rel_widths = c(.5, .5, .6)), 
                  plot_grid(spdimplot_tumor, img_hires, cnvevent_p, ncol=3, rel_widths = c(.5, .5, .4)), nrow=2)
  p_out <- plot_grid(p1, p2, nrow=2)
  if (!is.null(out_path)) {
    pdf(file = out_path, width=18, height=36)
    print(p_out)
    dev.off()
  }
  p_out
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

grp_dbsc_otl <- function (data, group, eps, minPts) {
  data_out <- data.frame(group, dbsc_group=NA) %>% set_rownames(rownames(data))
  group <- as.character(group)
  grp_unq <- unique(group)
  for (i in 1:length(grp_unq)) {
    grp_i <- grp_unq[i]
    data_i <- data[group==grp_i, ]
    dbsc_i <- dbscan::dbscan(data_i, eps=eps, minPts=minPts)
    data_out$dbsc_group[group==grp_i] <- dbsc_i$cluster
  }
  data_out
}

