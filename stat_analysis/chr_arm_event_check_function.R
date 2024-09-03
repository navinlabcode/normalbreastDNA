library(RColorBrewer)

###############################   colormaker   #################################
colormaker <- function(x, pattern = "Set2"){
  n_color = length(unique(x)); 
  if(n_color>8){
    minor_col<-setNames(colorRampPalette(brewer.pal(8, pattern))(n_color),sort(unique(x),decreasing = F)) 
  }else{minor_col<-setNames(brewer.pal(n_color, pattern),sort(unique(x),decreasing = F))}
  return(minor_col[1:n_color])
}

#################################################### plotHeatmap2  ####################################################
plotHeatmap_chrarm <- function (scCNA, assay = "segment_ratios", n_threads = 40, 
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
  
  #main heatmap data
  seg_data <- t(SummarizedExperiment::assay(scCNA, assay))
  
  #top chromosome annotation
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
  
  #bottom chr_arm event annotation
  chr_ranges$chr_arm = paste0(chr_ranges$seqnames,chr_ranges$arm)
  chr_ranges$chr_arm = factor(chr_ranges$chr_arm,levels = unique(chr_ranges$chr_arm))
  chr_arm_lengths <- rle(as.numeric(chr_ranges$chr_arm))$lengths #chr arm length
  col_split = factor(rep.int(x = unique(chr_ranges$chr_arm), times = chr_arm_lengths),levels = unique(chr_ranges$chr_arm))

  #get events count
  meta = colData(scCNA) %>% data.frame() 
  chr_event_sum = meta[,grepl('^chr\\d{1,2}|^chrX',colnames(meta))] %>% colSums()
  names(chr_event_sum) = gsub('chr23','chrX',names(chr_event_sum))
  chr_anno_df = data.frame(chr_arm = levels(chr_ranges$chr_arm),Gain = NA, Loss = NA)
  for(i in names(chr_event_sum)){
    num = chr_event_sum[i]
    index = str_split(i,'_') %>% unlist()
    chr_anno_df[chr_anno_df$chr_arm == index[1],index[2]] =  num
  }
  head(chr_anno_df)
  chr_anno_df[is.na(chr_anno_df)] =0
  #add color
  ampcol = data.frame(Gain = 0:10, ampcol = colorRampPalette(c('white','red'))(11))
  delcol = data.frame(Loss = 0:10, delcol = colorRampPalette(c('white','blue'))(11))
  chr_anno_df = left_join(chr_anno_df,ampcol);chr_anno_df$ampcol[is.na(chr_anno_df$ampcol)] = ampcol$ampcol[11]
  chr_anno_df = left_join(chr_anno_df,delcol);chr_anno_df$delcol[is.na(chr_anno_df$delcol)] = delcol$delcol[11]
  bha = HeatmapAnnotation(gain = anno_block(gp = gpar(fill = chr_anno_df$ampcol), labels = chr_anno_df$Gain),
                    loss = anno_block(gp = gpar(fill = chr_anno_df$delcol), labels = chr_anno_df$Loss)) 

  
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
                         bottom_annotation = bha,
                         column_split = col_split, column_gap = unit(0, "mm"),
                         show_column_names = FALSE, show_row_names = FALSE, show_heatmap_legend = TRUE)
  }else{
    complex_args <- list(use_raster = T, column_title = "genomic coordinates", raster_quality = 5, 
                         column_title_gp = grid::gpar(fontsize = 18), column_title_side = "bottom", 
                         row_title = paste0(nrow(seg_data_ordered), " cells"),
                         # row_title_rot = 0,
                         row_title_gp = grid::gpar(fontsize = 18), top_annotation = chr_bar, 
                         cluster_rows = FALSE, border = TRUE, cluster_columns = FALSE, 
                         col= color_heat,
                         bottom_annotation = bha,
                         column_split = col_split, column_gap = unit(0, "mm"),
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

