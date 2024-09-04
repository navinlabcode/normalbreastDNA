plotHeatmap_state <- function(df=NULL,
                           chr_ranges = NULL,
                           raster_quality = 30) {
  # browser()
  color_state_pal <- structure(
    c(
      "#6E3562",
      "#DFB160",
      "#F0ECEB",
      "#29A0B1",
      "#FF8976"
    ),
    names = c("A", "AAB", "AB", "B", "ABB")
  )
  
  #obtaining data
  seg_data_ordered <- t(df)
  chr_lengths <- rle(as.numeric(chr_ranges$seqnames))$lengths
  
  if (any(chr_ranges$seqnames == "24") ||
      any(chr_ranges$seqnames == "Y") ||
      any(chr_ranges$seqnames == "chrY")) {
    chr_binary <- rep(c(2, 1), length(chr_lengths) / 2)
  } else {
    chr_binary <- c(rep(c(2, 1), (length(chr_lengths) / 2)), 2)
  }
  
  chr <-
    data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))
  
  # getting lengths for chr numbers annotation
  chr_rl_c <- c(1, cumsum(chr_lengths))
  
  # creating a data frame to calculate rowMeans
  chr_df <-
    data.frame(
      a = chr_rl_c[1:length(chr_rl_c) - 1],
      b = chr_rl_c[2:length(chr_rl_c)]
    )
  chr_l_means <- round(rowMeans(chr_df))
  
  chrom.names <- c(1:22, "X", "Y")
  
  # creating the vector for chr number annotations
  v <- vector(length = sum(chr_lengths), mode = "character")
  suppressWarnings(v[chr_l_means] <- chrom.names)
  v[is.na(v)] <- ""
  
  # chr bar with the chr names
  chr_bar <-
    ComplexHeatmap::HeatmapAnnotation(
      chr_text = ComplexHeatmap::anno_text(v[1:ncol(seg_data_ordered)],
                                           gp = grid::gpar(fontsize = 14)
      ),
      df = as.character(chr[1:nrow(chr), ]),
      show_legend = FALSE,
      show_annotation_name = FALSE,
      which = "column",
      col = list(df = c("1" = "grey88", "2" = "black"))
    )
  
  # list with arguments for complexheatmap
  show_heatmap_legend=TRUE
  complex_args <- list(
    use_raster = TRUE,
    raster_quality = raster_quality,
    column_title = "genomic coordinates",
    column_title_gp = grid::gpar(fontsize = 18),
    column_title_side = "bottom",
    top_annotation = chr_bar,
    cluster_rows = FALSE,
    border = TRUE,
    cluster_columns = FALSE,
    show_column_names = FALSE,
    show_row_names = FALSE,
    row_names_side = "left",
    show_heatmap_legend = show_heatmap_legend
  )
  
  message("Plotting Heatmap.")
  
    do.call(ComplexHeatmap::Heatmap, c(
      list(
        matrix = seg_data_ordered,
        heatmap_legend_param = list(title = "Allelic imbalance"),
        col = color_state_pal
      ),
      complex_args
    ))
}
