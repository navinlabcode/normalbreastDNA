library(ArchR)
library(tidyverse)
library(copykit)
library(magrittr)
library(glue)

dblt_rate <- 4

samp_metadata <- readxl::read_xlsx("epi_enrich_coassay_path_v2_4sra.xlsx")
inputFiles <- samp_metadata$fragment_file
names(inputFiles) <- str_extract(inputFiles, pattern = "BCMHBCA[0-9]+[L|R]")

run = 'integration_19samples_0201'

## set up working space and parameters
setwd(glue("Projects/HBCA/coassay/atac/{run}"))
addArchRGenome("hg19")
addArchRThreads(threads = 100) 
marker_tbl <- read.csv("Projects/HBCA/allcell_markers.csv")
markerGenes_tbl <- marker_tbl%>% dplyr::filter(p_val_adj==0,avg_logFC>1.25)
seed=123
set.seed(seed)

## create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000,
  maxFrags = 100000,
  excludeChr = c("chrM", "chrY","KI270728.1", "KI270727.1", "GL000009.2", "GL000194.1",
                 "GL000205.2", "GL000195.1", "GL000219.1", "KI270734.1", "GL000213.1",
                 "GL000218.1", "KI270731.1", "KI270721.1", "KI270726.1", "KI270711.1",
                 "KI270713.1"),
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

## remove the cells from mda-mb231
proj63r <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "archr",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
metadf <- readxl::read_xlsx("epi_enrich_coassay_path_v2_4sra.xlsx")

samplesheet <- read_csv(glue("BCMHBCA63R/oldcellname_atac/GGGTTGAA/SampleSheet.csv"), skip = 5)

samplesheet <- samplesheet %>%
  mutate(index2_r = str_replace(index2, "GGGTTGAA", "")) %>%
  mutate(cellname = paste0("BCMHBCA63R","#",index, "+",index2_r)) %>%
  dplyr::select(Sample_ID,cellname)

# for(run in unique(merged_ckobj_filt$Patient_ID)[-which(unique(merged_ckobj_filt$Patient_ID)=="BCMHBCA63R")]){
for(id in unique(metadf$sample)[-which(unique(metadf$sample)%in%c("BCMHBCA63R"))]){
  p <- metadf$SampleSheet[metadf$sample==id]
  x <- read_csv(p, skip = 5)
  x <- x %>%
    mutate(index2_r = str_replace(index2, "GGGTTGAA", "")) %>%
    mutate(cellname = paste0(id,"#",index, "+",index2_r)) %>%
    dplyr::select(Sample_ID,cellname)
  samplesheet <- rbind(samplesheet,x)
}
matchname = left_join(as_tibble(proj63r$cellNames), samplesheet, by=c("value"="cellname"))

proj63r <- proj63r[matchname$value[!grepl("^C2_|^C50", matchname$Sample_ID)],]
saveArchRProject(ArchRProj = proj63r, outputDirectory = glue("rmmda"), dropCells = T, load = FALSE, overwrite = T)

ArrowFiles <- list.files("./rmmda/ArrowFiles/", pattern = "*.arrow")

## calculate doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

## create archr project and do sample qc
proj63r <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "archr",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)


## filter out cells possible doublets and cells with low TSSenrichment
idx <- proj63r$cellNames[which(proj63r$TSSEnrichment >= 7 & proj63r$TSSEnrichment <= 30 & proj63r$nFrags > 1000 & proj63r$nFrags < 40000)]


proj63r_rmdp <- filterDoublets(proj63r[idx,], filterRatio = dblt_rate)

df <- getCellColData(proj63r_rmdp, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 1)),
  ylim = c(0, quantile(df[,2], probs = 1))
) + geom_hline(yintercept = 7, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p

## dimension reduction and seurat clustering
proj63r_rmdp <- addIterativeLSI(
  ArchRProj = proj63r_rmdp,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 4, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.1,0.2,0.4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:30,
  force = T
)

proj63r_rmdp <- addClusters(
  input = proj63r_rmdp,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Seurat_clusters",
  resolution = 0.8,
  force = T
)

proj63r_rmdp <- addUMAP(
  ArchRProj = proj63r_rmdp, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = T,
  seed = seed
)

p1 <- plotEmbedding(ArchRProj = proj63r_rmdp, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj63r_rmdp, colorBy = "cellColData", name = "Seurat_clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = proj63r_rmdp, addDOC = FALSE, width = 5, height = 5)

projHeme2 <- addHarmony(
  ArchRProj = proj63r_rmdp,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)
projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Harmony_Seurat_clusters",
  resolution = 0.3,
  force = T
)

projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "Harmony",
  method = "scran",
  name = "Harmony_scran_clusters",
  knnAssign = 30,
  force = T
)

projHeme2 <- addUMAP(
  ArchRProj = projHeme2, 
  reducedDims = "Harmony", 
  name = "UMAPharmony", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)



p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Sample", embedding = "UMAPharmony")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Harmony_Seurat_clusters", embedding = "UMAPharmony")
# p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Harmony_hdbscan_clusters", embedding = "UMAPharmony") 
p3 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "Harmony_scran_clusters", embedding = "UMAPharmony") 

ggAlignPlots(p1, p2, p3, type = "h")
plotPDF(p3,p1,p2, name = "Plot-UMAP-Sample-Clusters-Harmony.pdf", ArchRProj = proj63r_rmdp, addDOC = FALSE, width = 5, height = 5)
p1 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "DoubletEnrichment", embedding = "UMAPharmony")
p2 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "DoubletScore", embedding = "UMAPharmony")
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-UMAP-Filtered-Doublet-Harmony.pdf", ArchRProj = proj63r_rmdp, addDOC = FALSE, width = 5, height = 5)

## get markers for each cluster
markersGS <- getMarkerFeatures(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Harmony_scran_clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes_tbl$gene,
  transpose = TRUE
)



## plot marker genes in a heatmap
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-scran-harmony", width = 12, height = 6, ArchRProj = projHeme2, addDOC = FALSE)

## plot marker genes in UMAP space
projHeme2 <- addImputeWeights(projHeme2)
#### choose marker genes to plot
markerGenes <- c('SAA2',"PI3", "PIGR", "S100A9",##lumsec
                 "CREM","RUNX3", "CD3D", "CD52", ##t cell
                 "TFF1","TFF3","MS4A7","AGR2", ##lumHR
                 "KRT14","KRT17", "KRT5", "CNN1",  ##basal
                 "COL3A1", "COL1A2", ##fibro
                 "MS4A6A", "GEM" ##myeloid
)

p <- plotEmbedding(
  ArchRProj = projHeme2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAPharmony",
  imputeWeights = getImputeWeights(projHeme2)
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = projHeme2, 
        addDOC = FALSE, width = 5, height = 5)
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
cowplot::ggsave2("archr/Plots/umap_gs_weighted_w_mkgene.pdf", width=8, height=6)



rnaobj <- readRDS("/volumes/USR2/rumwei/Analysis/HBCA/HBCA_new/hbca_sc_int_sub_20k_230105.rds")

projHeme2 <- addGeneIntegrationMatrix(
  ArchRProj = projHeme2, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Harmony",
  seRNA = rnaobj,
  addToArrow = FALSE,
  groupRNA = "celltype",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)


cM <- as.matrix(confusionMatrix(projHeme2$Harmony_scran_clusters, projHeme2$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) 

p5 <- plotEmbedding(ArchRProj = projHeme2, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAPharmony") 
plotPDF(p5, name = "Plot-UMAP-RNA-Mapped-Celltype-Clusters.pdf", ArchRProj = proj63r_rmdp, addDOC = FALSE, width = 5, height = 5)
saveArchRProject(ArchRProj = projHeme2, outputDirectory = glue("Save-Proj_dbr_{dblt_rate}"), dropCells = T, load = FALSE)

write_tsv(as_tibble(cM)%>%rownames_to_column("cluster"), glue("Save-Proj_dbr_{dblt_rate}/rna_inte_confusionmat.tsv"))

