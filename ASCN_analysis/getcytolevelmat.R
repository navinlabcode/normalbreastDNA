getcytolevelmat <- function(origin.gr, 
                            cytoband.gr=cytoband.gr,
                            origin.mat
                            ){
  resl <- BiocParallel::bplapply(1:length(cytoband.gr), function(i){
    .cytoperline(origin.gr=origin.gr,
                 cytoband.perline.gr=cytoband.gr[i,],
                 origin.mat=origin.mat)
  })
 res <- as.data.frame(do.call(rbind, resl))
}


.cytoperline <- function(origin.gr,
                         cytoband.perline.gr,
                         origin.mat
                         ){
  hits <- findOverlaps(cytoband.perline.gr, origin.gr)
  c(length(hits), apply(origin.mat[subjectHits(hits),], 2, Mode))
}
