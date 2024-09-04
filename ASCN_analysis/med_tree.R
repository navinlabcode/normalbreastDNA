library(readxl)
library(ape)
library(ggtree)
library(glue)

samp="BCMHBCA44L"

ascn_armlevel_a <- read_xlsx(glue("{samp}/ascn_armlevel..xlsx"), sheet = "alleleA")
ascn_armlevel_long <- cbind(ascn_armlevel_a, 
                              read_tsv("hg19_ck_cytoband.txt") %>% 
                                dplyr::select(chr, start, end)) %>%
  dplyr::select(-chrarm) %>%
  gather(key="sample_id",value="cn_a", -chr, -start, -end)

ascn_armlevel_long$cn_b <- read_xlsx(glue("{samp}/ascn_armlevel..xlsx"), sheet = "alleleB") %>% 
  gather(key="sample_id",value="cn_b", -chrarm) %>%
  dplyr::pull(cn_b)

colnames(ascn_armlevel_long)[1] <- "chrom"
write_tsv(ascn_armlevel_long, glue("{samp}/medicc/input.tsv"))

tree <- read.tree(glue("{samp}/medicc/output/input_final_tree.new"))
list_samples <-
  split(colnames(ascn_armlevel_a)[-1], c("1","2","3","1","4","4","4","4"))
tree <- ggtree::groupOTU(tree, list_samples)
treeplt<-ggtree::ggtree(ape::ladderize(tree),
                        ladderize = F,
                        size = .3) +
  ggtree::geom_tippoint(size=3, aes(color=group))+
  ggtree::geom_tiplab(size=5.5, aes(color=group), hjust=-0.05,alpha=1)+
  geom_text(aes(x=branch, label=round(branch.length, 0), vjust=-.5), size = 3) +
  # scale_colour_manual(values = subclones_pal(),breaks = names(subclones_pal())) +
  theme(legend.position = "none") +
  geom_treescale(x = 10)

cowplot::save_plot(plot = treeplt, filename = glue("{samp}/medicc/output/treeplt.pdf"), base_width = 4, base_height = 4)

