dge_results_down_genes_GSE7305 <- read_excel("up_down_genes_endo/deg_GSE7305_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE7307 <- read_excel("up_down_genes_endo/deg_GSE7307_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE11691 <- read_excel("up_down_genes_endo/deg_GSE11691_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE23339 <- read_excel("up_down_genes_endo/deg_GSE23339_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE25628 <- read_excel("up_down_genes_endo/deg_GSE25628_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE51981 <- read_excel("up_down_genes_endo/deg_GSE51981_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE87809 <- read_excel("up_down_genes_endo/deg_GSE87809_down_genes_endometriosis.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE105764 <- read_excel("up_down_genes_endo/deg_GSE105764_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE105765 <- read_excel("up_down_genes_endo/deg_GSE105765_down_genes_endometriosis.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)

dge_results_down_genes_GSE120103 <- read_excel("up_down_genes_endo/deg_GSE120103_down_genes.xlsx") %>%
  filter(logFC < -1 & adj.P.Val < 0.05)










down_endo_lists <- list(
  GSE7305 = dge_results_down_genes_GSE7305$GeneSymbol,
  GSE7307 = dge_results_down_genes_GSE7307$GeneSymbol,
  GSE11691 = dge_results_down_genes_GSE11691$GeneSymbol,
  GSE23339 = dge_results_down_genes_GSE23339$GeneSymbol,
  GSE25628 = dge_results_down_genes_GSE25628$GeneSymbol,
  GSE51981 = dge_results_down_genes_GSE51981$GeneSymbol,
  GSE87809 = dge_results_down_genes_GSE87809$GeneSymbol,
  GSE105764 = dge_results_down_genes_GSE105764$GeneSymbol,
  GSE105765 = dge_results_down_genes_GSE105765$GeneSymbol,
  GSE120103 = dge_results_down_genes_GSE120103$GeneSymbol
)


down_lists <- lapply(down_lists, function(vec) {
  toupper(trimws(vec))
})


common_down_genes_all <- Reduce(intersect, down_lists)
length(common_down_genes_all)

gene_counts <- table(unlist(down_lists))


common_down_genes_top5 <- names(gene_counts[gene_counts > 5])

sapply(down_lists, length)

length(common_up_genes_top5)
head(common_up_genes_top5)

gene_counts_df <- as.data.frame(gene_counts)

gene_counts_sorted <- gene_counts_df[order(-gene_counts_df$Freq), ]


top100_down_genes <- gene_counts_df %>%
  arrange(desc(Freq)) %>%
  head(100)

write.table(top100_down_genes$Var1, 
            file = "top100_common_down_genes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
write_xlsx(gene_counts_sorted, "total_down_endo_genes_sorted.xlsx")



