dge_results_up_genes_GSE7305 <- readxl::read_excel("up_down_genes_endo/deg_GSE7305_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE7307 <- readxl::read_excel("up_down_genes_endo/deg_GSE7307_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE11691 <- readxl::read_excel("up_down_genes_endo/deg_GSE11691_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE23339 <- readxl::read_excel("up_down_genes_endo/deg_GSE23339_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE25628 <- readxl::read_excel("up_down_genes_endo/deg_GSE25628_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE51981 <- readxl::read_excel("up_down_genes_endo/deg_GSE51981_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE87809 <- readxl::read_excel("up_down_genes_endo/deg_GSE87809_up_genes_endometriosis.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE105764 <- readxl::read_excel("up_down_genes_endo/deg_GSE105764_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE105765 <- readxl::read_excel("up_down_genes_endo/deg_GSE105765_up_genes_endometriosis.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_GSE120103 <- readxl::read_excel("up_down_genes_endo/deg_GSE120103_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)


up_endo_lists <- list(
  GSE7305 = dge_results_up_genes_GSE7305$GeneSymbol,
  GSE7307 = dge_results_up_genes_GSE7307$GeneSymbol,
  GSE11691 = dge_results_up_genes_GSE11691$GeneSymbol,
  GSE23309 = dge_results_up_genes_GSE23339$GeneSymbol,
  GSE25628 = dge_results_up_genes_GSE25628$GeneSymbol,
  GSE51981 = dge_results_up_genes_GSE51981$GeneSymbol,
  GSE87809 = dge_results_up_genes_GSE87809$GeneSymbol,
  GSE105764 = dge_results_up_genes_GSE105764$GeneSymbol,
  GSE105765 = dge_results_up_genes_GSE105765$GeneSymbol,
  GSE120103 = dge_results_up_genes_GSE120103$GeneSymbol
)

up_lists <- lapply(up_lists, function(gene_vec) {
  toupper(trimws(gene_vec))
})

gene_counts <- table(unlist(up_lists))

# Örneğin en az 5 veri setinde geçen genleri al
common_up_genes_top5 <- names(gene_counts[gene_counts > 5])

sapply(up_lists, length)

length(common_up_genes_top5)
head(common_up_genes_top5)

gene_counts_df <- as.data.frame(gene_counts)

gene_counts_sorted <- gene_counts_df[order(-gene_counts_df$Freq), ]

writeLines(common_up_genes_top5, "common_up_genes_top5.txt")

gene_counts_df <- as.data.frame(gene_counts)

writexl::write_xlsx(gene_counts_sorted, "total_up_endo_genes_sorted.xlsx")

top532_genes <- gene_counts_df %>%
  arrange(desc(Freq)) %>%
  head(532)

write.table(top532_genes$Var1, 
            file = "top532_common_up_genes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)


