dge_results_up_genes_auto_GSE1919 <- read_excel("up_down_genes_auto/up/deg_GSE1919_RA_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE11501 <- read_excel("up_down_genes_auto/up/deg_GSE11501_CD_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE21942<- read_excel("up_down_genes_auto/up/deg_GSE21942_MS_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE40611 <- read_excel("up_down_genes_auto/up/deg_GSE40611_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE43591 <- read_excel("up_down_genes_auto/up/deg_GSE43591_MS_up_genes_05.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE50772 <- read_excel("up_down_genes_auto/up/deg_GSE50772_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE51092 <- read_excel("up_down_genes_auto/up/deg_GSE51092_SS_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE55235 <- read_excel("up_down_genes_auto/up/deg_GSE55235_RA_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE55457 <- read_excel("up_down_genes_auto/up/deg_GSE55457_RA_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE61635 <- read_excel("up_down_genes_auto/up/deg_GSE61635_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE66795 <- read_excel("up_down_genes_auto/up/deg_GSE66795_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE72326 <- read_excel("up_down_genes_auto/up/deg_GSE72326_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE77298 <- read_excel("up_down_genes_auto/up/deg_GSE77298_RA_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE81622 <- read_excel("up_down_genes_auto/up/deg_GSE81622_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE84844 <- read_excel("up_down_genes_auto/up/deg_GSE84844_SS_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE87466 <- read_excel("up_down_genes_auto/up/deg_GSE87466_UC_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE87473 <- read_excel("up_down_genes_auto/up/deg_GSE87473_UC__up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)

dge_results_up_genes_auto_GSE138198 <- read_excel("up_down_genes_auto/up/HT_deg_GSE138198_up_genes.xlsx") %>%
  filter(logFC > 1 & adj.P.Val < 0.05)






up_auto_lists <- list(
  GSE1919  = dge_results_up_genes_auto_GSE1919$gene_symbol,
  GSE11501 = dge_results_up_genes_auto_GSE11501$GeneSymbol,
  GSE21942 = dge_results_up_genes_auto_GSE21942$GeneSymbol,
  GSE40611 = dge_results_up_genes_auto_GSE40611$gene_symbol,
  GSE43591 = dge_results_up_genes_auto_GSE43591$GeneSymbol,
  GSE50772 = dge_results_up_genes_auto_GSE50772$GeneSymbol,
  GSE51092 = dge_results_up_genes_auto_GSE51092$GeneSymbol,
  GSE55235 = dge_results_up_genes_auto_GSE55235$gene_symbol,
  GSE55457 = dge_results_up_genes_auto_GSE55457$gene_symbol,
  GSE61635 = dge_results_up_genes_auto_GSE61635$GeneSymbol,
  GSE66795 = dge_results_up_genes_auto_GSE66795$GeneSymbol,
  GSE72326 = dge_results_up_genes_auto_GSE72326$GeneSymbol,
  GSE77298 = dge_results_up_genes_auto_GSE77298$gene_symbol,
  GSE81622 = dge_results_up_genes_auto_GSE81622$GeneSymbol,
  GSE84844 = dge_results_up_genes_auto_GSE84844$gene_symbol,
  GSE87466 = dge_results_up_genes_auto_GSE87466$GeneSymbol,
  GSE87473 = dge_results_up_genes_auto_GSE87473$GeneSymbol,
  GSE138198= dge_results_up_genes_auto_GSE138198$gene_symbol
)


up_auto_lists <- lapply(up_auto_lists, function(genes) {
  toupper(trimws(genes))
})



up_auto_counts <- table(unlist(up_auto_lists))
up_auto_df     <- as.data.frame(up_auto_counts)
colnames(up_auto_df) <- c("GeneSymbol", "Frequency")

up_auto_sorted <- up_auto_df[order(-up_auto_df$Frequency), ]

write_xlsx(up_auto_sorted, "total_up_autoimmune_genes_sorted.xlsx")

auto_top150_up_genes <- up_auto_df %>%
  arrange(desc(Frequency)) %>%
  head(150)

write.table(auto_top150_up_genes$GeneSymbol, 
            file = "auto_top150_common_up_genes.txt", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)

