read_all_sheets <- function(file_path) {
  sheets <- excel_sheets(file_path)
  all_data <- lapply(sheets, function(sheet) {
    read_xlsx(file_path, sheet = sheet)
  })
  bind_rows(all_data)
}


autoimmune_dgea <- read_all_sheets("autoimmune_dgea_all.xlsx")
endo_dgea <- read_all_sheets("endo_dgea_all.xlsx")


autoimmune_dgea <- autoimmune_dgea %>%
  mutate(GeneSymbol = toupper(as.character(GeneSymbol)))

endo_dgea <- endo_dgea %>%
  mutate(GeneSymbol = toupper(as.character(GeneSymbol)))


auto_up <- autoimmune_dgea %>% filter(logFC > 1, adj.P.Val < 0.05)
auto_down <- autoimmune_dgea %>% filter(logFC < -1, adj.P.Val < 0.05)
endo_up <- endo_dgea %>% filter(logFC > 1, adj.P.Val < 0.05)
endo_down <- endo_dgea %>% filter(logFC < -1, adj.P.Val < 0.05)

common_up_genes <- intersect(auto_up$GeneSymbol, endo_up$GeneSymbol)
common_down_genes <- intersect(auto_down$GeneSymbol, endo_down$GeneSymbol)


select_top <- function(df, gene_list) {
  df %>%
    filter(GeneSymbol %in% gene_list) %>%
    arrange(adj.P.Val) %>%
    distinct(GeneSymbol, .keep_all = TRUE)
}

auto_up_filtered <- select_top(auto_up, common_up_genes)
endo_up_filtered <- select_top(endo_up, common_up_genes)
auto_down_filtered <- select_top(auto_down, common_down_genes)
endo_down_filtered <- select_top(endo_down, common_down_genes)


up_merged <- inner_join(auto_up_filtered, endo_up_filtered, by = "GeneSymbol", suffix = c("_auto", "_endo"))
down_merged <- inner_join(auto_down_filtered, endo_down_filtered, by = "GeneSymbol", suffix = c("_auto", "_endo"))


up_merged$Regulation <- "Upregulated"
down_merged$Regulation <- "Downregulated"


final_up <- up_merged %>%
  select(GeneSymbol, logFC_auto, adj.P.Val_auto, logFC_endo, adj.P.Val_endo, Regulation)

final_down <- down_merged %>%
  select(GeneSymbol, logFC_auto, adj.P.Val_auto, logFC_endo, adj.P.Val_endo, Regulation)

# 9. Birleştir ve dışa aktar
final_combined <- bind_rows(final_up, final_down)


write_xlsx(final_combined, "ortak_gen_istatistikleri_tum_sayfalar.xlsx")

# Bilgi mesajı
cat("İşlem tamamlandı. 'ortak_gen_istatistikleri_tum_sayfalar.xlsx' adlı dosya oluşturuldu.\n")


final_up_filtered <- final_up %>% filter(logFC_auto > 1.5 & logFC_endo > 1.5)
final_down_filtered <- final_down %>% filter(logFC_auto < -1.5 & logFC_endo < -1.5)

top_up <- final_up %>% arrange(desc(logFC_auto + logFC_endo)) %>% slice_head(n = 300)
top_down <- final_down %>% arrange(logFC_auto + logFC_endo) %>% slice_head(n = 300)



gen_data <- read_xlsx("ortak_gen_istatistikleri_tum_sayfalar.xlsx")

up_genes <- gen_data %>%
  filter(Regulation == "Upregulated") %>%
  distinct(GeneSymbol) %>%
  slice_head(n = 150)

down_genes <- gen_data %>%
  filter(Regulation == "Downregulated") %>%
  distinct(GeneSymbol) %>%
  slice_head(n = 150)

write.table(up_genes, "cmap_up_top150.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

write.table(down_genes, "cmap_down_top150.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)








library(clusterProfiler)
library(org.Hs.eg.db)

gene_symbols <- unique(c(up_genes$GeneSymbol, down_genes$GeneSymbol))

# ENTREZ ID’ye çevir
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL",
                   toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO analizleri (BP, MF, CC)
go_bp <- enrichGO(entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
go_mf <- enrichGO(entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
go_cc <- enrichGO(entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)

# KEGG ve Reactome
kegg <- enrichKEGG(entrez_ids$ENTREZID, organism = "hsa")
reactome <- enrichPathway(entrez_ids$ENTREZID, organism = "human")

go_bp_df <- as.data.frame(go_bp)
go_mf_df <- as.data.frame(go_mf)
go_cc_df <- as.data.frame(go_cc)

write_xlsx(go_bp_df, "common_go_bp.xlsx")
write_xlsx(go_mf_df, "common_go_mf.xlsx")
write_xlsx(go_cc_df, "common_go_cc.xlsx")

kegg_df <- as.data.frame(kegg)

reactome_df <- as.data.frame(reactome)


write_xlsx(kegg_df, "common_kegg.xlsx")
write_xlsx(reactome_df, "common_reactome.xlsx")

barplot(go_bp, showCategory = 20, title = "GO Biological Processes")
barplot(go_mf, showCategory = 20, title = "GO Molecular Functions")
barplot(go_cc, showCategory = 20, title = "GO Cellular Components")



dotplot(kegg, showCategory = 15, title = "KEGG Pathway Enrichment")

dotplot(reactome, showCategory = 15)






