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

final_combined <- bind_rows(final_up, final_down)

unique_list_data <- unique(final_combined$GeneSymbol)
print(unique_list_data)

writexl::write_xlsx(final_combined, "ortak_gen_istatistikleri.xlsx")



final_up_filtered <- final_up %>% filter(logFC_auto > 1.5 & logFC_endo > 1.5)
final_down_filtered <- final_down %>% filter(logFC_auto < -1.5 & logFC_endo < -1.5)

top_up <- final_up %>% arrange(desc(logFC_auto + logFC_endo)) %>% slice_head(n = 300)
top_down <- final_down %>% arrange(logFC_auto + logFC_endo) %>% slice_head(n = 300)

writexl::write_xlsx(final_up_filtered, "final_up_filtered.xlsx")
writexl::write_xlsx(final_down_filtered, "final_down_filtered.xlsx")

final_combined_cleaned <- bind_rows(final_up_filtered, final_down_filtered)


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

















library(xCell)
library(readxl)
library(dplyr)

df1 <- read_xlsx("xcell_endo/deg_GSE7305_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE7305")

df2 <- read_xlsx("xcell_endo/deg_GSE7307_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE7307")

df3 <- read_xlsx("xcell_endo/deg_GSE11691_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE11691")

df4 <- read_xlsx("xcell_endo/deg_GSE23339_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE23339")

df5 <- read_xlsx("xcell_endo/deg_GSE25628_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE25628")

df6 <- read_xlsx("xcell_endo/deg_GSE51981_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE51981")

df7 <- read_xlsx("xcell_endo/deg_GSE105764_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE105764")


endo_all <- bind_rows(df1, df2, df3, df4, df5, df6, df7)


auto1 <- read_xlsx("xcell_auto/deg_GSE11501_CD_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE11501")

auto2 <- read_xlsx("xcell_auto/deg_GSE1919_RA_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE1919")

auto3 <- read_xlsx("xcell_auto/deg_GSE21942_MS_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE21942")

auto4 <- read_xlsx("xcell_auto/deg_GSE40611_SS_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE40611")

auto5 <- read_xlsx("xcell_auto/deg_GSE43591_MS_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE43591")

auto6 <- read_xlsx("xcell_auto/deg_GSE50772_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE50772")

auto7 <- read_xlsx("xcell_auto/deg_GSE51092_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE51092")

auto8 <- read_xlsx("xcell_auto/deg_GSE55235_RA_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE55235")

auto9 <- read_xlsx("xcell_auto/deg_GSE55457_RA_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE55457")

auto10 <- read_xlsx("xcell_auto/deg_GSE61635_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE61635")

auto11 <- read_xlsx("xcell_auto/deg_GSE66795_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE66795")

auto12 <- read_xlsx("xcell_auto/deg_GSE72326_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE72326")

auto13 <- read_xlsx("xcell_auto/deg_GSE77298_RA_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE77298")

auto14 <- read_xlsx("xcell_auto/deg_GSE81622_SLE_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE81622")

auto15 <- read_xlsx("xcell_auto/deg_GSE84844_SS_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE84844")

auto16 <- read_xlsx("xcell_auto/deg_GSE87466_UC_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE87466")

auto17 <- read_xlsx("xcell_auto/deg_GSE87473_UC_xCell_Results.xlsx") %>%
  mutate(Dataset = "GSE87473")

auto18 <- read_xlsx("xcell_auto/xcell_results_deg_GSE138198_HT.xlsx") %>%
  mutate(Dataset = "GSE138198")

auto_all <- bind_rows(auto1, auto2, auto3, auto4, auto5, auto6, auto7, auto8, auto9, auto10, auto11, auto12, auto13, auto14, auto15, auto16, auto17, auto18) 


library(tidyr)
library(dplyr)

# Endometriozis verisi uzun formata
endo_long <- endo_all %>%
  pivot_longer(cols = -c(Cell_Type, Dataset), 
               names_to = "Sample", 
               values_to = "Score") %>%
  mutate(Group = "Endometriosis")

# Otoimmün verisi uzun formata
auto_long <- auto_all %>%
  pivot_longer(cols = -c(Cell_Type, Dataset), 
               names_to = "Sample", 
               values_to = "Score") %>%
  mutate(Group = "Autoimmune")
endo_long_clean <- endo_long %>%
  drop_na(Score)

auto_long_clean <- auto_long %>%
  drop_na(Score)



merged_long <- bind_rows(endo_long_clean, auto_long_clean)


stats <- merged_long %>%
  group_by(Cell_Type) %>%
  summarise(
    p_value = tryCatch(wilcox.test(Score ~ Group)$p.value, error = function(e) NA)
  ) %>%
  mutate(adj_p_value = p.adjust(p_value, method = "bonferroni")) %>%
  arrange(p_value)


sig_cells <- stats %>% filter(adj_p_value < 0.05)


library(ggplot2)

top20 <- sig_cells$Cell_Type[1:20]

plot_data <- merged_long %>% filter(Cell_Type %in% top20)

ggplot(plot_data, aes(x = Group, y = Score, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ Cell_Type, scales = "free_y") +
  theme_minimal() +
  labs(title = "xCell Skor Karşılaştırması – En Anlamlı Hucre Tipleri")



library(pheatmap)

heat_df <- merged_long %>%
  filter(Cell_Type %in% sig_cells$Cell_Type) %>%
  group_by(Cell_Type, Group) %>%
  summarise(mean_score = mean(Score, na.rm = TRUE)) %>%
  pivot_wider(names_from = Group, values_from = mean_score)

heat_matrix <- as.matrix(heat_df[,-1])
rownames(heat_matrix) <- heat_df$Cell_Type

pheatmap(heat_matrix, cluster_rows = TRUE, cluster_cols = FALSE,
         main = "xCell Ort. Skorları – Anlamlı Hücre Tipleri",
         display_numbers = TRUE)




