
all_up_endo <- unique(unlist(up_endo_lists))
all_up_auto <- unique(unlist(up_auto_lists))


common_up_genes <- intersect(all_up_endo, all_up_auto)

all_down_endo <- unique(unlist(down_endo_lists))
all_down_auto <- unique(unlist(down_auto_lists))

common_down_genes <- intersect(all_down_endo, all_down_auto)


length(common_up_genes)   # Ortak up-regüle gen sayısı
length(common_down_genes) # Ortak down-regüle gen sayısı

library(writexl)

write_xlsx(list(
  Common_Up_Genes = data.frame(Gene = common_up_genes),
  Common_Down_Genes = data.frame(Gene = common_down_genes)
), path = "common_up_down_genes.xlsx")


write_xlsx(list(
  Common_Up_Genes = data.frame(Gene = common_up_genes)
  ), 
  path = "common_up_genes.xlsx")



write_xlsx(list(
  Common_Down_Genes = data.frame(Gene = common_down_genes)
), 
path = "common_down_genes.xlsx")














venn.diagram(
  x = list(
    Endometriosis = all_up_endo,
    Autoimmune    = all_up_auto
  ),
  category.names = c("Endometriosis", "Autoimmune"),
  filename = "venn_up_genes.png",
  output = TRUE,
  imagetype = "png",
  height = 4000,      # piksel cinsinden yüksekliği ayarladık
  width = 8000,       # genişlik
  resolution = 300,   # yüksek çözünürlük
  compression = "lzw",
  fill = c("tomato", "skyblue"),
  alpha = 0.5,
  cex = 2,            # sayılar için yazı boyutu
  cat.cex = 2,        # kategori isimleri için yazı boyutu
  cat.col = c("tomato", "skyblue"),
  main = "Upregulated Genes Overlap"
)


venn.diagram(
  x = list(
    Endometriosis = all_down_endo,
    Autoimmune    = all_down_auto
  ),
  category.names = c("Endometriosis", "Autoimmune"),
  filename = "venn_down_genes.png",
  output = TRUE,
  imagetype = "png",
  height = 4000,
  width = 8000,
  resolution = 300,
  compression = "lzw",
  fill = c("darkgreen", "gold"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.col = c("darkgreen", "gold"),
  main = "Downregulated Genes Overlap"
)




