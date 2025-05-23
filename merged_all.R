
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
