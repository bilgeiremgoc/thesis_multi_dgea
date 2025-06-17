library(readr)
library(dplyr)


cmap_df <- read_tsv("cmap.gct", skip = 2)

cmap_clean <- cmap_df %>%
  select(desc, pert_type, cell_iname, norm_cs) %>%
  mutate(norm_cs = as.numeric(norm_cs))



# En güçlü terapötik etkiler (negatif)
top_negative <- cmap_clean %>% arrange(norm_cs) %>% slice_head(n = 20)

# En güçlü benzer etkiler (pozitif)
top_positive <- cmap_clean %>% arrange(desc(norm_cs)) %>% slice_head(n = 20)

write_xlsx(top_negative, "cmap_top_negative_tau.xlsx")
write_xlsx(top_positive, "cmap_top_positive_tau.xlsx")

