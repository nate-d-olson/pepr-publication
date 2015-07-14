library(Hmisc)
library(magrittr)
library(ggplot2)
library(peprr)
library(dplyr)
library(tidyr)
library(knitr)
source("rm_metadata.R")
peprDB <- dplyr::src_sqlite(db_path)

.get_chrom_names <- function(db_con){
    dplyr::tbl(db_con, "pur_join") %>%
        dplyr::select(CHROM) %>%
        dplyr::collect() %>%
        .$CHROM %>% unique() %>%
        grep("novel", ., invert = TRUE, value = TRUE)
}
chrom_names <- .get_chrom_names(peprDB)


pur_dat <- dplyr::tbl(peprDB, "pur_join") %>%
    dplyr::filter(CHROM %in% chrom_names) %>%
    dplyr::collect() 
pur_dat_id <- pur_dat %>% 
    mutate(plat1_group = ifelse(plat1 < 0.99, "MiSeq-Low", "MiSeq-High")) %>% 
    mutate(plat2_group = ifelse(plat2 < 0.99, "PGM-Low", "PGM-High")) %>% 
    mutate(pur_group = paste(plat1_group, plat2_group, sep = " "))

pur_dat_id_filt <- pur_dat_id %>%
    dplyr::filter((plat1 < 0.99 | plat2 < 0.99))

min_val <- min(c(pur_dat_id_filt$plat1, pur_dat_id_filt$plat2)) 

## values for text
pos_pur_gt99_both <- pur_dat_id %>% filter(pur_group == "MiSeq-High PGM-High") %>% nrow()
pur_dat_max <- pur_dat_id %>% rowwise() %>% mutate(max_pur = max(plat1, plat2))
pos_pur_gt99_one <- pur_dat_max %>% filter(max_pur > 0.99) %>% nrow()
pos_pur_gt97_one <- pur_dat_max %>% filter(max_pur > 0.97) %>% nrow()
min_max_pur_position  <- paste0(as.character(round(min(pur_dat_max$max_pur), 2)*100),"%") # minimum purity for at least one platform (min, max), return %

