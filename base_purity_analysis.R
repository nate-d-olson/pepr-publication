library(Hmisc)
library(magrittr)
library(ggplot2)
library(peprr)
library(dplyr)
library(tidyr)
library(knitr)

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
