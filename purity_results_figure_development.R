## Revisions to purity analysis plots
library(peprr)
library(dplyr)
library(knitr)
library(magrittr)

source("rm_metadata.R")
peprDB <- dplyr::src_sqlite(db_path)
.get_chrom_names <- function(db_con){
    dplyr::tbl(db_con, "pur_join") %>%
        dplyr::select(CHROM) %>%
        dplyr::collect() %>%
        .$CHROM %>% unique() %>%
        grep("novel", ., invert = TRUE, value = TRUE)
}

purity_scatter_plot <- function (db_con,
                                 plat1_name = "Miseq",
                                 plat2_name = "PGM") {
    chrom_names <- .get_chrom_names(db_con)
    dplyr::tbl(db_con, "pur_join") %>%
        dplyr::filter(CHROM %in% chrom_names,
                      (plat1 < 0.99 | plat2 < 0.99)) %>%
        dplyr::collect() %>%
        ggplot2::ggplot() +
        ggplot2::geom_point(ggplot2::aes(x = plat1, y = plat2),
                            alpha = 0.5) +
        ggplot2::labs(x = plat1_name, y = plat2_name) + ggplot2::theme_bw()
}

db_con <- peprDB

chrom_names <- .get_chrom_names(db_con)
purity_dat_filt <- dplyr::tbl(db_con, "pur_join") %>%
    dplyr::filter(CHROM %in% chrom_names,
                  (plat1 < 0.99 | plat2 < 0.99)) %>%
    dplyr::collect()

purity_dat <- dplyr::tbl(db_con, "pur_join") %>%
    dplyr::filter(CHROM %in% chrom_names) %>%
    dplyr::collect()

min_val <- min(c(purity_dat_filt$plat1, purity_dat_filt$plat2)) %>% round(digits = 1)

purity_plot <- ggplot2::ggplot(purity_dat, ggplot2::aes(x = plat1, y = plat2)) +
    ggplot2::geom_point(data = purity_dat_filt, ggplot2::aes(x = plat1, y = plat2),
                        alpha = 0.5) +
    ggplot2::labs(x = "MiSeq", y = "PGM") + ggplot2::theme_bw() + 
    ggplot2::xlim(min_val, 1) + ggplot2::ylim(min_val, 1)
    

ggExtra::ggMarginal(purity_plot, type = "histogram")


set_pur_group <- function(plat1, plat2){
    if(plat1 < 0.99){
        if(plat2 < 0.99){
            return("both_low")
        }else{
            return("plat1_low")
        }
    }else if(plat2 < 0.99){
        return("plat2_low")
    }else{
        return("both_high")
    }
}
pur_dat_id <- purity_dat %>% mutate(plat1_group = ifelse(plat1 < 0.99, "low","high"),
                                    plat2_group = ifelse(plat2 < 0.99, "low","high"),
                                    pur_group = paste(plat1_group, plat2_group, sep = "_"))
pur_dat_id_filt <- pur_dat_id %>%
    dplyr::filter((plat1 < 0.99 | plat2 < 0.99))

ggplot2::ggplot(pur_dat_id, ggplot2::aes(x = plat1, y = plat2)) +
    ggplot2::geom_point(data = pur_dat_id_filt, ggplot2::aes(x = plat1, y = plat2, color = pur_group),
                        alpha = 0.5) +
    ggplot2::labs(x = "MiSeq", y = "PGM") + ggplot2::theme_bw() + 
    ggplot2::xlim(min_val, 1) + ggplot2::ylim(min_val, 1)

ggplot(pur_dat_id) + geom_bar(aes(x = pur_group)) + theme_bw() + scale_y_log10()
## A table might be a better way to present the data


########################## Purity Group table ##################################
chrom_names <- .get_chrom_names(db_con)
purity_dat_id <- dplyr::tbl(db_con, "pur_join") %>%
    dplyr::filter(CHROM %in% chrom_names) %>%
    dplyr::collect() %>% 
    mutate(plat1_group = ifelse(plat1 < 0.99, "low","high"),
            plat2_group = ifelse(plat2 < 0.99, "low","high"),
            pur_group = paste(plat1_group, plat2_group, sep = "_")) 
pur_dat_id %>% group_by(pur_group) %>% summarize(count = n()) %>% kable()


## using a table eliminates the need for the marginal plot
#ggExtra::ggMarginal(purity_plot, type = "histogram")

########################## Purity Metrics ##################################
## VCF header 
##ALT=<ID=X,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of high-quality bases">
##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of high-quality non-reference bases">
##FORMAT=<ID=DPR,Number=R,Type=Integer,Description="Number of high-quality bases observed for each allele">
##FORMAT=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-fwd, ref-reverse, alt-fwd and alt-reverse bases">
##FORMAT=<ID=SP,Number=1,Type=Integer,Description="Phred-scaled strand bias P-value">

db_con <-  dplyr::src_sqlite("/Volumes/NDO_PROJECT/current_projects/micro_rm/pepr-data/MG002/MG002.sqlite")
chrom_names <- .get_chrom_names(db_con)
datasets <- dplyr::tbl(db_con, "cb_miseq") %>% 
                select(SAMPLE) %>% collect() %>% 
                unique() %>% .$SAMPLE

cb_miseq <- dplyr::tbl(db_con, "cb_miseq") %>%
    dplyr::filter(CHROM %in% chrom_names, SAMPLE == datasets[1]) %>%         
    dplyr::select(CHROM, POS, RPB, MQB, BQB, MQSB, MQ0F) %>%
    collect() 
cb_miseq <- cb_miseq %>% mutate(PLAT = "miseq")

datasets <- dplyr::tbl(db_con, "cb_pgm") %>% 
    select(SAMPLE) %>% collect() %>% 
    unique() %>% .$SAMPLE

cb_pgm <- dplyr::tbl(db_con, "cb_pgm") %>%
    dplyr::filter(CHROM %in% chrom_names, SAMPLE == datasets[1]) %>%         
    dplyr::select(CHROM, POS, RPB, MQB, BQB, MQSB, MQ0F) %>%
    dplyr::mutate(PLAT = "pgm") %>% collect()

cb_metric_long <- cb_pgm %>%
    tidyr::gather("metric","value",3:7) 

cb_metric_summary <- cb_metric_long %>% 
                        group_by(PLAT, metric) %>% 
                        summarise(na_count = sum(is.na(value)),
                                  min_val = min(value, na.rm = TRUE),
                                  max_val = max(value, na.rm = TRUE), 
                                  med_val = median(value, na.rm = TRUE), 
                                  q_9 = quantile(value, 0.9, na.rm = TRUE), 
                                  q_99 = quantile(value, 0.99, na.rm = TRUE))
### metric distributions
library(Hmisc)
library(ggplot2)
cb_metric_long <- cb_pgm %>%
    tidyr::gather("metric","value",3:7) %>% 
    group_by(metric)  %>% 
    filter(!is.na(value)) %>% 
    mutate(bin = ntile(value, n = 10))

## distribution of values in bins for metrics
ggplot(cb_metric_long) + geom_boxplot(aes(x = as.factor(bin), y = value, color = metric))

cb_pgm_metric <- cb_pgm %>%
    tidyr::gather("metric","value",3:7) %>% 
    group_by(metric) cb_pgm %>%
    tidyr::gather("metric","value",3:7) %>% 
    group_by(metric)  %>% 
    filter(!is.na(value)) %>% 
    mutate(bin = ntile(value, n = 10)) %>% 
    filter(!is.na(value)) %>% 
    mutate(bin = ntile(value, n = 10))

##  * filter cb_miseq and cb_pgm


##  * combine cb filtered tables
##  * cut metrics at log10 intervals
##  * plot scatter - filtering on different cut levels
