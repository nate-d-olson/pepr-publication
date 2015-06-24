## Fixing pilon load code ...
load_pilon <- function(pilon_dir, db_con){
    if(.check_db_table("pilon_changes", db_con)){
        return()
    }
    # added miseq to file search need to change to analyze multiple platforms without hard coding platform
    changes_file <- list.files(pilon_dir, pattern = "*miseq.changes",
                               full.names = TRUE, recursive = TRUE)
    if(length(as.vector(changes_file))> 1){
        warning("more than one pilon changes file in directory skipping loading into database")
    }else if(length(as.vector(changes_file))==0){
        warning("no changes file in pilon directory")
    }else {
        changes_df <- readr::read_table(changes_file, 
                                        col_names = c("chrom_ref", 
                                                      "chrom_pilon",
                                                      "seq_ref",
                                                      "seq_pilon"))
        if(nrow(changes_df)== 0){
            warning("no changes in pilon changes file")
            return()
        }
        changes_df <- changes_df %>%
            tidyr::separate(col = chrom_ref, 
                            into = c("chrom_ref","coord_ref"), sep = ":") %>%
            tidyr::separate(col = chrom_pilon, 
                            into = c("chrom_pilon","coord_pilon"), sep = ":") %>% 
            dplyr::copy_to(dest = db_con,df = changes_df,
                       name = "pilon_changes", temporary = FALSE)
    }
}

changes_file <- "../pepr-data/MG002/CFSAN030013_pilon/CFSAN030013_miseq.changes"
changes_col_names <- c("chrom_ref", "chrom_pilon","seq_ref","seq_pilon")
#try using fill = TRUE to fix issue with unequal number of columns in rows
changes_df <- read.table(file = changes_file, header = FALSE, sep = " ", fill = TRUE,
                         #col.names = changes_col_names,
                         stringsAsFactors = FALSE)
library(dplyr)

changes_readr <- read_table(changes_file, col_names = changes_col_names)
library(ggplot2)
ggplot(changes_df) + 
    geom_bar(aes(x = as.numeric(coord_pilon))) + 
    xlim(0,rm_genome_size)

changes_file <- "../pepr-data/MG002/CFSAN030013_pilon/CFSAN030013_pgm.changes"
ggplot(changes_df) + 
    geom_bar(aes(x = as.numeric(coord_pilon))) + 
    xlim(0,rm_genome_size)

## changes are all due to an insertion at 1902661

