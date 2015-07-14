pgm_med_read_length <- 0
pos_pur_gt99_both <- 0
pos_pur_gt99_one <- 0
pos_pur_gt95_one <- 0
min_max_pur_position  <- paste0(as.character(round(0, 2)*100),"%") # minimum purity for at least one platform (min, max), return %
max_contam <- paste0(as.character(round(1, 2)*100),"%") # maximum contamination per dataset