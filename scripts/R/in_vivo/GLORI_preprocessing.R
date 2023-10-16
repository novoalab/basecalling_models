################################################
## Preprocess and annotate GLORI-seq sites   ##
###############################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


### Task 2: Import cleaned and filter data ----

all_sites_hek <- read_tsv(here("known_sites","human_m6A", "raw","GLORI-seq_m6A_sites.tsv"))

### Task 3: Bin data and sample sites 

bin <- function(df){
  
  set.seed(42)
  
  temp_df <- df %>% 
    # Bin by stoichiometry and reproducibility
    mutate(m6A_group = case_when(
      df$m6A_level_rep1 >= 0.75 & df$m6A_level_rep2 >= 0.75 ~ "HIGH",
      df$m6A_level_rep1 < 0.75 & df$m6A_level_rep1 >= 0.5 & df$m6A_level_rep2 < 0.75 & df$m6A_level_rep2 >= 0.5  ~ "MEDIUM-HIGH",
      df$m6A_level_rep1 < 0.5 & df$m6A_level_rep1 > 0.25 & df$m6A_level_rep2 < 0.5 & df$m6A_level_rep2 > 0.25  ~ "MEDIUM-LOW",
      df$m6A_level_rep1 <= 0.25 & df$m6A_level_rep2 <= 0.25  ~ "LOW"
    )) %>% 
    # remove NA which correspond to values that are not reproducible in same bin
    drop_na()
  
  return(temp_df)
} 

### Task 4: Bring data into bed format compatible with eligos2

format_bed <- function(df){
  
    temp_df <- df %>% 
    # change chr annotation to fit ours
    mutate(Chr = sub("chr", "", Chr)) %>% 
    # make 0-based coordinated 1-based
    mutate(Sites_2 = Sites-1) %>% 
    mutate(fill_0 = 0) %>% 
    mutate(fill_1 = 1) %>% 
    select(Chr,Sites_2,Sites,Gene,fill_0,Strand,Sites,Sites_2,Cluster_info,m6A_group,m6A_level_rep1, m6A_level_rep2) %>% 
    arrange(Chr)
  
    return(temp_df)
}

### Task 5: Execute functions and save as bed

GLORI_sites_final <- bin(all_sites_hek) %>% 
  format_bed()

write_tsv(x = GLORI_sites_final,
          file = here("known_sites","human_m6A","GLORI_sites_all.bed"),
          col_names = F)













