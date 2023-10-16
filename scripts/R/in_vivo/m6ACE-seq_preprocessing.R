###################################
## Preprocess m6ACE-seq sites   ##
##################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


### Task 2: Import cleaned and filter data ----

ace <- read_csv(here("known_sites","human_m6A", "raw","m6ACE-Seq.csv"))

### Task 3: Bring data into bed format compatible with eligos2 output

format_bed <- function(df){
  
  temp_df <- df %>% 
    # change chr annotation to fit ours
    mutate(Chr = sub("chr", "", Chr)) %>% 
    # make 0-based coordinated 1-based
    mutate(fill_0 = 0) %>% 
    select(Chr,Start,End,Gene,fill_0,`padj-WT`,Motif) %>% 
    arrange(Chr)
  
  return(temp_df)
}


### Task 5: Execute functions and save as bed

ace_sites_final <- format_bed(ace)

write_tsv(x = ace_sites_final,
          file = here("known_sites","human_m6A", "filtered","m6ACE-Seq.bed"),
          col_names = F)
