#####################################################################################
## Homogenize miCLIP coordinates and check sites surrounding reported miCLIP sites ##
#####################################################################################

# Load library

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr")

lapply(pkgs, library, character.only = TRUE)

### Task 1: Import data

miclip <- read_xlsx(here('known_sites','human_m6A','raw', 'm6A_miCLIP_sites_HEK293T.xlsx'), col_names = F)

### Task 2: Create expand single position into 9mer and save

miclip_minusstrand <- miclip %>% 
  filter(...6 == "-") %>% 
  # modify 0-based coordinates to make them 1-based fitting remaining data sets
  mutate(...3 = ...3 +1,
         ...2 = ...2 +1,
  ) %>% 
  #keep infromation on kmer
  mutate(...4 = str_remove(...4, ".*_"))


miclip_plusstrand <- miclip %>% 
  filter(...6 == "+") %>% 
  # modify 0-based coordinates to make them 1-based fitting remaining data sets
  mutate(...3 = ...3 -1,
         ...2 = ...2 -1,
  ) %>% 
  #keep infromation on kmer
  mutate(...4 = str_remove(...4, ".*_"))


miclip_final <- rbind(miclip_minusstrand, miclip_plusstrand)

write_tsv(x = miclip_final,
          file = here("known_sites","human_m6A", "filtered","m6A_miCLIP_sites_HEK293T.bed"),
          col_names = F)










