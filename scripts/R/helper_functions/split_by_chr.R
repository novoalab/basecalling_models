###############################################################
### Split by chromosome and save each files one as a new df ###
###############################################################


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


### Task 2: Import cleaned and filter data ----

raw <- read.table(file = "/Users/gregordiensthuber/cluster/gdiensthuber/references/misc/human_m6A_sites/protein_coding_split_by_chr/Homo_sapiens.GRCh38.107_protein_coding_only.bed", sep = "\t")

list_of_chr <- raw %>% 
  group_by(V1) %>% 
  group_split()

### Task 3: Save each file as a for loop:


for (i in 1:length(list_of_chr)){
  
  write.table(list_of_chr[[i]], paste0("/Users/gregordiensthuber/cluster/gdiensthuber/references/misc/human_m6A_sites/protein_coding_split_by_chr/chr_",list_of_chr[[i]][[1]][1],".bed"),quote = F, col.names = F, row.names = F, sep = "\t")
  
}
