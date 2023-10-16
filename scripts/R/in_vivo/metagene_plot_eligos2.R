############################################
## Metagene plots on eligos2 intersection ##
############################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Read in bed files

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("Guitar")

library(Guitar)


stBedFiles <- list(here("data","eligos_2_human","hac","overlap_min_one_method_metagene.bed"),here("data","eligos_2_human","ivt","overlap_min_one_method_metagene.bed"))


txdb <- makeTxDbFromGFF(file = "/Users/gregordiensthuber/cluster/gdiensthuber/references/H_sapiens/Homo_sapiens.GRCh38.107.gtf")

guitarTxdb <- makeGuitarTxdb(txdb = txdb, txPrimaryOnly = FALSE)

GuitarPlot(txTxd = txdb,
           stBedFiles = stBedFiles,
           miscOutFilePrefix = "default_and_ivt_sites_min_1_method_support",
           headOrtail = FALSE,
           enableCI = FALSE,
           mapFilterTranscript = TRUE,
           pltTxType = c("mrna"))

