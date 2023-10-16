################################################
## Count kmers per bin from intersections    ##
###############################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


### Task 2: Import overlapped files generated with bedtools ----

overlap_hac <- read_tsv(here("data","eligos_2_human","hac","overlap_min_one_method.bed"), col_names = F)

non_overlap_hac <- read_tsv(here("data","eligos_2_human","hac","non_overlap_min_one_method.bed"), col_names = F)

overlap_ivt <- read_tsv(here("data","eligos_2_human","ivt","overlap_min_one_method.bed"), col_names = F)

non_overlap_ivt <- read_tsv(here("data","eligos_2_human","ivt","non_overlap_min_one_method.bed"), col_names = F)


### Task 3: Calculte number of kmers and extract top 10

summarize_data <- function(df,col_model,col_overlap){
  
  temp_df <- df %>% 
    count(X7) %>% 
    arrange(-n) %>% 
    #top_n(10) %>% 
    # add column that states whehter kmer is DRACH or no
    mutate(is_DRACH = factor(ifelse(str_detect(X7, "^[AGT][GA][A][C][ACT]$"), 1, 0), levels = c(1,0))) %>% 
    mutate(is_DRACH = case_when(
      is_DRACH == 1 ~ "DRACH",
      is_DRACH == 0 ~ "non-DRACH"
    ))
  
  temp_df$model <- col_model
  
  temp_df$overlap <- col_overlap
  
  
  return(temp_df)
  
}


overlap_hac_final <- summarize_data(df = overlap_hac, col_model = "default", col_overlap = "overlap")

non_overlap_hac_final <- summarize_data(df = non_overlap_hac, col_model = "default", col_overlap = "non-overlap")

overlap_ivt_final <- summarize_data(df = overlap_ivt, col_model = "ivt", col_overlap = "overlap")

## Remove final row as it is the eleventh value (equal n to postion 10)

overlap_ivt_final <- dplyr::slice(overlap_ivt_final, 1:(n() - 1)) 

non_overlap_ivt_final <- summarize_data(df = non_overlap_ivt, col_model = "ivt", col_overlap = "non-overlap")

final_df <- rbind(overlap_hac_final,non_overlap_hac_final,overlap_ivt_final,non_overlap_ivt_final)

final_df$model <- factor(final_df$model, levels = c("default","ivt"))

final_df_top5_non_overlap <- final_df %>% 
  filter(X7 %in% c("ACAGA","CCAGA","TCAGA","AAAGA","GCAGA","GGACT", "GGACA")) %>% 
  mutate(order = case_when(
    X7 == "ACAGA" ~ 1,
    X7 == "GGACT" ~ 2,
    X7 == "TCAGA" ~ 3,
    X7 == "CCAGA" ~ 4,
    X7 == "GGACA" ~ 5,
    X7 == "GCAGA" ~ 6,
    X7 == "AAAGA" ~ 7,
  ))
  

## Filter manually for sites previously determined by top_n() to include bo

### Task 4: Plot the distribution of sites for overlap and non_overlap

### Task 3: Plot relative ration of different stohciometries

# Set custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 20),
  legend.position = "top",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  axis.text.x = element_text(size = 15),
  strip.text.x = element_text(size = 15)
)  

# Plot barplots showing overall coutns

final_df_top5_non_overlap %>%
  filter(overlap == "non-overlap") %>% 
  #filter(X7 %in% c("GGACT","GGACA","GGACC","AGACT","TGACT")) %>% # These are top 5 overlapping
  ggplot(aes(x = reorder(X7, -order), fill = model, y=n)) +
  geom_bar(position = position_dodge2(preserve = "single"),stat = "identity") +
  geom_text(
    aes(x = X7, y = n, label = n, group = model), 
    size = 5, hjust = -0.5,
    position = position_dodge2(width = 1),
    inherit.aes = TRUE) +
  labs(x = "Top 5 kmers reported by each model",
       y = "number of reported sites") +
  scale_y_continuous(expand = c(0.02,0)) +
  scale_fill_manual(values=c('default' = '#E64B35FF', 'ivt' = "#00A087FF")) +
  theme_pubr() + t +
  coord_flip()

ggsave(filename = "non_overlapping_sites_top5.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","human_in_vivo", "kmer_analysis"),
       width = 12,
       height = 6,
       units = "in")

# Number of sites that are DRACH vs non

summarize_data <- function(df,col_model,col_overlap){
  
  temp_df <- df %>% 
    count(X7) %>% 
    arrange(-n) %>% 
    # add column that states whehter kmer is DRACH or no
    mutate(is_DRACH = factor(ifelse(str_detect(X7, "^[AGT][GA][A][C][ACT]$"), 1, 0), levels = c(1,0))) %>% 
    mutate(is_DRACH = case_when(
      is_DRACH == 1 ~ "DRACH",
      is_DRACH == 0 ~ "non-DRACH"
    ))
  
  temp_df$model <- col_model
  
  temp_df$overlap <- col_overlap
  
  return(temp_df)
  
}
  

test <- non_overlap_ivt %>% 
  count(X7) %>% 
  arrange(-n) %>% 
  # add column that states whehter kmer is DRACH or no
  mutate(is_DRACH = factor(ifelse(str_detect(X7, "^[AGT][GA][A][C][ACT]$"), 1, 0), levels = c(1,0))) %>% 
  mutate(is_DRACH = case_when(
    is_DRACH == 1 ~ "DRACH",
    is_DRACH == 0 ~ "non-DRACH"
  )) %>% 
  group_by(is_DRACH) %>% 
  summarize(sum = sum(n))
  


