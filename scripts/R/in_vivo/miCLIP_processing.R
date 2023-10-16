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

miclip_new <- miclip %>% 
  # modify 0-based coordinates to make them 1-based fitting remaining data sets
  mutate(...2 = ...3,
         ...3 = ...3 +1) %>% 
  #keep infromation on kmer
  mutate(...4 = str_remove(...4, ".*_"))

write_tsv(x = miclip_new,
          file = here("known_sites","human_m6A", "filtered","m6A_miCLIP_sites_HEK293T.bed"),
          col_names = F)

### Task 2.1. Save miCLIP files spanning different kmer lengths

miclip_1 <- miclip_new %>% 
  mutate(X2 = X2-1,
         X3 = X3 +1)

write_tsv(x = miclip_1,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_3mer.bed"),
          col_names = F)

miclip_2 <- miclip_new %>% 
  mutate(X2 = X2-2,
         X3 = X3 +2)

write_tsv(x = miclip_2,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_5mer.bed"),
          col_names = F)

miclip_3 <- miclip_new %>% 
  mutate(X2 = X2-3,
         X3 = X3 +3)

write_tsv(x = miclip_3,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_7mer.bed"),
          col_names = F)

miclip_4 <- miclip_new %>% 
  mutate(X2 = X2-4,
         X3 = X3 +4)

write_tsv(x = miclip_4,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_9mer.bed"),
          col_names = F)

miclip_5 <- miclip_new %>% 
  mutate(X2 = X2-5,
         X3 = X3 +5)

write_tsv(x = miclip_5,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_11mer.bed"),
          col_names = F)

miclip_10 <- miclip_new %>% 
  mutate(X2 = X2-10,
         X3 = X3 +10)

write_tsv(x = miclip_10,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_21mer.bed"),
          col_names = F)

miclip_50 <- miclip_new %>% 
  mutate(X2 = X2-50,
         X3 = X3 +50)

write_tsv(x = miclip_50,
          file = here("known_sites","human_m6A","hg38-m6A_sites_HEK293T_miCLIP_protein_coding_correct_start_101mer.bed"),
          col_names = F)



### Task 2: Visualize the relative improvement of ivt over hac at different kmers

### Values were manually extracted from performing bedtools intersect

df_kmer <- data.frame(model = c("default","ivt"),
                 "plus_minus_0" = c(1253,1392),
                 "plus_minus_1" = c(1280,1412),
                 "plus_minus_2" = c(2506,2815),
                 "plus_minus_3" = c(2643,2915),
                 "plus_minus_4" = c(3001,3098),
                 "plus_minus_5" = c(3101,3187),
                 "plus_minus_10" = c(3206,3264),
                 "plus_minus_50" = c(3699,3750))

df_plot <- df_kmer %>% 
  pivot_longer(cols = starts_with("plus_"),names_to = "kmer", values_to = "counts_overlap")

df_plot$model <- factor(df_plot$model, levels = c("default","ivt"))
df_plot$kmer <- factor(df_plot$kmer, levels = c("plus_minus_0","plus_minus_1","plus_minus_2","plus_minus_3","plus_minus_4","plus_minus_5","plus_minus_10","plus_minus_50"))

### Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  strip.text.x = element_text(size = 15)
)

### Plotting

ggplot(df_plot,aes(x = kmer, y = counts_overlap, fill = model)) +
  geom_col( position = position_dodge()) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(labels= c("+/-0","+/-1","+/-2","+/-3","+/-4","+/-5","+/-10","+/-50")) +
  scale_fill_manual(values = c('#E64B35FF',"#00A087FF")) +
  labs(x = paste("region considered for overlap", "surrounding reported miCLIP-site", sep = "\n"),
       y = "number of overlapping sites") +
  theme_pubr() + t

ggsave(filename = "overlap_with_neighboring_sites.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","miCLIP_analysis"),
       width = 6,
       height = 8,
       units = "in")

### Task 3: Convert this difference to ratios to additionally show that where the overall sites increase the most, the difference is also largest

calc_ratio <- function(df){
  
  temp_df <- df %>% 
    group_by(model,kmer) %>% 
    pivot_wider(names_from = model, values_from = counts_overlap) %>% 
    mutate(ratio_default = default/default * 100,
           ratio_ivt = ivt/default * 100) %>%
    pivot_longer(cols = starts_with("ratio"),
                 names_to = "model", values_to = "ratio") %>% 
    mutate(model = case_when(
      model == "ratio_default" ~ "default",
      model == "ratio_ivt" ~ "ivt"
    ))
  
  temp_df$kmer <- factor(temp_df$kmer, levels = c("plus_minus_0","plus_minus_1","plus_minus_2","plus_minus_3","plus_minus_4","plus_minus_5","plus_minus_10","plus_minus_50"))
  
  return(temp_df)
  
}

ratio_df <- calc_ratio(df_plot)

ggplot(ratio_df,aes(x = kmer, y = ratio, fill = model)) +
  geom_col( position = position_dodge()) +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent_format(scale = 1)) +
  scale_x_discrete(labels= c("+/-0","+/-1","+/-2","+/-3","+/-4","+/-5","+/-10","+/-50")) +
  scale_fill_manual(values = c('#E64B35FF',"#00A087FF")) +
  coord_cartesian(ylim=c(90,115)) +
  labs(x = paste("region considered for overlap", "surrounding reported miCLIP-site", sep = "\n"),
       y = "relative increase in recovered sites") +
  theme_pubr() + t

ggsave(filename = "overlap_with_neighboring_sites_relative.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","miCLIP_analysis"),
       width = 6,
       height = 8,
       units = "in")









