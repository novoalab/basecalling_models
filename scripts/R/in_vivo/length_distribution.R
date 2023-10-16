##############################################
### Model specific improvements on mapping ###
##############################################

######################
### Load Packages ####
######################

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)


##########################
#### 1: Data Import ######
##########################

# Task 1.1: Saving the file directories to 5mers


dir <- setwd(here("data","readlengths_human"))

file_list <- list.files(path = dir)


################################################################
# Task 1.4: Function for reading in all files from a directory #
################################################################

# Create a loop to read in every file of the directory and append it to the initialized data.frame plus add a new column that contains the name of the pool
# Comment: Could be sped up with fread and data.tables

read_dir <- function(file_list, work_dir){
  
  setwd(work_dir)
  
  dataset <- data.frame()
  
  for (i in 1:length(file_list)){
    temp_data <- read_table(file_list[i], col_names = F) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$species <-  gsub("\\.txt", "", file_list[i])#clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}


raw <- read_dir(file_list = file_list, work_dir = dir)

# Remove unnecessary column
raw <-  raw[,2:3]

############################################################################
# Task 1.5: Create function to split then bin summarizes and calculates n #
###########################################################################

bin_data <- function(df){
  
  temp_df <- df %>% 
    separate(col = species, into = c("model","replicate"), sep = "_") %>% 
    mutate(readlength_bin = case_when(
      X2 < 750 ~ "< 750",
      X2 >= 750 & X2 < 1500 ~ "750bp-1.5kb",
      X2 >= 1500 ~ " >= 1.5kb"
    )) %>% 
    group_by(model,replicate, readlength_bin) %>% 
    summarise(n = n())
  
  temp_df$readlength_bin <- factor(temp_df$readlength_bin, levels = c("< 750","750bp-1.5kb"," >= 1.5kb"))
  
  return(temp_df)
  
}

raw_binned <- bin_data(raw)

#####################################################################################
# Task 1.6: Create function to  calculate ratios and one to perform statistical test#
#####################################################################################

calc_ratio <- function(df){
  
  temp_df <- df %>% 
    group_by(replicate,readlength_bin) %>% 
    pivot_wider(names_from = model, values_from = n) %>% 
    mutate(ratio_default = default/default * 100,
           ratio_ivt = ivt/default * 100,
           ratio_sup = sup/default * 100) %>%
    pivot_longer(cols = starts_with("ratio"),
                 names_to = "model", values_to = "ratio") %>% 
    mutate(model = case_when(
      model == "ratio_default" ~ "default",
      model == "ratio_ivt" ~ "ivt",
      model == "ratio_sup" ~ "sup"
    ))
  
  temp_df$readlength_bin <- factor(temp_df$readlength_bin, levels = c("< 750","750bp-1.5kb"," >= 1.5kb"))
  
  return(temp_df)
  
}

ratios_df <- calc_ratio(raw_binned)


stat_test <- function(df){
  
  temp.stat.test <- df %>% 
    group_by(readlength_bin) %>% 
    t_test(ratio ~ model, ref.group = "default") %>% 
    adjust_pvalue(method = "hochberg") %>%
    add_significance("p.adj")
  
  temp.stat.test <- temp.stat.test %>%
    add_xy_position(x = "model", dodge = 0.8)
  
  temp.stat.test$y.position <- temp.stat.test$y.position + 0.05
  
  return(temp.stat.test)
  
}


stat.test  <- stat_test(ratios_df)

#######################################
# Task 1.7: Plot raw data and ratios #
######################################

# Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.text.x = element_text(angle=35, vjust = 0.65, hjust = 0.75),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  strip.text.x = element_text(size = 15)
)

# Task 1.7.1 raw data 
library(scales)

ggplot(raw_binned, aes(x = replicate, y = n, fill = model)) +
  geom_col(position = position_dodge()) +
  scale_y_log10(labels = label_comma())+
  #scale_y_continuous(labels = scales::comma, expand = c(0.025,0.025)) +
  labs(x="",y="number of reads") +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF")) +
  facet_wrap(~readlength_bin, nrow =1) +
  theme_pubr() + t

ggsave("mapping_length_log.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo","3_models"),
       width = 10,
       height = 5,
       units = "in")

# Task 1.7.2 ratios

ggbarplot(ratios_df, x="model", y="ratio", fill="model",
          add = "mean_sd", 
          facet.by = "readlength_bin", nrow = 1) +
  labs(x = '', y = "% reads mapped") +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = .025, size = 5, hide.ns = F) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0,0)) +
  scale_fill_manual(values = c('#E64B35FF', "#00A087FF",'#4DBBD5FF')) + theme_pubr() + t

ggsave("ratios_mapping_length_two_sided_padj_and_ns.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo","3_models"),
       width = 10,
       height = 5,
       units = "in")

###################################################
# Calculate global ratios for mRNA and mRNA-short #
###################################################

mRNA_per_rep <- raw_binned %>% 
  group_by(model, replicate) %>% 
  summarize(sum = sum(n))

mRNA_per_rep$method <- "mRNA"
mRNA_per_rep$mapping_params <- "minimap2 -ax map-ont -k14" 


mRNA_summary <- raw %>% 
  separate(col = species, into = c("model","replicate"), sep = "_") %>% 
  group_by(model, replicate) %>% 
  summarize(median = median(X2))

## Values obtained from MCF10 trna-seq runs

#trna_per_rep <- data.frame(model  = c("default", "sup", "default" ,"sup"),
#                          replicate  = c("rep1", "rep1", "rep2", "rep2"),
#                          sum  = c(1012220, 1077806, 722498, 760853),
#                          method  = c("tRNA", "tRNA", "tRNA","tRNA"))

#trna_per_rep$mapping_params <- "bwa mem xont2d W9-k5-T10"

unms <- data.frame(model  = c("default", "sup","default", "sup","default", "sup"),
         replicate  = c("rep1", "rep1","rep2", "rep2","rep3", "rep3"),
         sum  = c(40966, 50849, 137896, 226757, 60658, 91985),
         method  = c("mRNA-short", "mRNA-short","mRNA-short", "mRNA-short","mRNA-short", "mRNA-short"))

unms$mapping_params <- "minimap2 -ax map-ont -k14"

global_per_rep <- rbind(mRNA_per_rep, unms)

global_ratios <- global_per_rep %>% 
  filter(model != "ivt") %>% 
  group_by(replicate, method) %>% 
  pivot_wider(names_from = model, values_from = sum) %>% 
  mutate(ratio_default = (default/default * 100),
         ratio_sup = (sup/default * 100)) %>%
  pivot_longer(cols = starts_with("ratio"),
               names_to = "model", values_to = "ratio") %>% 
  mutate(model = case_when(
    model == "ratio_default" ~ "default",
    model == "ratio_sup" ~ "sup"
  ))


stat_test_global <- function(df){
  
  temp.stat.test <- df %>% 
    group_by(method) %>% 
    t_test(ratio ~ model, ref.group = "default", alternative = "less") %>% 
    adjust_pvalue(method = "hochberg") %>%
    add_significance("p.adj")
  
  temp.stat.test <- temp.stat.test %>%
    add_xy_position(x = "model", dodge = 0.8)
  
  temp.stat.test$y.position <- temp.stat.test$y.position + 0.05
  
  return(temp.stat.test)
  
}

stat.global <- stat_test_global(df = global_ratios)

###########################################
# Plot ratios of mRNA-long and mRNA short #
###########################################

# Task 1.7.2 ratios

ggbarplot(global_ratios, x="model", y="ratio", fill="model",
          add = c("mean_sd","jitter"), 
          facet.by = "method", nrow = 1) +
  labs(x = '', y = "% reads mapped") +
  stat_pvalue_manual(stat.global, label = "p.adj.signif", tip.length = .025, size = 5, hide.ns = F) +
  scale_y_continuous(labels = scales::percent_format(scale = 1),
                     expand = c(0,0)) +
  coord_cartesian(ylim = c(95, 175)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF')) + 
  theme_pubr() + t +
  theme(panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey90"))

ggsave("methods_mapping_improvement.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo","3_models"),
       width = 14,
       height = 7,
       units = "in")


### Readlgenths that generate above plot FASTQ


dir <- setwd(here("data","readlengths_UNM-S", "fastq"))

file_list <- list.files(path = dir)

raw_fq_unms <- read_dir(file_list = file_list, work_dir = dir)

process_fq_unms <- raw_fq_unms %>% 
  separate(species, into = c("model", "replicate"))

process_fq_unms$method <- "mRNA-short" 

## Get median

summary <- process_fq_unms %>% 
  group_by(replicate, model) %>% 
  summarize(median = median(X1))


########################################################
## Plot size distribution of mRNA and mRNA-short runs ##
########################################################

process_fq_unms %>% 
  ggplot(aes(x = replicate, y = X2, fill = model)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.1,
               position = position_dodge(width =0.4)) +
  # Central boxplot
  geom_boxplot(aes(fill = model),
               width=0.25, color="grey20",position = position_dodge(width =0.4),
               alpha = 1, lwd = 0.4, outlier.shape = NA) +
  # Distribution as half plot
  scale_y_continuous(limits = c(0,500)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF')) +
  labs(x="",y="read length [nt]") +
  #stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) 

ggsave("UNM-S_fastq_length.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo","3_models"),
       width = 10,
       height = 7,
       units = "in")


### Readlgenths that generate above plot BAM


dir <- setwd(here("data","readlengths_UNM-S", "bam_filtered"))

file_list <- list.files(path = dir)

raw_fq_unms <- read_dir(file_list = file_list, work_dir = dir)

process_fq_unms <- raw_fq_unms %>% 
  separate(species, into = c("model", "replicate"))

process_fq_unms$method <- "mRNA-short" 

process_fq_unms$model <- factor(process_fq_unms$model, level = c("sup","default"))

## Get median

summary <- process_fq_unms %>% 
  group_by(replicate, model) %>% 
  summarize(median = median(X1))

###############################################
## Plot size distribution of mRNA-short runs ##
###############################################

process_fq_unms %>% 
  ggplot(aes(x = replicate, y = X1, fill = model)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.2,
               position = position_dodge(width =0.85)) +
  # Central boxplot
  geom_boxplot(aes(fill = model),
               width=0.6, color="grey20",position = position_dodge(width =0.85),
               alpha = 1, lwd = 0.4, outlier.shape = NA) +
  # Distribution as half plot
  scale_y_continuous(limits = c(0,500)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF')) +
  labs(x="",y="read length [nt]") +
  #stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) +
  theme(axis.text.x = element_text(angle=0, vjust = 0, hjust = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor.x = element_line(color = "grey90"),
        ) +
  coord_flip()

ggsave("UNM-S_bam_filtered_length.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo","3_models"),
       width = 7,
       height = 7,
       units = "in")

#########################################
## Plot size distribution of mRNA runs ##
#########################################

human_lengths_processed <- raw %>% 
  separate(species, into = c("model", "replicate"))

human_lengths_processed$model <- factor(human_lengths_processed$model, levels = c("sup","default"))

## Get median

summary <- human_lengths_processed %>% 
  group_by(replicate, model) %>% 
  summarize(median = median(X2))

human_lengths_processed %>% 
  filter(model != "ivt") %>% 
  ggplot(aes(x = replicate, y = X2, fill = model)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.2,
               position = position_dodge(width =0.75)) +
  # Central boxplot
  geom_boxplot(aes(fill = model),
               width=0.6, color="grey20",position = position_dodge(width =0.75),
               alpha = 1, lwd = 0.4, outlier.shape = NA) +
  # Distribution as half plot
  scale_y_continuous(limits = c(0,3000)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF')) +
  labs(x="",y="read length [nt]") +
  #stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(angle=0, vjust = 0, hjust = 0),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor.x = element_line(color = "grey90")) +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) +
  coord_flip()

ggsave("human_bam_filtered_length.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo","3_models"),
       width = 7,
       height = 7,
       units = "in")






