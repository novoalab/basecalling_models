############################
## Distribution of errors ##
############################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


#### 1: Data Import ----

# Task 1.1: Saving the file directories to 5mers


dir_default <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_hac","per_position"))
dir_ivt <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_ivt_hac","per_position"))
dir_sup <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_sup","per_position"))
dir_ivt_sup <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_ivt_sup","per_position"))



# Task 1.2: save all file names in specified directory as a character vector

file_list_default <- list.files(path = dir_default)
file_list_ivt <- list.files(path = dir_ivt)
file_list_sup <- list.files(path = dir_sup)
file_list_ivt_sup <- list.files(path = dir_ivt_sup)

# Task 1.4: Function for reading in all files from a directory
# Create a loop to read in every file of the directory and append it to the initialized data.frame plus add a new column that contains the name of the pool
# Comment: Could be sped up with fread and data.tables

read_dir <- function(file_list, work_dir){
  
  setwd(work_dir)
  
  dataset <- data.frame()
  
  for (i in 1:length(file_list)){
    temp_data <- read_csv(file_list[i]) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$sample <-  gsub("\\.csv", "", file_list[i])#clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}


default_raw <- read_dir(file_list = file_list_default, work_dir = dir_default)
ivt_raw <- read_dir(file_list = file_list_ivt, work_dir = dir_ivt)
sup_raw <- read_dir(file_list = file_list_sup, work_dir = dir_sup)
ivt_sup_raw <- read_dir(file_list = file_list_ivt_sup, work_dir = dir_ivt_sup)




### 2: Filtering ------

filter_data <- function(df){
  
  df_clean <- df %>% 
    separate(col = sample, into = c("sample","model"), sep = "-") %>% 
    separate(col = sample, into = c("sample","modification"), sep = "_") %>% # Separate into three columns to contain modifications
    filter(!(sample %in% c("RNAAB089716","RNAAB090763","RNA310520191","RNA234958"))) %>%  # filter out replicates (lower # of overall reads)
    mutate(modification = recode(modification, "pseudoU" = "Y",
                                 "unmodified" = "UNM",
                                 "5hmC" = "hm5C")) %>% 
    mutate(model = recode(model, "ivt" = "ivt+default")) %>% 
    
    # Separate ACGT into individual bases
    
    separate(col= ACGT, into = c("A","C","G","T"), sep = ":") %>% 
    
    # Change them into numeric 
    
    mutate(A = as.numeric(A),
           C = as.numeric(C),
           G = as.numeric(G),
           `T` = as.numeric(`T`)) %>% 
    
    # remove first an last 20nt from each kmer -> This is done by epinano_to_kmer usually
    group_by(`#Ref`,model,modification) %>% 
    # This filtering was manually checked -> WORKS
    filter(pos > 20 & pos < (max(pos)-20)) %>% 
    ungroup() %>% 
    # Recalcualte coverage just to make sure it matches up
    mutate(cov = rowSums(cbind(as.numeric(A), as.numeric(C), as.numeric(G), as.numeric(`T`)))) %>% 
    # Get the sum of all bases and 
    group_by(model,modification,base) %>% 
    summarize(A = sum(`A`),
              C = sum(`C`),
              G = sum(`G`),
              `T` = sum(`T`)) %>% 
    pivot_longer(cols = `A`:`T`,
                 names_to = "ACGT",
                 values_to = "cov") %>% 
    ungroup() %>% 
    #grouping by model is useless here but keeps the variable for later
    group_by(model,modification,base) %>% 
    mutate(percentage = (cov/sum(cov))) %>% 
    ungroup() %>% 
    mutate(ACGT = case_when(
      base == "A" & ACGT == "A" ~ "Ref",
      base == "T" & ACGT == "T" ~ "Ref",
      base == "G" & ACGT == "G" ~ "Ref",
      base == "C" & ACGT == "C" ~ "Ref",
      TRUE ~ ACGT
    ))
    
  return(df_clean)
}


default_sum <- filter_data(default_raw)
ivt_sum <- filter_data(ivt_raw)
sup_sum <- filter_data(sup_raw)
ivt_sup_sum <- filter_data(ivt_sup_raw)

all_clean <- rbind(default_sum,ivt_sum,sup_sum, ivt_sup_sum)

all_clean$modification <- factor(all_clean$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
all_clean$model <- factor(all_clean$model, levels = c("default","ivt+default","sup"))
all_clean$base <- factor(all_clean$base, levels = c("A","T","G","C"))
all_clean$ACGT <- factor(all_clean$ACGT, levels = c("A","T","C","G","Ref"))

rm(default_raw)
rm(ivt_raw)
rm(sup_raw)
rm(ivt_sup_raw)

### Fitler out ivt+sup

all_clean <- all_clean %>% 
  filter(model != "ivt+sup")

#### 3: DataViz ----

### Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 20),
  legend.position = "top",
  axis.text = element_text(size=20),
  axis.title=element_text(size=20),
  axis.text.x = element_text(size = 20,face = "bold"),
  axis.title.x = element_blank(),
  strip.text.x = element_text(size = 15, face = "bold"),
  axis.title.y = element_blank()
)  


get_mismatch_plots <- function(list_of_plots, df){
  
  
  for (i in 1:length(list_of_plots)){
    
    df %>% 
      filter(modification == list_of_plots[[i]][1] & model == list_of_plots[[i]][2]) %>% 
      ggplot(aes(x = base, y = percentage, fill = ACGT)) + 
      geom_bar(position="stack", stat="identity") +
      scale_y_continuous(labels = scales::percent,
                         expand = c(0,0)) +
      geom_text(aes(label=paste0(sprintf("%1.1f", percentage*100))),
                position=position_fill(vjust=0.5), colour="white") +
      scale_fill_manual(values=c("#1fab89","#eb4d55","#1e56a0", "#f0cf85", "#888888")) +
      theme_pubr() + t 
      #coord_cartesian(ylim=c(0.3,1))
    
    ggsave(filename = paste0(list_of_plots[[i]][1],"_",list_of_plots[[i]][2],".pdf"),
           plot = last_plot(),
           path = here("results_Gregor","plots","epinanoRMS", "per_position", "mismatch_dir","0-100"),
           width = 8,
           height = 6,
           units = "in")
    
    
  }
}


list_of_plots <- list(c("m6A","default"),c("m5C","default"),c("hm5C","default"),c("ac4C","default"),c("Y","default"),c("m1Y","default"),c("m5U","default"),c("UNM","default"),
                      c("m6A","ivt+default"),c("m5C","ivt+default"),c("hm5C","ivt+default"),c("ac4C","ivt+default"),c("Y","ivt+default"),c("m1Y","ivt+default"),c("m5U","ivt+default"),c("UNM","ivt+default"),
                      c("m6A","sup"),c("m5C","sup"),c("hm5C","sup"),c("ac4C","sup"),c("Y","sup"),c("m1Y","sup"),c("m5U","sup"),c("UNM","sup"),
                      c("m6A","ivt+sup"),c("m5C","ivt+sup"),c("hm5C","ivt+sup"),c("ac4C","ivt+sup"),c("Y","ivt+sup"),c("m1Y","ivt+sup"),c("m5U","ivt+sup"),c("UNM","ivt+sup"))

get_mismatch_plots(list_of_plots = list_of_plots, df = all_clean)


### Plot only central position for each modification

get_mismatch_plots_mod_site <- function(list_of_plots_mod_site, df){
  
  
  for (i in 1:length(list_of_plots)){
    
    df %>% 
      filter(modification == list_of_plots[[i]][1] & base == list_of_plots[[i]][2]) %>% 
      ggplot(aes(x = model, y = percentage, fill = ACGT)) + 
      geom_bar(position="stack", stat="identity") +
      scale_y_continuous(labels = scales::percent,
                         expand = c(0,0)) +
      geom_text(aes(label=paste0(sprintf("%1.1f", percentage*100))),
                position=position_fill(vjust=0.5), colour="white") +
      scale_fill_manual(values=c( "A" = "#1fab89","T"="#eb4d55","C"="#1e56a0", "G"="#f0cf85", "Ref"="#888888")) +
      theme_pubr() + t
    
    ggsave(filename = paste0(list_of_plots[[i]][1],"_",list_of_plots[[i]][2],".pdf"),
           plot = last_plot(),
           path = here("results_Gregor","plots","epinanoRMS", "per_position", "mismatch_dir","at_modified_site"),
           width = 8,
           height = 6,
           units = "in")
    
    
  }
}

list_of_plots <- list(c("m6A","A"),c("m5C","C"),c("hm5C","C"),c("ac4C","C"),c("Y","T"),c("m1Y","T"),c("m5U","T"))

get_mismatch_plots_mod_site(list_of_plots = list_of_plots, df = all_clean)

#### GGALLUVIAL ####
library(ggalluvial)

all_clean_m1Y <- all_clean %>%
  filter(modification == "m1Y" & model %in% c("default","sup") & base == "T") %>% 
  select(model,ACGT,percentage)

is_alluvia_form(all_clean_Y, axes = 1:3, silent = TRUE)

ggplot(all_clean_m1Y,
       aes(x = model, stratum = ACGT, alluvium = ACGT,
           y = percentage,
           fill = ACGT, 
           label = percentage)) +
  scale_x_discrete(expand = c(.05, .05)) +
  scale_y_continuous(label = scales::percent_format(),
                     expand = c(.001, .001)) +
  scale_fill_manual(values=c( "A" = "#1fab89","T"="#eb4d55","C"="#1e56a0", "G"="#f0cf85", "Ref"="#888888")) +
  geom_flow(alpha = .7,  decreasing = TRUE) +
  geom_stratum( alpha = 1, decreasing = TRUE) +
  geom_text(aes(label = paste0(scales::percent(..count.., accuracy = .1))), stat = "stratum", size = 4, decreasing = TRUE) +
  theme_pubr() + t + 
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.x = element_blank())


ggsave(filename = "m1Y_default_vs_sup.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS", "per_position", "mismatch_dir","alluvial_plots"),
       width = 8,
       height = 6,
       units = "in")


### PLotting the central base only with associated 


all_clean %>% 
  filter(modification == "Y" & base == "T") %>% 
  ggplot(aes(x = model, y = percentage, fill = ACGT)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(labels = scales::percent,
                     expand = c(0,0)) +
  geom_text(aes(label=paste0(sprintf("%1.1f", percentage*100))),
            position=position_fill(vjust=0.5), colour="white") +
  scale_fill_manual(values=c( "A" = "#1fab89","T"="#eb4d55","C"="#1e56a0", "G"="#f0cf85", "Ref"="#888888")) +
  theme_pubr() + t +
  coord_cartesian(ylim=c(0.3,1))

