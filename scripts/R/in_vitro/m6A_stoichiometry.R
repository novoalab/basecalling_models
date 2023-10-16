##########################
### Data Preprocessing ###
##########################

######################
### Load Packages ####
######################

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr","rstatix","EnvStats")

lapply(pkgs, library, character.only = TRUE)


##########################
#### 1: Data Import ######
##########################

# Task 1.1: Saving the file directories to 5mers


dir_default <- setwd(here("data","epinano_m6A_stoich","default"))
dir_ivt <- setwd(here("data","epinano_m6A_stoich","ivt"))
dir_sup <- setwd(here("data","epinano_m6A_stoich","sup"))


##############################################################################
# Task 1.2: save all file names in specified directory as a character vector #
##############################################################################

file_list_default <- list.files(path = dir_default)
file_list_ivt <- list.files(path = dir_ivt)
file_list_sup <- list.files(path = dir_sup)

################################################################
# Task 1.4: Function for reading in all files from a directory #
################################################################

# Create a loop to read in every file of the directory and append it to the initialized data.frame plus add a new column that contains the name of the pool
# Comment: Could be sped up with fread and data.tables

read_dir <- function(file_list, work_dir){
  
  setwd(work_dir)
  
  dataset <- data.frame()
  
  for (i in 1:length(file_list)){
    temp_data <- read_tsv(file_list[i]) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$sample <-  gsub("\\.tsv", "", file_list[i])#clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}


default_raw <- read_dir(file_list = file_list_default, work_dir = dir_default)
ivt_raw <- read_dir(file_list = file_list_ivt, work_dir = dir_ivt)
sup_raw <- read_dir(file_list = file_list_sup, work_dir = dir_sup)

############################
#### 2: Data Wrangling #####
############################

####################################################################################
# Task 2.1: Function to generate new columns, remove replicates and rename samples #
####################################################################################

filter_data <- function(df){
  
  df_clean <- df %>% 
    separate(col = sample, into = c("expected_m6A","model"), sep = "-") %>% 
    mutate(expected_m6A = recode(expected_m6A,
                                 "UNM" = 0,
                                 "bc_1" = 12.5,
                                 "bc_2" = 75,
                                 "bc_3" = 50,
                                 "bc_4" = 25,
                                 "MOD" = 100)) %>% 
    # Calculate rowsums by splitting each column into 5 
    
    separate(mis, into = c("mis1", "mis2", "mis3", "mis4", "mis5"), sep = ",", remove = F, convert = T) %>% 
    separate(ins, into = c("ins1", "ins2", "ins3", "ins4", "ins5"), sep = ",", remove = F, convert = T) %>% 
    separate(del, into = c("del1", "del2", "del3", "del4", "del5"), sep = ",", remove = F, convert = T) %>% 
    separate(mean_q, into = c("mq1", "mq2", "mq3", "mq4", "mq5"), sep = ",", remove = F, convert = T) %>% 
    
    mutate(
      
      # Calculate corrected mismatch frequency so that sum of mis-,ins-,del- and match frequency = 1
      
      c_mis1 = mis1 * (1-(del1 + ins1)),
      c_mis2 = mis2 * (1-(del2 + ins2)),
      c_mis3 = mis3 * (1-(del3 + ins3)),
      c_mis4 = mis4 * (1-(del4 + ins4)),
      c_mis5 = mis5 * (1-(del5 + ins5)),
      
      # Calculate the Summed error per position and kmer; as.numeric is no longer needed as convert = T in separate
      
      SumErr_pos_minus_2 = rowSums(cbind(as.numeric(c_mis1), as.numeric(ins1), as.numeric(del1))),
      SumErr_pos_minus_1 = rowSums(cbind(as.numeric(c_mis2), as.numeric(ins2), as.numeric(del2))),
      SumErr_pos_0 = rowSums(cbind(as.numeric(c_mis3), as.numeric(ins3), as.numeric(del3))),
      SumErr_pos_plus_1 = rowSums(cbind(as.numeric(c_mis4), as.numeric(ins4), as.numeric(del4))),
      SumErr_pos_plus_2 = rowSums(cbind(as.numeric(c_mis5), as.numeric(ins5), as.numeric(del5))),
      SumErr_5mer = rowSums(cbind(as.numeric( SumErr_pos_minus_2), as.numeric(SumErr_pos_minus_1), as.numeric(SumErr_pos_0),
                                  as.numeric(SumErr_pos_plus_1), as.numeric(SumErr_pos_plus_2))),
      
      # Summed Error for each error feature and 5mer
      SumMis_5mer = rowSums(cbind(as.numeric(c_mis1), as.numeric(c_mis2), as.numeric(c_mis3), as.numeric(c_mis4), as.numeric(c_mis5))),
      
      SumIns_5mer = rowSums(cbind(as.numeric(ins1), as.numeric(ins2), as.numeric(ins3), as.numeric(ins4), as.numeric(ins5))),
      
      SumDel_5mer = rowSums(cbind(as.numeric(del1), as.numeric(del2), as.numeric(del3), as.numeric(del4), as.numeric(del5))),
      
      SumQ_5mer = rowSums(cbind(as.numeric(mq1), as.numeric(mq2), as.numeric(mq3), as.numeric(mq4), as.numeric(mq5))),
      
      log_SumMis_5mer = log10(SumMis_5mer),
      log_SumIns_5mer = log10(SumIns_5mer),
      log_SumDel_5mer  = log10(SumDel_5mer),
      log_SumErr_5mer = log10(SumErr_5mer),
      log_SumQ_5mer = log10(SumQ_5mer))
  
  
  return(df_clean)
}

default_clean <- filter_data(default_raw)
ivt_clean <- filter_data(ivt_raw)
sup_clean <- filter_data(sup_raw)

rm(default_raw)
rm(ivt_raw)
rm(sup_raw)

all_clean <- rbind(default_clean,ivt_clean,sup_clean)

###################################################################################################
# Task 2.2: Function to calculate delta SumErr and normalized delta SumErr for center mod of each #
###################################################################################################

filter_central <- function(df){
  
  temp_df <- df
  temp_df <- temp_df %>% 
    # Filter for modified central bases for each modification
    filter((ref_base == "A" )) %>% 
    # Annotate DRACH
    mutate(is_DRACH = factor(ifelse(str_detect(base, "^[AGT][GA][A][C][ACT]$"), 1, 0), levels = c(1,0))) %>% 
    # group-by X.Ref,base,model,ref_pos
    group_by(X.Ref,base,model,ref_pos) %>% 
    # Calculate delta per postion and on Summed values
    mutate(delta_SumErr = SumErr_5mer - SumErr_5mer[expected_m6A == 0],
           delta_SumMis = SumMis_5mer - SumMis_5mer[expected_m6A == 0],
           delta_SumIns = SumIns_5mer - SumIns_5mer[expected_m6A == 0],
           delta_SumDel = SumDel_5mer - SumDel_5mer[expected_m6A == 0],
           delta_SumQ = SumQ_5mer - SumQ_5mer[expected_m6A == 0],
           
           # Per position delta
           d_SumErr_pos_minus_2 = SumErr_pos_minus_2 - SumErr_pos_minus_2[expected_m6A == 0],
           d_SumErr_pos_minus_1 = SumErr_pos_minus_1 - SumErr_pos_minus_1[expected_m6A == 0],
           d_SumErr_pos_0 = SumErr_pos_0 - SumErr_pos_0[expected_m6A == 0 ],
           d_SumErr_pos_plus_1 = SumErr_pos_plus_1 - SumErr_pos_plus_1[expected_m6A == 0 ],
           d_SumErr_pos_plus_2 = SumErr_pos_plus_2 - SumErr_pos_plus_2[expected_m6A == 0 ],
           
           d_Mis_pos_minus_2 = as.numeric(mis1) - as.numeric(mis1)[expected_m6A == 0 ],
           d_Mis_pos_minus_1 = as.numeric(mis2) - as.numeric(mis2)[expected_m6A == 0 ],
           d_Mis_pos_0 = as.numeric(mis3) - as.numeric(mis3)[expected_m6A == 0 ],
           d_Mis_pos_plus_1 = as.numeric(mis4) - as.numeric(mis4)[expected_m6A == 0 ],
           d_Mis_pos_plus_2 = as.numeric(mis5) - as.numeric(mis5)[expected_m6A == 0 ],
           
           d_Ins_pos_minus_2 = as.numeric(ins1) - as.numeric(ins1)[expected_m6A == 0 ],
           d_Ins_pos_minus_1 = as.numeric(ins2) - as.numeric(ins2)[expected_m6A == 0 ],
           d_Ins_pos_0 = as.numeric(ins3) - as.numeric(ins3)[expected_m6A == 0 ],
           d_Ins_pos_plus_1 = as.numeric(ins4) - as.numeric(ins4)[expected_m6A == 0 ],
           d_Ins_pos_plus_2 = as.numeric(ins5) - as.numeric(ins5)[expected_m6A == 0 ],
           
           d_Del_pos_minus_2 = as.numeric(del1) - as.numeric(del1)[expected_m6A == 0 ],
           d_Del_pos_minus_1 = as.numeric(del2) - as.numeric(del2)[expected_m6A == 0 ],
           d_Del_pos_0 = as.numeric(del3) - as.numeric(del3)[expected_m6A == 0 ],
           d_Del_pos_plus_1 = as.numeric(del4) - as.numeric(del4)[expected_m6A == 0 ],
           d_Del_pos_plus_2 = as.numeric(del5) - as.numeric(del5)[expected_m6A == 0 ],
           
           d_Q_pos_minus_2 = as.numeric(mq1) - as.numeric(mq1)[expected_m6A == 0 ],
           d_Q_pos_minus_1 = as.numeric(mq2) - as.numeric(mq2)[expected_m6A == 0 ],
           d_Q_pos_0 = as.numeric(mq3) - as.numeric(mq3)[expected_m6A == 0 ],
           d_Q_pos_plus_1 = as.numeric(mq4) - as.numeric(mq4)[expected_m6A == 0 ],
           d_Q_pos_plus_2 = as.numeric(mq5) - as.numeric(mq5)[expected_m6A == 0 ],
           
           # Per position Signal to Noise with added pseudocount
           s_to_n_pos_minus_2 = ((SumErr_pos_minus_2 + 0.01) - (SumErr_pos_minus_2[expected_m6A == 0 ] + 0.01)) / (SumErr_pos_minus_2[expected_m6A == 0] + 0.01),
           s_to_n_pos_minus_1 = ((SumErr_pos_minus_1 + 0.01) - (SumErr_pos_minus_1[expected_m6A == 0 ] + 0.01)) / (SumErr_pos_minus_1[expected_m6A == 0] + 0.01),
           s_to_n_pos_0 = ((SumErr_pos_0 + 0.01) - (SumErr_pos_0[expected_m6A == 0 ] + 0.01)) / (SumErr_pos_0[expected_m6A == 0] + 0.01),
           s_to_n_pos_plus_1 = ((SumErr_pos_plus_1 + 0.01) - (SumErr_pos_plus_1[expected_m6A == 0 ] + 0.01)) / (SumErr_pos_plus_1[expected_m6A == 0] + 0.01),
           s_to_n_pos_plus_2 = ((SumErr_pos_plus_2 + 0.01) - (SumErr_pos_plus_2[expected_m6A == 0 ] + 0.01)) / (SumErr_pos_plus_2[expected_m6A == 0 ] + 0.01),
    ) %>% 
    dplyr::select(-mis1, -mis2, -mis3, -mis4, -mis5, -ins1, -ins2, -ins3, -ins4, -ins5, -del1, -del2, -del3, -del4, -del5,
                  -mq1, -mq2, -mq3, -mq4, -mq5)
  
}


all_clean_filtered <- all_clean %>% 
  filter_central()


### Plot Expected vs observed

### Custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "right",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15)
)  


### Point plot with fitted line GGACT

all_clean_filtered_GGACT <- all_clean_filtered %>% 
  filter(base == "GGACT",
         d_SumErr_pos_0 >= 0)

default_model <- all_clean_filtered_GGACT %>% 
  filter(model == "default")

default_model <- loess(d_SumErr_pos_0 ~ expected_m6A, data = default_model)

rmse_default <- sqrt(mean(resid(default_model)^2))

ivt_model <- all_clean_filtered_GGACT %>% 
  filter(model == "ivt",
         d_SumErr_pos_0 >= 0)

ivt_model <- loess(d_SumErr_pos_0 ~ expected_m6A, data = ivt_model)

rmse_ivt <- sqrt(mean(resid(ivt_model)^2))

sup_model <- all_clean_filtered_GGACT %>% 
  filter(model == "sup",
         d_SumErr_pos_0 >= 0)

sup_model <- loess(d_SumErr_pos_0 ~ expected_m6A, data = sup_model)

rmse_sup <- sqrt(mean(resid(sup_model)^2))

### Print rmse to stdout

rmse_default
rmse_ivt
rmse_sup

### Plot results

all_clean_filtered %>% 
  filter(base == "TGACT",
         d_SumErr_pos_0 >= 0) %>% 
  ggplot(aes(x = expected_m6A, y = d_SumErr_pos_0,color = model)) + 
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(aes(fill = model),method = "loess") +
  scale_color_manual(values = c('#E64B35FF',"#00A087FF",'#4DBBD5FF')) +
  labs(x = "m6ATP per IVT [%]",
       y = bquote(Delta*~.("SumErr"))) +
  #facet_wrap(~is_DRACH) +
  theme_pubr() + t

ggsave("TGACT_loess.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","m6A_stoichiometry","model_fitting"),
       width = 7,
       height = 5,
       units = "in")


### Point plot with fitted line DRACH vs non

### Default RMSE

default_model_DRACH <- all_clean_filtered %>% 
  filter(model == "default" & is_DRACH == 1, d_SumErr_pos_0 >= 0)

default_model_non_DRACH <- all_clean_filtered %>% 
  filter(model == "default" & is_DRACH == 0, d_SumErr_pos_0 >= 0)

default_model_DRACH <- loess(d_SumErr_pos_0 ~ expected_m6A, data = default_model_DRACH)
default_model_non_DRACH <- loess(d_SumErr_pos_0 ~ expected_m6A, data = default_model_non_DRACH)

rmse_default_DRACH <- sqrt(mean(resid(default_model_DRACH)^2))
rmse_default_non_DRACH <- sqrt(mean(resid(default_model_non_DRACH)^2))

### IVT RMSE

ivt_model_DRACH <- all_clean_filtered %>% 
  filter(model == "ivt" & is_DRACH == 1, d_SumErr_pos_0 >= 0)

ivt_model_non_DRACH <- all_clean_filtered %>% 
  filter(model == "ivt" & is_DRACH == 0, d_SumErr_pos_0 >= 0)

ivt_model_DRACH <- loess(d_SumErr_pos_0 ~ expected_m6A, data = ivt_model_DRACH)
ivt_model_non_DRACH <- loess(d_SumErr_pos_0 ~ expected_m6A, data = ivt_model_non_DRACH)

rmse_ivt_DRACH <- sqrt(mean(resid(ivt_model_DRACH)^2))
rmse_ivt_non_DRACH <- sqrt(mean(resid(ivt_model_non_DRACH)^2))

### SUP RMSE

sup_model_DRACH <- all_clean_filtered %>% 
  filter(model == "sup" & is_DRACH == 1, d_SumErr_pos_0 >= 0)

sup_model_non_DRACH <- all_clean_filtered %>% 
  filter(model == "sup" & is_DRACH == 0, d_SumErr_pos_0 >= 0)

sup_model_DRACH <- loess(d_SumErr_pos_0 ~ expected_m6A, data = sup_model_DRACH)
sup_model_non_DRACH <- loess(d_SumErr_pos_0 ~ expected_m6A, data = sup_model_non_DRACH)

rmse_sup_DRACH <- sqrt(mean(resid(sup_model_DRACH)^2))
rmse_sup_non_DRACH <- sqrt(mean(resid(sup_model_non_DRACH)^2))

### Print to standard out

rmse_default_DRACH
rmse_default_non_DRACH

rmse_ivt_DRACH
rmse_ivt_non_DRACH

rmse_sup_DRACH
rmse_sup_non_DRACH

### Plot

all_clean_filtered %>% 
  filter(d_SumErr_pos_0 >= 0) %>% 
  ggplot(aes(x = expected_m6A, y = d_SumErr_pos_0,color = model)) + 
  #geom_point(size = 2,alpha = 0.4) +
  geom_smooth(aes(fill = model),method = "loess") +
  scale_color_manual(values = c('#E64B35FF',"#00A087FF",'#4DBBD5FF')) +
  labs(x = "m6ATP per IVT [%]",
       y = bquote(Delta*~.("SumErr"))) +
  facet_wrap(~is_DRACH, labeller = as_labeller(c(`1` = "DRACH",`0` = "non-DRACH"))) +
  theme_pubr() + t

ggsave("DRACH_vs_non_DRACH_loess_no_points_long.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","m6A_stoichiometry","model_fitting"),
       #device = "png",
       #dpi = 400,
       width = 12,
       height = 8,
       units = "in")


### Prefilter data for boxplot

all_clean_filtered_bplt <- all_clean_filtered %>% 
  filter(expected_m6A != 0) %>% 
  unite(model_expct, c(model, expected_m6A),sep = "_", remove = F)

### Perform stat-test

# Perform paired t-test and correct for multiple hypothesis testing 

stat.test <- all_clean_filtered_bplt %>% 
  group_by(expected_m6A,is_DRACH) %>% 
  t_test(d_SumErr_pos_0 ~ model, ref.group = "default", alternative = "less") %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test

stat.test_DRACH <- stat.test %>% 
  filter(is_DRACH == 1)

stat.test_non_DRACH <- stat.test %>% 
  filter(is_DRACH == 0)


# add x_y_position

stat.test_DRACH <- stat.test_DRACH %>%
  add_xy_position(x = "model", dodge = 0.5)
stat.test_DRACH$y.position <- stat.test_DRACH$y.position - 0.025

stat.test_non_DRACH <- stat.test_non_DRACH %>%
  add_xy_position(x = "model", dodge = 0.5)
stat.test_non_DRACH$y.position <- stat.test_non_DRACH$y.position - 0.025

## Custom theme bxplt

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "right",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15)
) 


all_clean_filtered_bplt %>% 
  filter(is_DRACH == 1) %>% 
  ggplot(aes(x = as.factor(model), y = d_SumErr_pos_0)) + 
  geom_boxplot(aes(fill = model),
               width=0.25, color="grey20",position = position_dodge(width =0.85),
               notch = F, alpha = 0.9, lwd = 0.7) +
  geom_violin(aes(fill = model), position=position_dodge(width =0.85),
              alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.7) +
  facet_wrap(~expected_m6A, nrow = 1) +
  stat_n_text(size = 5) +
  stat_pvalue_manual(stat.test_DRACH, label = "p.adj.signif", tip.length = .01, size = 5) +
  scale_fill_manual(values = c('#E64B35FF',"#00A087FF",'#4DBBD5FF')) +
  labs(x="",y=expression(Delta*"SumErr")) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15,face = "bold"),
        legend.text = element_text( size = 15))

ggsave("DRACH_per_stoich_wide.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","m6A_stoichiometry","boxplots_pos0_vs_expected"),
       width = 15,
       height = 6,
       units = "in")








