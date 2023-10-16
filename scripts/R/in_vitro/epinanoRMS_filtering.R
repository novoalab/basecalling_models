##########################
### Data Preprocessing ###
##########################

######################
### Load Packages ####
######################

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr")

lapply(pkgs, library, character.only = TRUE)


##########################
#### 1: Data Import ######
##########################

# Task 1.1: Saving the file directories to 5mers


dir_default <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_hac","5mer"))
dir_ivt <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_ivt_hac","5mer"))
dir_sup <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_sup","5mer"))
dir_ivt_sup <- setwd(here("data","epinanoRMS","rna_r9.4.1_70bps_ivt_sup","5mer"))


##############################################################################
# Task 1.2: save all file names in specified directory as a character vector #
##############################################################################

file_list_default <- list.files(path = dir_default)
file_list_ivt <- list.files(path = dir_ivt)
file_list_sup <- list.files(path = dir_sup)
file_list_ivt_sup <- list.files(path = dir_ivt_sup)

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
ivt_sup_raw <- read_dir(file_list = file_list_ivt_sup, work_dir = dir_ivt_sup)

############################
#### 2: Data Wrangling #####
############################

####################################################################################
# Task 2.1: Function to generate new columns, remove replicates and rename samples #
####################################################################################

filter_data <- function(df){
  
  df_clean <- df %>% 
    
    separate(col = sample, into = c("sample","model"), sep = "-") %>% 
    separate(col = sample, into = c("sample","modification"), sep = "_") %>% # Separate into three columns to contain modifications
    filter(!(sample %in% c("RNAAB089716","RNAAB090763","RNA310520191","RNA234958"))) %>%  # filter out replicates (lower # of overall reads)
    mutate(modification = recode(modification, "pseudoU" = "Y",
                                 "unmodified" = "UNM",
                                 "5hmC" = "hm5C")) %>% 
    mutate(model = recode(model, "ivt" = "ivt+default")) %>% 
    
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
ivt_sup_clean <- filter_data(ivt_sup_raw)

rm(default_raw)
rm(ivt_raw)
rm(sup_raw)
rm(ivt_sup_raw)

all_clean <- rbind(default_clean,ivt_clean,sup_clean, ivt_sup_clean)

###################################################################################################
# Task 2.2: Function to calculate delta SumErr and normalized delta SumErr for center mod of each #
###################################################################################################

filter_central <- function(df){
  
  temp_df <- df
  temp_df <- temp_df %>% 
    # Filter for modified central bases for each modification
    filter((modification == "m6A" & ref_base == "A" |
              modification == "m5C" & ref_base == "C" |
              modification == "hm5C" & ref_base == "C" |
              modification == "ac4C" & ref_base == "C" |
              modification == "Y" & ref_base == "T" |
              modification == "m1Y" & ref_base == "T" |
              modification == "m5U" & ref_base == "T" |
              modification == "UNM")) %>% 
    # group-by X.Ref,base,model,ref_pos
    group_by(X.Ref,base,model,ref_pos) %>% 
    # Calculate delta per postion and on Summed values
    mutate(delta_SumErr = SumErr_5mer - SumErr_5mer[modification == "UNM"],
           delta_SumMis = SumMis_5mer - SumMis_5mer[modification == "UNM"],
           delta_SumIns = SumIns_5mer - SumIns_5mer[modification == "UNM"],
           delta_SumDel = SumDel_5mer - SumDel_5mer[modification == "UNM"],
           delta_SumQ = SumQ_5mer - SumQ_5mer[modification == "UNM"],
           
           # Per position delta
           d_SumErr_pos_minus_2 = SumErr_pos_minus_2 - SumErr_pos_minus_2[modification == "UNM" ],
           d_SumErr_pos_minus_1 = SumErr_pos_minus_1 - SumErr_pos_minus_1[modification == "UNM" ],
           d_SumErr_pos_0 = SumErr_pos_0 - SumErr_pos_0[modification == "UNM" ],
           d_SumErr_pos_plus_1 = SumErr_pos_plus_1 - SumErr_pos_plus_1[modification == "UNM" ],
           d_SumErr_pos_plus_2 = SumErr_pos_plus_2 - SumErr_pos_plus_2[modification == "UNM" ],
           
           d_Mis_pos_minus_2 = as.numeric(mis1) - as.numeric(mis1)[modification == "UNM" ],
           d_Mis_pos_minus_1 = as.numeric(mis2) - as.numeric(mis2)[modification == "UNM" ],
           d_Mis_pos_0 = as.numeric(mis3) - as.numeric(mis3)[modification == "UNM" ],
           d_Mis_pos_plus_1 = as.numeric(mis4) - as.numeric(mis4)[modification == "UNM" ],
           d_Mis_pos_plus_2 = as.numeric(mis5) - as.numeric(mis5)[modification == "UNM" ],
           
           d_Ins_pos_minus_2 = as.numeric(ins1) - as.numeric(ins1)[modification == "UNM" ],
           d_Ins_pos_minus_1 = as.numeric(ins2) - as.numeric(ins2)[modification == "UNM" ],
           d_Ins_pos_0 = as.numeric(ins3) - as.numeric(ins3)[modification == "UNM" ],
           d_Ins_pos_plus_1 = as.numeric(ins4) - as.numeric(ins4)[modification == "UNM" ],
           d_Ins_pos_plus_2 = as.numeric(ins5) - as.numeric(ins5)[modification == "UNM" ],
           
           d_Del_pos_minus_2 = as.numeric(del1) - as.numeric(del1)[modification == "UNM" ],
           d_Del_pos_minus_1 = as.numeric(del2) - as.numeric(del2)[modification == "UNM" ],
           d_Del_pos_0 = as.numeric(del3) - as.numeric(del3)[modification == "UNM" ],
           d_Del_pos_plus_1 = as.numeric(del4) - as.numeric(del4)[modification == "UNM" ],
           d_Del_pos_plus_2 = as.numeric(del5) - as.numeric(del5)[modification == "UNM" ],
           
           d_Q_pos_minus_2 = as.numeric(mq1) - as.numeric(mq1)[modification == "UNM" ],
           d_Q_pos_minus_1 = as.numeric(mq2) - as.numeric(mq2)[modification == "UNM" ],
           d_Q_pos_0 = as.numeric(mq3) - as.numeric(mq3)[modification == "UNM" ],
           d_Q_pos_plus_1 = as.numeric(mq4) - as.numeric(mq4)[modification == "UNM" ],
           d_Q_pos_plus_2 = as.numeric(mq5) - as.numeric(mq5)[modification == "UNM" ],
           
           # Per position Signal to Noise with added pseudocount
           s_to_n_pos_minus_2 = ((SumErr_pos_minus_2 + 0.01) - (SumErr_pos_minus_2[modification == "UNM" ] + 0.01)) / (SumErr_pos_minus_2[modification == "UNM"] + 0.01),
           s_to_n_pos_minus_1 = ((SumErr_pos_minus_1 + 0.01) - (SumErr_pos_minus_1[modification == "UNM" ] + 0.01)) / (SumErr_pos_minus_1[modification == "UNM"] + 0.01),
           s_to_n_pos_0 = ((SumErr_pos_0 + 0.01) - (SumErr_pos_0[modification == "UNM" ] + 0.01)) / (SumErr_pos_0[modification == "UNM"] + 0.01),
           s_to_n_pos_plus_1 = ((SumErr_pos_plus_1 + 0.01) - (SumErr_pos_plus_1[modification == "UNM" ] + 0.01)) / (SumErr_pos_plus_1[modification == "UNM"] + 0.01),
           s_to_n_pos_plus_2 = ((SumErr_pos_plus_2 + 0.01) - (SumErr_pos_plus_2[modification == "UNM" ] + 0.01)) / (SumErr_pos_plus_2[modification == "UNM" ] + 0.01),
    ) %>% 
    dplyr::select(-mis1, -mis2, -mis3, -mis4, -mis5, -ins1, -ins2, -ins3, -ins4, -ins5, -del1, -del2, -del3, -del4, -del5,
                  -mq1, -mq2, -mq3, -mq4, -mq5)
  
}

########################################
# Task 2.3 Subset tables with 4 models #
########################################

all_clean_central_4_models <- all_clean %>% 
  # Remove problematic 5mer that caused uenven number of rows for UNM and modified leading to errors -> check whether this difference already observed in EpinanoRMS ouput
  #filter(!(ref_pos == 24 & base == "TATAG")) %>% 
  filter_central()

# Sanity Check -> are calculations correct
test <- all_clean %>% 
  filter(base == "ATCCG" & ref_pos == 33) 

test1 <- test %>% 
  filter_central()

# Extract motifs

filter_m6A_DRACH <- function(df){
  
  temp_df <- df
  temp_df <- temp_df %>% 
    # Filter for modified central bases for each modification
    filter((modification == "m6A" |
            modification == "UNM")) %>% 
    filter(str_detect(base, "^[AGT][GA][A][C][ACT]$"))
}

m6A_DRACH_only_4_models <- filter_m6A_DRACH(all_clean_central_4_models)

filter_Y_Pus7 <- function(df){
  
  temp_df <- df
  temp_df <- temp_df %>% 
    # Filter for modified central bases for each modification
    filter((modification == "Y" |
              modification == "UNM")) %>% 
    filter(str_detect(base, "^[T][G][T][A][G]$"))
}

Y_Pus7_only_4_models <- filter_Y_Pus7(all_clean_central_4_models)

#########################################
# Task 2.4: Subset tables with 3 models #
#########################################

all_clean_3_models <- all_clean %>% 
  filter(model != 'ivt+sup')

all_clean_central_3_models <- all_clean_3_models %>% 
  filter_central()

m6A_DRACH_only_3_models <- filter_m6A_DRACH(all_clean_central_3_models)

Y_Pus7_only_3_models <- filter_Y_Pus7(all_clean_central_3_models)


##################################################
# Task 2.5: Write Table to processed data folder #
##################################################

######################
### 4 MODEL TABLES ###
######################

# all bases combined df -> 4 models

write.table(all_clean,
            file = here("data","epinanoRMS","data_processed","all_pos_no_filtering_4_models.txt"))

# central bases modified -> all modifications -> 4 models

write.table(all_clean_central_4_models,
            file = here("data","epinanoRMS","data_processed","central_mod_d_per_pos_4_models.txt"))


# m6A modified DRACH motif -> 4 models 

write.table(m6A_DRACH_only_4_models,
            file = here("data","epinanoRMS","data_processed","m6A_DRACH_d_per_pos_4_models.txt"))

# Y Pus7 consensus motif -> 4 models

write.table(Y_Pus7_only_4_models,
            file = here("data","epinanoRMS","data_processed","Y_Pus7_d_per_pos_4_models.txt"))


######################
### 3 MODEL TABLES ###
######################

#central bases modified -> 3 models

write.table(all_clean_3_models,
            file = here("data","epinanoRMS","data_processed","all_pos_no_filtering_3_models.txt"))


# central bases modified -> all modifications -> 4 models

write.table(all_clean_central_3_models,
            file = here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))


# m6A modified DRACH motif -> 4 models 

write.table(m6A_DRACH_only_3_models,
            file = here("data","epinanoRMS","data_processed","m6A_DRACH_d_per_pos_3_models.txt"))

# Y Pus7 consensus motif -> 4 models

write.table(Y_Pus7_only_3_models,
            file = here("data","epinanoRMS","data_processed","Y_Pus7_d_per_pos_3_models.txt"))





