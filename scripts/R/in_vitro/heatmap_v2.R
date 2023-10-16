################################################
### Heatmaps for Per Position Error Patterns ###
################################################


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import cleaned data ----

central_mod_4models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_4_models.txt"))

central_m6A_DRACH_4models <- read.table(here("data","epinanoRMS","data_processed","m6A_DRACH_d_per_pos_4_models.txt"))

central_mod_3models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))

central_m6A_DRACH_3models <- read.table(here("data","epinanoRMS","data_processed","m6A_DRACH_d_per_pos_3_models.txt"))
  
### Task 2.1: Set order for plotting and generate additional df ----

central_mod_4models$modification <- factor(central_mod_4models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_4models$model <- factor(central_mod_4models$model, levels = c("default","sup","ivt+default","ivt+sup"))

central_m6A_DRACH_4models$modification <- factor(central_m6A_DRACH_4models$modification, levels = c("UNM","m6A"))
central_m6A_DRACH_4models$model <- factor(central_m6A_DRACH_4models$model, levels = c("default","sup","ivt+default","ivt+sup"))

central_mod_3models$modification <- factor(central_mod_3models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_3models$model <- factor(central_mod_3models$model, levels = c("default","ivt+default","sup"))

central_m6A_DRACH_3models$modification <- factor(central_m6A_DRACH_3models$modification, levels = c("UNM","m6A"))
central_m6A_DRACH_3models$model <- factor(central_m6A_DRACH_3models$model, levels = c("default","ivt+default","sup"))

### Subset dataframe to default only

central_mod_4models_default <- central_mod_4models %>% 
  filter(model == "default")

central_mod_3models_default <- central_mod_3models %>% 
  filter(model == "default")

#### Task 3: Levelplot using ComplexHeatmap ----

# BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(tidyHeatmap)
library(viridis)
library(viridisLite)

list_of_plots  <- list(c("delta_SumErr","d_SumErr","SumErr"),c("delta_Mis","d_Mis","MisFreq"),c("delta_Ins","d_Ins", "InsFreq"), c("delta_Del","d_Del", "DelFreq"),
                       c("delta_Q","d_Q", "QScore"), c("Signal_to_Noise", "s_to_n", "SNR"))

colours <- list('model' = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt+default' = "#00A087FF",  'ivt+sup' = "#3C5488FF"),
                'modification' = c('m6A' = '#381a61', 'm5C' = '#f9d14a', 'hm5C' = "#ab3329", 'ac4C' = "#88a0dc",'Y' = '#e78429', 'm1Y' = '#ed968c','m5U' = '#7c4b73'))

get_heatmap_delta <- function(list_of_plots, df, colours, folder){
  for (i in 1:length(list_of_plots)){
    
    ### Step 1: Generate statistics rename columns and filter
    
    temp_df <- df %>% 
      
      select(model,modification,starts_with(paste0(list_of_plots[[i]][2]))) %>% 
      group_by(model,modification) %>% 
      summarise(med_pos_minus_2 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_minus_2"))),
                med_pos_minus_1 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_minus_1"))),
                med_pos_0 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_0"))),
                med_pos_plus_1 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_plus_1"))),
                med_pos_plus_2 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_plus_2")))) %>% 
      filter(modification != "UNM") %>% 
      # ungroup otherwise select doesn't allow removing this column
      ungroup() %>% 
      rename( "-2" = 3, "-1" = 4 , "0" = 5, "+1" = 6, "+2" = 7) %>% 
      # Finalize order for plotting
      arrange(modification, model) %>% 
      # Create new column with combined information and remove others
      mutate(mod_model = paste0(modification, '_', model), .before = model) 
    
    
    ### Step 2: Extract column information for RowAnnotations
    
    ann <- data.frame(temp_df$model)
    ann$modification <- temp_df$modification
    colnames(ann) <- c('model','modification')
    #colours <- list('model' = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF"),
    #'modification' = c('m6A' = '#381a61', 'm5C' = '#f9d14a', 'hm5C' = "#ab3329", 'ac4C' = "#88a0dc",'Y' = '#e78429', 'm1Y' = '#ed968c','m5U' = '#7c4b73'))
    leftAnn <- HeatmapAnnotation(df = ann,
                                 which = 'row',
                                 col = colours,
                                 annotation_width = unit(c(1, 4), 'cm'),
                                 gap = unit(1, 'mm'),
                                 show_annotation_name = FALSE,
                                 annotation_legend_param = list(
                                   labels_gp = grid::gpar(fontsize = 15),
                                   title_gp = grid::gpar(fontsize = 17)))   
    
    
    ### Step 3: Convert df into matrix with named columns
    
    
    temp_mx <- temp_df %>%  
      select(-model,-modification) %>%  
      # convert column to rownames and finally convert to matrix needed for corrplot
      column_to_rownames(var = "mod_model") %>% 
      data.matrix()
    
    ### Step 4: Plot heatmaps iteratively
    
    ht <- Heatmap(temp_mx, 
                  cluster_rows = F, cluster_columns = F,
                  column_names_rot = 0,
                  col = viridis(100),
                  column_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                  show_row_names = F,
                  left_annotation = leftAnn,
                  row_split = ann$modification,
                  heatmap_legend_param = list(
                    title = bquote(Delta*~.(list_of_plots[[i]][3])),
                    labels_gp = grid::gpar(fontsize = 10),
                    title_gp = grid::gpar(fontsize = 15)))
    
    tidyHeatmap::save_pdf(.heatmap = ht,
                          filename = here("results_Gregor","plots","epinanoRMS","per_kmer","heatmaps_5mer_delta",folder,paste0(list_of_plots[[i]][1],"_central_mod_no_scaling.pdf")),
                          width = 8,
                          height = 6,
                          units = "in")
    
  }
}


#################################
### 4 Models Heatmap unscaled ###
#################################

get_heatmap_delta(list_of_plots = list_of_plots, df = central_mod_4models, colours = colours, folder = "4_models/all_4_models")

get_heatmap_delta(list_of_plots = list_of_plots, df = central_mod_4models_default, colours = colours, folder = "4_models/default_only")

get_heatmap_delta(list_of_plots = list_of_plots, df = central_m6A_DRACH_4models, colours = colours, folder = "4_models/all_4_models_DRACH")

#################################
### 3 Models Heatmap unscaled ###
#################################

get_heatmap_delta(list_of_plots = list_of_plots, df = central_mod_3models, colours = colours, folder = "3_models/all_3_models")

get_heatmap_delta(list_of_plots = list_of_plots, df = central_mod_3models_default, colours = colours, folder = "3_models/default_only")

get_heatmap_delta(list_of_plots = list_of_plots, df = central_m6A_DRACH_3models, colours = colours, folder = "3_models/all_3_models_DRACH")

### Scale by row to show where the error occurs

get_heatmap_delta_scaled <- function(list_of_plots, df, colours, folder){
  for (i in 1:length(list_of_plots)){
    
    ### Step 1: Generate statistics rename columns and filter
    
    temp_df <- df %>% 
      
      select(model,modification,starts_with(paste0(list_of_plots[[i]][2]))) %>% 
      group_by(model,modification) %>% 
      summarise(med_pos_minus_2 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_minus_2"))),
                med_pos_minus_1 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_minus_1"))),
                med_pos_0 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_0"))),
                med_pos_plus_1 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_plus_1"))),
                med_pos_plus_2 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_plus_2")))) %>% 
      filter(modification != "UNM") %>% 
      # ungroup otherwise select doesn't allow removing this column
      ungroup() %>% 
      rename( "-2" = 3, "-1" = 4 , "0" = 5, "+1" = 6, "+2" = 7) %>% 
      # Finalize order for plotting
      arrange(modification, model) %>% 
      # Create new column with combined information and remove others
      mutate(mod_model = paste0(modification, '_', model), .before = model) 
    
    
    ### Step 2: Extract column information for RowAnnotations
    
    ann <- data.frame(temp_df$model)
    ann$modification <- temp_df$modification
    colnames(ann) <- c('model','modification')
    #colours <- list('model' = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF"),
    #'modification' = c('m6A' = '#381a61', 'm5C' = '#f9d14a', 'hm5C' = "#ab3329", 'ac4C' = "#88a0dc",'Y' = '#e78429', 'm1Y' = '#ed968c','m5U' = '#7c4b73'))
    rightAnn <- HeatmapAnnotation(df = ann,
                                  which = 'row',
                                  col = colours,
                                  annotation_width = unit(c(1, 4), 'cm'),
                                  gap = unit(1, 'mm'),
                                  show_annotation_name = FALSE,
                                  annotation_legend_param = list(
                                    labels_gp = grid::gpar(fontsize = 15),
                                    title_gp = grid::gpar(fontsize = 17)))   
    
    
    ### Step 3: Convert df into matrix with named columns
    
    
    temp_mx <- temp_df %>%  
      select(-model,-modification) %>%  
      # convert column to rownames and finally convert to matrix needed for corrplot
      column_to_rownames(var = "mod_model") %>% 
      data.matrix()
    ### This is ugly but works -> only part that differs between the above functions
    temp_mx_t <- t(temp_mx)
    temp_mx_t <- scale(temp_mx_t)
    temp_mx <- t(temp_mx_t)
    
    ### Step 4: Plot heatmaps iteratively
    
    ht <- Heatmap(temp_mx, 
                  cluster_rows = F, cluster_columns = F,
                  column_names_rot = 0,
                  col = viridis(100),
                  column_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                  show_row_names = F,
                  right_annotation = rightAnn,
                  heatmap_legend_param = list(
                    title = bquote(Delta*~.(list_of_plots[[i]][3])),
                    labels_gp = grid::gpar(fontsize = 10),
                    title_gp = grid::gpar(fontsize = 15)))
    
    tidyHeatmap::save_pdf(.heatmap = ht,
                          filename = here("results_Gregor","plots","epinanoRMS","per_kmer","heatmaps_5mer_delta",folder,paste0(list_of_plots[[i]][1],"_central_mod_no_scaling.pdf")),
                          width = 8,
                          height = 6,
                          units = "in")
    
  }
}

################################
### 4 Models Heatmap scaled ###
###############################

get_heatmap_delta_scaled(list_of_plots = list_of_plots, df = central_mod_4models, colours = colours, folder = "4_models/scaled_by_row_all_4_models")

get_heatmap_delta_scaled(list_of_plots = list_of_plots, df = central_mod_4models_default, colours = colours, folder = "4_models/scaled_by_row_default_only")

###############################
### 3 Models Heatmap scaled ###
###############################

get_heatmap_delta_scaled(list_of_plots = list_of_plots, df = central_mod_3models, colours = colours, folder = "3_models/scaled_by_row_all_3_models")

get_heatmap_delta_scaled(list_of_plots = list_of_plots, df = central_mod_3models_default, colours = colours, folder = "3_models/scaled_by_row_default_only")



### Scale by modification to show which model produces the strongest error term


get_heatmap_delta_scaled_mod <- function(list_of_plots, df, colours, folder){
  
  for (i in 1:length(list_of_plots)){
    
    ### Step 1: Generate statistics rename columns and filter
    
    temp_df <- df %>% 
      
      select(model,modification,starts_with(paste0(list_of_plots[[i]][2]))) %>% 
      group_by(model,modification) %>% 
      summarise(med_pos_minus_2 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_minus_2"))),
                med_pos_minus_1 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_minus_1"))),
                med_pos_0 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_0"))),
                med_pos_plus_1 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_plus_1"))),
                med_pos_plus_2 = median(!!sym(paste0(list_of_plots[[i]][2],"_pos_plus_2")))) %>% 
      filter(modification != "UNM") %>% 
      # ungroup otherwise select doesn't allow removing this column
      ungroup() %>% 
      rename( "-2" = 3, "-1" = 4 , "0" = 5, "+1" = 6, "+2" = 7) %>% 
      # Finalize order for plotting
      arrange(modification, model) %>% 
      # Create new column with combined information and remove others
      mutate(mod_model = paste0(modification, '_', model), .before = model) 
    
    
    ### Step 2: Extract column information for RowAnnotations
    
    ann <- data.frame(temp_df$model)
    ann$modification <- temp_df$modification
    colnames(ann) <- c('model','modification')
    #colours <- list('model' = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF"),
    #'modification' = c('m6A' = '#381a61', 'm5C' = '#f9d14a', 'hm5C' = "#ab3329", 'ac4C' = "#88a0dc",'Y' = '#e78429', 'm1Y' = '#ed968c','m5U' = '#7c4b73'))
    leftAnn <- HeatmapAnnotation(df = ann,
                                 which = 'row',
                                 col = colours,
                                 annotation_width = unit(c(1, 4), 'cm'),
                                 gap = unit(1, 'mm'),
                                 show_annotation_name = FALSE,
                                 annotation_legend_param = list(
                                   labels_gp = grid::gpar(fontsize = 15),
                                   title_gp = grid::gpar(fontsize = 17)))   
    
    
    ### Step 3: Convert df into matrix with named columns
    
    ### This part is ugly 
    temp_mx <- temp_df %>%  
      # pivot_longer to be able to apply scaling on modification groups
      pivot_longer(cols = `-2`:`+2`,
                   names_to = "position",
                   values_to = "sumerr") %>% 
      #group and scale
      group_by(modification) %>% 
      #vector conversion required otherwise df returned
      mutate(sumerr_scaled = as.vector(scale(sumerr))) %>% 
      ungroup() %>% 
      select(-sumerr) %>% 
      #repivot to wide
      pivot_wider(names_from = "position",
                  values_from = "sumerr_scaled") %>% 
      select(-model,-modification) %>%  
      # convert column to rownames and finally convert to matrix needed for corrplot
      column_to_rownames(var = "mod_model") %>% 
      data.matrix()
    
    ### Step 4: Plot heatmaps iteratively
    
    ht <- Heatmap(temp_mx, 
                  cluster_rows = F, cluster_columns = F,
                  column_names_rot = 0,
                  col = viridis(100),
                  column_names_gp = grid::gpar(fontsize = 20, fontface = "bold"),
                  show_row_names = F,
                  left_annotation = leftAnn,
                  row_split = ann$modification,
                  heatmap_legend_param = list(
                    title = "z-score",
                    at = c(-4, 0, +4),
                    labels_gp = grid::gpar(fontsize = 10),
                    title_gp = grid::gpar(fontsize = 15)))
    
    tidyHeatmap::save_pdf(.heatmap = ht,
                          filename = here("results_Gregor","plots","epinanoRMS","per_kmer","heatmaps_5mer_delta",folder,paste0(list_of_plots[[i]][1],"_scaled_by_mod.pdf")),
                          width = 8,
                          height = 6,
                          units = "in")
    
  }
}


###################################
### 3 Models Heatmap mod scaled ###
###################################

get_heatmap_delta_scaled_mod(list_of_plots = list_of_plots, df = central_mod_3models, colours = colours, folder = "3_models/test")

