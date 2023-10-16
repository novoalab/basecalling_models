
pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix","EnvStats")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import cleaned data ----

central_mod_3models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))

### Task 2.1: Set factor levels for desired order in plots

central_mod_3models$modification <- factor(central_mod_3models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_3models$model <- factor(central_mod_3models$model, levels = c("default","ivt+default","sup"))

### Task 3: Transform df and keep only pos_0

central_mod_3models <- central_mod_3models %>%
  pivot_longer(cols = starts_with("d_SumErr_pos_"),
               names_to = "position",
               values_to = "delta_per_pos") %>% 
  filter(position == "d_SumErr_pos_0")

### Task 4: Calculate pval and correct for multiple hypothesis testing

stat.test <- central_mod_3models %>% 
  group_by(modification) %>% 
  t_test(delta_per_pos ~ model, ref.group = "default") %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test


### Task 5: Plot per pos information and add significance levels ----

### Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "right",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)  

index_5mer <- list(c("m6A"),c("m5C"),c("hm5C"), c("ac4C"),
                   c("Y"),c("m1Y"),c("m5U"))


get_5mer_plots_delta <- function(list_of_comp, df, stats, path){
  
  for (i in 1:length(list_of_comp)){
    
    # Subset for correct modifications
    
    stat.mod <- stats %>% 
      filter(modification == list_of_comp[[i]][1])
    
    # add x_y_pos
    
    stat.mod <- stat.mod %>%
      add_xy_position(x = "model", step.increase = 0.075)
    stat.mod$y.position <- stat.mod$y.position
    
    
     df %>% 
      filter(modification == list_of_comp[[i]][1]) %>% 
      ggplot(aes(x = model, y = delta_per_pos)) +
      geom_boxplot(aes(fill = model),
                   width=0.35, color="grey20",position = position_dodge(width =0.95),
                   notch = T, alpha = 0.9, lwd = 0.4) +
      geom_violin(aes(fill = model), position=position_dodge(width =0.95),
                  alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.85) +
      stat_n_text(size = 5) +
      stat_pvalue_manual(stat.mod,  label = "p.adj.signif", tip.length = .02, size = 5) +
      scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt+default' = "#00A087FF")) +
      labs(x="",y=expression(Delta*"SumErr")) +
      theme_pubr() + t +
      theme(axis.text.x = element_text(size = 15,face = "bold"),
            legend.text = element_text( size = 15))
    
    ggsave(filename = paste0(list_of_comp[[i]][1],"_central_mod_deltaSumErr.pdf"),
           plot = last_plot(),
           path = path,
           width = 5,
           height = 7,
           units = "in")
    
    
  }
}

## For three models all mods

get_5mer_plots_delta(list_of_comp = index_5mer, df = central_mod_3models, stats = stat.test,
                     path = here("results_Gregor","plots","epinanoRMS", "per_kmer", "boxplots_pos0_pval_n","3_models"))








