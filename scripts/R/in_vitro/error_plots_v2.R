#####################################################
##         Summary Error Plots                     ##
#####################################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


### Task 2: Import cleaned and filter data ----

central_mod <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))

### Task 2.1: Set order for plotting and generate additional df ----

central_mod$modification <- factor(central_mod$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod$model <- factor(central_mod$model, levels = c("default","ivt+default","sup"))


#### 3: DataViz ----

### Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 25),
  legend.position = "top",
  
  axis.text = element_text(size=20),
  axis.title=element_text(size=20),
  axis.text.x = element_text(size = 20,face = "bold"),
  
  panel.grid.major.y = element_line(color = "grey90"),

  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)  

### For error plots keep subset UNM also by central base to make comparison fair

index <- list(c("SumErr_5mer","NULL"),c("SumMis_5mer","NULL"),c("SumDel_5mer","NULL"), c("SumIns_5mer","NULL"),
               c("log_SumQ_5mer","c(0,175)"))



get_error_plots_PER_MOD <- function(list_of_comp, df){
  
  for (i in 1:length(list_of_comp)){
  
    df %>% 
      select(modification,model,ref_base, list_of_comp[[i]][1]) %>% 
      ggplot(aes(x = modification, y = !!sym(list_of_comp[[i]][1]), fill = modification)) +
      geom_boxplot(aes(fill = modification),
                   width=0.4, color="grey20",position = position_dodge(width =0.85),
                   notch = TRUE, alpha = 0.9, lwd = 0.4,
                   outlier.size = 0.25
      ) +
      geom_violin(aes(fill = modification), position=position_dodge(width =0.85),
                  alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.7) +
      scale_y_continuous() +
      scale_fill_manual("modification", 
                        values = c('UNM' = 'grey50','m6A' = '#381a61', 'm5C' = '#f9d14a',
                                   'hm5C' = "#ab3329", 'ac4C' = "#88a0dc",'Y' = '#e78429',
                                   'm1Y' = '#ed968c','m5U' = '#7c4b73')) +
      labs(x="",y="") +
      
      # Use strip.text as y-axis label and rename default names
      facet_wrap(~model, ncol = 4) + 
      #strip.position = "left",
      #labeller = as_labeller(c(default = "Quality Score", sup = "", `ivt+default` = "",  `ivt+sup` = ""))) +
      
      # Theme settings such as strip text size 
      theme_pubr() + t +
      theme(axis.text = element_text(size=15),
            axis.title=element_text(size=20),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "none",
            legend.text = element_text( size = 20),
            # Change strip label appearance to make the look like y-axis
            strip.text = element_text(size = 20),
            strip.background = element_blank(),
            strip.placement = "outside",
            # Change panel spacing to prevent overlap
            panel.spacing = unit(1.5, "lines")) 
    
    ggsave(filename = paste0(list_of_comp[[i]][1],"_PER_MOD_for_paper.pdf"),
           plot = last_plot(),
           path = here("results_Gregor","plots","epinanoRMS", "per_kmer", "boxplots_5mer_error_plots","3_models"),
           width = 10,
           height = 3,
           units = "in")
    

  }
}

get_error_plots_PER_MOD(list_of_comp = index, df = central_mod)
  

