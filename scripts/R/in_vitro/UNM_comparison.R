#####################################################
##         UNM Error Signature                     ##
#####################################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","rstatix",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)


### Task 2: Import cleaned and filter data ----

central_mod <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))

### Task 2.1: Set order for plotting and generate additional df ----

central_mod$modification <- factor(central_mod$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod$model <- factor(central_mod$model, levels = c("default","ivt+default","sup"))

### Task 2.2: Perform statistical test ----

stat.test <- central_mod %>% 
  filter(modification == "UNM") %>% 
  wilcox_test(SumErr_5mer ~ model, ref.group = "default") %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test

# add x_y_position

stat.test <- stat.test %>%
  add_xy_position(x = "model", dodge = 0.8)
stat.test$y.position <- stat.test$y.position + 0.05

### Task 3: Create theme ----

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

### Task 5: Generate boxplots

library(EnvStats)
library(ggrepel)

central_mod %>% 
  filter(modification == "UNM") %>% 
  ggplot(aes(x=model, y =SumErr_5mer, fill = model)) +
  geom_boxplot(aes(fill = model),
               width=0.35, color="grey20",position = position_dodge(width =0.95),
               notch = T, alpha = 0.9, lwd = 0.4) +
  geom_violin(aes(fill = model), position=position_dodge(width =0.95),
              alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.85) +
  stat_n_text(size = 5) +
  stat_compare_means(comparisons = list(c("default","ivt+default"),c("default","sup")), tip.length = .02, size = 5, label = "p.signif") +
  scale_fill_manual(values = c('#E64B35FF',"#00A087FF",'#4DBBD5FF')) +
  labs(x="",y=expression("SumErr")) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15,face = "bold"),
        legend.text = element_text( size = 15))

ggsave(filename = "UNM_SumErr-3models.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","UNM_SumErr_stats"),
       width =5,
       height = 7,
       units = "in")

central_mod %>% 
  filter(modification == "UNM") %>% 
  group_by(model) %>% 
  summarize(median = median(SumErr_5mer))
  
  
  
