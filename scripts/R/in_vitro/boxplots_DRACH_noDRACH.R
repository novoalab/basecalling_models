################################################
### Boxplots showing DRACH vs no-DRACH       ###
################################################


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import cleaned data ----

central_mod_4models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_4_models.txt"))

### Task 2.1: Set order for plotting and generate additional df ----

central_mod_4models$modification <- factor(central_mod_4models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_4models$model <- factor(central_mod_4models$model, levels = c("default","sup","ivt+default","ivt+sup"))

### Task 3: Filter for m6A and tag DRACH vs. no-DRACH

# Subset for correct modification, central base and model
m6A_df_4models <- central_mod_4models %>%
  filter(modification == "m6A") %>% 
  # Calculate d_SumErr_5mer by summing across the respective per position rows
  mutate(d_SumErr_5mer = rowSums(select(.,d_SumErr_pos_minus_2:d_SumErr_pos_plus_2)),
         d_SumDel_5mer = rowSums(select(.,d_Del_pos_minus_2:d_Del_pos_plus_2)),
         d_SumIns_5mer = rowSums(select(.,d_Ins_pos_minus_2:d_Ins_pos_plus_2)),
         d_SumMis_5mer = rowSums(select(.,d_Mis_pos_minus_2:d_Mis_pos_plus_2))) %>% 
  # add column that states whehter kmer is DRACH or no
  mutate(is_DRACH = factor(ifelse(str_detect(base, "^[AGT][GA][A][C][ACT]$"), 1, 0), levels = c(1,0))) %>% 
  select(X.Ref,modification,model,base,ref_base,is_DRACH,d_SumErr_5mer, d_SumDel_5mer, d_SumIns_5mer, d_SumMis_5mer, d_SumErr_pos_0)

### Task 3.1: Take only 3 models and reorder
m6A_df_3models <- m6A_df_4models %>%
  filter(model != "ivt+sup")

m6A_df_3models$model <- factor(m6A_df_3models$model, levels = c("default","ivt+default","sup"))

### Task 3.2: Calculate pval

stat.test <- m6A_df_3models %>% 
  group_by(is_DRACH) %>% 
  wilcox_test(d_SumErr_pos_0 ~ model, ref.group = "default") %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test

# add x_y_position

stat.test <- stat.test %>%
  add_xy_position(x = "model", dodge = 0.8)
stat.test$y.position <- stat.test$y.position + 0.05


### Task 4: Set custom theme

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

m6A_df_3models %>% 
  # non-DRACH == 0 ; DRACH ==1
  #filter(is_DRACH == 1) %>% 
  ggplot(aes(x=model, y =d_SumErr_pos_0, fill = model)) +
  geom_boxplot(aes(fill = model),
               width=0.25, color="grey20",position = position_dodge(width =0.85),
               notch = TRUE, alpha = 0.9, lwd = 0.7) +
  geom_violin(aes(fill = model), position=position_dodge(width =0.85),
              alpha = 0.5, show.legend = F, color = NA, scale = "width", width = 0.7) +
  facet_wrap(~is_DRACH, labeller = as_labeller(c(`1` = "DRACH",`0` = "non-DRACH"))) +
  stat_n_text(size = 5) +
  stat_compare_means(comparisons = list(c("default","ivt+default"),c("default","sup")), tip.length = .02, size = 5, label = "p.signif") +
  scale_fill_manual(values = c('#E64B35FF',"#00A087FF",'#4DBBD5FF')) +
  labs(x="",y=expression(Delta*"SumErr")) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15,face = "bold"),
        legend.text = element_text( size = 15))

ggsave(filename = "DRACH_deltaSumErr_POS_0_facetted_3models_long.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","boxplots_DRACH","3_models"),
       width =8,
       height = 6,
       units = "in")



