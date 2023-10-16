########################################
### Paired analysis on DRACH motifs ###
#######################################


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import cleaned data ----

central_mod_4models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_4_models.txt"))

### Task 2.1: Set order for plotting and generate additional df ----

central_mod_4models$modification <- factor(central_mod_4models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_4models$model <- factor(central_mod_4models$model, levels = c("default","sup","ivt+default","ivt+sup"))

### Task 2.2: Subset for 3 models and change order 

central_mod_3models <- central_mod_4models %>% 
  filter(model != "ivt+sup")

central_mod_3models$model <- factor(central_mod_3models$model, levels = c("default","ivt+default", "sup"))

### Task 3: Function for filtering based on DRACH motif and getting summary stats

annotate_DRACH <- function(df){
  
  df_clean <- df %>% 
    filter(modification == "m6A") %>% 
    # add column that states whehter kmer is DRACH or no
    mutate(is_DRACH = factor(ifelse(str_detect(base, "^[AGT][GA][A][C][ACT]$"), 1, 0), levels = c(1,0))) %>% 
    filter(is_DRACH == 1) %>%
    mutate(base_and_model = str_c(base, "_",model)) %>% 
    select(X.Ref,modification,model,base,ref_base,is_DRACH,delta_SumErr, delta_SumDel, delta_SumIns, delta_SumMis, d_SumErr_pos_0)
  
}

m6A_df_3models <- annotate_DRACH(central_mod_3models)
m6A_df_4models <- annotate_DRACH(central_mod_4models)



annotate_DRACH_single_A <- function(df){
  
  df_clean <- df %>% 
    filter(modification == "m6A") %>% 
    # add column that states whehter kmer is DRACH or no
    mutate(is_DRACH = factor(ifelse(str_detect(base, "^[GT][G][A][C][CT]$"), 1, 0), levels = c(1,0))) %>% 
    mutate(base_and_model = str_c(base, "_",model)) %>% 
    select(X.Ref,modification,model,base,ref_base,is_DRACH,delta_SumErr, delta_SumDel, delta_SumIns, delta_SumMis, d_SumErr_pos_0)
  
}

m6A_df_3mod_singleA <- annotate_DRACH_single_A(central_mod_3models)
m6A_df_4mod_singleA <- annotate_DRACH_single_A(central_mod_4models)





### Task 3: Custom theme

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


########################
### DEFAULT vs. IVT ###
#######################

### Visualize the comparison per 
### Include adjusted pvalue

library(rstatix)
library(EnvStats)
library(ggpubr)

df_comb <-  m6A_df_3mod_singleA %>%  
  filter(model %in% c("default","ivt+default")) %>% 
  # Filter for single A, DRACH containing and select 4 random non-DRACH
  filter(is_DRACH == 1 | base %in% c("CCACG","GGATG","TCACT","TTAGC"))


# Perform paired t-test and correct for multiple hypothesis testing 

stat.test <- df_comb %>% 
  group_by(base) %>% 
  t_test(d_SumErr_pos_0 ~ model, paired = TRUE) %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test

# add x_y_position

stat.test <- stat.test %>%
  add_xy_position(x = "model", dodge = 0.8)
stat.test$y.position <- stat.test$y.position + 0.05

stat.test_plot1 <- stat.test %>% 
  filter(base %in% c("GGACC","GGACT","TGACC","TGACT"))

stat.test_plot2 <- stat.test %>% 
  filter(base %in% c("CCACG","GGATG","TCACT","TTAGC"))


# Arrange plotting order so that row 1 = DRACH, row 2 = non_DRACH and internatlly alphabetic order

#df_comb$base <- factor(df_comb$base, levels = c("GGACC","GGACT","TGACC","TGACT","CCACG","GGATG","TCACT","TTAGC"))

plot1 <- df_comb %>% filter(is_DRACH == 1)
plot2 <- df_comb %>% filter(is_DRACH == 0)



ggpaired(plot1, x = "model", y = "d_SumErr_pos_0",
         color = "model", line.color = "gray",
         width = 0.25, line.size = 1, point.size = 2) +
  stat_pvalue_manual(stat.test_plot1,  label = "p.adj.signif", tip.length = .05, size = 5) +
  scale_color_manual(values = c('#E64B35FF',"#00A087FF")) +
  labs(x="",y=expression(Delta*"SumErr")) +
  facet_wrap(.~base, nrow = 1) +
  stat_n_text(size = 5) +
  # add some padding 
  scale_y_continuous(expand = c(0.09,0.09)) +
  t

ggsave(filename = "d_SumErr_pos_0_pairwise_padj_BH_corr_n_scale_fixed_DRACH_fixed_y_LOL.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","pairwise_comp_DRACH","default_vs_ivt+default"),
       width = 12,
       height = 3,
       units = "in")



df_comb %>% 
  ggplot(aes(x = model , y = d_SumErr_pos_0, color = model)) +
  geom_boxplot() +
  facet_wrap(.~base, nrow = 2, scales = "free_y")


########################
### DEFAULT vs. SUP ###
#######################


df_comb_sup <-  m6A_df_3mod_singleA %>%  
  filter(model %in% c("default","sup")) %>% 
  # Filter for single A, DRACH containing and select 4 random non-DRACH
  filter(is_DRACH == 1 | base %in% c("CCACG","GGATG","TCACT","TTAGC"))


# Perform paired t-test and correct for multiple hypothesis testing 

stat.test <- df_comb_sup %>% 
  group_by(base) %>% 
  t_test(d_SumErr_pos_0 ~ model, paired = TRUE) %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test

# add x_y_position

stat.test <- stat.test %>%
  add_xy_position(x = "model")

# Somehow position was set to three here messing up the plot

stat.test$xmax <- 2
stat.test$y.position <- stat.test$y.position + 0.05


stat.test_plot1 <- stat.test %>% 
  filter(base %in% c("GGACC","GGACT","TGACC","TGACT"))

stat.test_plot2 <- stat.test %>% 
  filter(base %in% c("CCACG","GGATG","TCACT","TTAGC"))


# Arrange plotting order so that row 1 = DRACH, row 2 = non_DRACH and internatlly alphabetic order

#df_comb$base <- factor(df_comb$base, levels = c("GGACC","GGACT","TGACC","TGACT","CCACG","GGATG","TCACT","TTAGC"))

plot1 <- df_comb_sup %>% filter(is_DRACH == 1)
plot2 <- df_comb_sup %>% filter(is_DRACH == 0)



ggpaired(plot2, x = "model", y = "d_SumErr_pos_0",
         color = "model", line.color = "gray",
         width = 0.25, line.size = 1, point.size = 2) +
  stat_pvalue_manual(stat.test_plot2,  label = "p.adj.signif", tip.length = .02, size = 5) +
  scale_color_manual(values = c('#E64B35FF',"#4DBBD5FF")) +
  labs(x="",y=expression(Delta*"SumErr")) +
  facet_wrap(.~base, nrow = 1) +
  stat_n_text(size = 5) +
  # add some padding 
  scale_y_continuous(expand = c(0.09,0.09)) +
  t

ggsave(filename = "d_SumErr_pos_0_pairwise_padj_BH_corr_n_scale_fixed_non_DRACH_fixed_y.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","pairwise_comp_DRACH","default_vs_sup"),
       width = 12,
       height = 3,
       units = "in")



df_comb %>% 
  ggplot(aes(x = model , y = d_SumErr_pos_0, color = model)) +
  geom_boxplot() +
  facet_wrap(.~base, nrow = 2, scales = "free_y")


















