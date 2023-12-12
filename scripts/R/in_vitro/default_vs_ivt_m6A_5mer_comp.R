#################################################
### delta/delta SumErr to compare ivt-default ###
#################################################


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import cleaned data ---- 

central_mod_3models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))

### Task 2.1: Set order for plotting and generate additional df ----

central_mod_3models$modification <- factor(central_mod_3models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_3models$model <- factor(central_mod_3models$model, levels = c("default","ivt+default", "sup"))

### Task 3: Annotate DRACH motifs,calculate delta-delta SumErr 


filter_fun <- function(df){
  
  temp_df <- df
  temp_df <- df %>% 
    # Filter for modified central bases for each modification
    filter(modification == "m6A" &
           model %in% c("default","ivt+default")) %>% 
    # group-by X.Ref,base,model,ref_pos
    group_by(X.Ref,base,ref_pos) %>% 
    # Calculate double-delta per position and generate per 5mer value
    mutate(dd_SumErr_pos_minus_2 = d_SumErr_pos_minus_2 - d_SumErr_pos_minus_2[model == "default" ],
           dd_SumErr_pos_minus_1 = d_SumErr_pos_minus_1 - d_SumErr_pos_minus_1[model == "default" ],
           dd_SumErr_pos_0 = d_SumErr_pos_0 - d_SumErr_pos_0[model == "default" ],
           dd_SumErr_pos_plus_1 = d_SumErr_pos_plus_1 - d_SumErr_pos_plus_1[model == "default" ],
           dd_SumErr_pos_plus_2 = d_SumErr_pos_plus_2 - d_SumErr_pos_plus_2[model == "default" ]) %>% 
    mutate(dd_SumErr = rowSums(cbind(as.numeric(dd_SumErr_pos_minus_2), as.numeric(dd_SumErr_pos_minus_1), as.numeric(dd_SumErr_pos_0), as.numeric(dd_SumErr_pos_plus_1), as.numeric(dd_SumErr_pos_plus_2)))) %>%
    # Filer for only ivt+default as this contains relevant dd values also reduced runtime on case_when
    filter(model == "ivt+default") %>% 
    # Annotate DRACH_motifs
      mutate(is_DRACH = case_when(
        # Non-DRACH Single A
        str_detect(base, "^[GT][G][A][C][CT]$") == F ~ "other",
        # Single A DRACH
        str_detect(base, "^[GT][G][A][C][CT]$") == T ~ "DRACH",
      )) %>% 
  # Remove everything with more than one A
  filter(str_detect(base, "^[GTC][GTC][A][GTC][GTC]$"))
  
}


plot_frame_dd <- filter_fun(central_mod_3models)

### Task 4: Calculate test statistics based on m6A kmers

calcualte_pval <- function(df){
  
  temp_df <- df
  temp_df <- df %>% 
    # Filter for modified central bases for each modification
    filter(modification == "m6A" &
             model %in% c("default","ivt+default")) %>% 
    # Remove everything with more than one A
    filter(str_detect(base, "^[GTC][GTC][A][GTC][GTC]$")) %>% 
    # group-by base
    group_by(base) %>% 
    # Calculate adjusted pval for the paired data
    wilcox_test(d_SumErr_pos_0 ~ model, paired = TRUE, alternative = "less") %>% 
    adjust_pvalue(method = "hochberg") %>%
    add_significance("p.adj") %>% 
    mutate(p.adj_rev = 1-p.adj) %>% 
    # Annotate DRACH_motifs
    mutate(is_DRACH = case_when(
      # Non-DRACH Single A
      str_detect(base, "^[GT][G][A][C][CT]$") == F ~ "other",
      # Single A DRACH
      str_detect(base, "^[GT][G][A][C][CT]$") == T ~ "DRACH",
    ))
  
}


plot_frame_pval <- calcualte_pval(central_mod_3models)

### Task 5: Custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text(size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=10),
  axis.title=element_text(size=20),
  axis.text.y = element_text(size = 15),
  axis.text.x = element_text(family = "Courier", size = 8, angle = 35, hjust = 1),
  strip.text.x = element_text(size = 15, face = "bold")
)  

### Task 6: Plot ddSumErr for default vs. ivt+default

ggbarplot(plot_frame_dd, x = "base", y = "dd_SumErr_pos_0", fill = "is_DRACH",
          add = c("mean_se","point"), 
          sort.val = "asc",sort.by.groups = F) + 
  labs(x="",y=expression(Delta*Delta*"SumErr")) +
  scale_fill_manual(values = c("DRACH" = "#FC766AFF", "other" = "#5B84B1FF")) +
  t + coord_flip()


ggsave(filename = "ddSumErr_default_vs_ivt_singleA.pdf",
         plot = last_plot(),
         path = here("results_Gregor","plots","epinanoRMS","per_kmer","m6A_default_vs_ivt","delta_delta_barplot"),
         width = 10,
         height = 8,
         units = "in")

### Task 7: 

ggbarplot(plot_frame_pval, x = "base", y = "p.adj_rev", fill = "is_DRACH",
          sort.val = "desc",sort.by.groups = F, width = 0.75) + 
  labs(x="",y=expression(1-"p.adj")) +
  scale_fill_manual(values = c("DRACH" = "#E4CE81", "other" = "#B1B3B3")) +
  scale_y_continuous(expand = c(0.005,0.005)) +
  #pval cutoff refernce line
  geom_hline(yintercept = 0.95, color = "grey60", linetype = "dashed") +
  t


ggsave(filename = "default_vs_ivt_singleA_padj_upright_one_sided_wilcox.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS","per_kmer","m6A_default_vs_ivt","pval_plot"),
       width = 18,
       height = 6,
       units = "in")


  
