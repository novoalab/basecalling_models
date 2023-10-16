################################################################################
### determine the distribution of error terms between models for VVUVV kmers ###
################################################################################


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","ggdist",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import cleaned data ----

central_mod_3models <- read.table(here("data","epinanoRMS","data_processed","central_mod_d_per_pos_3_models.txt"))


### Task 2.1: Set order for plotting and generate additional df ----


central_mod_3models$modification <- factor(central_mod_3models$modification, levels = c("UNM","m6A","m5C","hm5C","ac4C","Y","m1Y","m5U"))
central_mod_3models$model <- factor(central_mod_3models$model, levels = c("sup","ivt+default","default"))


### Task 3: Calculate U -> C mismatch frequency for pseudoU containing only.

get_c_freq <- function(df){
  
  final_df <- df %>%
    # remove unnecessay columns and rows
    select(X.Ref,base,modification, model, ACGT, c_mis3, ins, del) %>% 
    filter(modification %in% c("m1Y","UNM")) %>% 
    # Split mismatch error into individual positions and only keep central
    separate(ACGT, into = c("-2","-1","0","+1","+2"), ",") %>% 
    select(-`-2`,-`-1`,-`+1`,-`+2`) %>% 
    # Split position 0 mismatch error into individual bases
    separate(`0`, into = c("A","C","G","T"), ":",convert = T) %>% 
    # Calculate c_freq
    mutate(c_freq = `C`/(`C`+`G`+`A`+`T`) * 100) %>% 
    # Split remaining errors per size
    separate(del, into = c("del_-2","del_-1","del_0","del_+1","del_+2"), ",") %>% 
    select(-c("del_-2","del_-1","del_+1","del_+2")) %>% 
    separate(ins, into = c("ins_-2","ins_-1","ins_0","ins_+1","ins_+2"), ",") %>% 
    select(-c("ins_-2","ins_-1","ins_+1","ins_+2")) %>%
    #rename
    rename(c_mis = c_mis3,
           del = del_0,
           ins = ins_0)
  
  final_df$del <- as.numeric(final_df$del)
  final_df$ins <- as.numeric(final_df$ins)
  
  
  return(final_df)
  
}

cfreq_central_mod_3models <- get_c_freq(central_mod_3models)

### Task 4: Bin data based on c_freq

bin_data <- function(df){
  
  final_df <- df %>% 
    # Bin data
    mutate(bin = case_when(
      df$c_freq >= 75 ~ "HIGH",
      df$c_freq < 75 & df$c_freq >= 50 ~ "MEDIUM-HIGH",
      df$c_freq < 50 & df$c_freq > 25  ~ "MEDIUM-LOW",
      df$c_freq <= 25 ~ "LOW"
    ))
  
  final_df$bin <- factor(final_df$bin, levels = c("LOW", "MEDIUM-LOW","MEDIUM-HIGH","HIGH"))
  
  return(final_df)
}


cfreq_central_mod_3models <- bin_data(cfreq_central_mod_3models)


### Task 5: Annotate Pus7 and TRUB1

annotate_data <- function(df){
  
  final_df <- df %>% 
    mutate(motif = case_when(
      grepl("^[T][G][T][A][G]$", base) ~ "PUS7",
      grepl("^[G][T][T][C][ATGC]$", base) ~ "TRUB1",
      TRUE ~ "other"
    )) %>% 
    mutate(single_U_motif = case_when(
      grepl("^[GCA][GCA][T][GCA][GCA]$", base) ~ "single",
      TRUE ~ "not-single"
    ) )
  
}

cfreq_central_mod_3models <- annotate_data(cfreq_central_mod_3models)

### Visualize the error rate for VVUVVV kmers

plot_df_errortype <- cfreq_central_mod_3models %>% 
  filter(single_U_motif == "single") %>% 
  filter(modification == "m1Y") %>% 
  pivot_longer(cols = c("c_mis","del","ins"),
               names_to = "error_type",
               values_to = "error_value")


### Calculate statistics to add manually:

stat.test <- plot_df_errortype %>%
  group_by(error_type) %>% 
  wilcox_test(error_value ~ model, ref.group = "default") %>% 
  adjust_pvalue(method = "hochberg") %>%
  add_significance("p.adj")
stat.test

### Add x.y position for plotting

stat.test <- stat.test %>%
  add_xy_position(x = "model", step.increase = 0.075)
stat.test$y.position <- stat.test$y.position

stat.cmis <- stat.test %>% filter(error_type == "c_mis")
stat.del <- stat.test %>% filter(error_type == "del")
stat.ins <-  stat.test %>% filter(error_type == "ins")


### Flipped boxplots to show the distribition of errors


### Custom Theme for flipped

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  panel.grid.major.x = element_line(color = "grey90"),
  panel.grid.minor.x = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)


### Flipped boxplot of individual errors

library(ggdist)

plot_df_errortype %>%
  filter(error_type == "ins") %>% 
  ggplot(aes(x = model, y = error_value)) +
  # Stat Boxxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.1) +
  # Central boxplot
  geom_boxplot(aes(fill = model),
               width=0.25, color="grey20",position = position_dodge(width =0.95),
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  # Distribution as half plot
  stat_halfeye(aes(fill = model), normalize = "xy", slab_colour = "black", slab_size = 0.5,
               alpha = 0.75, adjust = .5, width = .3, .width = 0, justification = -.7, point_colour = NA, show.legend = F) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt+default' = "#00A087FF")) +
  labs(x="",y=expression("Mismatch Frequency [%]")) +
  stat_pvalue_manual(stat.cmis,  label = "p.adj.signif", tip.length = .02, size = 5, coord.flip = T) +
  #stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  coord_flip()

ggsave(filename = "ins_error_VVUVV.pdf",
       plot = last_plot(),
       path = here("~/Library/CloudStorage/Dropbox/Oz-Gregor/4.IVTmod/results_Gregor/plots/epinanoRMS/per_kmer/m1Y_pos0_CtoU_mismatch"),
       width = 6,
       height = 4,
       units = "in")

### Get Mismatch error signatures --------

plot_df_mis_error <- plot_df_errortype %>% 
  # Quick workaround to remove triplicates although it does not change outcome
  filter(error_type == "c_mis") %>% 
  group_by(model) %>% 
  summarize(A_count = sum(A),
            C_count = sum(C),
            G_count = sum(G),
            `T_count` = sum(`T`)) %>% 
  rename(A = A_count,
         C= C_count,
         G = G_count,
         `Ref` = `T_count`) %>% 
  pivot_longer(cols = c("A","C","G","Ref"),
               names_to = "base",
               values_to = "cov") %>% 
  group_by(model) %>% 
  mutate(percentage = (cov/sum(cov))) %>% 
  ungroup()

plot_df_mis_error$model <- factor(plot_df_mis_error$model, levels = c("default", "ivt+default", "sup"))

### Mismatch Theme

### Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 20),
  legend.position = "top",
  axis.text = element_text(size=20),
  axis.title=element_text(size=20),
  axis.text.x = element_text(size = 20),
  axis.title.x = element_blank(),
  strip.text.x = element_text(size = 15),
  axis.title.y = element_blank()
)  

### Plot Bars

plot_df_mis_error %>% 
  filter(model == "ivt+default") %>% 
  ggplot(aes(x = model, y = percentage, fill = base)) + 
  geom_bar(position="stack", stat="identity") +
  scale_y_continuous(labels = scales::percent,
                     expand = c(0,0)) +
  geom_text(aes(label=paste0(sprintf("%1.1f", percentage*100))),
            position=position_fill(vjust=0.5), colour="white") +
  scale_fill_manual(values=c("#1fab89", "#1e56a0", "#f0cf85", "#888888")) +
  theme_pubr() + t +
  theme(
    axis.text = element_blank()
  ) +
  coord_flip() 
#coord_cartesian(ylim=c(0.3,1))

ggsave(filename = "VVUVV_kmer_ivt+default.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS", "per_position", "mismatch_dir","VVUVV_kmers_m1Y"),
       width = 8,
       height = 2,
       units = "in")

### Restructure data for ternary plot

plot_df_ternary <- plot_df_errortype %>% 
  filter(error_type == "c_mis") %>% 
  mutate(freq_A = (A / (rowSums(cbind(as.numeric(`A`), as.numeric(`C`), as.numeric(`G`)))))*100,
         freq_C = (C / (rowSums(cbind(as.numeric(`A`), as.numeric(`C`), as.numeric(`G`)))))*100,
         freq_G = (G / (rowSums(cbind(as.numeric(`A`), as.numeric(`C`), as.numeric(`G`)))))*100) %>% 
  drop_na() %>% 
  select(X.Ref:model,starts_with("freq_"),starts_with("error_")) %>% 
  rename(A = freq_A,
         C = freq_C,
         G = freq_G)

plot_df_ternary$model <- factor(plot_df_ternary$model, levels = c("default","ivt+default","sup"))

# install.packages('ggtern')
library(ggtern)
library(rcartocolor)

ggtern(plot_df_ternary, aes(A,G,C)) +
  #geom_mask()+
  geom_point(size = 1,aes(colour= error_value), alpha = 0.6, shape = 16) +
  scale_color_carto_c(palette = "Magenta")+
  labs(color = "Mismatch Frequency") +
  ggtitle("Ternary Diagram for U positions in  Curlcakes") + 
  theme_bw() +
  facet_wrap(~model, nrow=3)

ggsave(filename = "VVUVV_ternary_plot.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","epinanoRMS", "per_position", "mismatch_dir","VVUVV_kmers_m1Y"),
       width = 8,
       height = 10,
       units = "in")







