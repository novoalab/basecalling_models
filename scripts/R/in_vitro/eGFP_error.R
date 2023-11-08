############################################################################
### m1Y based error characterisation of synthetic eGFP-based RNA vaccine ###
############################################################################

#Preprocessing

# mop_preprocess
# bam2select.py INBAM OUTBAM 20000 (with idensity conda_env)
# reads were matched between sup and default


### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","ggdist",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import data ----

default <- read_csv(here("data","eGFP_epinano","VAC_mU_fast5_s_20k_sup_match.per.site.baseFreq.csv"))

default$model <-  "default"

sup <- read_csv(here("data","eGFP_epinano","VAC_mU_fast5_s_20k.per.site.baseFreq.csv"))

sup$model <- "sup"

comb <- rbind(default, sup)

### Task 3: Calculate U -> C mismatch frequency for pseudoU containing only.

get_c_freq <- function(df){
  
  final_df <- df %>%
    # remove unnecessay columns and rows
    select(`#Ref`,pos,base, model, ACGT, mis, ins, del, model) %>% 
    # correct mismatch frequency to total 1
    mutate(c_mis = mis * (1-(del + ins))) %>%
    # Split mismatch error into individual positions and only keep central
    separate(ACGT, into = c("A","C","G","T"), ":",convert = T) %>% 
    # Calculate c_freq
    mutate(c_freq = `C`/(`C`+`G`+`A`+`T`) * 100) %>% 
  
  return(final_df)
  
}

eGFP_final <- get_c_freq(comb)

### Visualize the error rate for VVUVVV kmers

eGFP_final_error <- eGFP_final %>% 
  filter(base == "T") %>% 
  pivot_longer(cols = c("c_mis","del","ins"),
               names_to = "error_type",
               values_to = "error_value")

eGFP_final_error$model <- factor(eGFP_final_error$model, levels = c("sup","default"))

### Calculate statistics to add manually:

stat.test <- eGFP_final_error %>%
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

eGFP_final_error %>%
  filter(error_type == "del") %>% 
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
  stat_pvalue_manual(stat.del,  label = "p.adj.signif", tip.length = .02, size = 5, coord.flip = T) +
  #stat_n_text(size = 5) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  coord_flip()

ggsave(filename = "del_error_eGFP.pdf",
       plot = last_plot(),
       path = here("~/Library/CloudStorage/Dropbox/Oz-Gregor/4.IVTmod/results_Gregor/plots/eGFP_vaccine/boxplot_halfeye/"),
       width = 6,
       height = 5,
       units = "in")

### Showcase where the error occurs using a overlapping barplot

eGFP_final %>% 
  filter(base == "T") %>% 
  ggplot(aes(x = pos, y = c_freq, fill = model)) +
  geom_col(position = "identity", alpha = 0.45, width = 2.5) +
  scale_y_continuous(expand = c(0,0),
                     limits = c(0,100)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF')) +
  theme_pubr() + t +
  xlab("position [nt]") + ylab("U:C mismatch frequency [%]")+
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey90")
        )



ggsave(filename = "barplot_per_position.pdf",
       plot = last_plot(),
       path = here("~/Library/CloudStorage/Dropbox/Oz-Gregor/4.IVTmod/results_Gregor/plots/eGFP_vaccine/overlapping_barplot/"),
       width = 15,
       height = 5,
       units = "in")


### Global U:C mismatch freq

# Calculate median values for each category

median_values <- eGFP_final %>%
  filter(base == "T") %>% 
  group_by(model) %>%
  summarise(median_c_freq = median(c_freq))


eGFP_final %>% 
  filter(base == "T") %>% 
  ggplot(aes(x = c_freq)) +
  geom_density(aes(fill = model,
                   color = model), alpha = 0.6) +
  labs(x=expression("U:C mismatch frequency [%]"),y="") +
  scale_y_continuous(expand = c(0.02,0)) +
  scale_x_continuous(limits = c(0,100)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt+default' = "#00A087FF")) +
  theme_pubr() + t +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15),
        panel.grid.major.x = element_line(color = "grey90"),
        panel.grid.minor.x = element_line(color = "grey90"),
        panel.border = element_blank()) +
  geom_vline(data = median_values, aes(xintercept = median_c_freq, color = model),
             linetype = "dashed", show.legend = F) +
  geom_text_repel(data = median_values, aes(x = median_c_freq, y = 0.025, label = round(median_c_freq,2),
                                            color = model), box.padding = 0.5, min.segment.length = 0.1,show.legend = F) 


ggsave(filename = "barplot_per_position.pdf",
       plot = last_plot(),
       path = here("~/Library/CloudStorage/Dropbox/Oz-Gregor/4.IVTmod/results_Gregor/plots/eGFP_vaccine/density_plot/"),
       width = 7,
       height = 5,
       units = "in")

#######################################################
#### Compare 20 lowest U:C between default and sup ####
######################################################


sup_sites <- eGFP_final %>% 
  filter(base == "T", model == "sup") %>%
  select(pos, c_freq, model)
  
default_sites <- eGFP_final %>% 
  filter(base == "T", model == "default") %>%
  select(pos, c_freq, model)

final_sites <- sup_sites %>% 
  full_join(default_sites, by = "pos",
            suffix = c("_sup", "_default")) %>% 
  top_n(-20, c_freq_default) %>% 
  pivot_longer(cols = starts_with("c_freq"),
               names_to = "model",
               values_to = "c_freq")

final_sites$model <- factor(final_sites$model, levels = c("c_freq_default","c_freq_sup"))

summary_stats <- final_sites %>% 
  group_by(model) %>% 
  summarize(median = median(c_freq))

comparison <- list(c("c_freq_default","c_freq_sup"))

ggpaired(final_sites, x = "model", y = "c_freq", id = "pos",
         color = "model", line.color = "gray",
         width = 0.25, line.size = 1, point.size = 2) +
  scale_color_manual(values = c('#E64B35FF',"#4DBBD5FF")) +
  stat_compare_means(method = "wilcox.test", paired = T, label = "p.signif",
                     comparisons = comparison) +
  labs(x="",y="U:C mismatch frequency [%]") +
  theme_pubr() + t +
  theme(panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90"),
        panel.grid.minor.y = element_line(color = "grey90")
  )

ggsave(filename = "paired_low_20_position.pdf",
       plot = last_plot(),
       path = here("~/Library/CloudStorage/Dropbox/Oz-Gregor/4.IVTmod/results_Gregor/plots/eGFP_vaccine/paired/"),
       width = 6,
       height = 5,
       units = "in")


  











