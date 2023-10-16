
### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce","ggdist",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "rstatix")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import data ----

default_mod <- read_csv(here("data","eGFP_epinano","VAC_mU_fast5_s_20k_sup_match.per.site.baseFreq.csv"))

default_mod$model <-  "default"
default_mod$modification <-  "MOD"

sup_mod <- read_csv(here("data","eGFP_epinano","VAC_mU_fast5_s_20k.per.site.baseFreq.csv"))

sup_mod$model <- "sup"
sup_mod$modification <- "MOD"

### UNM ###

default_unmod <- read_csv(here("data","eGFP_epinano","VAC_S4_U_fast5_s_20k_sup_match.per.site.baseFreq.csv"))

default_unmod$model <-  "default"
default_unmod$modification <-  "UNM"

sup_unmod <- read_csv(here("data","eGFP_epinano","VAC_S4_U_fast5_s_20k.per.site.baseFreq.csv"))

sup_unmod$model <- "sup"
sup_unmod$modification <- "UNM"


comb <- rbind(default_mod, sup_mod, default_unmod, sup_unmod)

### Task 3: Calculate U -> C mismatch frequency for pseudoU containing only.

get_c_freq <- function(df){
  
  final_df <- df %>%
    # remove unnecessay columns and rows
    select(`#Ref`,pos,base, model, ACGT, mis, ins, del, model, modification) %>% 
    # correct mismatch frequency to total 1
    mutate(c_mis = mis * (1-(del + ins))) %>% 
    mutate(sum_err = rowSums(cbind(as.numeric(c_mis), as.numeric(ins), as.numeric(del)))) %>% 
    # Split mismatch error into individual positions and only keep central
    separate(ACGT, into = c("A","C","G","T"), ":",convert = T) %>% 
    # Calculate c_freq
    mutate(c_freq = `C`/(`C`+`G`+`A`+`T`) * 100) %>%
    unite("base_model_mod", c("base","model","modification"), remove = F)
    
    return(final_df)
  
}

eGFP_final <- get_c_freq(comb)

eGFP_final$model <- factor(eGFP_final$model, levels = c("sup","default"))
eGFP_final$modification <- factor(eGFP_final$modification, levels = c("UNM","MOD"))

### Boxlpots to determine the background error


### Custom Theme for flipped

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  axis.text.x = element_text(angle = 35, vjust = 0.5),
  panel.grid.major.y = element_line(color = "grey90"),
  panel.grid.minor.y = element_line(color = "grey90"),
  panel.border = element_rect(colour = "black", fill = NA),
  strip.text.x = element_text(size = 15, face = "bold")
)


### Flipped boxplot of individual errors

library(EnvStats)

eGFP_final %>%
  filter(base %in% c("A","C","G", "T"),
         # Filter out poly-A-tail that inflates error on A
         pos < 1061) %>% 
  ggplot(aes(x = base_model_mod, y = sum_err, fill = base)) +
  # Stat Boxxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.5) +
  # Central boxplot
  geom_boxplot(aes(fill = base),
               width=0.75, color="grey20",
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  scale_fill_manual(values=c(A = "#1fab89", `T` = "#eb4d55", C = "#1e56a0", G = "#f0cf85")) +
  labs(x="",y=expression("Summed Error")) +
  theme_pubr() + t +
  stat_n_text() +
  theme(axis.text.x = element_text(size = 15),
        legend.text = element_text( size = 15)) + 
  facet_wrap(~modification)

ggsave(filename = "ins_per_base_error_eGFP_boxplot.pdf",
       plot = last_plot(),
       path = here("~/Library/CloudStorage/Dropbox/Oz-Gregor/4.IVTmod/results_Gregor/plots/eGFP_vaccine/error_direction/"),
       width = 10,
       height = 5,
       units = "in")




