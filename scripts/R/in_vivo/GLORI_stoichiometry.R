####################################################################################
## Plot distribution of GLORI-sites for protein coding genes filterec by coverage ##
####################################################################################

### Task 1 : Load packages ----

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "forcats")

lapply(pkgs, library, character.only = TRUE)

### Task 2: Import detected sites for different model that overlap with GLORI-seq ----

default <- read_tsv(here('data','eligos_2_human','hac','GLORI_overlap_hac.bed'), col_names = F)

ivt <- read_tsv(here('data','eligos_2_human','ivt','GLORI_overlap_ivt.bed'), col_names = F)

### Task 3:  Create a shared column model and merge the two dataframes -----


default$model <- "default"
ivt$model <- "ivt"

merge_df <- rbind(default, ivt) %>% 
  group_by(model) %>% 
  count(X16) %>% 
  arrange(-n)

### Task 4: Visualize the distribution of sites ----

# Set custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 20),
  legend.position = "top",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  axis.text.x = element_text(size = 15),
  strip.text.x = element_text(size = 15)
)  


# Plot sites

merge_df %>% 
  ggplot(aes(x = reorder(X16, +n), fill = model, y=n)) +
  geom_bar(position = position_dodge2(preserve = "single"),stat = "identity") +
  geom_text(
    aes(x = X16, y = n, label = n, group = model), 
    size = 5, hjust = -0.5,
    position = position_dodge2(width = 1),
    inherit.aes = TRUE) +
  labs(x = "",
       y = "number of reported sites") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values=c('default' = '#E64B35FF', 'ivt' = "#00A087FF")) +
  theme_pubr() + t +
  coord_flip()

ggsave(filename = "default_ivt_glori_stoich.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","human_in_vivo", "GLORI_stoichiometry"),
       width = 12,
       height = 4,
       units = "in")
