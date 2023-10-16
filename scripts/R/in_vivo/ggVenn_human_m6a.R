##################################################
## Plot overlap between models and ground-truth ##
##################################################

# Load library

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr")

lapply(pkgs, library, character.only = TRUE)

### Task 1: Import data

glori <- read_tsv(here('known_sites','human_m6A', 'filtered','GLORI_min_cov20_3models_protein_coding_only.bed'), col_names = F)

miclip <- read_tsv(here('known_sites','human_m6A', 'filtered','m6A_miCLIP_sites_HEK293T_min_cov20_3models_protein_coding_only.bed'), col_names = F)

ace <- read_tsv(here('known_sites','human_m6A', 'filtered','m6ACE-Seq_min_cov20_3models_protein_coding_only.bed'), col_names = F)

default <- read_tsv(here('data','eligos_2_human','hac', 'hac_wt_vs_ko_merge_baseExt0.A.filtered_single_pos.bed'), col_names = F)
  
ivt <- read_tsv(here('data','eligos_2_human','ivt', 'merge_baseExt0.A.filtered_single_pos.bed'), col_names = F)

sup <- read_tsv(here('data','eligos_2_human','sup', 'merge_baseExt0.A.filtered.bed'), col_names = F)

### Task 3: Generate column based on chr/start/strand to overlap unambigously

glori <- glori %>% 
  mutate(overlap = str_c(X1, '_',X2,'_',X6))

miclip <- miclip %>% 
  mutate(overlap = str_c(X1, '_',X2,'_',X6))

ace <- ace %>% 
  mutate(overlap = str_c(X1, '_',X2,'_',X6))

default <- default %>% 
  mutate(overlap = str_c(X1, '_',X2,'_',X6))

ivt <- ivt %>% 
  mutate(overlap = str_c(X1, '_',X2,'_',X6)) %>% 
  #there is single pos NA that needs to be checked in previous lines of code!
  drop_na()

sup <- sup %>% 
  mutate(overlap = str_c(X1, '_',X2,'_',X6))


### Task 4: Prepare data for plotting

library(ggvenn)

plot_list <- list("default" = default$overlap,
                  "ivt" = ivt$overlap,
                  "sup" = sup$overlap,
                  "GLORI-seq" = glori$overlap,
                  "miCLIP" = miclip$overlap,
                  "m6ACE-seq" = ace$overlap)

### Task 5: Plot all overlaps of interset

### Task 5.1: default + ivt + glori
 
ggvenn(plot_list, c("default","ivt","GLORI-seq"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#E64B35FF', "#00A087FF", "#5B5EA6FF")
       )


ggsave(filename = "default_ivt_glori.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.2: default + ivt + miclip

ggvenn(plot_list, c("default","ivt","miCLIP"),
       show_percentage = F,
       fill_alpha = 0.6,
       stroke_color = NA,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#E64B35FF', "#00A087FF", "#F5CB5CFF")
)


ggsave(filename = "default_ivt_miclip.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.3: default + ivt + sup

ggvenn(plot_list, c("default","ivt","sup"),
       show_percentage = F,
       fill_alpha = 0.6,
       stroke_color = NA,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#E64B35FF', "#00A087FF", "#4DBBD5FF")
)


ggsave(filename = "default_ivt_sup.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.4: glori + miclip

ggvenn(plot_list, c("GLORI-seq","miCLIP"),
       show_percentage = F,
       stroke_color = NA,
       set_name_size = 8,
       fill_alpha = 0.7,
       text_size = 8,
       fill_color = c("#5B5EA6FF", "#F5CB5CFF")
)


ggsave(filename = "glori_miclip.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.4: glori + miclip + ace-seq

ggvenn(plot_list, c("GLORI-seq","miCLIP", "m6ACE-seq"),
       show_percentage = F,
       stroke_color = NA,
       set_name_size = 8,
       fill_alpha = 0.7,
       text_size = 8,
       fill_color = c("#5B5EA6FF", "#F5CB5CFF", "#D5DCF9FF")
)


ggsave(filename = "glori_miclip_ace.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Tasl 5.5: default + glori+ miclip

ggvenn(plot_list, c("default","GLORI-seq","miCLIP"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#E64B35FF', "#5B5EA6FF", "#F5CB5CFF")
)


ggsave(filename = "default_glori_miclip.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.6: ivt + glori+ miclip

ggvenn(plot_list, c("ivt","GLORI-seq","miCLIP"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#00A087FF', "#5B5EA6FF", "#F5CB5CFF")
)


ggsave(filename = "ivt_glori_miclip.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.7: sup + glori+ miclip

ggvenn(plot_list, c("sup","GLORI-seq","miCLIP"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c("#4DBBD5FF", "#5B5EA6FF", "#F5CB5CFF")
)


ggsave(filename = "sup_glori_miclip.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 8,
       height = 8,
       units = "in")

### Task 5.8: default + glori + miclip + ace

ggvenn(plot_list, c("default","GLORI-seq","miCLIP", "m6ACE-seq"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#E64B35FF', "#5B5EA6FF", "#F5CB5CFF", "#D5DCF9FF")
)


ggsave(filename = "default_glori_miclip_ace.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 9,
       height = 9,
       units = "in")

### Task 5.9: ivt + glori + miclip + ace

ggvenn(plot_list, c("ivt","GLORI-seq","miCLIP", "m6ACE-seq"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#00A087FF', "#5B5EA6FF", "#F5CB5CFF", "#D5DCF9FF")
)


ggsave(filename = "ivt_glori_miclip_ace.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 9,
       height = 9,
       units = "in")

### Task 5.9: sup + glori + miclip + ace

ggvenn(plot_list, c("sup","GLORI-seq","miCLIP", "m6ACE-seq"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#4DBBD5FF', "#5B5EA6FF", "#F5CB5CFF", "#D5DCF9FF")
)


ggsave(filename = "sup_glori_miclip_ace.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 9,
       height = 9,
       units = "in")

### Task 5.10: default + ace

ggvenn(plot_list, c("default", "m6ACE-seq"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#E64B35FF', "#D5DCF9FF")
)


ggsave(filename = "default_ace.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 9,
       height = 9,
       units = "in")

### Task 5.11: ivt + ace

ggvenn(plot_list, c("ivt", "m6ACE-seq"),
       show_percentage = F,
       stroke_color = NA,
       fill_alpha = 0.6,
       set_name_size = 8,
       text_size = 8,
       fill_color = c('#00A087FF', "#D5DCF9FF")
)


ggsave(filename = "ivt_ace.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo","venn_diagrams","sites_protein_coding_min_cov20"),
       width = 9,
       height = 9,
       units = "in")

### Task 5.8: Create a plot showing the overlap in sites supported by at least one technology

### Custom Theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "bottom",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  strip.text.x = element_text(size = 15)
)

### Plot 

overlapping_df <- data.frame(model = c("default","default","default","default","ivt","ivt","ivt","ivt"),
                             method = c("default","miCLIP","GLORI-seq","m6ACE-seq","ivt","miCLIP","GLORI-seq","m6ACE-seq"),
                             value = c(66.7, 29.9, 27.6, 55.9, 76.7, 30.0, 28.1, 56.7))

#overlapping_df$method <- factor(overlapping_df$method, levels = c("default","ivt","miCLIP","GLORI-seq"))
# rewrite this using map()

overlapping_df %>% 
  filter(model == "default") %>% 
ggplot(aes(x = reorder(method, -value), y = value, fill = method, label=scales::percent(value, scale = 1))) +
  geom_col(alpha = 0.8) +
  geom_text(
    aes(x = method, y = value),
    vjust = -0.5, size = 5) +
  scale_y_continuous(expand = c(0,0),
                     labels = scales::percent_format(scale = 1),
                     limits = c(0,100)) +
  scale_fill_manual(values = c("default" = "#E64B35FF","ivt" = '#00A087FF', "miCLIP" = "#F5CB5CFF", "GLORI-seq" = "#5B5EA6FF", "m6ACE-seq" = "#D5DCF9FF")) +
  labs(y = "sites supported by \u2265 2 methods",
       x ="") +
  theme_pubr() + t

ggsave(filename = "default_greater_two_v2.pdf",
       plot = last_plot(),
       path = here("results_Gregor", "plots","human_in_vivo", "barplot_greater_two_methods" , "sites_protein_coding_min_cov20"),
       width = 5,
       height = 6,
       units = "in")


### Task 5.9 Modify df to suit UpSetR input

library(UpSetR)

listInput <- list(default = default$overlap,
                  ivt = ivt$overlap,
                  sup = sup$overlap,
                  `miCLIP` = miclip$overlap,
                  `GLORI-seq` = glori$overlap,
                  `m6ACE-seq` = ace$overlap)

upset_plot <- upset(fromList(listInput), order.by = "degree", decreasing = F,
                    scale.sets = "log10",
                    scale.intersections = "log10")

upset_plot

pdf(file=here("results_Gregor", "plots","human_in_vivo","upset_plot", "sites_protein_coding_min_cov20","upset_plot_log.pdf"),
    width = 10,
    height = 5)

print(upset_plot)


dev.off()



# Potential colors for m6ACE-seq 
#7D5BA6FF -> Royal purple
#1F01B9FF -> Dark Blue
#C45BAAFF -> Sky magenta
#D5DCF9FF -> Lavender

################
## Test Zone ##
###############












