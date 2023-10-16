#####################################
### Model accuracy across species ###
#####################################

######################
### Load Packages ####
######################

pkgs <- c("backports","tidyverse","here","skimr","dplyr", "ggplot2", "ggsci","ggforce",
          "janitor","readxl","xlsx", "MetBrewer","ggrepel", "usethis", "ggpubr", "ggridges")

lapply(pkgs, library, character.only = TRUE)


##########################
#### 1: Data Import ######
##########################

# Task 1.1: Saving the file directories to 5mers


dir <- setwd(here("data","bc_accuracy"))


##############################################################################
# Task 1.2: save all file names in specified directory as a character vector #
##############################################################################

file_list <- list.files(path = dir)

################################################################
# Task 1.4: Function for reading in all files from a directory #
################################################################

# Create a loop to read in every file of the directory and append it to the initialized data.frame plus add a new column that contains the name of the pool
# Comment: Could be sped up with fread and data.tables

read_dir <- function(file_list, work_dir){
  
  setwd(work_dir)
  
  dataset <- data.frame()
  
  for (i in 1:length(file_list)){
    temp_data <- read_tsv(file_list[i]) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$species <-  gsub("\\.tsv", "", file_list[i])#clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}


raw <- read_dir(file_list = file_list, work_dir = dir)


###############################################################
# Task 1.5: Ridgeline plot of BC accuracy across model species#
###############################################################

### Set custom theme

t <- theme(
  legend.title = element_blank(),
  legend.text = element_text( size = 15),
  legend.position = "right",
  axis.text = element_text(size=15),
  axis.title=element_text(size=15),
  #panel.grid.major.y = element_line(color = "grey90"),
  #panel.grid.minor.y = element_line(color = "black"),
  #axis.text.y=element_blank(),  #remove y axis labels
  #axis.ticks.y=element_blank()  #remove y axis ticks
)  


processed <- raw %>% 
  unite("species_model", sep = "_", remove = F, c(species,File))

processed$species_model <- factor(processed$species_model,
                                  levels = c("a_thaliana_sup","a_thaliana_ivt","a_thaliana_default",
                                             "yeast_sup","yeast_ivt","yeast_default",
                                             "x_laevis_sup","x_laevis_ivt","x_laevis_default",
                                             "mouse_sup","mouse_ivt","mouse_default",
                                             "human_sup","human_ivt","human_default"))

processed %>% 
  ggplot(aes(x = Identity, y = species_model, color = File, fill = File)) +
  geom_density_ridges(panel_scaling = F, alpha = 0.7) +
  scale_x_continuous(expand = c(0, 0),
                     breaks= seq(0.75, 1, 0.05),
                     limits = c(0.725,1.025)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF")) +
  scale_color_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF")) +
  labs(x = "Identity to reference", y = "") +
  coord_cartesian(clip = "off") +
  theme_ridges() + t


ggsave(filename = "bc_accuracy_cross_species.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","accuracy_plots","transcriptome"),
       width = 7,
       height = 9,
       units = "in")

median_table <- processed %>% 
  group_by(species_model) %>% 
  summarize(median = round(median(Identity), digits = 3))
  

#############################################
### Figure 1B Human accuraccy with values ###
#############################################

processed %>% 
  filter(species == "human") %>% 
  ggplot(aes(x = Identity, color = File)) +
  geom_density(size=1.2) +
  scale_x_continuous(expand = c(0, 0),
                     breaks= seq(0.75, 1, 0.05),
                     limits = c(0.725,1.025)) +
  theme_bw() + t +
  scale_color_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF")) +
  labs(x = "Identity to reference", y = "") +
  coord_cartesian(clip = "off") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  
ggsave(filename = "bc_accuracy_h_sapiens.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","accuracy_plots","transcriptome"),
       width = 8,
       height = 5,
       units = "in")
  
  
  
  
  
  


