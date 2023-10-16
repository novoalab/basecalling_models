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


dir <- setwd(here("data","readlengths_all_species"))


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
    temp_data <- read_table(file_list[i], col_names = F) #each file will be read in, specify which columns you need read in to avoid any errors # specifying col_types is essential to see spike_ins
    temp_data$species <-  gsub("\\.txt", "", file_list[i]) #clean the data as needed, in this case I am creating a new column that indicates which file each row of data came from
    dataset <- rbind(dataset, temp_data) #for each iteration, bind the new data to the building dataset
  }
  
  rm(i)
  rm(temp_data)
  
  return(dataset)
}


raw <- read_dir(file_list = file_list, work_dir = dir)

processed <- raw %>% 
  separate(species, into = c("model","species_2"), remove = F)


##############################################################
# Task 1.5: Boxplot of read lengths across different species #
##############################################################

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



processed$species_2 <- factor(processed$species_2,
                                  levels = c("human","mouse","xlaevis","yeast","athaliana"))



processed %>% 
  ggplot(aes(x = X2, y = model, fill = model)) +
  # Stat Boxplot to add errorbar optic
  stat_boxplot(geom = "errorbar",
               lwd = 0.4,
               width = 0.3) +
  # Central boxplot
  geom_boxplot(width=0.25, color="grey20",
               alpha = 1, lwd = 0.4, notch = T, outlier.shape = NA) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c('default' = '#E64B35FF', 'sup' = '#4DBBD5FF', 'ivt' = "#00A087FF")) +
  labs(x="read length [nt]",y="") +
  scale_x_log10(expand = c(0, 0)) +
  facet_wrap(~species_2, ncol = 1) +
  theme_bw() + t
  


ggsave(filename = "read_length_species.pdf",
       plot = last_plot(),
       path = here("results_Gregor","plots","basecalling","in_vivo", "3_models"),
       width = 7,
       height = 12,
       units = "in")

median <- processed %>% 
  group_by(species) %>% 
  summarize(median = median(X2))

