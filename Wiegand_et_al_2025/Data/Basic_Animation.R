# Load necessary library
library(ggplot2)
library(grid)
library(png)
library(gganimate)
library(ggimage)
library(dplyr)



setwd("C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Local_Files/Wiegand_et_al_2025_Local/Data")

Aggregation_Metric_Exponent_12 <- read.csv("Aggregation_Metric_Exponent_12.csv")
Crowding_index_Num_Neighbors_12 <- read.csv("Crowding_index_Num_Neighbors_12.csv")
dispersal_centers <- read.csv("dispersal_centers.csv")
Time_Series_05 <- read.csv("Time_Series_05.csv")
Spatial_Locations_05 <- read.csv("Spatial_Locations_05.csv")


head(Spatial_Locations_05)


# Find the maximum timestep for each Index
max_timesteps <- Spatial_Locations_05 %>%
  group_by(Index) %>%
  summarize(max_timestep = max(timestep))

# Join the max_timesteps back to the original data
Spatial_Locations_05 <- Spatial_Locations_05 %>%
  left_join(max_timesteps, by = "Index")

# Create the Died column
Spatial_Locations_05 <- Spatial_Locations_05 %>%
  mutate(Died = ifelse(timestep == max_timestep, 1, 0))

# Remove the max_timestep column if not needed
Spatial_Locations_05 <- Spatial_Locations_05 %>%
  select(-max_timestep)


Spatial_Locations_05 <- subset(Spatial_Locations_05,timestep>1)



# Example usage
Tree_WD <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Images"
setwd(Tree_WD)
# Get a list of all files in the working directory
a <- list.files()

# Loop through each file and save it with the new name format
for (i in seq_along(files)) {
  file_name <- files[i]
  new_name <- paste0("Tree_Sp_", i)
  file.copy(file_name, new_name)
}


seed_img_path <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Seed_1.png"

