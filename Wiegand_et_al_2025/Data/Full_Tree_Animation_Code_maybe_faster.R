# Load necessary libraries
library(ggplot2)
library(grid)
library(png)
library(gganimate)
library(ggimage)
library(magick)
library(dplyr)

add_image_BAD <- function(plot, img_path, x, y, width, height) {
  img <- rasterGrob(png::readPNG(img_path), width = unit(width, "npc"), height = unit(height, "npc"))
  plot + annotation_custom(img, xmin = x , xmax = x , ymin = y , ymax = y )
}

add_image <- function(plot, img, x, y, width, height) {
  img_grob <- rasterGrob(img, width = unit(width, "npc"), height = unit(height, "npc"))
  plot + annotation_custom(img_grob, xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.25, ymax = y + 0.75)
}

add_image_b <- function(plot, img_path, x, y, width, height) {
  img <- rasterGrob(png::readPNG(img_path), width = unit(width, "npc"), height = unit(height, "npc"))
  plot + annotation_custom(img, xmin = x - 0.03, xmax = x + 0.03, ymin = y - 0.03, ymax = y + 0.03)
}

# Will take some pngs and combine them to gif (used for the X's for mortality)
create_gif_from_plots <- function(plot_list, gif_name = "combined_plots.gif", fps = 5, width = 800*3, height = 900*3) {
  # Create a temporary directory to save the plots
  temp_dir <- tempdir()
  
  # Save each plot in the list as a PNG file with specified dimensions
  plot_files <- lapply(seq_along(plot_list), function(i) {
    file_path <- file.path(temp_dir, paste0("plot", i, ".png"))
    ggsave(file_path, plot = plot_list[[i]] + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
           width = width / 100, height = height / 100, dpi = 100)
    file_path
  })
  
  # Read the images and combine them into a GIF
  images <- lapply(plot_files, image_read)
  
  
  # Repeat the last frame 10 times
  last_frame <- images[[length(images)]]
  for (i in 1:6) {
    images <- append(images, list(last_frame))
  }
  
  gif <- image_animate(image_join(images), fps = fps)
  
  # Save the GIF
  image_write(gif, gif_name)
  
  # Delete the temporary PNG files
  unlink(plot_files)
  
  # Return the path to the GIF
  #return(gif_name)
}


create_gif_from_plots <- function(plot_list, gif_name = "combined_plots.gif", fps = 5, width = 800*3, height = 900*3) {
  # Create a temporary directory to save the plots
  temp_dir <- tempdir()
  
  # Save each plot in the list as a PNG file with specified dimensions
  plot_files <- lapply(seq_along(plot_list), function(i) {
    file_path <- file.path(temp_dir, paste0("plot", i, ".png"))
    ggsave(file_path, plot = plot_list[[i]] + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")), 
           width = width / 300, height = height / 300, dpi = 300)
    file_path
  })
  
  # Read the images and combine them into a GIF
  images <- lapply(plot_files, image_read)
  
  # Repeat the last frame 10 times
  last_frame <- images[[length(images)]]
  for (i in 1:10) {
    images <- append(images, list(last_frame))
  }
  
  gif <- image_animate(image_join(images), fps = fps)
  
  # Save the GIF
  image_write(gif, gif_name)
  
  # Delete the temporary PNG files
  unlink(plot_files)
}


# Function will prompt creation of all animations / save them to specified working directory  
Animate_All <- function(Spatial_Locations_Data,TimeSteps,Frames_Seeds,Tree_Images,Tree_Images_list,Seed_Image,X_Image_Path,alphabetical_vector,WD_Save,Scale_Factor,Mx,My,DELTA){
  
  Spatial_Locations_Data$y <- Spatial_Locations_Data$y + DELTA
  Spatial_Locations_Data$y_parent <- Spatial_Locations_Data$y_parent + DELTA
  
  for(tt in 2:TimeSteps){
    
    Spatial_Locations_TS <- subset(Spatial_Locations_Data,timestep==tt)
    Spatial_Locations_TS <- Spatial_Locations_TS %>% arrange(desc(y))
    Species <- unique(sort(Spatial_Locations_TS$Species))
    
    Living_Adults_TS <- subset(Spatial_Locations_TS,Birth != tt)
    X_Y_Living       <- Living_Adults_TS %>% select(x, y, species)
    
    New_Adults       <- subset(Spatial_Locations_TS,Birth == tt)
    X_Y_New          <- New_Adults %>% select(x, y, species)
    X_Y_Parent       <- New_Adults %>% select(x_parent, y_parent)
    X_Y_Parent$y     <- X_Y_Parent$y 
    
    Deaths           <- subset(Spatial_Locations_TS,Died==1)
    X_Y_Deaths       <- Deaths %>% select(x, y, species)
    X_Y_Deaths$y     <- X_Y_Deaths$y 
    
    # make data for births 
    
    animation_data_list <- list()
    # Create a data frame for the animation
    
    for(jj in 1:nrow(New_Adults)){
      
      x_start <- X_Y_Parent$x[jj]
      y_start <- X_Y_Parent$y[jj]
      
      x_end   <- X_Y_New$x[jj]
      y_end   <- X_Y_New$y[jj]
      
      Total_Move <- (sqrt((x_end-x_start)^2 +(y_end-y_start)^2 ))/Scale_Factor
      
      Movement_Frames <- round(Frames_Seeds*(Total_Move/Frames_Seeds))  
      End_Frames      <- Frames_Seeds - Movement_Frames
      
      Total_Frames <- Movement_Frames+End_Frames
      
      # paths of seed
      x_path <- seq(x_start, x_end, length.out = Movement_Frames)
      y_path <- seq(y_start, y_end, length.out = Movement_Frames)
      
      # make dataframe of seed paths 
      animation_data_list[[jj]] <- data.frame(
        frame = 1:(Total_Frames),
        x   = c(x_path, rep(x_end, End_Frames)),
        y   = c(y_path, rep(y_end, End_Frames)),
        img = c(rep(Seed_Image, Total_Frames))
      )
      
    }
    
    # implement birth animations 
    Births_anim <- create_animation(X_Y_Living, animation_data_list, Tree_Images,Tree_Images_list,Mx,My)
    
    # Name births animation 
    Name_animation_Births <- paste0("TimeStep_",alphabetical_vector[tt],"_Births",".gif")
    
    # set working directory to save and save
    setwd(WD_Save)
    anim_save(Name_animation_Births, Births_anim)
    
    # take last frame of animation 
    last_frame <- image_read(Births_anim)[length(image_read(Births_anim))]
    last_frame_plot <- ggplot() +
      annotation_custom(rasterGrob(last_frame, width = unit(1, "npc"), height = unit(1, "npc")))
    
    # Create a new plot on top of the last frame
    new_plot <- last_frame_plot
    
    # Add the image at each coordinate
    List_Plots <- list()
    List_Plots[[1]] <- new_plot
    
    for(dd in 1:(nrow(X_Y_Deaths))) {
      List_Plots[[dd+1]] <- add_image_b(List_Plots[[dd]],X_Image_Path,X_Y_Deaths$x[dd]/(Mx), (X_Y_Deaths$y[dd])/(My),.5,.5) 
    }
    
    # Name deaths animation 
    Name_animation_Deaths <- paste0("TimeStep_",alphabetical_vector[tt],"_Deaths",".gif")
    
    # make deaths animation 
    create_gif_from_plots(List_Plots,Name_animation_Deaths,2,800,900)
    
  }
  
}

# Function to create the animation
create_animation <- function(Adult_Positions, animation_data, Tree_Images,Tree_Images_list,Mx,My) {
  # Create a data frame for the grid
  grid_data <- expand.grid(x = 0:Mx, y = 0:My)
  
  # Create the texture
  set.seed(48354)
  texture <- rasterGrob(matrix(runif(250000), nrow = 500, ncol = 500), 
                        width = unit(1, "npc"), height = unit(1, "npc"), 
                        interpolate = TRUE)
  
  # Create the base animation plot
  anim_plot <- ggplot(grid_data, aes(x, y)) +
    annotation_custom(texture, xmin = 0, xmax = Mx, ymin = 0, ymax = My) + # change these boundaries to make raster match 
    geom_tile(fill = "#034700", color = NA, alpha = .9, height = 1.0065, width = 1.0065, linetype = 1) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "grey20", color = "grey20"),
      plot.margin = unit(c(-.5/3, -1/3, -.5/3, -1/3), "cm"),
      panel.spacing = unit(c(0, 0, 0, 0), "cm"))
  
  # Loop through each dataframe in the list and add layers to the plot
  for (jj in seq_along(animation_data)) {
    anim_plot <- anim_plot +
      geom_point(data = animation_data[[jj]], aes(x, y), size = 0) +
      geom_line(data = animation_data[[jj]], aes(x=x, y=y), color = "yellow2", size = 1) +
      geom_image(data = animation_data[[jj]], aes(x, y, image = img, frame = frame), size = .025) + 
      transition_states(frame, transition_length = 1, state_length = 1) +
      enter_fade() + exit_fade() +  transition_reveal(frame)
  }
  
  # add all the adults that are already there
  for(pp in 1:nrow(Adult_Positions)){
    
    x1 <- Adult_Positions$x[pp]
    y1 <- Adult_Positions$y[pp]
    
    Species <- Adult_Positions$species[pp]
    
    IMAGE <- Tree_Images_list[[Species]]
    
    anim_plot <- add_image(anim_plot,IMAGE, x1, y1, 6, 6)
  }  
  
  # Make Animation 
  #Animation <- animate(anim_plot, nframes = n_frames + pause_frames+15, fps = 10,height = 900, width =1600,end_pause = 15,renderer = gifski_renderer())
  anim_plot_New <- animate(anim_plot, nframes = Frames_Seeds + 15, res = 75,fps = 7, height = 900, width = 800, end_pause = 15, renderer = gifski_renderer())
  return(anim_plot_New)
}

# Will combine gifs in folder into single gif
combine_gifs <- function(folder_path, output_gif = "combined.gif", fps = 10) {
  # List all GIF files in the folder
  gif_files <- list.files(folder_path, pattern = "\\.gif$", full.names = TRUE)
  
  # Read all GIFs
  gifs <- lapply(gif_files, image_read)
  
  # Combine all GIFs into one
  combined_gif <- image_join(gifs)
  
  # Animate the combined GIF
  animated_gif <- image_animate(combined_gif, fps = fps)
  
  # Save the combined GIF
  image_write(animated_gif, output_gif)
  
  # Return the path to the combined GIF
  return(output_gif)
}


# set working directory for simulation data
setwd("C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Local_Files/Wiegand_et_al_2025_Local/Data")

# load all simulation data
Aggregation_Metric_Exponent_12 <- read.csv("Aggregation_Metric_Exponent_12.csv")
Crowding_index_Num_Neighbors_12 <- read.csv("Crowding_index_Num_Neighbors_12.csv")
dispersal_centers <- read.csv("dispersal_centers.csv")
Time_Series_05 <- read.csv("Time_Series_05.csv")
Spatial_Locations_05 <- read.csv("Spatial_Locations_05.csv")

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

# Ignore 1st time step... 
Spatial_Locations_05 <- subset(Spatial_Locations_05,timestep>1)

# set WD for tree images 
Tree_WD <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Images"
Tree_WD <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Tree_Images_3"
setwd(Tree_WD)
# Get a list of all files in the working directory
Tree_Images <- list.files()

Tree_image_files <- list.files(Tree_WD, pattern = "\\.png$", full.names = TRUE)
#Tree_Images_list <- lapply(Tree_image_files, image_read)
Tree_Images_list      <- lapply(Tree_image_files, png::readPNG)



# number of timesteps considered
TimeSteps <- 2

# Generate all combinations of three letters
combinations <- expand.grid(letters, letters, letters)

# Combine the letters into a single string for each combination
alphabetical_vector <- apply(combinations, 1, paste0, collapse = "")

# Select the first n elements
alphabetical_vector <- sort(alphabetical_vector[1:TimeSteps])

# Scale factor sects the number of frames
Scale_Factor <- 10#5

# Dimensions 
Mx <- 80
My <- 90
DELTA <- 2.2

# number of frames per timestep  
Frames_Seeds <- round(sqrt((Mx)^2 +(My)^2) /Scale_Factor)

# set image paths 
Seed_Image <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Seed_1.png"
X_Image_Path <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Other_Images/Red_X.png"
X_Image <- image_read(X_Image_Path)

# where to save? 
WD_Save <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Animation_Testing_2"

Animate_All(Spatial_Locations_05,TimeSteps,Frames_Seeds,Tree_Images,Tree_Images_list,Seed_Image,X_Image_Path,alphabetical_vector,WD_Save,Scale_Factor,Mx,My,DELTA)


combine_gifs(WD_Save,"Combinedt10.gif",5)









Test <- Animate_All(Spatial_Locations_05,TimeSteps,Frames_Seeds,Tree_Images,Tree_WD,Seed_Image,X_Image,alphabetical_vector,WD_Save,Scale_Factor)

last_frame <- image_read(Test)[length(image_read(Test))]

last_frame <- image_read(Births_anim)[length(image_read(Births_anim))]



last_frame_plot <- ggplot() +
  annotation_custom(rasterGrob(last_frame, width = unit(1, "npc"), height = unit(1, "npc")))




# Read the image to be added
X_Image_Path <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Other_Images/Red_X.png"
X_Image <- image_read(X_Image_Path)






# Create a function to add the image at specified coordinates
add_image_2 <- function(plot, x, y, image) {
  plot + annotation_custom(rasterGrob(image, width = unit(10, "npc"), height = unit(10, "npc")), xmin = x, xmax = x, ymin = y, ymax = y)
}

# Create a new plot on top of the last frame
new_plot <- last_frame_plot

# Add the image at each coordinate
for (i in 1:nrow(X_Y_Deaths)) {
  new_plot <- add_image_2(new_plot, X_Y_Deaths$x[i], X_Y_Deaths$y[i], X_Image)
}


# Add the image at each coordinate
dd <- 2
List_Plots <- list()
List_Plots[[1]] <- new_plot
for (dd in 1:(nrow(X_Y_Deaths))) {
  List_Plots[[dd+1]] <- add_image_b( List_Plots[[dd]],X_Image_Path,(X_Y_Deaths$x[dd]+4)/(Mx+5), (X_Y_Deaths$y[dd]+5)/(My+10),.35,.35) #add_image_2(new_plot, X_Y_Deaths$x[dd], X_Y_Deaths$y[dd], X_Image)
}
List_Plots[[dd+1]]


new_plot <- add_image_2(new_plot, 1, 1, X_Image)

new_plot <- add_image_b(new_plot,X_Image_Path, .5,.5,.35,.35)

new_plot

create_gif_from_plots(List_Plots)



create_gif_from_plots <- function(plot_list, gif_name = "combined_plots.gif", fps = 5, width = 800*3, height = 900*3) {
  # Create a temporary directory to save the plots
  temp_dir <- tempdir()
  
  # Save each plot in the list as a PNG file with specified dimensions
  plot_files <- lapply(seq_along(plot_list), function(i) {
    file_path <- file.path(temp_dir, paste0("plot", i, ".png"))
    ggsave(file_path, plot = plot_list[[i]], width = width / 300, height = height / 300, dpi = 300)
    file_path
  })
  
  # Read the images and combine them into a GIF
  images <- lapply(plot_files, image_read)
  gif <- image_animate(image_join(images), fps = fps)
  
  # Save the GIF
  image_write(gif, gif_name)
  
  # Delete the temporary PNG files
  unlink(plot_files)
  
  # Return the path to the GIF
  return(gif_name)
}

















# Animate the new plot
new_anim <- new_plot +
  transition_states(states = 1:nrow(X_Y_Deaths), transition_length = 2, state_length = 1) +
  ease_aes('linear')

# Render the animation
animated_gif <- animate(new_anim, nframes = nrow(X_Y_Deaths), res = 150, fps = 10, height = 900*3, width = 1600*3/2, end_pause = 15, renderer = gifski_renderer())

# Save the new animation
anim_save("new_animation.gif", animated_gif)




















# Define the length of the vector
n <- 100  # Replace with your desired length

# Generate all combinations of three letters
combinations <- expand.grid(letters, letters, letters)

# Combine the letters into a single string for each combination
alphabetical_vector <- apply(combinations, 1, paste0, collapse = "")

# Select the first n elements
alphabetical_vector <- sort(alphabetical_vector[1:n])

# Print the vector
print(alphabetical_vector)











grid_data <- expand.grid(x = (0-5):(192+5), y = (0-5):(108+5))

# Create the texture
texture <- rasterGrob(matrix(runif(250000), nrow = 500, ncol = 500), 
                      width = unit(1, "npc"), height = unit(1, "npc"), 
                      interpolate = TRUE)

Base_Plot <- ggplot(grid_data, aes(x, y)) +
  annotation_custom(texture, xmin = .5-5, xmax = 192.45+5, ymin = .55-5, ymax = 108.5+5) + # change these boundaries to make raster match 
  geom_tile(fill = "#034700", color = NA,alpha=.9, height=1.0065,width=1.0065, linetype = 1) +
  #geom_tile(fill = "#7B3F00", color = NA,alpha=.9, height=1.0065,width=1.0065, linetype = 1) +
  
  theme_void()


Spatial_Locations_05_A <- subset(Spatial_Locations_05,timestep==2)
Spatial_Locations_05_B <- subset(Spatial_Locations_05,timestep==3)

Spatial_Locations_05_A <- Spatial_Locations_05_A %>% arrange(desc(y))
Spatial_Locations_05_A <- subset(Spatial_Locations_05_A,species==1)

for(jj in 1:nrow(Spatial_Locations_05_A)){
  
  if(Spatial_Locations_05_A$Birth[jj] != Spatial_Locations_05_A$timestep[1]){
    
    x1 <- Spatial_Locations_05_A$x[jj]
    y1 <- Spatial_Locations_05_A$y[jj]
    
    
    Species <- Spatial_Locations_05_A$species[jj]
    
    setwd(Tree_WD)
    IMAGE <-  paste0(Tree_WD,"/",files[Species])
    
    Base_Plot <- add_image(Base_Plot,IMAGE, x1, y1, 6, 6)
  }
}  
#Base_Plot  

plant_img_path <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images"
setwd(plant_img_path)
ggsave("my_plot2.png", plot = Base_Plot, width = 16*3, height = 9*3, units = "in")





# Example usage
plant_img_path <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Tree_1.png"
seed_img_path <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Seed_1.png"


# Example usage
Plant_Image <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Tree_1.png"
Seed <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Seed_1.png"

create_animation(5, 5, 16, 8, Plant_Image, Seed)

x1 <- 5
y1 <- 5
x2 <- 16
y2 <- 8
n_frames <- 10
pause_frames <- 10

Number_Trees <- 50

for(jj in 1:Number_Trees){
  
  x <- sample(1:16,1)  
  y <- sample(1:9,1)  
  
  anim_plot <- add_image(anim_plot,plant_img_path, x, y, 1, 1)
}

animate(anim_plot, nframes = n_frames + pause_frames+15, fps = 10,height = 900, width =1600,end_pause = 15,renderer = gifski_renderer())


anim_save("Trees.gif",anim_plot,nframes = n_frames + pause_frames+15, fps = 10,height = 900, width =1600,end_pause = 15,renderer = gifski_renderer())

