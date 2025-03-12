# Load necessary libraries
library(ggplot2)
library(grid)
library(png)
library(gganimate)
library(ggimage)
library(magick)

add_image <- function(plot, img_path, x, y, width, height) {
  img <- rasterGrob(png::readPNG(img_path), width = unit(width, "npc"), height = unit(height, "npc"))
  plot + annotation_custom(img, xmin = x , xmax = x , ymin = y , ymax = y )
}

add_image <- function(plot, img_path, x, y, width, height) {
  img <- rasterGrob(png::readPNG(img_path), width = unit(width, "npc"), height = unit(height, "npc"))
  plot + annotation_custom(img, xmin = x - 0.5, xmax = x + 0.5, ymin = y - 0.25, ymax = y + 0.75)
}


Animate_All <- function(Spatial_Locations_Data,TimeSteps,Frames_Seeds,Tree_Images,Tree_WD,Seed_Image,X_Image,alphabetical_vector,WD_Save,Scale_Factor){
  
  for(tt in 2:TimeSteps){
    
    Spatial_Locations_TS <- subset(Spatial_Locations_Data,timestep==tt)
    Spatial_Locations_TS <- Spatial_Locations_TS %>% arrange(desc(y))
    Species <- unique(sort(Spatial_Locations_TS$Species))
    
    Living_Adults_TS <- subset(Spatial_Locations_TS,Birth != tt)
    X_Y_Living       <- Living_Adults_TS %>% select(x, y, species)
    
    New_Adults       <- subset(Spatial_Locations_TS,Birth == tt)
    X_Y_New          <- New_Adults %>% select(x, y, species)
    X_Y_Parent       <- New_Adults %>% select(x_parent, y_parent)
    
    Deaths           <- subset(Spatial_Locations_TS,Died==1)
    X_Y_Deaths       <- Deaths %>% select(x, y, species)
    
    
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
    Births_anim <- create_animation(X_Y_Living, animation_data_list, Tree_Images,Tree_WD)
    
    Name_animation <- paste0("TimeStep_",alphabetical_vector[tt],".gif")
    
    setwd(WD_Save)
    anim_save(Name_animation, Births_anim)
    
    #anim_save(Name_animation,Births_anim,nframes = Frames_Seeds +15, fps = 10,height = 900, width =1600,end_pause = 15,renderer = gifski_renderer())

    # add deaths 
  #  last_frame <- Births_anim[ZFrames_Seeds]
  #  last_plot <- ggplot_build(last_frame)
    
    
  #  final_plot <- last_plot +
  #    geom_image(data = X_Y_Deaths, aes(x = x, y = y), image = X_Image, size = 0.5) 
    
    #anim_with_images <- final_plot +
  #    transition_reveal(along = x)  # Adjust the transition as needed
    
  #  final_plot_New <- animate(anim_with_images, nframes =20, res = 150,fps = 10, height = 900*3, width = 1600*3, end_pause = 15, renderer = gifski_renderer())
  #  Name_animation <- paste0("TimeStep_",alphabetical_vector[tt],"_2",".gif")
    
  #  anim_save(Name_animation, final_plot_New)
    
    # remove dead individuals 
  
  }
  
}
    

# Function to create the animation
create_animation <- function(Adult_Positions, animation_data, Tree_Images,Tree_WD) {
  # Create a data frame for the grid
  grid_data <- expand.grid(x = (0-5):(192+5), y = (0-5):(108+5))
  
  # Create the texture
  texture <- rasterGrob(matrix(runif(250000), nrow = 500, ncol = 500), 
                        width = unit(1, "npc"), height = unit(1, "npc"), 
                        interpolate = TRUE)
  
  # Create the base animation plot
  anim_plot <- ggplot(grid_data, aes(x, y)) +
    annotation_custom(texture, xmin = .5-5, xmax = 192.45+5, ymin = .55-5, ymax = 108.5+5) + # change these boundaries to make raster match 
    geom_tile(fill = "#034700", color = NA, alpha = .9, height = 1.0065, width = 1.0065, linetype = 1) +
    theme_void() +
    theme(
      panel.background = element_rect(fill = "grey20", color = "grey20"),
      plot.margin = unit(c(-.5, -1.25, -.5, -1.25), "cm"),
      panel.spacing = unit(c(0, 0, 0, 0), "cm")
    )
  
  # Loop through each dataframe in the list and add layers to the plot
  for (jj in seq_along(animation_data)) {
    anim_plot <- anim_plot +
      geom_point(data = animation_data[[jj]], aes(x, y), size = 0) +
      geom_line(data = animation_data[[jj]], aes(x=x, y=y), color = "red", size = .5) +
      geom_image(data = animation_data[[jj]], aes(x, y, image = img, frame = frame), size = .025) + 
      transition_states(frame, transition_length = 1, state_length = 1) +
      enter_fade() + exit_fade() +  transition_reveal(frame)
  }
  
  # add all the adults that are already there
  for(pp in 1:nrow(Adult_Positions)){
    
      x1 <- Adult_Positions$x[pp]
      y1 <- Adult_Positions$y[pp]
      
      Species <- Adult_Positions$species[pp]

            IMAGE <-  paste0(Tree_WD,"/",Tree_Images[Species])
      
      anim_plot <- add_image(anim_plot,IMAGE, x1, y1, 6, 6)
  }  
  
  # Make Animation 
  #Animation <- animate(anim_plot, nframes = n_frames + pause_frames+15, fps = 10,height = 900, width =1600,end_pause = 15,renderer = gifski_renderer())
  anim_plot_New <- animate(anim_plot, nframes = Frames_Seeds + 15, res = 150,fps = 10, height = 900*3, width = 1600*3, end_pause = 15, renderer = gifski_renderer())
  return(anim_plot_New)
}

# set WD for tree images 
Tree_WD <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Images"
setwd(Tree_WD)
# Get a list of all files in the working directory
Tree_Images <- list.files()

# number of timesteps considered
TimeSteps <- 2

# Generate all combinations of three letters
combinations <- expand.grid(letters, letters, letters)

# Combine the letters into a single string for each combination
alphabetical_vector <- apply(combinations, 1, paste0, collapse = "")

# Select the first n elements
alphabetical_vector <- sort(alphabetical_vector[1:TimeSteps])

# Scale factor sects the number of frames
Scale_Factor <- 16#5

# number of frames per timestep  
Frames_Seeds <- round(sqrt((192)^2 +(108)^2) /Scale_Factor)

# set seed image path 
Seed_Image <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images/Seed_1.png"

X_Image_Path <- "C:/Users/smith/OneDrive/Desktop/Videos and More/Paper_Highlights/Paper_Highlights/Wiegand_et_al_2025/Other_Images/Red_X.png"
X_Image <- image_read(X_Image_Path)

# where to save? 
WD_Save <- "C:/Users/smith/OneDrive/Desktop/Teaching_Related_Doucments/Course Doucments/Mathematical Biology Course/General_Animation_Code/Images"


Animate_All(Spatial_Locations_05,TimeSteps,Frames_Seeds,Tree_Images,Tree_WD,Seed_Image,X_Image,alphabetical_vector,WD_Save,Scale_Factor)




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

