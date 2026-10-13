require(nls2)
require(RColorBrewer)
require(lattice)
require(plyr)
require(mixtools)
require(scales)
library(RColorBrewer)
library(plotly)
library(lubridate)
library(tidyr)

split_to_list<-function(input,size){
  containers<-list()
  for (i in 1:(floor(length(input)/size)+1)){
    containers[[i]]<-list()
  }
  for (i in 1:length(input)){
    j<-ceiling(i/size)
    containers[[j]]<-append(containers[[j]],input[i])
  }
  containers
}


### ### ### ### ### ### ### ### ### ###
##### "read.batch.beads" function #####
### ### ### ### ### ### ### ### ### ###
read.batch.beads <- function (path,inc_date=F,inc_plate=F) {
  
  setwd(path)
  temp <- list.files(pattern="*.csv")
  plate_lab <- substr(temp,1,nchar(temp)-4)
  plate <- list()
  
  for(i in 1:length(temp)) { 
    
    plate_raw <- read.csv(temp[i])
    startrow <- grep("DataType",plate_raw[,1])[1]
    plate[[i]] <- read.csv(temp[i], skip=startrow+1)
    
    section1 <- grep("Count",plate[[i]]$Sample)[1]
    section2 <- grep("Avg Net MFI",plate[[i]]$Sample)
    
    plate[[i]] <- plate[[i]][(section1+2):(section2-2),-ncol(plate[[i]])]
    
    plate[[i]][,3:ncol(plate[[i]])] <- sapply(plate[[i]][,3:ncol(plate[[i]])],as.character)
    plate[[i]][,3:ncol(plate[[i]])] <- sapply(plate[[i]][,3:ncol(plate[[i]])],as.numeric)
    
    if (inc_date==T) {
      plate[[i]]$date <- plate_raw[grep("Date",plate_raw[,1]),2]
      plate[[i]]$date <- as.Date(plate[[i]]$date,"%m/%d/%Y")
    }
    
    if (inc_plate==T) {
      plate[[i]]$plate <- i
      
    }
  }
  
  names(plate) <- plate_lab
  
  return(plate)
  
}


### ### ### ### ### ### ### ### ##
##### "join.plates" function #####
### ### ### ### ### ### ### ### ##

# Binds all plates into a single object

join.plates <- function(plates,exc_blank=F) {
  
  output <- c()
  
  col_order <- colnames(plates[[1]])
  
  for (i in 1:length(plates)) {
    
    plates[[i]] <- as.data.frame(plates[[i]])
    plates[[i]] <- plates[[i]][col_order]
    output <- rbind(output,plates[[i]])
    
  }
  
  return(output)
  
}

### ### ### ### ### ### ### ### ### ### #
##### "plot_beads_by_well" function #####
### ### ### ### ### ### ### ### ### ### #
plot_beads_by_well <- function(plate_data, plate_name = NULL, threshold = input$bead_threshold) {
  # Standardize column name
  if (!"Well" %in% names(plate_data)) {
    if ("Location" %in% names(plate_data)) {
      plate_data <- plate_data %>% dplyr::rename(Well = Location)
    } else {
      stop("No 'Well' or 'Location' column found in plate_data")
    }
  }
  
  # Clean Location from 1(1,A1) -> 1
  plate_data$Well <- sub(".*,(.*)\\).*", "\\1", plate_data$Well)
  
  # Create vector of breaks (first well of each row)
  x_breaks <- unique(sub("(.)\\d+", "\\11", plate_data$Well))
  
  # Remove Total Events (if it is in the dataset)
  if("Total.Events" %in% names(plate_data)) {
    plate_data <- plate_data %>%
      select(-c(Total.Events))
  }
  
  ## Convert to long format to plot
  bead_data <- plate_data %>%
    pivot_longer(cols = -c(Sample, Well), names_to = "Antigen", values_to = "BeadCount")
  
  bead_data$below_threshold <- bead_data$BeadCount < threshold
  
  ## Plot
  p <- ggplot(bead_data, aes(x = Well, y = BeadCount, group = Antigen)) +
    geom_line() +
    geom_point(aes(color = below_threshold), size = 1) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), guide = "none") +
    geom_hline(yintercept = threshold, linetype = "dashed") +
    scale_x_discrete(breaks = x_breaks) +
    labs(title = ifelse(is.null(plate_name), "Bead counts by well", paste("Bead counts:", plate_name))) +
    theme_minimal() +
    theme(strip.text = element_text(size=20),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(hjust = 1, size=15), # angle = 45 (makes the label tilted if needing to fit for size)
          axis.title = element_text(size = 18),
          plot.title  = element_text(size = 20, face = "bold"),
          legend.position="none") + # Remove legend
    facet_wrap(~Antigen, scales = "free_y", ncol = 2)
  
  return(p)
  
}