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
library(proto)
library(nls2)
library(lattice)
library(plyr)
library(tidyverse)

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


### ### ### ### ### ### ### ### #
##### "read.batch" function #####
### ### ### ### ### ### ### ### #

# Reads raw csv file from xPONENT software
# Removes header rows
# Keeps only median MFI values
# Removes last column "Total.Events"
# Ensures all MFI values are in numeric format

read.batch <- function (path, inc_date=F, inc_plate=F, save_path=NULL) {

  setwd(path)
  temp <- list.files(pattern="*.csv")
  plate_lab <- substr(temp,1,nchar(temp)-4)
  plate <- list()

  for(i in 1:length(temp)) { 

    plate_raw <- read.csv(temp[i])
    startrow <- grep("DataType",plate_raw[,1])[1]
    plate[[i]] <- read.csv(temp[i], skip=startrow+1)
    
    section <- grep("DataType",plate[[i]]$Location)[1]
    plate[[i]] <- plate[[i]][1:(section-1),-ncol(plate[[i]])]
    
    plate[[i]][,3:ncol(plate[[i]])] <- sapply(plate[[i]][,3:ncol(plate[[i]])],as.character)
    plate[[i]][,3:ncol(plate[[i]])] <- sapply(plate[[i]][,3:ncol(plate[[i]])],as.numeric)

    if (inc_date==T) {
      # This line of code reads in the date, but uses the package lubridate to make the date format ok with multiple date types as long as it is m/d/y
      plate[[i]]$date <- lubridate::mdy(plate_raw[grep("Date", plate_raw[,1]), 2])
    }
    
    if (inc_plate==T) plate[[i]]$plate <- i
    
  }
  
  names(plate) <- plate_lab
  
  if (is.null(save_path)==F) {
   for (i in 1:length(temp)){
     write.csv(plate[[i]],file=paste0(save_path,names(plate)[i],"_clean",".csv"),row.names=F)
   }
  }
  
  return(plate)
  
}




### ### ### ### ### ### ### ### ### ###
##### "read.batch.std_v2" function #####
### ### ### ### ### ### ### ### ### ###
read.batch.std_v2 <- function(path, std_labels = NULL, inc_date = FALSE, inc_plate = FALSE) {
  
  setwd(path)
  temp <- list.files(pattern = "*.csv")
  plate_lab <- substr(temp, 1, nchar(temp) - 4)
  plate <- list()
  
  for (i in 1:length(temp)) {
    
    plate_raw <- read.csv(temp[i])
    startrow  <- grep("DataType", plate_raw[, 1])[1]
    raw       <- read.csv(temp[i], skip = startrow + 1)
    names(raw) <- trimws(names(raw))
    raw$Sample <- trimws(as.character(raw$Sample))
    
    # TEMPORARY DEBUG — remove once fixed
    message("=== Plate ", i, " unique Sample values ===")
    message(paste(unique(raw$Sample), collapse = "\n"))
    
    # --- locate the two blocks ---
    med_start  <- grep("^Median$",         raw$Sample, ignore.case = TRUE)[1]
    med_end    <- grep("Net MFI",          raw$Sample, ignore.case = TRUE)[1]
    dil_start  <- grep("Dilution Factor",  raw$Sample, ignore.case = TRUE)[1]
    dil_end    <- grep("Analysis Types",   raw$Sample, ignore.case = TRUE)[1]
    
    if (any(is.na(c(med_start, med_end, dil_start, dil_end)))) {
      stop(paste("Plate", i, ": could not locate Median or Dilution Factor sections"))
    }
    
    # --- extract MFI block (rows P1, P2 ...) ---
    mfi_block <- raw[(med_start + 2):(med_end - 2), ]
    mfi_block <- mfi_block[, !names(mfi_block) %in% "Total.Events"]
    ag_cols   <- names(mfi_block)[3:ncol(mfi_block)]        # antigen column names
    mfi_block[, ag_cols] <- lapply(mfi_block[, ag_cols], function(x) suppressWarnings(as.numeric(as.character(x))))
    
    # --- extract Dilution Factor block ---
    dil_block <- raw[(dil_start + 2):(dil_end - 2), ]
    dil_block <- dil_block[, c("Sample", ag_cols[ag_cols %in% names(dil_block)])]
    
    # keep only antigen columns that exist in both blocks
    shared_ags <- intersect(ag_cols, names(dil_block))
    dil_block[, shared_ags] <- lapply(dil_block[, shared_ags], function(x) suppressWarnings(as.numeric(as.character(x))))
    
    # rename antigen columns in dilution block so they don't clash
    names(dil_block)[names(dil_block) %in% shared_ags] <- paste0(shared_ags, "_dil")
    
    # --- optionally filter to user-supplied std labels (e.g. P1, P2 ...) ---
    if (!is.null(std_labels) && length(std_labels) > 0) {
      pattern  <- paste(trimws(std_labels), collapse = "|")
      mfi_block <- mfi_block[grepl(pattern, mfi_block$Sample, ignore.case = TRUE), ]
      dil_block <- dil_block[grepl(pattern, dil_block$Sample, ignore.case = TRUE), ]
    }
    
    # --- merge MFI + dilution by Sample name ---
    combined <- merge(mfi_block, dil_block, by = "Sample", all.x = TRUE)
    
    # for each antigen, use the *mean* dilution across wells as the x-axis value
    # (Dilution Factor block often has one value per standard point)
    for (ag in shared_ags) {
      dil_col <- paste0(ag, "_dil")
      if (dil_col %in% names(combined)) {
        combined[[paste0(ag, "_dilution")]] <- combined[[dil_col]]
        combined[[dil_col]] <- NULL
      }
    }
    
    plate[[i]] <- combined
    
    if (inc_date) {
      plate[[i]]$date <- as.Date(plate_raw[grep("Date", plate_raw[, 1]), 2], "%m/%d/%Y")
    }
    if (inc_plate) {
      plate[[i]]$plate <- i
    }
  }
  
  names(plate) <- plate_lab
  return(plate)
}
### ### ### ### ### ### ### ### ### #
##### "read.batch.std" function #####
### ### ### ### ### ### ### ### ### #

read.batch.std <- function(path, std_labels, inc_date = FALSE, inc_plate = FALSE) {
  
  setwd(path)
  temp <- list.files(pattern = "*.csv")
  plate_lab <- substr(temp, 1, nchar(temp) - 4)
  plate <- list()
  
  for(i in 1:length(temp)) { 
    
    plate_raw <- read.csv(temp[i])
    startrow <- grep("DataType",plate_raw[,1])[1]
    plate[[i]] <- read.csv(temp[i], skip=startrow+1) # original code
    #plate[[i]] <- read.csv(temp[i], skip=startrow)
    
    # Normalize column names
    names(plate[[i]]) <- trimws(names(plate[[i]]))
    
    sample_col <- grep("Sample", names(plate[[i]]), ignore.case = TRUE, value = TRUE)
    
    if (length(sample_col) == 0) {
      stop(paste("Plate", i, ": no Sample column found"))
    }
    
    sample_col <- sample_col[1]
    plate[[i]]$Sample <- trimws(as.character(plate[[i]][[sample_col]]))
    
    section1 <- grep("^\\s*Median\\b", trimws(plate[[i]]$Sample), ignore.case = TRUE)
    section1_end <- grep("Net MFI", plate[[i]]$Sample, ignore.case = TRUE)
    
    section2 <- grep("Dilution Factor", plate[[i]]$Sample, ignore.case = TRUE)
    section2_end <- grep("Analysis Types", plate[[i]]$Sample, ignore.case = TRUE)
    
    if (length(section1) == 0 || length(section1_end) == 0 ||
        length(section2) == 0 || length(section2_end) == 0) {
      stop(paste("Plate", i, ": standard curve sections not found"))
    }
    
    # ONLY NOW take the first index
    section1 <- section1[1]
    section1_end <- section1_end[1]
    section2 <- section2[1]
    section2_end <- section2_end[1]
    
    # Create Median and Dilution data blocks
    plate_i_Median <- plate[[i]][(section1 + 2):(section1_end - 2), -ncol(plate[[i]])]
    plate_i_dilution <- plate[[i]][(section2 + 2):(section2_end - 2), -ncol(plate[[i]])]
    
    # Create regex pattern from user input
    pattern <- paste(std_labels, collapse = "|")

    # Filter rows in both datasets to only contain manual std values
    plate_i_Median  <- plate_i_Median[grepl(pattern, plate_i_Median$Sample, ignore.case = TRUE), ]
    plate_i_dilution <- plate_i_dilution[grepl(pattern, plate_i_dilution$Sample, ignore.case = TRUE), ]


    # Make first 2 columns character and the rest numeric for both datasets
    plate_i_Median[, 1:2] <- lapply(plate_i_Median[, 1:2],as.character)
    plate_i_Median[,3:ncol(plate_i_Median)] <- lapply(plate_i_Median[,3:ncol(plate_i_Median)],as.numeric)

    plate_i_dilution[, 1:2] <- lapply(plate_i_dilution[, 1:2],as.character)
    plate_i_dilution[,3:ncol(plate_i_dilution)] <- lapply(plate_i_dilution[,3:ncol(plate_i_dilution)],as.numeric)
    
    # Combine datasets into main dataset
    plate_i_combined <- merge(
      plate_i_Median, 
      plate_i_dilution, 
      by = c("Sample"), 
      all.x = TRUE
    )
    
    # Save as plate[[i]]
    plate[[i]] <- plate_i_combined
    
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

### ### ### ### ### ### ### ### #
##### "get_clean_plate_map" #####
### ### ### ### ### ### ### ### #
get_clean_plate_map <- function(df, plate_name, plate_size) {
  # 1. Filter for the specific plate
  res <- df[df$Plate == plate_name, ]
  
  # 2. Define valid row letters
  rows <- if(plate_size == "96") LETTERS[1:8] else LETTERS[1:16]
  
  # 3. Keep only rows that fall within those letters
  # This ensures Well_Type and all other columns remain attached
  res <- res[res$Row_ID %in% rows, ]
  
  return(list(data = res, rows = rows))
}
# get_clean_plate_map <- function(data, plate_name, plate_size) {
#   # 1. Filter for the specific plate
#   df <- data[data$Plate == plate_name, ]
#   
#   # 2. STRICT REGEX FILTER
#   # This looks for the xPONENT well format: 1(1,A1)
#   # It effectively ignores the "Warnings/Errors" rows at the bottom
#   df <- df[grepl("\\([0-9]+,[A-P][0-9]+\\)", df$Location), ]
#   
#   # 3. Standardize Extraction
#   df$Well_ID <- sub(".*,([A-P][0-9]+)\\).*", "\\1", df$Location)
#   df$Row_ID  <- substr(df$Well_ID, 1, 1)
#   df$Col_ID  <- as.numeric(gsub("[A-P]", "", df$Well_ID))
#   
#   # 4. Final safety: Remove any row where Row_ID is NA
#   df <- df[!is.na(df$Row_ID), ]
#   
#   # 5. Define Row Limits based on size
#   valid_rows <- if(plate_size == "96") LETTERS[1:8] else LETTERS[1:16]
#   df <- df[df$Row_ID %in% valid_rows, ]
#   
#   return(df)
# }

### ### ### ### ### ### ### ##
##### "convert_dilution" #####
### ### ### ### ### ### ### ##
# Helper function for the function plot.std.curve3_interactive 
# Converts dilutions from fraction format to log numeric

convert_dilution <- function(dilution_str) {
  if (grepl("/", dilution_str)) {
    parts <- strsplit(dilution_str, "/")[[1]]
    return(as.numeric(parts[1]) / as.numeric(parts[2]))
  } else {
    return(as.numeric(dilution_str))
  }
}


### ### ### ### ### ### ### ### ### ### ### ### ##
##### "plot.std.curve3_interactive" function #####
### ### ### ### ### ### ### ### ### ### ### ### ##
plot.std.curve3_interactive <- function(data_list, antigen, plate_selected, dilutions, std_label = NULL) {
  
  # Convert dilution strings → numeric dilution factors
  dilution_numeric <- sapply(dilutions, convert_dilution)
  
  curves <- list()
  n_dil <- length(dilutions)
  
  # Extract standard curves for each plate
  for (i in seq_along(data_list)) {
    plate_data <- data_list[[i]]
    if (!(antigen %in% colnames(plate_data))) next
    std_curve <- get.standard(plate_data,dilutions,std_label = std_label)
    curves[[ names(data_list)[i]]] <- std_curve[[antigen]]
  }
  
  # Convert list → tidy dataframe
  df <- do.call(rbind, lapply(names(curves), function(plate) {
    data.frame(
      Plate = plate,
      Dilution = factor(dilutions, levels = dilutions),
      MFI = curves[[plate]]
    )
  }))
  
  # Build ggplot
  p <- ggplot(df, aes(x = Dilution, y = log(MFI), color = Plate, group = Plate)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    labs(
      title = paste("Standard Curve for", antigen),
      x = "Dilution",
      y = "Log MFI"
    ) +
    theme_minimal(base_size = 14)
  
  # Convert ggplot → plotly
  return(p)
}

# OTHER OLD CODE
# plot.std.curve3_interactive <- function(data, antigen, dilutions, std_label, negs=NULL, blanks=NULL, plate_labels=NULL) {
#   
#   curves <- list()
#   max.set <- c()
#   plate_labels_def <- c()
#   
#   # Split the user-defined standard labels
#   std_labels <- unlist(strsplit(std_label, ","))   # e.g., "P1|P2|P3|P4|P5|P6|P7|P8"
#   
#   for (j in 1:length(data)){
#     
#     if (!(antigen %in% colnames(data[[j]]))) {
#       print(paste0(" '", antigen, "' does not exist in dataset ", j, " - check spelling "))
#       next
#     }
#     
#     # Filter standard rows
#     std_rows <- data[[j]]$Sample %in% std_labels       # Subset only standard rows
#     standard_data <- data[[j]][std_rows, ]
#     
#     # Call function get.standard on subset
#     std_curve <- get.standard(standard_data, std_label = std_label, dilutions = dilutions, n_points = length(dilutions))
#     
#     if (!is.null(std_curve)) {
#       curves[[j]] <- std_curve[[antigen]]            # extract only antigen column
#       max.set <- c(max.set, curves[[j]])
#       plate_labels_def <- c(plate_labels_def, paste0("plate", j))
#     } else {
#       curves[[j]] <- NA
#     }
#     
#   }
#   
#   pnames <- names(data)
#   max.set <- max.set[is.finite(max.set)]
#   max <- max(max.set)
#   
#   n_dil <- length(dilutions)
#   
#   # Plot with plotly
#   fig <- plot_ly(x = 1:n_dil, y = curves[[1]], mode="lines", type="scatter", name=pnames[1], showlegend=F)
#   for (i in 2:length(data)){
#     fig <- fig %>% add_trace(x = 1:n_dil, y = curves[[i]], mode="lines", type="scatter", name=pnames[i], showlegend=F)
#   }
#   
#   fig <- fig %>% layout(
#     xaxis = list(
#       tickmode = "array",
#       ticktext = as.list(dilutions),
#       tickvals = as.list(1:n_dil)
#     ),
#     title = antigen
#   )
#   
#   return(fig)
# }

# ORIGINAL CODE
# plot.std.curve3_interactive <- function(data, antigen, dilutions, negs=NULL, blanks=NULL, plate_labels=NULL) {
#   curves <- list()
#   max.set <- c()
#   plate_labels_def <- c()
#   
#   for (j in 1:length(data)){
#     
#     if ((antigen %in% colnames(data[[j]]))==F) {
#       print(paste0(" '",antigen,"'"," does not exist in dataset ",j," - check spelling "))
#       next
#     } else {
#       curves[[j]] <- data[[j]][,which(names(data[[j]])==antigen)]
#       max.set <- c(max.set,curves[[j]])
#       plate_labels_def <- c(plate_labels_def,paste0("plate",j))
#     }
#     
#   }
#   pnames<-names(data)
#   max.set <- max.set[is.finite(max.set)]
#   max <- max(max.set)
#   
#   conc <- sapply(strsplit(dilutions, split = "/"),function(x) as.numeric(x[1]) / as.numeric(x[2]))
#   
#   n_dil <- length(dilutions)
#   fig <- plot_ly(x = 1:n_dil, y = curves[[1]], mode="lines", type="scatter", name=pnames[1], showlegend=F)
#   for (i in 2:length(data)){
#     fig <- fig %>% add_trace(x = 1:n_dil, y = curves[[i]], mode="lines", type="scatter", name=pnames[i], showlegend=F)
#   }
#   fig <- fig %>% layout(
#     xaxis = list(
#       tickmode = "array",
#       ticktext = as.list(dilutions),
#       tickvals = as.list(1:n_dil)
#     ),
#     title = antigen
#   )
#   
#   # fig<- plot_ly(x = 1:6, y= curves[[1]],mode="lines",type="scatter",name=pnames[1],showlegend = F)
#   # 
#   # for (i in 2:length(data)){
#   #   fig <- fig %>% add_trace(x = 1:6, y= curves[[i]],mode="lines",type="scatter",name=pnames[i],showlegend = F)
#   # }
#   # fig <- fig %>% layout( 
#   #   xaxis = list(
#   #     tickmode="array",
#   #     ticktext=as.list(dilutions[1:6]), 
#   #     tickvals=as.list(1:5)
#   #     ),
#   #     title = antigen
#   #   )
#   fig
#   
# }


### ### ### ### ### ### ### ### ###
##### "blank.adjust" function #####
### ### ### ### ### ### ### ### ###

# Calculates mean MFI values for "blank" samples for each antigen
# Subtracts mean blank value from samples MFI values
# Sets any negative MFI values after subtraction to 0

blank.adjust <- function (data, blank_label, bkg_select) {

  blank_mean <- list()
  blank <- data[grep(blank_label,data$Sample),]

  if (nrow(blank)==0) {

    print(paste0("Error: No samples with label"," '",blank_label,"' ","included in this dataset"))

  } else {

    blank_mean <- apply(blank[,3:ncol(data)], 2, mean)

    for (i in 3:ncol(data)) {

    # Subtract
      if(bkg_select() == "Subtract"){
      data[-grep(blank_label,data$Sample),i] <- data[-grep(blank_label,data$Sample),i] - blank_mean[i-2]

        for (j in 1:nrow(data)){
  
          if (is.na(data[j,i])==T) data[j,i] <- data[j,i]
          else if (data[j,i] < 0) data[j,i] <- 0
        }
      }
      
    # Divide
      if(bkg_select() == "Divide"){
        data[-grep(blank_label,data$Sample),i] <- data[-grep(blank_label,data$Sample),i] / blank_mean[i-2]
        
        for (j in 1:nrow(data)){
          
          if (is.na(data[j,i])==T) data[j,i] <- data[j,i]
          else if (data[j,i] < 0) data[j,i] <- 0
        }
      }
    }
  }

  return(data)

}

    # ORIGINAL FUNCTION
          # blank.adjust <- function (data, blank_label) {
          # 
          #   blank_mean <- list()
          #   blank <- data[grep(blank_label,data$Sample),]
          #   
          #   if (nrow(blank)==0) {
          #     
          #     print(paste0("Error: No samples with label"," '",blank_label,"' ","included in this dataset"))
          #     
          #   } else {
          #   
          #     blank_mean <- apply(blank[,3:ncol(data)], 2, mean)
          # 
          #     for (i in 3:ncol(data)) {  
          #     
          #       data[-grep(blank_label,data$Sample),i] <- data[-grep(blank_label,data$Sample),i] - blank_mean[i-2]
          #     
          #       for (j in 1:nrow(data)){
          #       
          #         if (is.na(data[j,i])==T) data[j,i] <- data[j,i]
          #         else if (data[j,i] < 0) data[j,i] <- 0
          #         
          #       }
          #     }
          #   }
          #   
          #   return(data)
          #   
          # }

blank.adjust.batch <- function(plates,blank_label) {

  blank <- list()
  blank_mean <- list()

  for (i in 1:length(plates)){
    
    blank[[i]] <- plates[[i]][grep(blank_label,plates[[i]]$Sample),]
    
      if (nrow(blank[[i]])==0) {
      
        print(paste0("Error: No samples with label"," '",blank_label,"' ","included in"," ",names(plates)[i]))
        
      } else {
    
        blank_mean[[i]] <- apply(blank[[i]][,3:ncol(blank[[i]])],2,mean)
    
      }
  }

  for (i in 1:length(plates)){
    
    if (nrow(blank[[i]])==0) {
      
      next
      
    } else {
    
      for (j in 3:ncol(plates[[i]])) {  
      
        plates[[i]][-grep(blank_label,plates[[i]]$Sample),j] <- plates[[i]][-grep(blank_label,plates[[i]]$Sample),j] - blank_mean[[i]][j-2]
      
        for (k in 1:nrow(plates[[i]])){
        
          if (is.na(plates[[i]][k,j])==T) plates[[i]][k,j] <- plates[[i]][k,j]
        
          else if (plates[[i]][k,j] < 1) plates[[i]][k,j] <- 1
        
        }
      }
    }
  }
  
  return(plates)
  
}

### ### ### ### ### ### ### ##
##### "exclude" function #####
### ### ### ### ### ### ### ##

## removes blank, neg or pos control entries from dataset

exclude <- function(data,blank_label=NULL,neg_label=NULL,std_label=NULL) {
  
  if(is.null(blank_label)==F) data <- data[-grep(blank_label,data$Sample),] else data <- data
  if(is.null(neg_label)==F) data <- data[-grep(neg_label,data$Sample),] else data <- data
  if(is.null(std_label)==F) data <- data[-grep(std_label,data$Sample),] else data <- data
        
  return(data)

}

### ### ### ### ### ### ### ### ###
##### "get.standard" function #####
### ### ### ### ### ### ### ### ###
# 1. Normalise dilution labels
# Examples: 1/100, 1:100 and 1 / 100 become 1/100

normalise_dilution <- function(x) {
  
  x <- gsub("\\s+", "", trimws(as.character(x)))
  
  vapply(strsplit(x, "[/:]", perl = TRUE), function(parts) {
    
    if (length(parts) == 2 &&
        all(grepl("^[0-9]+$", parts))) {
      
      paste0(
        as.integer(parts[1]), "/",
        as.integer(parts[2])
      )
      
    } else {
      NA_character_
    }
    
  }, character(1))
}


# 2. Identify standard names and their dilutions

identify.standards <- function(data, dilutions) {
  
  if (!"Sample" %in% names(data)) {
    stop("The plate data must contain a 'Sample' column.")
  }
  
  dilutions <- normalise_dilution(dilutions)
  dilutions <- unique(dilutions[!is.na(dilutions)])
  
  result <- data.frame(
    Sample = as.character(data$Sample),
    Standard = NA_character_,
    Dilution = NA_character_,
    IsStandardCandidate = FALSE,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(result$Sample)) {
    
    sample_name <- trimws(result$Sample[i])
    
    if (is.na(sample_name) || !nzchar(sample_name)) {
      next
    }
    
    for (dilution in dilutions) {
      
      # Match a dilution at the END of the sample name.
      # Accept / or : as the dilution separator,
      # and underscore, space, colon or hyphen before it.
      
      parts <- strsplit(dilution, "/", fixed = TRUE)[[1]]
      
      dilution_pattern <- paste0(
        parts[1], "[/:]", parts[2]
      )
      
      pattern <- paste0(
        "[_[:space:]:-]*",
        dilution_pattern,
        "$"
      )
      
      if (grepl(pattern, sample_name, ignore.case = TRUE)) {
        
        standard_name <- sub(
          pattern, "", sample_name,
          ignore.case = TRUE
        )
        
        standard_name <- trimws(standard_name)
        
        if (nzchar(standard_name)) {
          
          result$Standard[i] <- standard_name
          result$Dilution[i] <- dilution
          result$IsStandardCandidate[i] <- TRUE
          
        }
        
        break
      }
    }
  }
  
  result
}


# 3. Extract the selected standard curve

get.standard <- function(data, dilutions, std_label = NULL, n_points = NULL) {
  
  if (!"Sample" %in% names(data)) {
    stop("The plate data must contain a 'Sample' column.")
  }
  
  original_dilutions <- trimws(as.character(dilutions))
  dilutions <- normalise_dilution(dilutions)
  
  if (length(dilutions) == 0 || anyNA(dilutions)) {
    stop("Please enter valid dilutions, e.g. 1/50, 1/100, 1/200.")
  }
  
  dilutions <- unique(dilutions)
  
  identified <- identify.standards(data, dilutions)
  
  # Match the requested standard name, ignoring case
  if (is.null(std_label)) {
    keep <- identified$IsStandardCandidate
  } else {
    keep <- !is.na(identified$Standard) &
      identified$IsStandardCandidate &
      tolower(identified$Standard) == tolower(std_label)
  }
  
  if (!any(keep)) {
    stop(
      paste0(
        "No samples found for standard '", std_label,
        "' with the specified dilutions."
      )
    )
  }
  
  std <- data[keep, , drop = FALSE]
  
  std$Dilution <- identified$Dilution[keep]
  std$Standard <- identified$Standard[keep]
  
  # Keep the dilution order entered by the user.
  # Multiple rows per dilution are allowed for replicates.
  std$dil_no <- match(std$Dilution, dilutions)
  
  std <- std[order(std$dil_no), , drop = FALSE]
  
  # Warn about expected dilutions that were not found
  missing_dilutions <- setdiff(
    dilutions,
    unique(std$Dilution)
  )
  
  if (length(missing_dilutions) > 0) {
    warning(
      paste0(
        "Standard '", std_label,
        "' is missing dilution(s): ",
        paste(missing_dilutions, collapse = ", ")
      )
    )
  }
  
  rownames(std) <- NULL
  
  std
}
  
  # data$Sample <- as.character(data$Sample)
  # dilutions <- trimws(as.character(dilutions))
  # 
  # # Add a dilution number column
  # data$dil_no <- NA
  # for (i in seq_along(dilutions)) {
  #   hits <- grep(dilutions[i], data$Sample)
  #   if (length(hits) > 0) data$dil_no[hits] <- i
  # }
  
#   # Order by dilution
#   data <- data[order(data$dil_no, na.last = NA), ]
#   
#   # Numeric columns
#   data[, 3:ncol(data)] <- sapply(data[, 3:ncol(data)], as.numeric)
#   return(data)
# }

    # OTHER OLD CODE
      # get.standard <- function(data, std_label, dilutions, n_points) {
      #   
      #   # Make sure character vectors
      #   data$Sample <- as.character(data$Sample)
      #   std_label <- as.character(std_label)
      #   dilutions <- as.character(dilutions)
      #   
      #   # Instead of regex, match exact Sample names
      #   matches <- which(data$Sample %in% std_label)
      #   
      #   if(length(matches) == 0) {
      #     warning("No standard curve samples found for std_label: ", std_label)
      #     return(NULL)
      #   }
      #   
      #   if(length(matches) < n_points) {
      #     warning("Fewer standard samples than n_points. Using all available.")
      #     std <- data[matches, ]
      #   } else {
      #     std <- data[matches[n_points:1], ]
      #   }
      #   
      #   # Safely convert numeric columns
      #   num_cols <- sapply(std, is.numeric)
      #   std[,num_cols] <- lapply(std[,num_cols], as.numeric)
      #   
      #   # Assign dilution numbers
      #   std$dil_no <- NA
      #   for (i in seq_along(dilutions)) {
      #     hits <- grep(dilutions[i], std$Sample)
      #     if(length(hits) > 0){
      #       std$dil_no[hits] <- i
      #     }
      #   }
      #   
      #   # Order by dilution number
      #   std <- std[order(std$dil_no, na.last = NA), ]
      #   
      #   # Convert columns to numeric
      #   std[,3:ncol(data)] <- sapply(std[,3:ncol(data)], as.character)
      #   
      #   if ("date" %in% colnames(std)) {
      #     std[,3:(which(colnames(std)=="date")-1)] <- sapply(std[,3:(which(colnames(std)=="date")-1)], as.numeric)
      #   } else {
      #     std[,3:ncol(data)] <- sapply(std[,3:ncol(data)], as.numeric)
      #   }
      #   
      #   return(std)
      # }

          # ORIGINAL CODE
          # get.standard <- function(data,std_label = std_label, dilutions,n_points){
          #   
          #   # Convert std_label to regex that matches any substring
          #   pattern <- paste0("(?i)", std_label)  # (?i) for case-insensitive
          #   
          #   std <- data[grep(std_label, data$Sample)[n_points:1],]
          #   
          #   # for (i in 1:length(dilutions)) {
          #   #   std$dil_no[grep(dilutions[i], std$Sample)] <- i
          #   # }
          #   
          #   for (i in seq_along(dilutions)) {
          #     hits <- grep(dilutions[i], std$Sample)
          #     if(length(hits) > 0){
          #       std$dil_no[hits] <- i
          #     }
          #   }
          #   
          #   std <- std[order(std$dil_no),]
          #   
          #   std[,3:ncol(data)] <- sapply(std[,3:ncol(data)],as.character)
          #   
          #   if ("date" %in% colnames(std)) {
          #   
          #     std[,3:which(colnames(std)=="date")-1] <- sapply(std[,3:which(colnames(std)=="date")-1],as.numeric)
          #   
          #   } else std[,3:ncol(data)] <- sapply(std[,3:ncol(data)],as.numeric)
          #   
          #   return(std)
          #   
          # }  



### ### ### ### ### ### ### #
##### "coefv" function #####
### ### ### ### ### ### ### #

coefv<-function(data, antigens, std_label = "(?i)CP3|Std Curve", dilutions, target_dilution){
  sapply(antigens,function(a){
    d<-sapply(names(data),function(p){
      tmp<-get.standard(data[[p]],std_label = std_label,dilutions = dilutions,n_points = 6)
      idx < -which(grepl(std_label, tmp$Sample) & grepl(paste0(target_dilution,"$"),tmp$Sample))
      tmp[idx,a]
    })
    sd(d)/mean(d)
  })
}


### ### ### ### ### ### ### ### #
##### "plot.coefv" function #####
### ### ### ### ### ### ### ### #

plot.coefv<-function(data,antigens,std_label,dilutions,target_dilution){
  vals<-coefv(data,antigens,std_label,dilutions,target_dilution)
  # plot(vals,type="o",pch=19,ylim=c(0,1),xlab="Antigen",ylab="Coefficient of variation",bty='L',main="Coefficient of variation")
  # abline(h=0.2,col="red",lty=3,lwd=1.5)
  # abline(h=0.4,col="red",lty=2,lwd=1.5)
  plot_ly(x = 1:length(vals), y= vals,type="scatter",mode="lines+markers",name="coefv",showlegend = F) %>% 
  layout( 
    title=paste("Coefficient of variation",target_dilution),
    xaxis = list(
      tickmode="array",
      ticktext=names(vals), 
      tickvals=as.list(1:length(vals))
    ),
    yaxis = list(
      range=c(0,max(vals))
    )
  )

}


### ### ### ### ### ### ### ### ### ##
##### "view.plate.data" function #####
### ### ### ### ### ### ### ### ### ##

view.plate.data <- function(data,antigen) {
  
  if ((antigen %in% colnames(data))==F) {
    
    print("Antigen name entered does not exist in dataset - check spelling")
    
  } else {
    
  layout <- sample_ids <- matrix (nrow=8,ncol=12)
  colnames(layout) <- colnames(sample_ids) <- seq(1,12,1)
  rownames(layout) <- rownames(sample_ids) <- c("A","B","C","D","E","F","G","H")
  
  for (i in 1:length(rownames(layout))) {
    
    for (j in 1:length(colnames(layout))) {
      
      rowtemp <- rownames(layout)[i]
      coltemp <- colnames(layout)[j]
      
      if (length(grep(paste0(rowtemp,coltemp,")"),data$Location))==0) {
        layout[i,j] <- sample_ids[i,j] <- NA
      } else {
  
        layout[i,j] <- as.numeric(as.character(data[grep(paste0(rowtemp,coltemp,")"),data$Location),antigen]))
        sample_ids[i,j] <- as.character(data[grep(paste0(rowtemp,coltemp,")"),data$Location),"Sample"])
      
      }
    
    }
    
  }
  
  output <- list(layout,sample_ids)
  names(output) <- c("sample.data","sample.ids")
  
  return(output)
  
  }
  
}

### ### ### ### ### ### ### ### ### #
##### "normalise.conc" function #####
### ### ### ### ### ### ### ### ### #

norm.conc <- function(data,dilutions,standards,ag_list,plot=F,path=getwd(),filename="data_conc.csv") {

  sample_conc <- data.frame(matrix(nrow=nrow(data), ncol=length(ag_list)+2))
  colnames(sample_conc) <- c(colnames(data)[1:2],ag_list)
  sample_conc[,1] <- as.character(data[,1])
  sample_conc[,2] <- as.character(data[,2])

  fit.mat <- as.data.frame(matrix(ncol=length(ag_list)+1,nrow=100))
  fit.mat[,1] <- seq(1,100,1)
  colnames(fit.mat) <- c("conc_no",ag_list)
  
  # par(mfrow=c(3,3))

  for (i in 14:length(ag_list)){ 
    
    conc <- sapply(strsplit(dilutions, split = "/"),
                 function(x) as.numeric(x[1]) / as.numeric(x[2]))
    
    std <- standards[,ag_list[i]]
    
    if (all(is.na(std))==T) {
      
      print(paste0("No standard curve data for ",ag_list[i]))
      next
      
    }

    max <- max(std[which(is.finite(std))])  
  
    if (is.na(std[6])==T|is.na(std[5])==T) {
      
    } else if (std[6] < std[5]) {
        std <- std[-6]
        conc <- conc[-6]
    }
    
    len_std <- length(std)
    
    fit.data <- data.frame("MFI"=std,"conc"=log(conc))
    
    cc1 <- try( starting.values <- getInitial(MFI ~ SSlogis(conc,maxMFI,Xmid,Scale),data=fit.data), silent=T)
    cc2 <- try( nls(MFI ~ (maxMFI / (1 + exp((Xmid - conc)/Scale))),data=fit.data, start=starting.values), silent=T)
    
    if (is(cc1,"try-error")|is(cc2,"try-error")) {
      
      if (is.finite(std[len_std])==T) std_last <- std[len_std] else std_last <- std[len_std-1]
      
      starting.values <- expand.grid(maxMFI=seq(0, std_last*2.5, len = 20),
                                     Xmid=seq(-50, 0, len = 20),
                                     Scale=seq(0, 50, len = 40))
      
      fit <- nls2(MFI ~ (maxMFI / (1 + exp((Xmid - conc)/Scale))), algorithm="brute-force",data=fit.data, start=starting.values)
      
    } else {
      
      starting.values <- getInitial(MFI ~ SSlogis(conc,maxMFI,Xmid,Scale),data=fit.data)
      fit <- nls(MFI ~ (maxMFI / (1 + exp((Xmid - conc)/Scale))),data=fit.data, start=starting.values)
      
    }
    
    params <- as.list(coef(fit))
    
    scale <- params$Scale
    mid <- params$Xmid
    max <- params$maxMFI
    
    x1 <- seq(min(fit.data$conc),max(fit.data$conc),length=100)
    y1 <- predict(fit, newdata=list(conc=x1))
    fit.mat[,i+1] <- y1
    
    logk <- params$Xmid - params$Scale*(log(params$maxMFI/(params$maxMFI*0.5) - 1))
    km <- exp(logk)
    
    if (plot==T) {
    
      plot(log(conc),std, type="p", pch=19, col="dodgerblue4",
         ylab="MFI", xlab="Dilution factor",axes=F, 
         ylim=c(0,max+4000), xlim=c(min(log(conc)),max(log(conc))+1),
         main=paste0(ag_list[i]," Standard curve "), cex.main=0.8, cex.lab=0.8 )
      axis(1, at=c(log(conc),max(log(conc))+1) , labels=c(dilutions[1:len_std],""), cex.axis=0.8)
      axis(2, cex.axis=0.8)
      lines(x1,y1,type="l",lwd=2, col="dodgerblue4")
      abline(v=logk,lty=3)
    
    }
    
    for (j in 1:nrow(data)){
    
      if (data[j,i+2]<=0) mm_par <- (max/rnorm(1,mean=10,sd=3))-1
      else mm_par <- (max/data[j,i+2])-1

      par_temp <- mid-scale*(log(mm_par))
    
      sample_conc[j,i+2] <- exp(par_temp)
      
    }
    
  }
  
  sample_conc[,1:2] <- data[,1:2]

  write.csv(sample_conc,file=paste0(path,filename),row.names=F)
  
  output <- list(sample_conc,fit.mat)
  names(output) <- c("sample_conc","fit_mat")
  
  return(output)
  
}

### ### ### ### ### ### ### ### ### ### ###
##### "normalise.conc.batch" function #####
### ### ### ### ### ### ### ### ### ### ###

norm.conc.batch <- function(data_batch,standards_batch,ag_list,dilutions,plot=F,plate_label) {
  
  output1 <- list()
  output2 <- list()
  
  for (i in 1:length(data_batch)){
    
    aux.data <- norm.conc(data=data_batch[[i]],standards=standards_batch[[i]],ag_list=ag_list,dilutions=dilutions,plot=F)
   
    output1[[i]] <- aux.data$sample_conc
    output2[[i]] <- aux.data$fit_mat
    output1[[i]][,"plate"] <- plate_label[[i]]

  }
  
  output <- list(output1,output2)
  names(output) <- c("sample_conc","fit_mat")
  
  return(output)
  
}


### ### ### ### ### ### ### ### #
##### "plot.plate" function #####
### ### ### ### ### ### ### ### #

plot.plate <- function(data1,antigen,plotfile="plot1.pdf",plot=F,path=getwd()) {
  
  cols <- colorRampPalette(c("white","palevioletred3"))(250)
  
    plate <- levelplot(t(data1[8:1,]),col.regions=cols,main=paste0(antigen," MFI values"),xlab=NULL,ylab=NULL,border="white")
    
    if (plot==T) {
      pdf(paste0(path,plotfile))
      plot(plate)
      dev.off()
    }
    
    plot(plate)

}

### ### ### ### ### ### ### ### ### #
##### "levey.jennings" function #####
### ### ### ### ### ### ### ### ### #
# Jody original
levey.jennings <- function (data, std_label, blank_label, dil_high,dil_mid,dil_low,ag_list,by_var="date",subref=NULL,labels=F) {

  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  std_all <- data[grep(std_label, data$Sample),]
  blank_all <- data[grep(blank_label, data$Sample),]
  
  std_all_high <- std_all[grep(dil_high, std_all$Sample),]
  std_all_mid <- std_all[grep(dil_mid, std_all$Sample),]
  std_all_low <- std_all[grep(dil_low, std_all$Sample),]

  if (by_var=="date") {
  
    par(oma=c(3,1,1,1),mfrow=c(2,4),mar=c(5.1,4.1,4.1,2.1))
    
    for (i in 1:length(ag_list)) {
      
    # highest dilution
      
      plot (std_all_high$date, std_all_high[,ag_list[i]], pch=21,col="black",bg=alpha("darkblue",0.5),
            ylim=c(0,max(std_all_high[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="date",cex.axis=1,cex.main=1.2,cex.lab=1.4,cex=1.4,
            main=paste0(ag_list[i]," ",dil_high))
      
      if(is.null(subref)==T) {
        
        mean <- mean(std_all_high[,ag_list[i]],na.rm=T)
        sd <- sd(std_all_high[,ag_list[i]],na.rm=T) 
        
      } else if(is.null(subref)==F) {
      
        qc_set <- subref[grep(std_label, subref$Sample),]
        qc_set_high <- qc_set[grep(dil_high, qc_set$Sample),]
        
        mean <- mean(qc_set_high[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_high[,ag_list[i]],na.rm=T) 
      
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)

      if(labels==T) {text(std_all_high$date, std_all_high[,ag_list[i]],labels=std_all_high$plate,pos=3,cex=0.6)}
      
    # mid dilution
      
      plot (std_all_mid$date, std_all_mid[,ag_list[i]], pch=21,col="black",bg=alpha("darkred",0.5),
            ylim=c(0,max(std_all_high[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="date",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,cex=1.4,
            main=paste0(ag_list[i]," ",dil_mid))
      
      mean <- mean(std_all_mid[,ag_list[i]],na.rm=T)
      sd <- sd(std_all_mid[,ag_list[i]],na.rm=T)
      
      if(is.null(subref)==T) {
        
        mean <- mean(std_all_mid[,ag_list[i]],na.rm=T)
        sd <- sd(std_all_mid[,ag_list[i]],na.rm=T)
        
      } else if(is.null(subref)==F)  {
        
        qc_set <- subref[grep(std_label, subref$Sample),]
        qc_set_mid <- qc_set[grep(dil_mid, qc_set$Sample),]
        
        mean <- mean(qc_set_mid[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_mid[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      if(labels==T) {text(std_all_mid$date, std_all_mid[,ag_list[i]],labels=std_all_mid$plate,pos=3,cex=0.6)}
      
    # low dilution
      
      plot (std_all_low$date, std_all_low[,ag_list[i]], pch=21,col="black",bg=alpha("goldenrod2",0.5),
            ylim=c(0,max(std_all_high[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="date",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,cex=1.4,
            main=paste0(ag_list[i]," ",dil_low))
      
      mean <- mean(std_all_low[,ag_list[i]],na.rm=T)
      sd <- sd(std_all_low[,ag_list[i]],na.rm=T)
      
      if(is.null(subref)==T) {
        
        mean <- mean(std_all_low[,ag_list[i]],na.rm=T)
        sd <- sd(std_all_low[,ag_list[i]],na.rm=T)
        
      } else if(is.null(subref)==F)  {
        
        qc_set <- subref[grep(std_label, subref$Sample),]
        qc_set_low <- qc_set[grep(dil_low, qc_set$Sample),]
        
        mean <- mean(qc_set_low[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_low[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      if(labels==T) {text(std_all_low$date, std_all_low[,ag_list[i]],labels=std_all_low$plate,pos=3,cex=0.6)}
      
      # blanks
      
      plot (blank_all$date, blank_all[,ag_list[i]], pch=21,col="black",bg=alpha("black",0.3),
            ylim=c(0,max(blank_all[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="date",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,
            main=paste0(ag_list[i]," ","Background"))
      
      mean <- mean(blank_all[,ag_list[i]],na.rm=T)
      sd <- sd(blank_all[,ag_list[i]],na.rm=T)
      
      if(is.null(subref)==T) {
        
        mean <- mean(blank_all[,ag_list[i]],na.rm=T)
        sd <- sd(blank_all[,ag_list[i]],na.rm=T)
        
      } else if(is.null(subref)==F)  {
        
        qc_set_blank <- subref[grep(blank_label, subref$Sample),]
        
        mean <- mean(qc_set_blank[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_blank[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      text(blank_all$date, blank_all[,ag_list[i]],labels=blank_all$plate,pos=3,cex=0.8)
      
      if (is.wholenumber(i/3)) {
      par(fig= c(0,1,0,1), oma=c(0,0,0,0),mar=c(0.5,0.5,0.5,0.5),new=TRUE)
      plot(0, 0, type = "n", bty="n", xaxt = "n", yaxt = "n")
      legend("bottom", c("Mean", "+/- 2 SD", "+/- 1 SD", "Data"), xpd = TRUE, horiz = TRUE, 
             inset = c(0,0), bty="n",lty = c(1,2,3,NA), pch=c(NA,NA,NA,19),cex=1.4)
      
      par(oma=c(3,1,1,1),mfrow=c(2,4),mar=c(5.1,4.1,4.1,2.1))
        
      }
    }
  }
  
  if (by_var=="plate") {
    
    par(oma=c(3,1,1,1),mfrow=c(2,2),mar=c(5.1,4.1,4.1,2.1))
    
    for (i in 1:length(ag_list)) {
      
      # highest dilution
      
      plot (std_all_high$plate, std_all_high[,ag_list[i]], pch=21,col="black",bg=alpha("darkblue",0.5),
            ylim=c(0,max(std_all_high[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="plate",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,cex=1.4,
            main=paste0(ag_list[i]," ",dil_high))
      
      if(is.null(subref)==T) {
        
        mean <- mean(std_all_high[,ag_list[i]],na.rm=T)
        sd <- sd(std_all_high[,ag_list[i]],na.rm=T) 
        
      } else if(is.null(subref)==F) {
        
        qc_set <- subref[grep(std_label, subref$Sample),]
        qc_set_high <- qc_set[grep(dil_high, qc_set$Sample),]
        
        mean <- mean(qc_set_high[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_high[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      if(labels==T) {text(std_all_high$plate, std_all_high[,ag_list[i]],labels=std_all_high$plate,pos=3,cex=0.6)}
      
      # mid dilution
      
      plot (std_all_mid$plate, std_all_mid[,ag_list[i]], pch=21,col="black",bg=alpha("darkred",0.5),
            ylim=c(0,max(std_all_high[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="plate",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,cex=1.4,
            main=paste0(ag_list[i]," ",dil_mid))
      
      if(is.null(subref)==T) {
        
        mean <- mean(std_all_mid[,ag_list[i]],na.rm=T)
        sd <- sd(std_all_mid[,ag_list[i]],na.rm=T) 
        
      } else if(is.null(subref)==F) {
        
        qc_set <- subref[grep(std_label, subref$Sample),]
        qc_set_mid <- qc_set[grep(dil_mid, qc_set$Sample),]
        
        mean <- mean(qc_set_mid[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_mid[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      if(labels==T) {text(std_all_mid$plate, std_all_mid[,ag_list[i]],labels=std_all_mid$plate,pos=3,cex=0.6)}
      
      # low dilution
      
      plot (std_all_low$plate, std_all_low[,ag_list[i]], pch=21,col="black",bg=alpha("goldenrod2",0.5),
            ylim=c(0,max(std_all_high[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="plate",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,cex=1.4,
            main=paste0(ag_list[i]," ",dil_low))
      
      if(is.null(subref)==T) {
        
        mean <- mean(std_all_low[,ag_list[i]],na.rm=T)
        sd <- sd(std_all_low[,ag_list[i]],na.rm=T)
        
      } else if(is.null(subref)==F)  {
        
        qc_set <- subref[grep(std_label, subref$Sample),]
        qc_set_low <- qc_set[grep(dil_low, qc_set$Sample),]
        
        mean <- mean(qc_set_low[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_low[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      if(labels==T) { text(std_all_low$plate, std_all_low[,ag_list[i]],labels=std_all_low$plate,pos=3,cex=0.6)}
      
      # blanks
      
      plot (blank_all$plate, blank_all[,ag_list[i]], pch=21,col="black",bg=alpha("black",0.3),
            ylim=c(0,max(blank_all[,ag_list[i]],na.rm=T)*1.25),
            ylab="median MFI", xlab="plate",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,
            main=paste0(ag_list[i]," ","Background"))
      
      if(is.null(subref)==T) {
        
        mean <- mean(blank_all[,ag_list[i]],na.rm=T)
        sd <- sd(blank_all[,ag_list[i]],na.rm=T)
        
      } else if(is.null(subref)==F)  {
        
        qc_set_blank <- subref[grep(blank_label, subref$Sample),]
        
        mean <- mean(qc_set_blank[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_blank[,ag_list[i]],na.rm=T) 
        
      }
      
      if(is.null(subref)==T) {
        
        mean <- mean(blank_all[,ag_list[i]],na.rm=T)
        sd <- sd(blank_all[,ag_list[i]],na.rm=T)
        
      } else if(is.null(subref)==F)  {
        
        qc_set_blank <- subref[grep(blank_label, subref$Sample),]
        
        mean <- mean(qc_set_blank[,ag_list[i]],na.rm=T)
        sd <- sd(qc_set_blank[,ag_list[i]],na.rm=T) 
        
      }
      
      abline(h=mean)
      abline(h=mean+2*sd,lty=2)
      abline(h=mean-2*sd,lty=2)
      abline(h=mean+sd,lty=3)
      abline(h=mean-sd,lty=3)
      
      text(blank_all$plate, blank_all[,ag_list[i]],labels=blank_all$plate,pos=3,cex=0.8)
      
      if (is.wholenumber(i/3) | i==length(ag_list)) {
        par(fig= c(0,1,0,1), oma=c(0,0,0,0),mar=c(0.5,0.5,0.5,0.5),new=TRUE)
        plot(0, 0, type = "n", bty="n", xaxt = "n", yaxt = "n")
        legend("bottom", c("Mean", "+/- 2 SD", "+/- 1 SD", "Data"), xpd = TRUE, horiz = TRUE, 
               inset = c(0,0), bty="n",lty = c(1,2,3,NA), pch=c(NA,NA,NA,19),cex=1.4)
        
        par(oma=c(3,1,1,1),mfrow=c(2,4),mar=c(5.1,4.1,4.1,2.1))
      }
    }
  }
}


### ### ### ### ### ### ### ### ### ##
##### "normalise.loess" function #####
### ### ### ### ### ### ### ### ### ##

norm.loess <- function (data,dilutions,ag_list,ref_fit,plate_fit,plot_fit=F,lm=F,path=NULL,file_name=NULL,plate_label) {

  mean_mat <- as.data.frame(matrix(ncol=length(ag_list)+1,nrow=100))
  colnames(mean_mat) <- c("conc_no",ag_list)

  mean_mat$conc_no <- seq(1,100,1)

  for (a in 1:length(ag_list)){
   
    ag_mfi_all <- c()
    
    for (j in 1:length(ref_fit)) {
      
      ag_mfi_all <- cbind(ag_mfi_all,ref_fit[[j]][,ag_list[a]])
      
    }
    
    mean_mat[,ag_list[a]] <- apply(ag_mfi_all,1,mean,na.rm=T)
    
  }
  
  plates_norm <- list()
  
  aux.plate <- as.data.frame(data[,c("Location","Sample",ag_list)])
  aux.fit <- plate_fit
  
  plates_norm <- matrix(ncol=ncol(aux.plate),nrow=nrow(aux.plate))
  colnames(plates_norm) <- colnames(aux.plate)
  plates_norm[,1] <- as.character(aux.plate[,1])
  plates_norm[,2] <- as.character(aux.plate[,2])
  
  # par(mfrow=c(3,3))
    
  for (j in ag_list) {
      
    stnd <- c()
      
    stnd$mean <- mean_mat[,j]
    stnd$plate <- aux.fit[,j]
    stnd$delt <- stnd$plate-stnd$mean
    stnd <- as.data.frame(stnd)
    stnd <- stnd[which(is.na(stnd$plate)==F),]

    fit1 <- loess(delt~mean,stnd)
    fit_x <- seq(0,max(stnd$mean)*1.5,0.5)
    fit1_pred <- predict(fit1,fit_x,se=T)

    # lm fit

    fit2 <- lm(delt~mean,stnd)
    fit2_x <- as.data.frame(seq(0,max(stnd$mean)*1.5,0.5))
    colnames(fit2_x) <- "mean"
    fit2_x <- data.frame(mean = fit2_x$mean)
    fit2_pred <- predict.lm(fit2,fit2_x,se.fit=T)

    fit2_mat <- cbind(fit2_x,fit2_pred$fit,fit2_pred$se.fit)
    colnames(fit2_mat) <- c("x","y","se")
    fit2_mat$li <- fit2_mat$y - 2*fit2_mat$se
    fit2_mat$ui <- fit2_mat$y + 2*fit2_mat$se
      
    if (plot_fit==T) {
        
      ymax <- max(na.omit(fit1_pred$fit),fit2_pred$se.fit)
      ymin <- min(na.omit(fit1_pred$fit),fit2_pred$se.fit)
          
      plot(unlist(fit2_x),fit2_pred$fit,col="white",main=j,xlab="mean MFI (ref plates)",ylab="plate MFI - mean MFI (ref)",ylim=c(ymin,ymax))

      points(stnd$mean,stnd$delt,cex=0.8)
      lines(fit_x,fit1_pred$fit, col="red",lwd=2)
      
      if (lm==T) {
        lines(unlist(fit2_x),fit2_pred$fit, col="blue",lwd=2)
        lines(unlist(fit2_x),fit2_pred$fit + 2*fit2_pred$se.fit,lty=3,col="blue")
        lines(unlist(fit2_x),fit2_pred$fit - 2*fit2_pred$se.fit,lty=3,col="blue")
      }
      
    }

    if (lm==T) {
      
      for (k in 1:nrow(aux.plate)) {
        plates_norm[k,j] <- aux.plate[k,j]-fit2_mat[which(round(fit2_mat$x)==round(aux.plate[k,j]))[1],"y"]
        if (plates_norm[k,j]<0) plates_norm[k,j] <- rnorm(1,mean=10,sd=3)
        
      }
      
    }
      
    if (lm==F) {
      
      for (k in 1:nrow(aux.plate)) {
      
        min_adj <- stnd[which(stnd$mean==min(stnd$mean)),"delt"]/stnd[which(stnd$mean==min(stnd$mean)),"plate"]
        max_adj <- stnd[which(stnd$mean==max(stnd$mean)),"delt"]/stnd[which(stnd$mean==max(stnd$mean)),"plate"]

        if (abs(min_adj)>1) min_adj1 <- 0 else min_adj1 <- min_adj
        if (abs(max_adj)>1) max_adj1 <- 0 else max_adj1 <- max_adj
          
        if(is.na(aux.plate[k,j])==T) plates_norm[k,j] <- aux.plate[k,j]
        else if(is.na(aux.plate[k,j])==F & aux.plate[k,j]<=max(stnd$mean) & aux.plate[k,j]>=min(stnd$mean)) plates_norm[k,j] <- aux.plate[k,j]-predict(fit1,aux.plate[k,j])
        else if(is.na(aux.plate[k,j])==F & aux.plate[k,j]>max(stnd$mean)) plates_norm[k,j] <- aux.plate[k,j]-max_adj1*aux.plate[k,j]
        else if(is.na(aux.plate[k,j])==F & aux.plate[k,j]<min(stnd$mean)) plates_norm[k,j] <- aux.plate[k,j]-min_adj1*aux.plate[k,j]       

        if (plates_norm[k,j]<=0 & is.na(aux.plate[k,j])==F) plates_norm[k,j] <- rnorm(1,mean=10,sd=3)
          
      }
        
    }
      
  }
  
  if(is.null(path)==F&is.null(file_name)==T) write.csv(plates_norm,file=paste0(path,"norm_plate.csv"),row.names=F)
  else if(is.null(path)==F&is.null(file_name)==F) write.csv(plates_norm,file=paste0(path,file_name),row.names=F)
  
  
  return(plates_norm)    

}


### ### ### ### ### ### ### ### ### #
##### "normalise.frac" function #####
### ### ### ### ### ### ### ### ### #

norm.frac <- function (data,std_label,ag_list,path=NULL) {
 
  # create data objects for each plate to store normalised values
  
  sample_frac <- list()
  sample_PP <- list()
  
  for (i in 1:length(data)){
    sample_PP[[i]] <- data.frame(matrix(nrow=nrow(data[[i]]), ncol=ncol(data[[i]])))
    colnames(sample_PP[[i]]) <- colnames(data[[i]])
    sample_PP[[i]][,1] <- as.character(data[[i]][,1])
    sample_PP[[i]][,2] <- as.character(data[[i]][,2])
    
    sample_frac[[i]] <- data.frame(matrix(nrow=nrow(data[[i]]), ncol=ncol(data[[i]])))
    colnames(sample_frac[[i]]) <- colnames(data[[i]])
    sample_frac[[i]][,1] <- as.character(data[[i]][,1])
    sample_frac[[i]][,2] <- as.character(data[[i]][,2])
  }
  
  # create matrix of mean standard curve values across all plates
  
  std_all <- data[[1]][grep(std_label, data[[1]][,"Sample"]),]
  
  for (i in 2:length(data)){
    std_all <- rbind(std_all,data[[i]][grep(std_label, data[[i]][,"Sample"]),])
  }
  
  std_mid_all <- std_all[grep("1/250", std_all$Sample),]
  
  std_mid_avg <- vector()
  
  std <- list()
  
  for (i in 1:length(data)){
    
    std[[i]] <- data[[i]][grep("PP", data[[i]]$Sample)[6:1],]
    
    for (j in 1:length(ag_list)){ 
      
      std_mid_avg[j] <- mean(std_mid_all[,ag_list[j]])
      
      frac <- log(std[[i]][4,ag_list[j]])/log(std_mid_avg[j])
      
      for (k in 1:nrow(data[[i]])){
        
        if (is.na(data[[i]][k,ag_list[j]])==T){
          next
        }
        
        sample_PP[[i]][k, ag_list[j]] <- data[[i]][k,ag_list[j]]/std[[i]][4,ag_list[j]]
        sample_frac[[i]][k, ag_list[j]] <- data[[i]][k,ag_list[j]]/frac
        
      }
    }
    
    sample_PP[[i]][,1:2] <- sample_frac[[i]][,1:2] <- data[[i]][,1:2]
    
    if(is.null(path)==F) {
      write.csv(sample_PP[[i]],file=paste0(path,"norm_PP",i,".csv"),row.names=F)
      write.csv(sample_frac[[i]],file=paste0(path,"norm_frac",i,".csv"),row.names=F)
    }
  }
  
  output <- list(sample_PP,sample_frac)
  names(output) <- c("PP","FP")
  
  return(output)
  
}

### ### ### ### ### ### ### ### #
##### "plate2list" function #####
### ### ### ### ### ### ### ### #

# converting plate plans to importable sample id lists and vice versa
# appends sample IDs by column in descending rows (i.e., A1-H1, A2-H2 etc)

plate2list <- function(data,no_plates,names=NULL,path){
  
  id_list <- list()
  
  for (i in 1:no_plates) {
    
    start <- grep("A",data[,1])[i]
    
    for (j in 1:12) {
    
      id_col <- as.vector(data[start:(start+7),j+1])
      
      if (j==1) id_list[[i]] <- id_col else id_list[[i]] <- c(id_list[[i]],id_col)
    
    }
    
    if (is.null(names)==T) {
      write.table(id_list[[i]],file=paste0(path,"sample_ids",i,".csv"),sep=",",row.names=F,col.names=F)
      write.table(id_list[[i]],file=paste0(path,"sample_ids",i,".txt"),sep=",",row.names=F,col.names=F)
      
    } else {
      write.table(id_list[[i]],file=paste(path,names[[i]],".csv"),sep=",",row.names=F,col.names=F)
      write.table(id_list[[i]],file=paste(path,names[[i]],".txt"),sep=",",row.names=F,col.names=F)
    }
  }
  
}


### ### ### ### ### ### ### ### ### ##
##### "ag column print" function #####
### ### ### ### ### ### ### ### ### ##

ag.col.print <- function(data) {
  
  for (i in 1:length(data)) {
    
    print(names(data)[[i]])
    print(names(data[[i]]))
    
  }
  
}


### ### ### ### ### ### ### ### ### ### #
##### "ag column mismatch" function #####
### ### ### ### ### ### ### ### ### ### #

ag.match.check <- function(data, plateprint=F) {

  ag_all <- list()
  ag_unlist <- c()
  
  for (i in 1:length(data)) {
    
    ag_all[[i]] <- colnames(data[[i]])
    ag_unlist <- c(ag_unlist,colnames(data[[i]]))
    
  }
  
  ag_same <- Reduce(intersect, ag_all)
  
  ag_mismatch <- unique(ag_unlist[-which(ag_unlist %in% ag_same)])
  
  print("Antigens common across all plates:")
  print(ag_same)
  print("")
  
  if (length(ag_mismatch)>0) {
    
    print("Antigens NOT common across all plates. Check spelling.")
    print(ag_mismatch)

  } else print("All antigens (column headings) are matching in all plates and can be merged")
  
  if (length(ag_mismatch)>0&plateprint==T) {
    
    ag_miss <- list()
    
    print("")
        print("Missing antigen headings by plate:")

    for (i in 1:length(data)) {
      
      ag_miss[[i]] <- ag_mismatch[which(ag_mismatch %in% ag_all[[i]]==F)]
      
      print(c(names(data)[[i]],ag_miss[[i]]))
      
    }
    
  }
  
  return(ag_mismatch)

}


### ### ### ### ### ### ### ### ### 
##### "ag column add" function #####
### ### ### ### ### ### ### ### ### 

# add missing columns to plate after ag.match.check function is run

ag.column.add <- function(data,ag_mismatch) {
  
  for (i in 1:length(data)) {
    
    for (j in ag_mismatch) {
      
      if (j %in% names(data[[i]])) data[[i]][,j] <- data[[i]][,j] else data[[i]][,j] <- NA
      
    }
    
  }
  
  return(data)
  
}



### ### ### ### ### ### ### 
##### Unused functions #####
### ### ### ### ### ### ### 

bead.check <- function(data, min=30, ag_list) {

  # Checks if bead counts fall below a specified threshold by antigen
  # Prints list of samples per antigen
    
  for (i in ag_list) {
    
    bead_low <- which(data[,i]<min)
    
    if (length(bead_low)>0) {
      
      bead_low_table <- as.data.frame(cbind(as.character(data[bead_low,"Location"]),as.character(data[bead_low,"Sample"]),data[bead_low,i]))
      names(bead_low_table) <- c("well","sample_id","beads")
      
      print(i)
      print(bead_low_table)
      
    } else {
      print(i)
      print("All wells are above the minimum bead count for this antigen")
    }
    
  }
  
}


read.batch.test <- function (path, inc_date=F, inc_plate=F, save_path=NULL) {
  
  setwd(path)
  temp <- list.files(pattern="*.csv")
  plate_lab <- substr(temp,1,nchar(temp)-4)
  plate <- list()
  
  for(i in 1:length(temp)) { 
    
    plate_raw <- read.csv(temp[i])
    startrow <- grep("DataType",plate_raw[,1])[1]
    plate[[i]] <- read.csv(temp[i], skip=startrow+1)
    
  }
}



read.plate <- function (path,file_name,inc_date=F) {
  
  data_raw <-  read.csv(paste0(path,file_name))
  startrow <- grep("DataType",data_raw[,1])[1]
  data <- read.csv(paste0(path,file_name), skip=startrow+1)
  
  section <- grep("1,",data$Location, invert=T)[1]
  
  data <- data[1:(section-1),-ncol(data)]
  
  data[,3:ncol(data)] <- sapply(data[,3:ncol(data)],as.character)
  data[,3:ncol(data)] <- sapply(data[,3:ncol(data)],as.numeric)
  
  if (inc_date==T) {
    data$date <- data_raw[grep("Date",data_raw[,1]),2]
    data$date <- as.Date(data$date,"%m/%d/%Y")
  }
  
  return(data)
  
}

read.plate.beads <- function (path,file_name,inc_date=F) {
  
  data_raw <-  read.csv(paste0(path,file_name))
  startrow <- grep("DataType",data_raw[,1])[1]
  data <- read.csv(paste0(path,file_name), skip=startrow+1)
  
  section1 <- grep("Count",data$Sample)[1]
  section2 <- grep("Avg Net MFI",data$Sample)
  
  data <- data[(section1+2):(section2-2),-ncol(data)]
  
  data[,3:ncol(data)] <- sapply(data[,3:ncol(data)],as.character)
  data[,3:ncol(data)] <- sapply(data[,3:ncol(data)],as.numeric)
  
  if (inc_date==T) {
    data$date <- data_raw[grep("Date",data_raw[,1]),2]
    data$date <- as.Date(data$date,"%m/%d/%Y")
  }
  
  return(data)
  
}



##### "plot.std.curve" function
# Plots standard curve for 1 antigen
plot.std.curve1 <- function(data, plotfile="plot1.pdf", antigen, dilutions, ref_curve, save_file=F, path=getwd()) {
  
  ref_curve <- ref_curve[,which(colnames(ref_curve)==antigen)]
  
  if ((antigen %in% colnames(data))==F) {
    
    print("Antigen name entered does not exist in dataset - check spelling")
    
  } else {
    
    max.set <- data[,which(colnames(data)==antigen)]
    max.set <- cbind(max.set[is.finite(max.set)],ref_curve)
    max.set <- max.set[!is.na(max.set)]
    
    max <- max(max.set)
    
    conc <- sapply(strsplit(dilutions, split = "/"),
                   function(x) as.numeric(x[1]) / as.numeric(x[2]))
    
    # saves file #
    
    if (save_file==T) {
      
      pdf(paste0(path,plotfile))
      
      par(mfrow=c(1,1))
      
      plot(log(conc),data[,which(colnames(data)==antigen)],
           type="o",pch=19,col="red",axes=F,bty='L',
           ylab="MFI",xlab="Dilution factor",
           ylim=c(0,max*1.2),xlim=c(min(log(conc)),max(log(conc))+1),
           main=paste(colnames(data[which(colnames(data)==antigen)]),"Standard curve",sep='\n'))
      axis(1,at=c(log(conc),max(log(conc))+1),labels=c(dilutions,""),cex.axis=0.9)
      axis(2,cex.axis=0.9)
      lines(log(conc),ref_curve,type="o",col="black",pch=19)
      lines(log(conc),data[,which(colnames(data)==antigen)],type="o",col="red",pch=19)
      
      legend("topleft",cex=0.9,pch=19,bty="n",
             c("Reference curve","Plate curve"),col=c("black","red"))
      
      dev.off()
    }
    
    # plots for viewing in R #
    
    par(mfrow=c(1,1))
    
    plot(log(conc),data[,which(colnames(data)==antigen)],
         type="o",pch=19, col="red",axes=F,bty='L',
         ylim=c(0,max*1.2),xlim=c(min(log(conc)),max(log(conc))+1),
         ylab="MFI",xlab="Dilution factor",
         main=paste(colnames(data[which(colnames(data)==antigen)]),"Standard curve",sep='\n'))
    axis(1,at=c(log(conc),max(log(conc))+1),labels=c(dilutions,""),cex.axis=0.8)
    axis(2,cex.axis=0.8)
    lines(log(conc),ref_curve,type="o",col="black",pch=19)
    lines(log(conc),data[,which(colnames(data)==antigen)],type="o",col="red",pch=19)
    
    legend("topleft",cex=0.9,pch=19,bty="n",
           c("Reference curve","Plate curve"),col=c("black","red"))
    
  }
  
}



# Plots standard curve for multiple antigens in 3x3 matrix
plot.std.curve2 <- function(data,plotfile="plot1.pdf",antigen.list,dilutions,save_file=F,ref_curve,path=getwd(),ref_label="Reference",plate_label="Plate") {
  
  if (save_file==T) {
    
    pdf(paste0(path,plotfile))
    
    # par(mfrow=c(3,3))
    
    for (j in 1:length(antigen.list)){
      
      if ((antigen.list[j] %in% colnames(data))==F) {
        next
      } else {
        
        ref_curve1 <- ref_curve[,which(names(ref_curve)==antigen.list[j])]
        
        max.set <- data[,which(colnames(data)==antigen.list[j])]
        max.set <- cbind(max.set[is.finite(max.set)],ref_curve1)
        max.set <- max.set[!is.na(max.set)]
        
        max <- max(max.set)
        
        conc <- sapply(strsplit(dilutions, split = "/"),
                       function(x) as.numeric(x[1]) / as.numeric(x[2]))
        
        plot(log(conc),data[,which(colnames(data)==antigen.list[j])], type="o", pch=19, col="red",
             ylab="MFI", xlab="Dilution factor",axes=F,bty='L', 
             ylim=c(0,max*1.2), xlim=c(min(log(conc)),max(log(conc))+1),
             main=paste(colnames(data[which(colnames(data)==antigen.list[j])]),"Standard curve", sep='\n'),cex.main=0.9 )
        axis(1, at=c(log(conc),max(log(conc))+1) , labels=c(dilutions,""), cex.axis=0.8)
        axis(2, cex.axis=0.8)
        lines(log(conc),ref_curve1,type="o",col="black",pch=19)
        lines(log(conc),data[,which(colnames(data)==antigen.list[j])],type="o",col="red",pch=19)
        
        legend("topleft", cex=0.8, pch=19,bty="n",y.intersp=0.8,
               c(ref_label,plate_label),col=c("black","red"))
      }
      
    }
    
    dev.off()
    
  }
  
  # par(mfrow=c(3,3))
  
  for (i in 1:length(antigen.list)){
    
    if ((antigen.list[i] %in% colnames(data))==F) {
      print(paste0(" '",antigen.list[i],"'"," does not exist in dataset - check spelling "))
      next
    } else {
      
      ref_curve1 <- ref_curve[,which(colnames(ref_curve)==antigen.list[i])]
      
      max.set <- data[,which(colnames(data)==antigen.list[i])]
      max.set <- cbind(max.set[is.finite(max.set)],ref_curve1)
      max.set <- max.set[!is.na(max.set)]
      
      max <- max(max.set)
      
      conc <- sapply(strsplit(dilutions, split = "/"),
                     function(x) as.numeric(x[1]) / as.numeric(x[2]))
      
      # plots for viewing in R #
      
      plot(log(conc),data[,which(colnames(data)==antigen.list[i])],
           type="o",pch=19,col="red",axes=F,bty='L',cex.main=0.9,
           ylab="MFI", xlab="Dilution factor", 
           ylim=c(0,max*1.2), xlim=c(min(log(conc)),max(log(conc))+1),
           main=paste(colnames(data[which(colnames(data)==antigen.list[i])]),"Standard curve",sep='\n'))
      axis(1,at=c(log(conc),max(log(conc))+1),labels=c(dilutions,""),cex.axis=0.8)
      axis(2,cex.axis=0.8)
      lines(log(conc),ref_curve1,type="o",col="black",pch=19)
      lines(log(conc),data[,which(colnames(data)==antigen.list[i])],type="o",col="red",pch=19)
      
      legend("topleft", cex=0.9, pch=19,bty="n",
             c(ref_label,plate_label),col=c("black","red"))
    }
  }
}



# Plots standard curve for multiple plates (up to 5) and multiple antigens in 3x3 matrix
# data list of standards from multiple plates (maximum is 5 at a time)
plot.std.curve3 <- function(data,plotfile="plot1.pdf",antigen.list,dilutions,save_file=F,path=getwd(),negs=NULL,blanks=NULL,plate_labels=NULL) {
  
  if (save_file==T) {
    
    pdf(paste0(path,plotfile))
    
    # par(mfrow=c(3,2))
    
    for (i in 1:length(antigen.list)){
      
      curves <- list()
      max.set <- c()
      plate_labels_def <- c()
      
      for (j in 1:length(data)){
        
        if ((antigen.list[i] %in% colnames(data[[j]]))==F) {
          print(paste0(" '",antigen.list[i],"'"," does not exist in dataset ",j," - check spelling "))
          next
        } else {      
          curves[[j]] <- data[[j]][,which(names(data[[j]])==antigen.list[i])]
          max.set <- c(max.set,curves[[j]])
          plate_labels_def <- c(plate_labels_def,paste0("plate",j))
        }
      }
      
      max.set <- max.set[is.finite(max.set)]
      max <- max(max.set)
      
      conc <- sapply(strsplit(dilutions, split = "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
      cols1 <- c("blue2","red3","green4","purple3","goldenrod2")
      
      plot(log(conc),curves[[1]], type="o", pch=19, col=cols1[[1]],ylab="MFI", xlab="Dilution factor",axes=F,bty='L', 
           ylim=c(0,max*1.2), xlim=c(min(log(conc)),max(log(conc))+1),
           main=paste(antigen.list[i],"Standard curve", sep='\n'),cex.main=0.9 )
      axis(1, at=c(log(conc),max(log(conc))+1) , labels=c(dilutions,""), cex.axis=0.8)
      axis(2, cex.axis=0.8)
      
      for (j in 2:length(data)) { lines(log(conc),curves[[j]],type="o",col=cols1[[j]],pch=19) }
      
      if (length(plate_labels)==0) { 
        legend("topleft",cex=0.6,pch=19,bty="n",y.intersp=0.8,plate_labels_def,col=cols1[1:length(data)]) 
      } else { legend("topleft",cex=0.8,pch=19,bty="n",y.intersp=0.6,plate_labels,col=cols1[1:length(data)]) }
    }
    
    dev.off()
    
  }
  
  # par(mfrow=c(3,3))
  
  for (i in 1:length(antigen.list)){
    
    curves <- list()
    max.set <- c()
    plate_labels_def <- c()
    
    for (j in 1:length(data)){
      
      if ((antigen.list[i] %in% colnames(data[[j]]))==F) {
        print(paste0(" '",antigen.list[i],"'"," does not exist in dataset ",j," - check spelling "))
        next
      } else {
        curves[[j]] <- data[[j]][,which(names(data[[j]])==antigen.list[i])]
        max.set <- c(max.set,curves[[j]])
        plate_labels_def <- c(plate_labels_def,paste0("plate",j))
      }
      
    }
    
    max.set <- max.set[is.finite(max.set)]
    max <- max(max.set)
    
    conc <- sapply(strsplit(dilutions, split = "/"),function(x) as.numeric(x[1]) / as.numeric(x[2]))
    cols1 <- c("green4","blue2","red3","purple3","goldenrod2")
    
    # plots for viewing in R #
    
    plot(log(conc),curves[[1]], type="o", pch=19, col=cols1[[1]],
         ylab="MFI", xlab="Dilution factor",axes=F,bty='L', 
         ylim=c(0,max*1.2), xlim=c(min(log(conc)),max(log(conc))+1),
         main=paste(antigen.list[i],"Standard curve", sep='\n'),cex.main=0.9 )
    axis(1, at=c(log(conc),max(log(conc))+1) , labels=c(dilutions,""), cex.axis=0.8)
    axis(2, cex.axis=0.8)
    if (is.null(blanks)==F) abline(h=mean(blanks[,antigen.list[i]]),lty=3,lwd=1.5)
    if (is.null(negs)==F) abline(h=mean(negs[,antigen.list[i]]),col="red",lty=3,lwd=1.5)
    
    for (j in 2:length(data)) { lines(log(conc),curves[[j]],type="o",col=cols1[[j]],pch=19) }
    
    if (length(plate_labels)==0) { 
      legend("topleft",cex=0.8,pch=19,bty="n",y.intersp=0.4,plate_labels_def,col=cols1[1:length(data)]) 
    } else if(is.null(negs)==F&is.null(blanks)==F) {
      legend("topleft",cex=0.8,pch=c(19,19,NA,NA),lty=c(NA,NA,3,3),bty="n",y.intersp=0.4,
             c(plate_labels,"Neg","Blank"),col=c(cols1[1:length(data)],"red","black"),seg.len=1)
    } else {
      legend("topleft",cex=0.5,pch=19,bty="n",y.intersp=0.4,
             c(plate_labels),col=c(cols1[1:length(data)]),seg.len=1)   
    }
  }
}

plot.std.curve3_single_ag <- function(data,antigen,dilutions,negs=NULL,blanks=NULL,plate_labels=NULL) {
  curves <- list()
  max.set <- c()
  plate_labels_def <- c()
  
  for (j in 1:length(data)){
    
    if ((antigen %in% colnames(data[[j]]))==F) {
      print(paste0(" '",antigen,"'"," does not exist in dataset ",j," - check spelling "))
      next
    } else {
      curves[[j]] <- data[[j]][,which(names(data[[j]])==antigen)]
      max.set <- c(max.set,curves[[j]])
      plate_labels_def <- c(plate_labels_def,paste0("plate",j))
    }
    
  }
  
  max.set <- max.set[is.finite(max.set)]
  max <- max(max.set)
  
  conc <- sapply(strsplit(dilutions, split = "/"),function(x) as.numeric(x[1]) / as.numeric(x[2]))
  cols1 <- c("green4","blue2","red3","purple3","goldenrod2")
  # cols1 <- brewer.pal(length(data),"Paired")
  
  # plots for viewing in R #
  
  plot(log(conc),curves[[1]], type="o", pch=19, col=cols1[[1]],
       ylab="MFI", xlab="Dilution factor",axes=F,bty='L', 
       ylim=c(0,max*1.2), xlim=c(min(log(conc)),max(log(conc))+1),
       main=paste(antigen,"Standard curve", sep='\n'),cex.main=0.9 )
  axis(1, at=c(log(conc),max(log(conc))+1) , labels=c(dilutions,""), cex.axis=0.8)
  axis(2, cex.axis=0.8)
  if (is.null(blanks)==F) abline(h=mean(blanks[,antigen]),lty=3,lwd=1.5)
  if (is.null(negs)==F) abline(h=mean(negs[,antigen]),col="red",lty=3,lwd=1.5)
  
  for (j in 2:length(data)) { lines(log(conc),curves[[j]],type="o",col=cols1[[j]],pch=19) }
  
  if (length(plate_labels)==0) { 
    legend("topleft",cex=0.8,pch=19,bty="n",y.intersp=0.4,plate_labels_def,col=cols1[1:length(data)]) 
  } else if(is.null(negs)==F&is.null(blanks)==F) {
    legend("topleft",cex=0.8,pch=c(19,19,NA,NA),lty=c(NA,NA,3,3),bty="n",y.intersp=0.4,
           c(plate_labels,"Neg","Blank"),col=c(cols1[1:length(data)],"red","black"),seg.len=1)
  } else {
    legend("topleft",pch=19,bty="n",
           c(plate_labels),col=c(cols1[1:length(data)]))   
  }
}



# "cutoff" function
  # calculates seropositivity cutoff based on finite mixture model using both raw MFI and log transformed MFI using normalmixEM

cutoff <- function (data,breaks=20,ag_list,neg_label=NULL,std_label=NULL,blank_label=NULL,no_sd=2,print=F) {
  
  if(is.null(std_label)==F) data <- data[-grep(std_label,data$Sample),] else data <- data
  
  if(is.null(neg_label)==F) {
    negs <- data[grep(neg_label,data$Sample),]
    data <- data[-grep(neg_label,data$Sample),]
  }
  
  if(is.null(blank_label)==F) {
    blanks <- data[grep(blank_label,data$Sample),]
    data <- data[-grep(blank_label,data$Sample),]
  }
  
  # par(mfrow=c(3,4))
  
  len <- length(ag_list)
  
  cut_all <- vector(mode="numeric",length=len)
  names(cut_all) <- ag_list
  log_cut_all <- vector(mode="numeric",length=len)
  names(log_cut_all) <- ag_list
  
  for (i in ag_list){
    
    print("")
    print(paste0("Calculating cutoff for ",i))
    
    plot.data <- na.omit(data[,i])
    
    if(is.null(neg_label)==F) {
      neg.data <- negs[,i]
      neg.mean <- mean(neg.data,na.rm=T)
    }
    
    if(is.null(blank_label)==F) {
      blank.data <- blanks[,i]
      blank.mean <- mean(blank.data)
    }
    
    labels1 <- as.numeric(c(10,100,1000,10000,100000))
    
    print("MFI")
    fit <- normalmixEM(na.omit(plot.data),2)
    
    if(print==T) {
      print("FMM 2 component (MFI)")
      summary(fit)
    }
    
    mu1 <- fit$mu[1]
    mu2 <- fit$mu[2]
    sig1 <- fit$sigma[1]
    sig2 <- fit$sigma[2]
    
    min_comp1 <- which(fit$mu == min(fit$mu))
    
    cut <- cut_all[i] <- fit$mu[min_comp1]+no_sd*sqrt(fit$sigma[min_comp1])
    
    plot_x <- seq(0,max(plot.data,na.rm=T),1)
    
    gauss1a <- dnorm(plot_x,mu1,sig1)
    gauss2a <- dnorm(plot_x,mu2,sig2)
    
    print("logMFI")
    log_fit <- normalmixEM(na.omit(log(plot.data[plot.data>1.0001])),2)
    
    if(print==T) {
      print("FMM 2 component (log MFI)")
      summary(log_fit)
    }
    
    log_mu1 <- log_fit$mu[1]
    log_mu2 <- log_fit$mu[2]
    log_sig1 <- log_fit$sigma[1]
    log_sig2 <- log_fit$sigma[2]
    
    min_comp2 <- which(log_fit$mu == min(log_fit$mu))
    
    log_cut <- log_cut_all[i] <- log_fit$mu[min_comp2]+no_sd*sqrt(log_fit$sigma[min_comp2])
    
    gauss1b <- dnorm(log(plot_x),log_mu1,log_sig1)
    gauss2b <- dnorm(log(plot_x),log_mu2,log_sig2)
    
    hist(plot.data,breaks=breaks, main=i,xlab="MFI",freq = F,
         ylim=c(0,max(gauss1a,gauss2a)))
    
    if(is.null(neg_label)==F) abline(v=neg.mean,col="springgreen4",lwd=2)
    #abline(v=neg.data,col=alpha("springgreen4",0.5),lty=3,lwd=0.2)
    if(is.null(blank_label)==F) abline(v=blank.mean,col="darkgrey",lwd=1)
    #abline(v=blank.data,col=alpha("darkgrey",0.5),lty=3,lwd=0.2)
    abline(v=cut,col="red",lwd=2)
    lines(plot_x,gauss1a,col="dodgerblue4")
    lines(plot_x,gauss2a,col="dodgerblue4")
    
    if(is.null(neg_label)==F&is.null(blank_label)==F) {
      
      legend("topleft",c(paste0("cutoff (MFI): ",round(cut,0)),paste0("neg ctrls: ",round(neg.mean,0)),
                         paste0("background: ",round(blank.mean,0))),lty=c(1,1,1),
             col=c("red","springgreen4","darkgrey"),cex=0.9,bty="n",y.intersp=0.2,x.intersp=0.2,seg.len=0.5)
    } else {
      
      legend("topleft",paste0("cutoff (MFI): ",round(cut,0)),lty=1,col="red",cex=0.9,bty="n",y.intersp=0.2,x.intersp=0.2,seg.len=0.5,text.col="red")
      
    }
    
    hist(log(plot.data[plot.data>1.0001]),breaks=breaks, main=i,xaxt='n',xlab="MFI (log scale)",freq=F,
         ylim=c(0,max(gauss1b,gauss2b)))
    axis(side=1, at=log(labels1), labels=labels1)
    if(is.null(neg_label)==F) abline(v=log(neg.mean),col="springgreen4",lwd=1)
    #abline(v=log(neg.data),col=alpha("springgreen4",0.5),lty=3,lwd=0.2)
    if(is.null(blank_label)==F) abline(v=log(blank.mean),col="darkgrey",lwd=1)
    #abline(v=log(blank.data),col=alpha("darkgrey",0.5),lty=3,lwd=0.2)
    abline(v=log_cut,col="red",lwd=2)
    lines(log(plot_x),gauss1b,col="dodgerblue4")
    lines(log(plot_x),gauss2b,col="dodgerblue4")
    
    if(is.null(neg_label)==F&is.null(blank_label)==F) {
      
      legend("topleft",c(paste0("cutoff (MFI): ",round(exp(log_cut),0)),paste0("neg ctrls: ",round(neg.mean,0)),
                         paste0("background: ",round(blank.mean,0))),lty=c(1,1,1),
             col=c("red","springgreen4","darkgrey"),cex=0.9,bty="n",y.intersp=0.2,x.intersp=0.2,seg.len=0.5)
      
    } else {
      
      legend("topleft",paste0("cutoff (MFI): ",round(exp(log_cut),0)),lty=1,col="red",cex=0.9,bty="n",y.intersp=0.2,x.intersp=0.2,seg.len=0.5,text.col="red")
    }
  }
  
  log_cut_all <- exp(log_cut_all)
  
  output <- list(cut_all,log_cut_all)
  names(output) <- c("MFI","logMFI")
  
  return(output)
  
}


# "seropos" function
  # Using the cutoff values calculated in cutoff function
  # Assigns seropositivity values for all antigens and adds to an existing dataframe

seropos <- function(data,cutoffs,ag_list) {
  
  new_data <- data
  
  for (i in ag_list) {
    
    var_new1 <- paste0(i,"_pos") 
    #var_new2 <- paste0(i,"_poslog") 
    
    new_data[,var_new1] <- NA
    new_data[,var_new1][new_data[,i]>=cutoffs[[i]]] <- 1
    new_data[,var_new1][new_data[,i]<cutoffs[[i]]]  <- 0
    
    #new_data[,var_new2] <- NA
    #new_data[,var_new2][new_data[,i]>=cutoffs[[2]][[i]]]  <- 1
    #new_data[,var_new2][new_data[,i]<cutoffs[[2]][[i]]]  <- 0
    
  }
  
  return(new_data)
  
}


### New for Mika

## Recognise standard names that end in
identify.standards <- function(data, dilutions) {
  
  # Accept either one plate data frame or a list of plates
  if (is.data.frame(data)) {
    data <- list(plate1 = data)
  }
  
  dilutions <- trimws(as.character(dilutions))
  
  # Escape dilution strings before using them in regular expressions
  escape.regex <- function(x) {
    gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
  }
  
  dilution.pattern <- paste(
    escape.regex(dilutions),
    collapse = "|"
  )
  
  results <- list()
  
  for (plate.name in names(data)) {
    
    plate <- data[[plate.name]]
    
    if (!"Sample" %in% names(plate)) next
    
    sample.names <- trimws(as.character(plate$Sample))
    
    # Find a user-specified dilution at the end of each sample name
    matched.dilution <- regmatches(
      sample.names,
      regexpr(
        paste0("(", dilution.pattern, ")$"),
        sample.names,
        perl = TRUE
      )
    )
    
    has.dilution <- matched.dilution != ""
    
    # Remove the dilution and its preceding separator to obtain the name
    standard.name <- sample.names
    standard.name[has.dilution] <- sub(
      paste0("[_ -]?", dilution.pattern, "$"),
      "",
      sample.names[has.dilution],
      perl = TRUE
    )
    
    results[[plate.name]] <- data.frame(
      Plate = plate.name,
      Sample = sample.names,
      Standard = standard.name,
      Dilution = matched.dilution,
      IsStandardCandidate = has.dilution,
      stringsAsFactors = FALSE
    )
  }
  
  do.call(rbind, results)
}


### New functions

## ### ### ### ### ### ### ###
##### identify.standards #####
## ### ### ### ### ### ### ###
identify.standards <- function(data, dilutions) {

# Accept either one plate or a list of plates
if (is.data.frame(data)) {
  data <- list(plate1 = data)
}

if (is.null(names(data)) || any(names(data) == "")) {
  names(data) <- paste0("plate", seq_along(data))
}

dilutions <- unique(trimws(as.character(dilutions)))
dilutions <- dilutions[!is.na(dilutions) & dilutions != ""]

results <- lapply(names(data), function(plate_name) {
  
  plate <- data[[plate_name]]
  
  if (!"Sample" %in% names(plate)) {
    return(NULL)
  }
  
  sample_names <- trimws(as.character(plate$Sample))
  
  matched_dilution <- rep(NA_character_, length(sample_names))
  standard_name <- rep(NA_character_, length(sample_names))
  
  for (i in seq_along(sample_names)) {
    
    sample <- sample_names[i]
    
    if (is.na(sample) || sample == "") next
    
    # Find dilutions that occur at the end of this sample name
    matches <- dilutions[
      vapply(dilutions, function(x) endsWith(sample, x), logical(1))
    ]
    
    if (length(matches) == 0) next
    
    # Use the longest match if more than one dilution matches
    dilution <- matches[which.max(nchar(matches))]
    
    # Remove the dilution suffix and its preceding separator
    prefix <- substr(
      sample,
      1,
      nchar(sample) - nchar(dilution)
    )
    
    prefix <- sub("[_ -]+$", "", prefix)
    
    if (prefix == "") next
    
    matched_dilution[i] <- dilution
    standard_name[i] <- prefix
  }
  
  data.frame(
    Plate = plate_name,
    Sample = sample_names,
    Standard = standard_name,
    Dilution = matched_dilution,
    IsStandardCandidate = !is.na(matched_dilution),
    stringsAsFactors = FALSE
  )
})

results <- Filter(Negate(is.null), results)

if (length(results) == 0) {
  return(data.frame(
    Plate = character(),
    Sample = character(),
    Standard = character(),
    Dilution = character(),
    IsStandardCandidate = logical()
  ))
}

do.call(rbind, results)
}