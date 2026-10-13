# app.R

## Upload UI
  ### Background adjustment
    # # Add level 5 header for Background adjustment
    # h5(strong("Background adjustment"), style = "font-size:120%"),
    # 
    #   # Adds a paragraph of text under the heading
    #   p(
    #     "Check this box if you want to adjust for background. This takes",
    #     "the mean of all background readings for all plates for each antigen",
    #     "in turn, and subtracts this from all other readings for that antigen"
    #     , style = "font-size:120%"),
    #   
    #   checkboxInput("background_adjustment", "Adjust for background", value = FALSE), # subtraction. do the same with division
    #  
    # ),

## Dilutions

  # Defines a default set of dilutions
  dilutions <- c("1/31250", "1/6250", "1/1250", "1/1000", "1/250", "1/50", "1/10")
  
  output$dilutions <- renderText(paste("The default dilutions are set to", paste(dilutions,collapse=", "),". If you have your own dilutions, then you can manually update this here. Please separate the dilutions with a comma."))
  
  # Update dilutions interactively (flexibility)
  observeEvent(input$dilutionButton, {
  
      # Reads user input (dilutionInput) and parses it into a proper format
      dilutions <- parse_dilution_string(input$dilutionInput)
      output$dilutions <- renderText(paste("The dilutions are set to", paste(dilutions,collapse=", "),". You can update this here."))
  })


## Standard Curve
  # Reads user input (std_labelInput) and parses it into a proper format. Update std_label interactively (flexibility)
  observeEvent(input$std_labelButton, {
    std_label <- parse_dilution_string(input$std_labelInput)
    output$std_label <- renderText(paste("The standard curve labels are set to ", paste(std_label,collapse=", "),". You can update this here."))
    })
  
  # Data preparation - Standard curve
  std_curve_data <- reactive({
    req(input$std_labelInput)   # ensures input is not NULL

    # Parse and clean user input
    std_labels <- trimws(unlist(strsplit(input$std_labelInput, ",")))

    # Filter standard rows dynamically
    std_rows <- d()$combinedplates$Sample %in% std_labels
    standard_data <- d()$combinedplates[std_rows, ]

    # Pass filtered data to get.standard
    get.standard(
      data = standard_data,
      std_label = std_labels(),  # character vector
      dilutions = dilutions(),
      n_points = length(std_labels())
    )
  })
  
  # Output
    # ORIGINAL CODE
      output$std_curves <- renderUI({
          req(input$fileUpload)
          fluidRow(
              lapply(ag(),function(antigen){
                  id <- paste0("plot_std_", antigen)
                  plotOutput(outputId = id)
    
                  std_curves <- lapply(d()$platenames, function(x) {
    
                    # Generates interactive standard curve per antigen and plate
                    get.standard(
                          data = d()$plates[[x]],
                          std_label = std_label(), # standard curves are defined as those that contain CP3, the word Std Curve, or WHO (for WHO)
                          dilutions = dilutions,
                          n_points = length(dilutions))
                    })
    
                  names(std_curves) <- d()$platenames
    
                  column(6,
                      renderPlotly({
                          # Generates interactive standard curve plots per antigen and plate
                          plot.std.curve3_interactive(
                              data = std_curves,
                              antigen = antigen,
                              dilutions = dilutions,
                              plate_labels = input$std_curve_input )
                      }) )
              }) )
      })
      
      # Update plate selection dropdown
        observe({
          req(std_data())
          updateSelectInput(
            session,
            "plate_selection",
            choices = names(std_data()$plates),
            selected = names(std_data()$plates)[1]
          )
        })
  
  # Define std_label as default values or user assigned
  std_labels <- if (input$std_labelInput == "") {
    c("CP3", "Std Curve", "WHO")
  } else {
    trimws(unlist(strsplit(input$std_labelInput, ",")))
  }
  
  
  

# luminex_dx_functions_v11.R

## read.batch.std
    read.batch.std <- function(path, std_labels, inc_date = FALSE, inc_plate = FALSE) {
  
      setwd(path)
      temp <- list.files(pattern = "*.csv")
      plate_lab <- substr(temp, 1, nchar(temp) - 4)
      plate <- list()
  
      for(i in 1:length(temp)) {
  
        # Read full file
        plate_raw <- read.csv(temp[i])
  
        # Locate main "DataType" row
        startrow <- grep("DataType", plate_raw[,1])[1]
        plate[[i]] <- read.csv(temp[i], skip=startrow+1)
  
        # Identify Dilution + Median sections
          dt_rows <- grep("DataType", plate[[i]][,1])
          dil_row <- dt_rows[ grep("Dilution Factor", plate[[i]][dt_rows, 2], ignore.case = TRUE)[1] ]
          med_row <- dt_rows[ grep("Median",   plate[[i]][dt_rows, 2], ignore.case = TRUE)[1] ]
          # dilution_row <- grep("dilution",  plate[[i]][(startrow+1):nrow( plate[[i]]), 1], ignore.case = TRUE)[1]
          # median_row   <- grep("median",  plate[[i]][(startrow+1):nrow( plate[[i]]), 1], ignore.case = TRUE)[1]
  
        # Extract Dilution section
          next_dt_after_dil <- dt_rows[dt_rows > dil_row][1]
          dil_section <- dat[(dil_row + 1):(next_dt_after_dil - 1), ]
          # dilution_row <- dilution_row + startrow
          # median_row   <- median_row + startrow
  
        # Extract dilution and median sections
          next_dt_after_med <- dt_rows[dt_rows > med_row][1]
          med_section <- dat[(med_row + 1):(next_dt_after_med - 1), ]
          # dilution_section <-  plate[[i]][(dilution_row + 1):(median_row - 1), ]
          # median_section   <-  plate[[i]][(median_row + 1):nrow( plate[[i]]),   ]
  
        # Numeric columns
          dil_section[,3:ncol(dil_section)] <- sapply(dil_section[,3:ncol(dil_section)], as.numeric)
          med_section[,3:ncol(med_section)] <- sapply(med_section[,3:ncol(med_section)], as.numeric)
        # Reassign column names from the row immediately after section header
        colnames(dilution_section) <-  plate[[i]][dilution_row + 1, ]
        colnames(median_section)   <-  plate[[i]][median_row + 1, ]
  
        # Filter dilution rows based on std_labels
        # std_labels can be regex or exact names
        match_rows <- grepl(paste(std_labels, collapse = "|"), dilution_section[,2])
  
        dilution_filtered <- dilution_section[match_rows, ]
        if (nrow(dilution_filtered) == 0) {
          warning("No rows matching std_labels found in file: ", temp[i])
          next
        }
  
        # Keep only columns 3+ from median
        median_trim <- median_section[, c(1, 3:ncol(median_section)) ]
  
        # Merge on the first column
        merged_df <- merge(dilution_filtered, median_trim, by = colnames(dilution_filtered)[1])
  
        #Convert numeric columns
        num_cols <- names(merged_df)[sapply(merged_df, function(x) all(grepl("^[0-9.]+$", x)))]
        merged_df[num_cols] <- lapply(merged_df[num_cols], as.numeric)
  
  
        if (inc_date==T) {
          plate[[i]]$date <-  plate[[i]][grep("Date", plate[[i]][,1]),2]
          plate[[i]]$date <- as.Date(plate[[i]]$date,"%m/%d/%Y")
        }
  
        if (inc_plate==T) {
          plate[[i]]$plate <- i
  
        }
  
        plate[[i]] <- merged_df
      }
  
      # Assign plate names
      names(plate) <- plate_lab
      return(plate)
    }

        # Identify sections
        section1 <- grep("Dilution", plate[[i]]$Sample)
        section2 <- grep("Median", plate[[i]]$Sample)
        section3 <- grep("Avg Net MFI",plate[[i]]$Sample[1])
  
        # Subset rows in the dilution-to-median section
        plate[[i]] <- plate[[i]][(section1+2):(section2-2),-ncol(plate[[i]])]
  
        # Keep only standard curve rows
        plate[[i]] <- plate[[i]][plate[[i]]$Sample %in% std_label, ]
  
        # Convert numeric columns
        plate[[i]][, 3:ncol(plate[[i]])] <- sapply(plate[[i]][, 3:ncol(plate[[i]])], as.numeric)
  
        if (inc_date) {
          plate[[i]]$date <- as.Date(plate_raw[grep("Date", plate_raw[,1]),2], format="%m/%d/%Y")
        }
        if (inc_plate) {
          plate[[i]]$plate <- i
        }
  
      }
  
      names(plate) <- plate_lab
      return(plate)
    }

  # Potential code/old code
    read.batch.std <- function (path, std_label, inc_date=F, inc_plate=F) {
  
      setwd(path)
      temp <- list.files(pattern="*.csv")
      plate_lab <- substr(temp,1,nchar(temp)-4)
      plate <- list()
  
      # Split std_label input into vector (handles comma + space)
      std_labels <- trimws(unlist(strsplit(std_label, ",")))
  
      for(i in 1:length(temp)) {
  
        plate_raw <- read.csv(temp[i])
        startrow <- grep("DataType",plate_raw[,1])[1]
        plate[[i]] <- read.csv(temp[i], skip=startrow+1)
  
        section1 <- grep("Dilution Factor",plate[[i]]$Sample)
        section2 <- grep("Median",plate[[i]]$Sample)
  
        plate[[i]] <- plate[[i]][(section1+2):(section2-2),-ncol(plate[[i]])]
  
        plate[[i]] <- plate[[i]][plate[[i]]$Sample %in% std_labels, ]
  
        #plate[[i]][,3:ncol(plate[[i]])] <- sapply(plate[[i]][,3:ncol(plate[[i]])],as.character)
        plate[[i]][,3:ncol(plate[[i]])] <- sapply(plate[[i]][,3:ncol(plate[[i]])],as.numeric)
  
        if (inc_date==T) {
          plate[[i]]$date <- plate_raw[grep("Date",plate_raw[,1]),2]
          plate[[i]]$date <- as.Date(plate[[i]]$date,"%m/%d/%Y")
        }
  
        if (inc_plate==T) {
          plate[[i]]$plate <- i }
  
      }
  
      names(plate) <- plate_lab
  
      return(plate)
  
    }