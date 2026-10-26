library(shiny)
library(shinyFiles)
library(plotly)
library(glue)
library(rsconnect)
library(DT)

# Define UI ----
# Define UI for random distribution app ----


# TAB LAYOUT

intro_tab_ui <- function(){
  tabPanel(
    "Intro",
    sidebarLayout(
      
      sidebarPanel(
        "", width=3,
        
        h3(strong("Variable requirements")), 
        p("Please read the following variable input requirements before proceeding with the app", 
          style = "font-size:100%"),
        
        # Dilutions
        h3(tags$u("Dilutions")),
          tags$li("Dilutions should be formatted with 1/n. For example: 1/100, 1/50, 1/25..", style = "font-size:100%"),
          tags$li("Each dilution should be seperated by a comma and a space.", style = "font-size:100%"),
          tags$li("Dilutions should be listed from highest to lowest.", style = "font-size:100%"),
          tags$div(HTML('<span style="color:green; font-weight:bold;">Correct example format:</span> 1/50, 1/100, 1/200...'), style = "font-size:100%"),
          tags$div(HTML('<span style="color:red; font-weight:bold;">Incorrect example format:</span> 0.25, 100, 1:200'), style = "font-size:100%"),
          
        # Standard curve
        h3(tags$u("Standard Curve")),
          tags$li("Each standard curve should be seperated by a comma and a space.", style = "font-size:100%"),
          tags$li("Please ensure that the standard curve input is in the correct capitalisation.", style = "font-size:100%"),
          tags$li("If you have multiple standard curve controls on a plate, then only enter one set of standards at a time.", style = "font-size:100%"),
          tags$div(HTML('<span style="color:green; font-weight:bold;">Correct example format:</span> P1, P2, P3, P4, P5, P6, P7, P8'), style = "font-size:100%"),
        
        # Background samples
        h3(tags$u("Background samples"))
      ),
      
      mainPanel(
        fluidRow(
          column(12, align="center",
                 img(src='logo.png', width=200, height=200))
        ),
        
        fluidRow(
          column(12, align="center",
                 h1(strong("Luminex Data Analysis (LDA)"))
          )
        ),
        
        # Intro to the app
        p("Welcome to the Luminex Data Analysis app. This app is designed to help you pre-process your raw MBA Luminex dataset into an analysable format. Please read the variable input requirements on the left panel before proceeding with the app.",
          style = "font-size:100%"),
        
        hr(),
        
        # About
        h2(strong("About")),
        p("This R Shiny app was created by researchers in the Drakeley group at the London School of Hygiene and Tropical Medicine to standardise the quality control steps of analysing raw multiplex datasets from a Luminex machine.", 
          style = "font-size:100%"),
        
        div(
          h3("Composition"),
          tags$ul(
            tags$li(strong("Bead count:"), " Visualise wells with low bead counts on each plate. Default threshold is 30 beads."),
            tags$li(strong("Concordance:"), " A way to visualise variation across plates."),
            tags$li(strong("Standard curve:"), " Plots MFI values against log dilution values."),
            tags$li(strong("Levey-Jennings plots:"), " Tracks control samples over time."),
            tags$li(strong("Coefficient of variation plots:"), " Measures variation between samples in a plate."),
            tags$li(strong("Normalisation:"), " Seeks differences between units or batches."),
            tags$li(strong("Loess:"), " Visualise differences between raw and normalised plate."),
            tags$li(strong("Download:"), " Plates can be merged, cleaned, and downloaded for downstream analysis.")
          ), style = "font-size:100%"
        ),
        
        p("No coding experience is necessary to interpret and use this app.", 
          style = "font-size:100%"),
        
        hr(),
        
        # Data Storage
        h2(strong("Data Storage")),
        
        hr(),
        
        # Contact information
        h2(strong("Contact"))
      )
    )
  )
}

upload_tab_ui <- function() {
  tabPanel(
    "Upload",
    sidebarLayout(
      # --- SIDEBAR: ALL CONTROLS ---
      sidebarPanel(
        width = 3,
        h3(strong("File upload"), style = "font-size:130%"),
        p("Upload your data using the file upload button below.", style = "font-size:120%"),
        fileInput(
          inputId = "fileUpload",
          label = "",
          multiple = TRUE,
          buttonLabel = "Browse...",
          placeholder = "No file selected"
        ),
        
        hr(),
        
        # Plate size (96 well plate or 384 well)
        h3(strong("Select Plate Size"), style = "font-size:130%"),
        selectInput(
          inputId = "plate_size",
          label = "",
          choices = c("96 Well" = 96, "384 Well" = 384),
          selected = 96
        ),
        
        hr(),
        
        # Dilutions
        h3(strong("Dilutions"), style = "font-size:130%"),
        textInput(inputId = "dilutionInput", label = ""),
        actionButton("dilutionButton", "Submit dilutions"),
        hr(),
        
        # Standard curve
        h3(strong("Standard Curve"), style = "font-size:130%"),
        textInput(inputId = "std_labelInput", label = ""),
        actionButton("std_labelButton", "Submit standard curves"),
        hr(),
        
        # Background
        h3(strong("Background samples"), style = "font-size:130%"),
        textInput(inputId = "bkg_labelInput", label = ""),
        actionButton("bkg_labelButton", "Submit background samples"),
        hr(),
        
        # Controls
        h3(strong("Control samples (optional)"), style = "font-size:130%"),
        textInput(inputId = "control_labelInput", label = ""),
        actionButton("control_labelButton", "Submit control samples")
      ),
      
      # Main panel
      mainPanel(
        tabsetPanel(
          id = "upload_subtabs",
          
          # SUB-TAB 1: Plate Map
          tabPanel(
            "Plate Map",
            fluidPage(
              br(),
              fluidRow(
                # Dropdown for plate selection
                column(6, uiOutput("plate_selection_control"))
              ),
              hr(),
              h4("Plate Layout Verification", style = "font-size: 24px", align = "left"),
              # Plot output with the fixed height
              plotOutput("plate_layout_map", height = "700px")
            )
          ),
          
          # Sub tab 2: Upload summary
          tabPanel(
            "Upload Summary",
            fluidPage(
              fluidRow(
                column(12, align = "center", h1("Luminex Data Analysis"))
              ),
              p("Please upload your data on the left. Once uploaded, you can use the 'Plate Map' tab to visualize the layout.", 
                style = "font-size:120%"),
              hr(),
              h5(strong("Data Overview"), style = "font-size:150%", align = "center"),
              div(style = "font-size:120%;", tableOutput("upload_summary")),
              br(),
              
              fluidRow(
                column(6, tableOutput("plate_table")),
                column(6, tableOutput("antigen_table"))
              ),
              hr(),
              h4("Sample Names and Standard Detection"),
              p("Review the detected standard names and dilution values."),
              DT::dataTableOutput("sample_names_table")
            )
          )
        )
      )
    )
  )
}

bead_count_ui <- function() {
  
  tabPanel("Bead Count",
           sidebarLayout(
             sidebarPanel(width=3,
                          
                          # Add level 5 header for Bead threshold
                          h5(strong("Bead count threshold"), style = "font-size:120%"),
                          
                          p(
                            "Low bead counts may not ensure statistically valid results.",
                            "Please add the minimum number beads required to proceed with analysis (usually 30 or 50).", style = "font-size:120%"),
                          
                          hr(),
                          
                          # Introduce the bead colour meaning
                          p(),
                          p(
                            "Each dot is the beads per antigen per well. ",
                            "A black dot represents the bead count of the antigen in the well being over the specified threshold. ",
                            "A red dot represents the bead count of the antigen in the well being under the specified threshold.", style = "font-size:120%"),
                          
                          hr(),
                          
                          # User input for bead threshold
                          numericInput(inputId = "bead_threshold",
                                       label = "Please enter a minimum bead count: ",
                                       value = 30, # default
                                       min = 0),
                          
                          # Plate selection
                          selectInput(inputId = "plate_selection",
                                      label = "Select Plate:",
                                      choices = NULL),
                          
                          hr(),
                          
                          h5(strong("Plates with low bead counts")),
                          
                          p("The following plates have at least one well with at an antigen bead count below your threshold."),
                          
                          uiOutput("low_bead_plates_ui")

             ),
             mainPanel(uiOutput("bead_count_ui_dynamic")))
           )
}

concordance_plots_ui <- function(){
    
    tabPanel("Concordance plots",
             sidebarLayout(
               sidebarPanel(width=3,
                            
                            # Add level 5 header for Bead threshold
                            h5(strong("Concordance"), style = "font-size:110%"),
                            
                            p("Concordance is a way to view the technical variation across plates. This page will compare MFI values across two plates.", style = "font-size:100%"),
                            
                            hr(),
                            
                            # Plate selection (x axis)
                            selectInput(
                              inputId = "plate_selection_x",
                              label = "Select Plate (x axis):",
                              choices = NULL
                            ),
                            
                            # Plate selection (y axis)
                            selectInput(
                              inputId = "plate_selection_y",
                              label = "Select Plate (y axis):",
                              choices = NULL
                            ),
                            
                            uiOutput("antigen_filter_ui")
               ),
                            
               mainPanel(
                 h4("MFI Concordance Between Plates"),
                 p("Each point represents a matched standard and dilution."),
                 uiOutput("concordance_plot_ui"),
                 hr(),
                 h4("Concordance Statistics"),
                 DT::dataTableOutput("concordance_table")
               )
    ))
  }

std_curve_plots_ui <- function(){
  tabPanel("Standard Curves",
           sidebarLayout(
             sidebarPanel(width = 3,
                          p("Standard curves show MFI plotted against the log of the dilution.", style="font-size:110%"),
                          hr(),
                          
                          div(
                            style = "background-color: #f8f9fa; padding: 12px; border-radius: 6px; border: 1px solid #dee2e6;",
                            h4(strong("Standard Verification"), style = "font-size: 110%; margin-top: 0; color: #2c3e50;"),
                            p("Active standard label from Upload page:", style = "font-size: 95%; margin-bottom: 5px;"),
                            verbatimTextOutput("verified_std_label")
                          ),
                          hr(),
                          
                          p(strong("Select plates to compare (max 5):"), style="font-size:110%"),
                          uiOutput("std_curve_plate_selection"),
                          hr(),
                          
                          p(strong("Select standard series:"), style = "font-size:110%"),
                          uiOutput("std_curve_standard_selection"),
                          hr(),
                          
                          p(strong("Select antigen:"), style = "font-size:110%"),
                          uiOutput("std_curve_antigen_selection")
                          # p(strong("Select antigen:"), style="font-size:120%") #,
                          # uiOutput("std_curve_antigen_selection")
             ),
             mainPanel(
               h4("Standard Curve Plot", align = "center"),
               plotlyOutput("std_curve_plot", height = "600px"),
               hr(),
               h4("Standard Curve Data", align = "center"),
               DT::dataTableOutput("std_table")
             )
           )
  )
}

levey_jennings_plots_ui <- function(){

    tabPanel("Plate variation", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "Levey-Jennings plots tracks how control samples behave over time across different plates.",
                    " The Levey-Jennings plot is a plot shows the median MFI for each antigen and",
                    " plate. The lines indicate 1 and 2 standard deviations from the mean. ",
                    "This plot can be used to identify antigens that are not performing well."
                ),
                
                hr(),
                h5("Levey-Jennings controls"),
                uiOutput("levey_jennings_select_points"),
                
                hr(),
                
                p(
                    "The coefficient of variation is a measure of the variation between samples",
                    "in a plate. It is calculated as the standard deviation divided by the mean."
                ),
            
            ),
            
                
            mainPanel(
                h4("Levey-Jennings plots"),
                uiOutput("levey_jennings_plots"),
                
                hr(),
                
                h4("Coefficient of variation plots"),
                uiOutput("coefv_plots")
                
            )
        )
    )
}

normalisation_ui <- function(){
    tabPanel("Normalisation", 
        sidebarLayout(
            sidebarPanel("", width=3,
                
                         p("A curve is fitted to the control samples for each antigen. The curve",
                           "is then used to estimate 100 MFI values at 100 incremental dilutions",
                           "within curve's original range.", style = "font-size:130%"),
                         
                         p("Different ways of fitting the curve have been implemented. Please",
                           "choose one and wait a few seconds for the visualisation to appear:", style = "font-size:130%"),
                         
                         p(strong("poly method"),": fits a cubic polynomial curve to the data", style = "font-size:130%"),
                         p(strong("nls method"), ": uses non-linear least squares fitting using the nls2 library", style = "font-size:130%"),
                         p(strong("kernsmooth method"), ": estimates a regression function using local polynomials using the KernSmooth library", style = "font-size:130%"),
                         p(strong("cgam method"), ": fits a generalised additive model using the cgam library", style = "font-size:130%"),
                         uiOutput("normControls")
                         ),
            mainPanel(
                uiOutput("normalisation_plots"),
                # uiOutput("normalised_ag_plot")
            )
        )
    )
}

loess_ui <- function(){
    tabPanel("Loess", 
        sidebarLayout(
            sidebarPanel("", width=3,
                p(
                    "MFIs are subtracted between the chosen reference and each",
                    "target plate curve at 100 dilutions for each antigen,",
                    "producing Δ MFI values for the dilutions. A loess",
                    "curve is then fitted to the Δ MFI values for each",
                    "plate/antigen. The loess curve is then used to adjust the",
                    "raw MFI values for each plate/antigen.",
                    ),
                p(paste(
                    "Select a plate to use as a reference and use the tabs to",
                    "visualise the outputs"
                )),
                uiOutput("loessControls1"),
                uiOutput("loessControls2"),
                uiOutput("download_button")
            ),
            mainPanel(
                tabsetPanel(
                    id="loess_tabs",
                    type = "tabs",
                    tabPanel("loess_plots",
                        br(),
                        p(paste(
                            "These plots visualise the difference between",
                            "the reference raw data and the normalised plate"
                        )),
                        uiOutput("loess_plots")
                    ),
                    tabPanel("raw_comparison",
                        br(),
                        p(paste(
                            "These plots visualise the loess adjusted MFI",
                            "value for each plate and each antigen. Each",
                            "line represents an individual plate"
                        )),
                        uiOutput("loess_ag_plot"),
                        uiOutput("loess_data")
                    )
                )

            )
        )
    )
}


download_data_ui <- function() {
  
  tabPanel("Download data",
           sidebarLayout(
             sidebarPanel(width=3,
                         
                          # Add description that this page is to adjust for data if happy. After all the steps if it looks good then download
                          
                          # Adjusting low bead count
                          h5(strong("Adjusting low bead count"), style = "font-size:120%"),
                            
                            # User input for bead threshold
                            numericInput(inputId = "bead_threshold",
                                         label = "Please enter a minimum bead count: ",
                                         value = 30, # default
                                         min = 0),
                            
                            # Adds a paragraph of text under the heading
                            p("Check this box if you want to remove samples with a low bead count for the cleaned dataset",
                              style = "font-size:120%"),
                            
                            checkboxInput("bead_adjustment", "Remove low beads", value = FALSE),
                          
                          # Inserts a horizontal line (a divider in the UI)
                          hr(),
                          
                          # Background adjustment
                          h5(strong("Remove background"), style = "font-size:120%"),
                          
                            # Adds text
                            p("Please select how you would like to remove the background MFI", style = "font-size:120%"),
                            
                            # Adds a paragraph of text under the heading
                            p(
                              radioButtons(
                                inputId = "backgroundremoval",
                                label = "Choose one option:",
                                choices = c("Subtract", "Divide"),
                                selected = "Subtract"  # default choice
                              )
                            )
                          
             ),
             mainPanel()
             )
  )
}


ui <- fluidPage(
  
  # Increase the size of the tab panel
    tags$head(tags$style(HTML("
      /* Tab headers */
      .nav-tabs > li > a {
        font-size: 18px !important;       /* text size */
        padding: 15px 20px !important;    /* increase height & horizontal padding */
      }
  
      /* Active tab */
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus,
      .nav-tabs > li.active > a:hover {
        font-size: 18px !important;
        padding: 15px 20px !important;
      }
  
      /* Space between tabs and content */
      .tab-content {
        padding-top: 20px !important;
      }
    "))),
  
  # Tabset panel
    tabsetPanel(
      id = "tabs_id",
      type = "tabs",
      intro_tab_ui(),
      upload_tab_ui(),
      bead_count_ui(),
      concordance_plots_ui(),
      std_curve_plots_ui(),
      levey_jennings_plots_ui(),
      normalisation_ui(),
      loess_ui(),
      download_data_ui()
    )
)


# Define server logic ----

parse_dilution_string <- function(input){
    trimws(unlist(strsplit(trimws(input),",")))
}
  
server <- function(input, output, session) {
  
  # Read in functions
    source("luminex_dx_functions_v11.R")
    source("normalisation_functions.R")
    source("bead_plot_function.R")
  
  
  # Set up dilutions
    # Defines a default set of dilutions
    dilutions <- c("1/31250","1/6250","1/1250","1/250","1/50")
    
    # Add text to home button on what the value the dilutions are set to
    output$dilutions <- renderText(paste("The default dilutions are set to", paste(dilutions,collapse=", "),". If you have your own dilutions, then you can manually update this here. Please separate the dilutions with a comma."))
  
    # Reactive dilutions (default or user input)
    user_dilutions <- reactiveVal(dilutions)
    
    observeEvent(input$dilutionButton, {
      parsed <- parse_dilution_string(input$dilutionInput)
      if(length(parsed) > 0){
        user_dilutions(parsed)
        output$dilutions <- renderText(paste("The dilutions are set to", paste(parsed, collapse=", "),". You can update this here."))
      }
    })
    
  
  # Set up standard curves
    output$std_labels <- renderText({paste("Please manually list standard curve values here. Each value should be seperated with a comma.")})
    
    std_labels <- reactive({
      parse_dilution_string(input$std_labelInput)
    })
    
    # Update std_labels interactively (flexibility)
    observeEvent(input$std_labelButton, {
      s_labels <- parse_dilution_string(input$std_labelInput)
      
      if(length(s_labels) > 0){
        output$std_labels <- renderText(
          paste("The standard curve labels are set to ", paste(s_labels, collapse=", "), ". You can update this here.")
        )
      }
    })
    
    # On standard curve page
      # Capture the standard label when the submit button is clicked
      submitted_std_label <- eventReactive(input$std_labelButton, {
        req(input$std_labelInput)
        input$std_labelInput
      }, ignoreNULL = FALSE)
      
      # Render the verification text in the sidebar panel
      output$verified_std_label <- renderText({
        lbl <- submitted_std_label()
        if (is.null(lbl) || lbl == "") {
          "Warning: No standard label submitted yet. Please submit one on the Upload tab."
        } else {
          paste0("✓ Active Standard: ", lbl)
        }
      })
    
    
  # Set up background values 
    bkg_labels <- ("(?i)(blank|background)") # contains any cell that has blank or background
    
    output$bkg_labels <- renderUI({
      HTML(paste0(
        "The default background samples are set to variables that contain ",
        "<i>blank</i> or <i>background</i>. ",
        "If your background samples are labelled differently, then you can manually update this here. ",
        "Please separate each background sample name with a comma."
      ))
    })
    
    # Update bkg_labels interactively based on user input
    observeEvent(input$bkg_labelButton, {
      
      # Parse user input: split by comma, trim whitespace
      labels <- trimws(unlist(strsplit(input$bkg_labelInput, ",")))
      
      if(length(labels) > 0){
        # Combine into regex pattern
        regex_pattern <- paste(labels, collapse = ",")
        bkg_labels(regex_pattern)  # update reactiveVal
        
        # Update displayed text
        output$bkg_labels <- renderText({
          paste("The background samples are now set to: ", paste(labels, collapse = ", "))
        })
      }
    })
    
    # Set up controls
      output$control_labels <- renderText({paste("Please manually list standard curve values here. Each value should be seperated with a comma.")})
      
      control_labels <- reactive({
        parse_dilution_string(input$control_labelInput)
      })
      
      # Update control_labels interactively (flexibility)
      observeEvent(input$control_labelButton, {
        s_labels <- parse_dilution_string(input$control_labelInput)
        
        if(length(s_labels) > 0){
          output$control_labels <- renderText(
            paste("The controls are set to ", paste(s_labels, collapse=", "), ". You can update this here.")
          )
        }
      })
    
  # File upload and data preparation
    d <- reactive({
      
      # Make sure a file has been uploaded
      req(input$fileUpload)
      
      # Extract the folder path and the file label, drop .csv extension
      raw_data_path <- dirname(input$fileUpload[1,"datapath"])
      plate_lab <- substr(input$fileUpload$name,1,nchar(input$fileUpload$name)-4)
      
      # Reading uploaded data
      d = list()
      d$rawdatapath <- raw_data_path
      d$uploadedfilenumber <- nrow(input$fileUpload)
      
        ## Load batch data from the uploaded files
        d$plates <- read.batch(path=raw_data_path)
        
        ## Assign user-friendly names to the plates
        names(d$plates)<-plate_lab
      
      # Join plates
        ## Store metadata: plate names, number of plates
        d$platenames <- names(d$plates)
        d$numplates <- length(d$platenames)
      
        ## Combine plates into one dataset
        d$combinedplates <- join.plates(d$plates)
        
          ### Skip NAs
          d$combinedplates <- na.omit(d$combinedplates)
          for (i in seq_along(d$plates)){
            d$plates[[i]] <- na.omit(d$plates[[i]])
          }
        
        ## Extracts antigen names (columns 3 onward = data columns)
        d$ags <- colnames(d$combinedplates)[3:ncol(d$combinedplates)]
        
        ## Extracts Sample names
          d$sample_type <- colnames(d$combinedplates)[2]
          
          # Replace any cell containing "Unknown" (case-insensitive) with "Sample"
          d$sample_type <- ifelse(grepl("Unknown", d$sample_type, ignore.case = TRUE), "Sample", d$sample_type)[2]
      
        ## Metadata with dates - include date & plate info
        d$datesplates <- read.batch(path=raw_data_path,inc_date=T,inc_plate=T)
        d$combineddatesplates <- join.plates(d$datesplates)
        d
    })
    
  # Create dataset for plate mapping
  td <- reactive({
      req(input$fileUpload)
      files <- input$fileUpload
      
      plate_list <- list()
      all_data <- data.frame()
      
      for (i in 1:nrow(files)) {
        # skip=41 based on your file structure
        temp_df <- read.csv(files$datapath[i], skip = 41, stringsAsFactors = FALSE)
        
        # CLEANING THE LOCATION: "1(1,A1)" -> "A1"
        # This regex looks for the characters after the comma and before the closing bracket
        temp_df$Well_ID <- sub(".*,([A-H][0-9]+)\\).*", "\\1", temp_df$Location)
        temp_df$Row_ID  <- substr(temp_df$Well_ID, 1, 1)
        temp_df$Col_ID  <- as.numeric(gsub("[A-H]", "", temp_df$Well_ID))
        
        plate_name <- sub(".csv", "", files$name[i])
        temp_df$Plate <- plate_name
        plate_list[[plate_name]] <- temp_df
        all_data <- rbind(all_data, temp_df)
      }
      
      # Identify Antigens: based on your CSV, it's usually columns 3 to 15
      # Adjust index if necessary
      ags <- names(temp_df)[3:15] 
      
      return(list(
        plates = plate_list,
        combinedplates = all_data,
        ags = ags,
        platenames = names(plate_list),
        numplates = length(plate_list),
        uploadedfilenumber = nrow(files)
      ))
    })
    
  # Reactive to "tag" wells based on user inputs
  tagged_data <- reactive({
        req(td())
        df <- td()$combinedplates
        
        # Standardize the Well Extraction
        # This regex handles formats like "1(1,A1)" or "1(1,P24)"
        # It captures the letter (A-P) and all following digits (96 and 384 well)
        df$Well_ID <- sub(".*,([A-P][0-9]+)\\).*", "\\1", df$Location)
        df$Row_ID  <- substr(df$Well_ID, 1, 1)
        df$Col_ID  <- as.numeric(gsub("[A-P]", "", df$Well_ID))
        
        ## Remove NA row
        df <- df[!is.na(df$Row_ID) & df$Row_ID != "" & !is.na(df$Col_ID), ]
        
        # 3. Create Search Patterns from UI Inputs
        # Helper to handle empty inputs or comma-separated lists
        get_pattern <- function(input_val, default = "___NONE___") {
          if (is.null(input_val) || input_val == "") return(default)
          # Split by comma, trim whitespace, and create an OR regex (e.g. "pos|control")
          vals <- trimws(strsplit(input_val, ",")[[1]])
          paste0("(?i)(", paste(vals, collapse = "|"), ")")
        }
        
        std_pattern  <- get_pattern(input$std_labelInput)
        ctrl_pattern <- get_pattern(input$control_labelInput)
        # Default background includes "blank" or "background" even if input is empty
        bkg_pattern  <- get_pattern(input$bkg_labelInput, default = "(?i)(blank|background)")
        
        # 4. Identify Well Types
        # We use the 'Sample' column (usually column 2)
        sample_names <- df$Sample # Or df[[2]] if the column name varies
        
        df$Well_Type <- dplyr::case_when(
          grepl(std_pattern,  sample_names) ~ "Standard",
          grepl(ctrl_pattern, sample_names) ~ "Control",
          grepl(bkg_pattern,  sample_names) ~ "Background",
          TRUE ~ "Sample"
        )
        
        return(df)
    })
    
  # Select plate for plate map
  output$plate_selection_control <- renderUI({
      req(td()) # Ensure data is uploaded
      selectInput(
        inputId = "plate_select_map", 
        label = "Select Plate to View:", 
        choices = td()$platenames
      )
    })
  
  # Identify all standard curves from uploaded plate data
  standard_info <- reactive({
    
    req(d(), user_dilutions())
    
    all_standards <- lapply(d()$plates, function(plate_data) {
      identify.standards(plate_data, user_dilutions())
    })
    
    unique_standards <- unique(unlist(lapply(
      all_standards,
      function(x) {
        x$Standard[!is.na(x$Standard) & nzchar(x$Standard)]
      }
    )))
    
    sort(unique_standards)
  })
  
  # ---- CONCORDANCE BETWEEN TWO PLATES ----
  # Update the plate dropdowns when files have been uploaded
  observe({
    req(d())
    plates <- d()$platenames
    
    updateSelectInput(
      session, "plate_selection_x",
      choices = plates,
      selected = plates[1]
    )
    
    updateSelectInput(
      session, "plate_selection_y",
      choices = plates,
      selected = if (length(plates) >= 2) plates[2] else plates[1]
    )
  })
  
  # Antigen selector: select one or more antigens, or all antigens
  output$antigen_filter_ui <- renderUI({
    req(d())
    
    selectizeInput(
      inputId = "concordance_antigens",
      label = "Select antigens:",
      choices = d()$ags,
      selected = d()$ags,
      multiple = TRUE,
      options = list(plugins = list("remove_button"))
    )
  })
  
  
  # ---- Prepare matchd concordance across plates ----
  # Prepare matched standard-control MFI measurements
  concordance_data <- reactive({
    
    req(d())
    req(input$plate_selection_x, input$plate_selection_y)
    
    # Read standard labels from Upload page
    labels <- trimws(as.character(submitted_standards()))
    labels <- unique(labels[!is.na(labels) & nzchar(labels)])
    
    # Message if no standards entered on Upload page
    validate(
      need(
        length(labels) > 0,
        paste(
          "Please enter the full standard sample labels on the Upload page",
          "and click 'Submit standard curves'."
        )
      )
    )
    
    # Message if select same plate
    validate(
      need(
        input$plate_selection_x != input$plate_selection_y,
        "Please select two different plates."
      )
    )
    
    plates <- d()$plates
    # Convert the selected plate values to single character strings
    plate_x_name <- as.character(
      unlist(input$plate_selection_x, use.names = FALSE)
    )
    plate_y_name <- as.character(
      unlist(input$plate_selection_y, use.names = FALSE)
    )
    
    validate(
      need(
        length(plate_x_name) == 1 && nzchar(plate_x_name),
        "Please select a plate for the x-axis."
      ),
      need(
        length(plate_y_name) == 1 && nzchar(plate_y_name),
        "Please select a plate for the y-axis."
      )
    )
    
    plate_x_name <- plate_x_name[1]
    plate_y_name <- plate_y_name[1]
    
    # Find the plate positions by name
    plate_x_index <- match(plate_x_name, names(plates))
    plate_y_index <- match(plate_y_name, names(plates))
    
    validate(
      need(!is.na(plate_x_index), "The selected x-axis plate was not found."),
      need(!is.na(plate_y_index), "The selected y-axis plate was not found."),
      need(plate_x_name != plate_y_name, "Please select two different plates.")
    )
    
    # Use numeric indices to retrieve the plate data
    plate_x <- plates[[plate_x_index]]
    plate_y <- plates[[plate_y_index]]
    
    antigens <- intersect(d()$ags, 
                          intersect(
                            names(plate_x),
                            names(plate_y)
                            ))
    
    validate(
      need(length(antigens) > 0,
           "The selected plates have no antigens in common.")
    )
    
    # Normalise only letter case and leading/trailing whitespace.
    # Keep punctuation such as "/" and ":" unchanged.
    label_keys <- tolower(labels)
    
    # Extract measurements for the exact submitted sample labels
    extract_plate <- function(pd) {
      
      sample_names <- trimws(as.character(pd$Sample))
      sample_keys <- tolower(sample_names)
      
      keep <- !is.na(sample_keys) &
        sample_keys %in% label_keys
      
      if (!any(keep)) {
        return(NULL)
      }
      
      selected <- pd[keep, , drop = FALSE]
      selected_keys <- sample_keys[keep]
      
      # Convert the selected samples into one row per
      # Sample x Antigen combination
      long_data <- do.call(
        rbind,
        lapply(antigens, function(antigen) {
          
          data.frame(
            SampleKey = selected_keys,
            Antigen = antigen,
            MFI = suppressWarnings(
              as.numeric(as.character(selected[[antigen]]))
            ),
            stringsAsFactors = FALSE
          )
        })
      )
      
      # Remove missing or non-numeric MFI values
      long_data <- long_data[
        is.finite(long_data$MFI),
        ,
        drop = FALSE
      ]
      
      if (nrow(long_data) == 0) {
        return(NULL)
      }
      
      # Average replicate wells for the same sample and antigen
      aggregate(
        MFI ~ SampleKey + Antigen,
        data = long_data,
        FUN = mean
      )
    }
    
    x <- extract_plate(plate_x)
    y <- extract_plate(plate_y)
    
    validate(
      need(
        !is.null(x) && nrow(x) > 0,
        paste(
          "None of the submitted sample labels matched Plate X:",
          input$plate_selection_x
        )
      ),
      need(
        !is.null(y) && nrow(y) > 0,
        paste(
          "None of the submitted sample labels matched Plate Y:",
          input$plate_selection_y
        )
      )
    )
    
    # Rename MFI columns, then match by sample identifier and antigen
    names(x)[names(x) == "MFI"] <- "MFI_X"
    names(y)[names(y) == "MFI"] <- "MFI_Y"
    
    matched <- merge(
      x,
      y,
      by = c("SampleKey", "Antigen"),
      all = FALSE
    )
    
    validate(
      need(
        nrow(matched) > 0,
        paste(
          "Samples were found on both plates, but no matching",
          "sample and antigen pairs could be formed."
        )
      )
    )
    
    # Display the original submitted sample label
    matched$Sample <- labels[
      match(matched$SampleKey, label_keys)
    ]
    
    matched
  })
  
  # ------ Size concordance plots like with the bead-count plots ------ #
  output$concordance_plot_ui <- renderUI({
    
    req(concordance_data(), input$concordance_antigens)
    
    selected_antigens <- intersect(
      input$concordance_antigens,
      unique(concordance_data()$Antigen)
    )
    
    req(length(selected_antigens) > 0)
    
    n_ags <- length(selected_antigens)
    n_rows <- ceiling(length(selected_antigens) / 3)
    plot_height <- max(600, n_rows * 320)
    
    plotlyOutput("concordance_plot", height = paste0(plot_height, "px"))
  })
  
  # output$concordance_plot_ui <- renderUI({
  #   
  #   req(concordance_data(), input$concordance_antigens)
  #   
  #   selected_antigens <- intersect(
  #     input$concordance_antigens,
  #     unique(concordance_data()$Antigen)
  #   )
  #   
  #   req(length(selected_antigens) > 0)
  #   
  #   # Three antigen plots per row
  #   n_rows <- ceiling(length(selected_antigens) / 3)
  #   
  #   # Approximate square panels
  #   plot_height <- max(400, n_rows * 300)
  #   
  #   plotlyOutput(
  #     "concordance_plot",
  #     height = paste0(plot_height, "px")
  #   )
  # })
  
  # ----- Plot one concordance scatterplot per antigen ----- #
  output$concordance_plot <- plotly::renderPlotly({
    
    req(concordance_data(), input$concordance_antigens)
    
    df <- concordance_data()
    selected_antigens <- input$concordance_antigens
    
    df <- df[df$Antigen %in% selected_antigens, , drop = FALSE]
    
    validate(need(nrow(df) > 0, "No matched measurements for the selected antigens."))
    
    # Remove missing or non-finite measurements
    df <- df[is.finite(df$MFI_X) & is.finite(df$MFI_Y), ,drop = FALSE]
    
    validate(need(nrow(df) > 0, "No valid paired MFI measurements available."))
    
    # Order in antigen selection
    df$Antigen <- factor(df$Antigen, levels = selected_antigens)
    
    n_ags <- length(selected_antigens)
    n_rows <- ceiling(length(selected_antigens) / 3)
    plot_height <- max(600, n_rows * 320)
    
    # Create informative hover labels for every point
    df$HoverLabel <- paste0(
      "<b>Standard sample:</b> ", df$Sample,
      "<br><b>Antigen:</b> ", df$Antigen,
      "<br><b>Plate X MFI:</b> ",
      format(round(df$MFI_X, 1), big.mark = ","),
      "<br><b>Plate Y MFI:</b> ",
      format(round(df$MFI_Y, 1), big.mark = ","),
      "<br><b>Difference (Y - X):</b> ",
      format(round(df$MFI_Y - df$MFI_X, 1), big.mark = ",")
    )
    
    # Plot
    p <- ggplot(df, aes(x = MFI_X, y = MFI_Y)) +
      geom_abline(
        intercept = 0, slope = 1,
        linetype = "dashed", colour = "cornflowerblue") +
      geom_point(aes(text = HoverLabel), size = 2.5) +
      facet_wrap(~ Antigen, scales = "free", ncol = 3) +
      labs(
        title = paste("Plate Concordance:", input$plate_selection_x, "vs", input$plate_selection_y),
        x = paste(input$plate_selection_x, " MFI"),
        y = paste(input$plate_selection_y, " MFI")
      ) +
        theme_minimal(base_size = 13) +
        theme(
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1,
          plot.margin = margin(8, 8, 8, 8)
        )
    
    # Convert the complete faceted plot to Plotly
    n_rows <- ceiling(length(unique(df$Antigen)) / 3)
    plot_height <- max(600, n_rows * 320)
      
    # Convert to an interactive Plotly chart
    plotly::ggplotly(p, height = plot_height, tooltip = "text") |>
      plotly::layout(
        autosize = TRUE,
        margin = list(l = 65, r = 25, t = 80, b = 65)
      ) |>
      plotly::config(responsive = TRUE)
   })
  #   
  #   # Arrange the separate plots into three columns
  #   n_rows <- ceiling(length(plot_list) / 3)
  #   
  #   combined_plot <- plotly::subplot(
  #     plot_list,
  #     nrows = n_rows,
  #     shareX = FALSE, shareY = FALSE,
  #     titleX = TRUE, titleY = TRUE, margin = 0.06
  #   )
  #   
  #   # Add an overall title and responsive layout
  #   combined_plot |>
  #     plotly::layout(
  #       title = list(
  #         text = paste("Plate Concordance:", input$plate_selection_x, "vs", input$plate_selection_y),
  #         x = 0.5),
  #       autosize = TRUE,
  #       margin = list(l = 65, r = 25, t = 85, b = 65)
  #     ) |>
  #     plotly::config(
  #       responsive = TRUE
  #     )
  # })
  
  # ---- Concordance stats table by antigen ----- #
  # Calculate concordance statistics separately for each antigen
  output$concordance_table <- DT::renderDataTable({
  
    req(concordance_data())
    req(input$concordance_antigens)
    df <- concordance_data()
    df <- df[df$Antigen %in% input$concordance_antigens, , drop = FALSE]
    
    validate(
      need(nrow(df) > 0, "No matched measurements available.")
    )
    
    # Split measurements by antigen
    antigen_groups <- split(df, as.character(df$Antigen))
    
    results_list <- lapply(names(antigen_groups), function(antigen_name) {
      
      z <- antigen_groups[[antigen_name]]
      x <- z$MFI_X
      y <- z$MFI_Y
      # Remove missing or non-finite paired measurements
      keep <- is.finite(x) & is.finite(y)
      x <- x[keep]
      y <- y[keep]
      
      n <- length(x)
      
      pearson_r <- if (
        n >= 2 && stats::sd(x) > 0 && stats::sd(y) > 0
      ) {
        stats::cor(x, y, method = "pearson")
      } else {
        NA_real_
      }
      
      data.frame(
        Antigen = antigen_name,
        # 'Matched points' = n,
        'Pearson (r)' = pearson_r,
        'MFI difference (Plate Y minus Plate X)' = mean(y - x),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
    
    result <- do.call(rbind, results_list)
    rownames(result) <- NULL
    
    DT::datatable(
      result,
      rownames = FALSE,
      options = list(pageLength = 25, scrollX = TRUE)
    ) |>
      DT::formatRound(
        columns = c(
          "Pearson (r)", "MFI difference (Plate Y minus Plate X)"
        ),
        digits = 3
      )
  })
  
  
  # Show the detected standard curves in a dropdown
  output$std_curve_standard_selection <- renderUI({
    
    req(standard_info())
    
    selectInput(
      inputId = "std_curve_standard",
      label = "Standard curve:",
      choices = standard_info(),
      selected = standard_info()[1]
    )
  })
    
    
  # Render the Plate Map Plot
  output$plate_layout_map <- renderPlot({
      req(tagged_data())
      req(input$plate_select_map)
      req(input$plate_size)
      
      # Filter for the selected plate
      df_plot <- tagged_data() %>%
        filter(Plate == input$plate_select_map) %>%
        filter(!is.na(Row_ID) & Row_ID != "" & !is.na(Col_ID))
      
      # Determine Grid Dimensions based on UI selection
      if (input$plate_size == "96") {
        rows_to_show <- LETTERS[1:8] # c("A","B","C","D","E","F","G","H")
        cols_to_show <- 1:12
      } else {
        rows_to_show <- LETTERS[1:16] # A through P
        cols_to_show <- 1:24
      }
      
      # Clean data for plate map plot
      if (input$plate_size == "96") {
        df_plot <- tagged_data() %>%
        filter(Plate == input$plate_select_map) %>%
        filter(Row_ID %in% rows_to_show) %>%   # Strictly A:H
        filter(Col_ID %in% cols_to_show) %>%   # Strictly 1:12
        mutate(Row_ID = factor(Row_ID, levels = rev(rows_to_show)))
      }
      
      # Plate layout map
      ggplot(df_plot %>%
        filter(!is.na(Row_ID), Row_ID != "NA"), aes(x = Col_ID, y = factor(Row_ID, levels = rev(rows_to_show)))) +
        geom_point(aes(color = Well_Type), size = ifelse(input$plate_size == "96", 18, 7)) +
        geom_text(aes(label = Well_ID), size = ifelse(input$plate_size == "96", 6, 3)) +
        scale_x_continuous(breaks = cols_to_show, limits = range(cols_to_show), position = "top") +
        scale_color_manual(values = c(
          "Standard"   = "chocolate3", 
          "Control"    = "darkseagreen", 
          "Background" = "deepskyblue3", 
          "Sample"     = "linen"
        )) +
        theme_minimal() +
        labs(title = paste("Plate Layout:", input$plate_select_map), x = "Column", y = "Row", color = "Sample type") +
        theme(panel.grid.minor = element_blank(),
              legend.position = "bottom",
              axis.title = element_text(size = 17),
              axis.text = element_text(size = 15),
              plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    })
    
  # Data preparation beads - Reactive for bead count
  bead_data <- reactive({
        
        # Make sure a file has been uploaded
        req(input$fileUpload)
        
        # Extract the folder path and the file label, drop .csv extension
        raw_data_path <- dirname(input$fileUpload[1,"datapath"])
        plate_lab <- substr(input$fileUpload$name,1,nchar(input$fileUpload$name)-4)
        
        # Reading uploaded data
        bead_data = list()
        bead_data$rawdatapath <- raw_data_path
        bead_data$uploadedfilenumber <- nrow(input$fileUpload)
        
          ## Load batch data from the uploaded files
          bead_data$plates <- read.batch.beads(path = raw_data_path)
        
          ## Assign user-friendly names to the plates
          names(bead_data$plates) <- plate_lab
        
        # Join plates
          ## Store metadata: plate names, number of plates
          bead_data$platenames <- names(bead_data$plates)
          bead_data$numplates <- length(bead_data$platenames)
          
          ## Combine plates into one dataset
          bead_data$combinedplates <- join.plates(bead_data$plates)
          
          ## Skip NAs
          bead_data$combinedplates <- na.omit(bead_data$combinedplates)
          for (i in seq_along(bead_data$plates)){
            bead_data$plates[[i]] <- na.omit(bead_data$plates[[i]])
        }
        
        ## Extracts antigen names (columns 3 onward = data columns)
        bead_data$ags <- colnames(bead_data$combinedplates)[3:ncol(bead_data$combinedplates)]
        
        ## Metadata with dates - include date & plate info
        bead_data$datesplates <- read.batch(path=raw_data_path,inc_date=T,inc_plate=T)
        bead_data$combineddatesplates <- join.plates(bead_data$datesplates)
        bead_data
        })
  
  # Bead threshold adjustment
    ## Set default beads
    bead_control <- reactiveValues(threshold = 30)
    
    ## Update bead threshold from 30 if the user defines an alternative amount
    observeEvent(input$bead_threshold, {
      # Only update if numeric and not NA
      if (!is.null(input$bead_threshold) && !is.na(as.numeric(input$bead_threshold))) {
        bead_control$threshold <- as.numeric(input$bead_threshold)
      }
    })
    
    ## Update plate selection dropdown
    observe({
      req(bead_data())
      updateSelectInput(
        session,
        "plate_selection",
        choices = names(bead_data()$plates),
        selected = names(bead_data()$plates)[1]
      )
    })
  
    
  # Identify plates with bead counts below the selected threshold
    low_bead_plates <- reactive({
      
      req(bead_data())
      req(input$bead_threshold)
      
      plates <- bead_data()$plates
      threshold <- as.numeric(input$bead_threshold)
      
      low_counts <- vapply(plates, function(plate) {
        
        # Antigen bead-count columns start at column 3
        if (ncol(plate) < 3) {
          return(0L)
        }
        
        bead_values <- unlist(
          lapply(plate[, 3:ncol(plate), drop = FALSE], function(x) {
            suppressWarnings(as.numeric(as.character(x)))
          }),
          use.names = FALSE
        )
        
        sum(
          !is.na(bead_values) &
            is.finite(bead_values) &
            bead_values < threshold
        )
        
      }, integer(1))
      
      data.frame(
        Plate = names(low_counts),
        Low_bead_count_results = as.integer(low_counts),
        stringsAsFactors = FALSE
      )
    })
    
  # Show all plates with low bead counts in the sidebar
  output$low_bead_plates_ui <- renderUI({
      
    flagged <- low_bead_plates()
    flagged_plates <- flagged[flagged$Low_bead_count_results > 0, , drop = FALSE]
      
      if (nrow(flagged) == 0) {
        return(
          p(
            "No plates have bead counts below the current threshold.",
            style = "color: green;"
          )
        )
      }
      
    # Display plate names and the number of low bead-count results
    tags$div(
      lapply(seq_len(nrow(flagged_plates)), function(i) {
        
        tags$div(
          style = paste(
            "padding: 8px 10px;",
            "margin-bottom: 6px;",
            "border-left: 3px solid #dc3545;",
            "color: red;"
          ),
          
          tags$strong(
            paste0("\u26A0 ", flagged_plates$Plate[i])
          ),
          
          tags$p(
            paste(
              flagged_plates$Low_bead_count_results[i],
              "bead-count result(s) below threshold"
            ),
            style = "margin: 4px 0 0 0; font-size: 12px;"
          )
        )
      })
    )
  })
    
    
  # Data preparation - Reactive for standard curve
    std_data <- reactive({
      req(d())
      
      list(
        plates = d()$plates,
        platenames = d()$platenames
      )
    })
    
    # Reactive expression to find all unique matching standard patterns in the uploaded plates
    available_standards <- reactive({
      req(std_data())
      
      identified <- identify.standards(
        data = std_data()$plates,
        dilutions = user_dilutions()
      )
      
      standards <- unique(
        identified$Standard[
          identified$IsStandardCandidate &
            !is.na(identified$Standard) &
            identified$Standard != ""
        ]
      )
      
      sort(standards)
    })
      
  
  # Create a reactive that extracts antigen column names from the uploaded data for app use
  ag <- reactive({
        colnames(d()$plates[[1]])[3:ncol(d()$plates[[1]])]
    })

  # Normalisation
  nd <- reactive({
        # ag_list <- colnames(d()$plates[[1]])[3:ncol(d()$plates[[1]])]
        progress <- Progress$new(session, min=1, max=d()$numplates)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')
        
        # Use normalise_plates function
        normalise_plates(
            all_plates = d()$plates,
            dilutions = user_dilutions(),
            ag_list = ag(),
            fit_type = input$fitting_function,
            progress = progress )
    })

  # Loess adjustment - apply loess regression across plates to smooth out systematic plate-to-plate variation
  ld <- reactive({
        
        # reference plate selected by the user
        req(input$loess_ref_plate)
      
        progress <- Progress$new(session, min=1, max=d()$numplates)
        on.exit(progress$close())

        progress$set(message = 'Calculation in progress',
                    detail = 'This may take a while...')

        d <- loess_adjustment(
            all_plates = d()$plates,
            fitted_data = nd()$data,
            ref_plate = input$loess_ref_plate,
            ag_list = ag(),
            progress = progress
        )
        d$data <- lapply(d$data, as.data.frame)
        d
    })
    
  

  # Upload summary: number of files, number of plates, antigens
        output$upload_summary<- renderText({
              if (!is.null(input$fileUpload)){
                  paste(d()$uploadedfilenumber," files were uploaded, leading to ", d()$numplates," plates being loaded with the following plates and antigens:")
              } else {
                  "Please upload a file using the upload button on the left sidebar."
              }
          })
    
        # Tables of antigen names
        output$antigen_table<-renderTable({
            req(input$fileUpload)
            data.frame(Antigens=d()$ags)
        })
          
        # Tables of Sample names
        output$sample_type_table<-renderTable({
          req(input$fileUpload)
          data.frame(Sample=d()$sample_type)
        })
          
        # Tables of plate names
        output$plate_table<-renderTable({
            req(input$fileUpload)
            data.frame(Plates=d()$platenames)
        })
        
        # Show the sample names and automatically detected standards
        output$sample_names_table <- DT::renderDataTable({
          
          req(std_data())
          
          identified <- identify.standards(
            data = std_data()$plates,
            dilutions = user_dilutions()
          )
          
          DT::datatable(
            identified,
            options = list(
              pageLength = 15,
              scrollX = TRUE
            ),
            rownames = FALSE
          )
        })
          
    output$download_button <- renderUI({
        req(ld())
        downloadButton("download_combined", "Download normalised dataset")
    })

    # After normalisation and adjustment, allows the user to download the final dataset as normdata.csv
    output$download_combined <- downloadHandler(
        
        filename = function() {
        paste0("normdata.csv")
        },
        content = function(file) {
        write.csv(ld()$data, file)
        }
    )

    output$raw_data_files <- renderTable({
        req(input$fileUpload)
        as.data.frame(names(d()$plates))

    })


    output$moreControls <- renderUI({
        req(input$fileUpload)
        tagList(
            selectInput("ref_plate", "Reference plate", d()$platenames),
            selectInput("standard_plate", "Reference plate", d()$platenames),
            selectInput("selected_ag", "Select antigen", d()$ags) )
    })

    output$normControls <- renderUI({
        req(input$fileUpload)
        tagList(
            selectInput("fitting_function","Select fitting function",c("poly","nls","kernsmooth","cgam")),
            renderText("fitting_function_text"),
            selectInput("norm_plot_plate", "Select plate to visualise", d()$platenames), )
    })

    output$fitting_function_text <- renderText({
        req(input$fitting_function)
        if (input$fitting_function=="poly"){
            "The poly method fits a cubic polynomial curve to the data"
        } else if (input$fitting_function=="nls"){
            "The nls method uses non-linear least squares fitting using the nls2 library"
        } else if (input$fitting_function=="kernsmooth"){
            "The kernsmooth method estimates a regression function using local polynomials using the KernSmooth library"
        } else if (input$fitting_function=="cgam"){
            "The cgam method fits a generalised additive model using the cgam library"
        }
    })
  
    output$loessControls1 <- renderUI({
        req(input$fileUpload)   
        tagList(
            selectizeInput(
                inputId = "loess_ref_plate",
                label = "Select reference plate",
                choices = d()$platenames,
                options = list(
                    onInitialize = I('function() { this.setValue(""); }') )
            ) )
    })

    output$loessControls2 <- renderUI({
        req(input$loess_ref_plate)   
        tagList(
            selectizeInput(
                inputId = "loess_plot_plate",
                label = "Select plate to visualise",
                choices = d()$platenames) )
    })

    # Standard Curves
    output$std_table <- DT::renderDataTable({
      
      req(std_data())
      req(user_dilutions())
      
      plates <- std_data()$plates
      dilutions <- user_dilutions()
      
      # Get the standard names submitted on the Upload page
      labels <- submitted_standards()
      
      # Identify matching standards on every plate
      standard_tables <- lapply(names(plates), function(plate_name) {
        
        plate_data <- plates[[plate_name]]
        
        identified <- identify.standards(
          plate_data,
          dilutions
        )
        
        # Retain only identified standard samples
        keep <- identified$IsStandardCandidate &
          !is.na(identified$Standard)
        
        # Restrict to submitted standard names
        if (length(labels) > 0) {
          keep <- keep &
            tolower(identified$Standard) %in% tolower(labels)
        }
        
        if (!any(keep)) {
          return(NULL)
        }
        
        result <- plate_data[keep, , drop = FALSE]
        
        result$Plate <- plate_name
        result$Standard <- identified$Standard[keep]
        result$MatchedDilution <- identified$Dilution[keep]
        
        result$DilutionOrder <- match(
          result$MatchedDilution,
          normalise_dilution(dilutions)
        )
        
        result <- result[
          order(result$Standard, result$DilutionOrder),
          ,
          drop = FALSE
        ]
        
        result
      })
      
      # Remove plates without matching standards
      standard_tables <- Filter(
        Negate(is.null),
        standard_tables
      )
      
      validate(need(
        length(standard_tables) > 0,
        "No standard-curve samples matched your submitted standard names and dilutions."
      ))
      
      # Combine matching samples across all plates
      combined <- do.call(rbind, standard_tables)
      rownames(combined) <- NULL
      
      DT::datatable(
        combined,
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    submitted_standards <- reactiveVal(character(0))
    
    observeEvent(input$std_labelButton, {
      
      labels <- trimws(
        unlist(strsplit(input$std_labelInput, ",", fixed = TRUE))
      )
      
      labels <- unique(labels[nzchar(labels)])
      
      submitted_standards(labels)
      
      output$std_labels <- renderText({
        if (length(labels) == 0) {
          "No standard curves submitted."
        } else {
          paste(
            "Submitted standard curves:",
            paste(labels, collapse = ", ")
          )
        }
      })
    })
    
      
      # Allow user to select up to 5 plates for comparison (compare standard curves visually) 
    output$std_curve_plate_selection <- renderUI({
      req(d())
      selectizeInput(
        inputId = "std_curve_input",
        label = "",
        choices = d()$platenames,
        selected = d()$platenames[1],   # ← defaults to first plate so plot renders immediately
        multiple = TRUE,
        options = list(maxItems = 5)
      )
    })
    
    # Parse standard labels entered by the user (comma or newline separated)
    parsed_standards <- reactive({
      req(input$std_labelInput)
      labs <- trimws(unlist(strsplit(input$std_labelInput, "[,\\n]+")))
      labs[labs != ""]
    })
    
    # Update the dropdown to show the parsed standards instead of antigens
    output$std_curve_standard_selection <- renderUI({
      
      standards <- available_standards()
      
      validate(
        need(
          length(standards) > 0,
          paste(
            "No standard series detected.",
            "Check the sample names and submitted dilutions on the Upload tab."
          )
        )
      )
      
      selectInput(
        inputId = "std_curve_standard_choice",
        label = "",
        choices = standards,
        selected = standards[1]
      )
    })
    
    # Dropdown for antigens
    output$std_curve_antigen_selection <- renderUI({
      
      req(d())
      
      selectInput(
        inputId = "std_curve_antigen",
        label = "",
        choices = d()$ags,
        selected = d()$ags[1]
      )
    })
    
    output$std_curve_plot <- renderPlotly({
      
      req(std_data())
      req(input$std_curve_input)
      req(input$std_curve_antigen)
      req(input$std_curve_standard_choice)
      
      selected_std <- input$std_curve_standard_choice
      selected_plates <- input$std_curve_input
      antigen <- input$std_curve_antigen
      dilutions <- normalise_dilution(user_dilutions())
      
      req(length(dilutions) > 0, !anyNA(dilutions))
      
      # Extract the denominator: 1/50 -> 50, 1/100 -> 100
      dilution_factor <- function(x) {
        vapply(strsplit(x, "/", fixed = TRUE), function(parts) {
          if (length(parts) == 2) {
            as.numeric(parts[2]) / as.numeric(parts[1])
          } else {
            NA_real_
          }
        }, numeric(1))
      }
      
      # Process each selected plate
      df_list <- lapply(selected_plates, function(plate_name) {
        
        pd <- std_data()$plates[[plate_name]]
        
        if (is.null(pd) || !"Sample" %in% names(pd)) {
          return(NULL)
        }
        
        if (!antigen %in% names(pd)) {
          return(NULL)
        }
        
        # Identify standards using the user's dilution list
        identified <- identify.standards(pd, dilutions)
        
          print(data.frame(
            Sample = pd$Sample,
            Standard = identified$Standard,
            Dilution = identified$Dilution,
            Candidate = identified$IsStandardCandidate
          ))
          
          print(paste("Selected standard:", selected_std))
          print(paste("Entered dilutions:", paste(dilutions, collapse = ", ")))
        
        # Keep ONLY the selected standard and requested dilutions
          print(unique(identified$Standard))
          print(unique(identified$Dilution))
          print(selected_std)
          
          keep <- !is.na(identified$Standard) &
            tolower(trimws(identified$Standard)) ==
            tolower(trimws(selected_std)) &
            identified$Dilution %in% dilutions
          
          print(table(keep, useNA = "ifany"))
        
        if (!any(keep)) {
          return(NULL)
        }
        
        matched_dilutions <- identified$Dilution[keep]
        
        data.frame(
          Plate = plate_name,
          Sample = as.character(pd$Sample[keep]),
          Standard = identified$Standard[keep],
          Dilution = matched_dilutions,
          DilutionNumeric = dilution_factor(matched_dilutions),
          MFI = suppressWarnings(
            as.numeric(as.character(pd[[antigen]][keep]))
          ),
          stringsAsFactors = FALSE
        )
      })
      
      # Remove plates without matching samples
      df_list <- Filter(
        function(x) !is.null(x) && nrow(x) > 0,
        df_list
      )
      
      validate(need(
        length(df_list) > 0,
        paste0(
          "No samples matched standard '", selected_std,
          "' and the entered dilutions: ",
          paste(dilutions, collapse = ", ")
        )
      ))
      
      df <- do.call(rbind, df_list)
      
      # Retain valid positive measurements and dilution factors
      df <- df[
        !is.na(df$MFI) &
          is.finite(df$MFI) &
          df$MFI > 0 &
          !is.na(df$DilutionNumeric) &
          is.finite(df$DilutionNumeric) &
          df$DilutionNumeric > 0,
        ,
        drop = FALSE
      ]
      
      validate(need(
        nrow(df) > 0,
        "No valid positive MFI and dilution values are available to plot."
      ))
      
      # Order points according to the user's entered dilutions
      df$DilutionOrder <- match(df$Dilution, dilutions)
      
      df <- df[
        order(df$Plate, df$DilutionOrder),
        ,
        drop = FALSE
      ]
      
      df$LogDilution <- log10(df$DilutionNumeric)
      df$LogMFI <- log10(df$MFI)
      
      p <- ggplot(
        df,
        aes(
          x = LogDilution,
          y = LogMFI,
          color = Plate,
          group = Plate,
          text = paste0(
            "Plate: ", Plate,
            "<br>Sample: ", Sample,
            "<br>Standard: ", Standard,
            "<br>Dilution: ", Dilution,
            "<br>MFI: ", round(MFI, 1)
          )
        )
      ) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        labs(
          title = paste(
            "Standard Curve -", selected_std, "-", antigen
          ),
          x = "Log10 Dilution Factor",
          y = "Log10 MFI",
          color = "Plate"
        ) +
        theme_minimal(base_size = 13) +
        theme(
          axis.text.x = element_text(angle = 35, hjust = 1)
        )
      
      ggplotly(p, tooltip = "text", height = 550)
    })
          

    # QC
    ## Levey-Jennings controls for tracking high, mid, and low dilution points across plates
    output$levey_jennings_select_points <- renderUI({
        tagList(
            selectInput("levey_high", "Select high dilution", user_dilutions(), selected="1/10"),
            selectInput("levey_mid", "Select mid dilution", user_dilutions(), selected="1/50"),
            selectInput("levey_low", "Select low dilution", user_dilutions(),selected="1/250") )
    })

    ## Plots coefficient of variation across dilutions to see consistency of replicates
    output$coefv_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            column(12,
                renderPlotly({
                    plot.coefv(
                        data = d()$plates,
                        ag(),
                        std_label = std_labels(),
                        dilutions=user_dilutions(),
                        target_dilution = "1/250" ) 
                })  ),
            column(12,
                renderPlotly({
                    plot.coefv(
                        data = d()$plates,
                        ag(),
                        std_label = std_labels(),
                        dilutions=user_dilutions(),
                        target_dilution = "1/1250" )
                }) )
        ) })

    ## Generate Levey-Jennings plots for each antigen
    output$levey_jennings_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(antigen) {
                column(6,
                    renderPlot({
                        levey.jennings(
                            data = d()$combineddatesplates,
                            std_label = std_labels(),
                            blank_label = bkg_labels, # "Background0",
                            dil_high = input$levey_high,
                            dil_mid = input$levey_mid,
                            dil_low = input$levey_low,
                            by_var = "plate",
                            ag_list = c(antigen) )
                }) )
            }) )
    })
    
  # Generate low bead plot for each antigen
    output$bead_count_ui_dynamic <- renderUI({
      req(bead_data())
      req(input$plate_selection)
      
      n_ags <- length(bead_data()$ags)
      n_rows <- ceiling(n_ags / 2)
      plot_height <- max(600, n_rows * 220)
      
      plotlyOutput("bead_count_plots", height = paste0(plot_height, "px"))
    })
    
    output$bead_count_plots <- renderPlotly({
      req(bead_data())
      req(input$plate_selection)
      
      plate_selected <- bead_data()$plates[[input$plate_selection]]
      
      n_ags <- length(bead_data()$ags)
      n_rows <- ceiling(n_ags / 2)
      plot_height <- max(400, n_rows * 220)
      
      p <- plot_beads_by_well(
        plate_selected,
        plate_name = input$plate_selection,
        threshold = bead_control$threshold
      )
      
      ggplotly(p, height = plot_height, tooltip = c("Well", "BeadCount"))
    })
    
  # Normalised Data & Plots
    ## Shows normalised matrix
    output$normalised_matrix<-renderTable({
        nd()$data[["MP01"]]
    })
    
    ## Produces normalisation diagnostic plots per antigen and plate
    output$normalisation_plots <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(a) {
                column(4,
                    renderPlot({
                        nd()$plots[[input$norm_plot_plate]][[a]]

                    }) )
            }) )
    })

    # Loess Adjustment Visualisation - Compares raw vs loess-adjusted plots for each antigen and plate
    output$loess_ag_plot <- renderUI({
        req(input$fileUpload)
        fluidRow(
            lapply(ag(), function(a) {
                column(6,
                    renderPlot({
                        par(mfrow=c(1,2))

                        # Compares raw vs loess-adjusted antigen plots
                        plot.ag.plates.raw(data=d()$plates, ag=a, dilutions=user_dilutions())
                        plot.ag.plates.raw(data=ld()$data, ag=a, dilutions=user_dilutions())
                        par(mfrow=c(1,1))
                    }) )
            }) )
    })

    # Loess-adjusted diagnostic plots
    output$loess_plots <- renderUI({
        req(input$fileUpload)
        req(input$loess_ref_plate)
        fluidRow(
            lapply(ag(), function(a) {
                column(6,
                    renderPlot({
                        ld()$plots[[input$loess_plot_plate]][[a]]

                    }) )
            }) )
    })

    # Display loess-adjusted dataset as a table
    output$loess_data <- renderTable({
        ld()$data[[1]]
    })}


    # # Downloading cleaned data
    # output$cleaned_data <- renderDataTable({
    #   
    #   # Remove background
    #   
    #     ## Create variable that is the selected background removal
    #     bkg_select <- reactive({
    #       input$backgroundremoval
    #     })

          ## Create cleaned dataset

    # })
  


# Run the app ----
shinyApp(ui = ui, server = server)

# rsconnect::deployApp(
#   appName = "Luminex_Data_Analysis"
# )