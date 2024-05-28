# UI logic for lipid heatmap #3
# This section of code only describes the Lipid Heatmap functions in MetaboLink, and do not work independently from all the other code done in MetaboLink. 
# Meaning if testing of Lipid Heatmap is wanted, go to the test server linked in the B.Sc. ‘Methods and Materials’ section.



# Lipid Heatmap 
tabPanel("Lipid Heatmap",
         useShinyjs(), # Enable Shiny JavaScript operations for enhanced UI interactions.
         tabsetPanel(
           # Heatmap Visualization tab: Provides UI elements for processing and visualizing heatmap data.
           tabPanel("Heatmap Visualization",
                    actionButton("run_process", "Run Data Processing"), # Button to trigger data processing.
                    actionButton("show_help", "Show User Guide"), # Button to display a user guide or help information.
                    radioButtons("selected_dataset", "Select data frame:", # Radio buttons for selecting the dataset type.
                                 choices = c("Original Data" = "original", "Merged Data" = "merged"),
                                 selected = "original"),
                    
                    
                    textOutput("upload_message"), # Message area for displaying upload status or instructions.
                    column(width = 4,  
                           # Group selection
                           uiOutput("select_group_ui_heatmap"), 
                           tableOutput("numerator_table"),
                           tableOutput("denominator_table"),
                    ),
                    column(width = 4,
                           uiOutput("select_lipid_ui"), # lipid selection
                    ),
                    column(width = 4,  
                    
                           uiOutput("p_value_max_ui"), # P-value threshold
                           uiOutput("logFC_input_ui"), # logFC-value threshold
                    ),
                    uiOutput("lipid_display_message"), # Not in use, for future development. Will display is users selected lipids are not in between threshold.
                    
                    column(width = 12,
                           plotOutput("heatmapPlot", width = "100%", height = "650px") # display lipid heatmap
                    ),
                    
                    # Table in its own row, below the plot
                    column(width = 12,
                           dataTableOutput("pValueTable")
                    ), # Lines below in # are for testing are making sure of corrcet data. 
                    uiOutput("select_group_ui"), # Dynamic UI for selecting groups, populated server-side.
                    tableOutput("selected_group_table"), # Displays selected groups for verification.
           ),
           
           # Data of groups in Heatmap tab: Displays the data behind the groups used in the heatmap.
           tabPanel("Data of groups in Heatmap",
                    uiOutput("table_message"), # Dynamic message about the data table, populated server-side.
                    conditionalPanel(
                      condition = "input.run_process > 0", # Only display the following if data processing has been triggered.
                      tableOutput("groups_table"), # Shows table of groups.
                      textOutput("limitation_notice"), # Notice about any limitations or considerations.
                      textOutput("rows_removed_text"), # Information on data rows removed during processing.
                      tableOutput("raw_data_table") # Displays the raw data table.
                    ),
                    verbatimTextOutput("error_message") # Area to display any error messages.
           )
         )
),
