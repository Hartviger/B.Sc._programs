library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
#library(shinyalert)
library(shinyWidgets)
library(spsComps)
library(DT)
library(dplyr)
library(plotly)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(impute)
library(randomForest)
library(writexl)
library(igraph)
library(stringi)
library(BiocManager)
library(shinycssloaders)
library(jsonlite)
library(reshape2)
library("lipidomeR") # new
library(ggh4x) # new


options(repos = BiocManager::repositories())
source("functions.R")
source("plotfunctions.R")
source("analysis_functions.R")
source("normalization.R")

shinyUI(dashboardPage(
  dashboardHeader(
    title = tags$div(
      tags$img(  src = "/Users/jakobhartvig/Desktop/uni/7_semester/ISA_kandidat/September/MetaboLink_test_of_test_server/logo.png",
                 style = "height: auto; width: auto; max-height: 30px; max-width: 53px; vertical-align: middle; margin-right: 10px;"
      ), # Adjust path and styling
      "MetaboLink TEST SERVER"
    ),
    titleWidth = 400,
    dropdownMenu(type = "notifications",
                 icon = icon("question-circle"),
                 badgeStatus = NULL,
                 headerText = "Help",
                 notificationItem("User manual", icon = icon("book"),
                                  href = "https://github.com/anitamnd/jlspec_2_0/wiki"),
                 notificationItem("Source code and installation", icon = icon("file"),
                                  href = "https://github.com/anitamnd/jlspec_2_0"),
                 notificationItem("Institution", icon = icon("university"),
                                  href = "https://www.sdu.dk/en")
    )
  ),
  
  
  # Sidebar 
  
  dashboardSidebar(
    width = "400",
    useShinyjs(),
    tags$style(HTML(".panel-primary {color: #000000;}")),
    tags$head(tags$script(src="CallShiny.js")),
    extendShinyjs(script="CallShiny.js", functions=c("retrieve_results","send_message","run_button")),
    fluidPage(
      fluidRow(
        selectizeInput("selectDataset", "Active dataset",
                       choices = NULL, width = "100%",
                       options = list(placeholder = "Please upload a file to start")
        )
      ),
      bsCollapse(
        id = "menu", multiple = FALSE, open = "Data input",
        bsCollapsePanel("Data input",
                        style = "primary",
                        fluidRow(style = "padding: 0px;",
                                 selectInput("fileType", "Select file format",
                                             choices = c("Samples in columns", "Samples in rows"),
                                             selected = "Samples in columns",
                                             width = "100%"
                                 )
                        ),
                        fluidRow(style = "padding: 0px;",
                                 fileInput("inputFile", "Upload file (.txt or .csv)",
                                           accept = c("txt/csv", "text/comma-seperated-values, text/plain", ".csv"),
                                           width = "100%"
                                 )
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, style = "padding-left:0px;",
                                 bsButton("submit", "Submit", width = "100%")),
                          column(6, style = "padding-left:0px;",
                                 bsButton("example", "Load example", width = "100%"))
                        )
        ),
        bsCollapsePanel("Blank filtration",
                        style = "primary",
                        fluidRow(
                          style = "padding: 0px;",
                          column(12,
                                 sliderInput("signalStrength", "Signal strength above blank", 0, 10, 5, step = 0.1, width = "100%"),
                                 style = "padding: 0px"
                          )
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(12,
                                 checkboxInput("discardBlank", "Discard blank", value = T, width = "100%"),
                                 style = "padding: 0px; margin-top: -25px; margin-bottom: -15px; margin-left: 5px;"
                          ),
                          column(12,
                                 checkboxInput("keepIS", "Keep internal standards", value = T, width = "100%"),
                                 style = "padding: 0px; margin-top: -15px; margin-bottom: -15px; margin-left: 5px;"
                          ),
                          column(12,
                                 checkboxInput("newFileBF", "Save as new file", value = T, width = "100%"),
                                 style = "padding: 0px; margin-top: -15px; margin-bottom: -15px; margin-left: 5px;"
                          )
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("blankFiltrate", "Blank filtrate", width = "100%"),
                                 style = "padding-left:0px; margin-top: 10px;"
                          ),
                          column(6, bsButton("saveBF", "Save", width = "100%"),
                                 style = "padding-left:0px; margin-top: 10px;"
                          )
                        )
        ),
        bsCollapsePanel("Missing value filtration",
                        style = "primary",
                        fluidRow(
                          style = "padding: 0px;",
                          column(12,
                                 sliderInput("cutoffNAs", "Minimum percentage of messured values", 0, 100, 80, step = 5, width = "100%"),
                                 style = "padding: 0px"
                          )
                        ),
                        fluidRow(
                          column(6,
                                 prettyCheckboxGroup("filterNAmethod", "", choices = c("in QC", "in group", "entire data"))
                          )),
                        fluidRow(
                          column(6,
                                 checkboxInput("mvf_newsave", "Save as new file", value = T, width = "100%")
                          )),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("runFilterNA", "Run", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("saveFilterNA", "Save", width = "100%"), style = "padding-left:0px;")
                        )
        ),
        bsCollapsePanel("Imputation",
                        style = "primary",
                        fluidRow(selectInput("imputationMethod", "Imputation method", choices = c("KNN", "Min/X", "Median"), width = "100%")),
                        fluidRow(hidden(div(id = "imp_remaining_hide", selectInput("remainingNAs", "Remaining missing values", choices = c("zero", "Min/X", "Median"), width = "100%")))),
                        fluidRow(hidden(div(id = "imp_minx_hide", sliderInput("imputationMinX", "Divide min by", min = 1, max = 10, value = 1, step = 1, width = "100%")))),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, prettyCheckbox("imp_onlyQC", "Only imputate QC")),
                          column(6)
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6,
                                 checkboxInput("newFileImp", "Save as new file", value = T, width = "100%"), 
                                 style = "padding: 0px; margin-top: -10px; margin-left: 10px; margin-right: -10px;")
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("runImputation", "Run", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("saveImputation", "Save", width = "100%"), style = "padding-left:0px;")
                        )
        ),
        bsCollapsePanel("IS normalization",
                        style = "primary",
                        fluidRow(
                          selectInput("isMethod", "Method", choices = c("Nearest RT", "Same lipid structure"), selected = "Nearest RT", width = "100%")
                        ),
                        fluidRow(
                          checkboxGroupInput("isChoose", NULL, choices = NULL, selected = NULL, inline = FALSE)
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, checkboxInput("normalizeQC", "Normalize QC", value = T, width = "100%"), style = "padding: 0px; margin-top: 0px; margin-left: 10px; margin-right: -10px;"),
                          column(6, checkboxInput("newFileIS", "Save as new file", value = T, width = "100%"), style = "padding: 0px; margin-top: 0px; margin-left: 10px; margin-right: -10px;")
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("normalizeIS", "Normalize", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("optimizeIS", "Optimize", width = "100%"), style = "padding-left:0px;")
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("removeIS", "Remove IS", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("saveIS", "Save", width = "100%"), style = "padding-left:0px;")
                        )
        ),
        bsCollapsePanel("Drift correction",
                        style = "primary",
                        fluidRow(selectInput("driftMethod", "Signal correction method", choices = c("QC-RFSC (random forrest)", "QC-RLSC (robust LOESS)"), width = "100%")),
                        fluidRow(div(id = "dc_ntree_hide", sliderInput("driftTrees", "ntree", min = 100, max = 1000, value = 500, step = 100, width = "100%"))),
                        fluidRow(hidden(div(id = "dc_qcspan_hide", sliderInput("driftQCspan", "QCspan", min = 0.2, max = 0.75, value = 0.7, step = 0.05, width = "100%")))),
                        fluidRow(hidden(div(id = "dc_degree_hide", sliderInput("driftDegree", "degree", min = 0, max = 2, value = 2, step = 1, width = "100%")))),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(12, checkboxInput("newFileDrift", "Save as new file", value = T, width = "100%"), style = "padding: 0px; margin-top: -10px; margin-left: 10px; margin-right: -10px;"),
                          column(6, bsButton("runDrift", "Run", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("saveDrift", "Save", width = "100%"), style = "padding-left:0px;")
                        )
        ) %>% bsTooltip("For correct usage, please refer to the user manual.", placement = "bottom", trigger = "hover"),
        bsCollapsePanel("Normalization",
                        style = "primary",
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(12, selectInput("normMethod", "Function", choices = c("QC (PQN)", "Sum", "Median", "Sample amount"), width = "100%"), style = "padding-left:0px;")
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(12, checkboxInput("newFileNorm", "Save as new file", value = F, width = "100%"), style = "padding: 0px; margin-top: -10px; margin-left: 10px; margin-right: -10px;"),
                          column(6, bsButton("normalize", "Run", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("saveNormalization", "Save", width = "100%"), style = "padding-left:0px;")
                        )
        ),
        bsCollapsePanel("Log transform and scaling",
                        style = "primary",
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, style = "padding-left:0px;",
                                 selectInput("logTransform", "Log transform", choices = c("None", "log2", "log10", "ln"), width = "100%") %>%
                                   bsTooltip("Do not use log transformation for negative values.", placement = "top", trigger = "hover")
                          ),
                          column(6, selectInput("scaling", "Data scaling", choices = c("None", "Mean center", "Auto scale"), width = "100%"), style = "padding-left:0px;")
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("transform", "Run", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("saveTransform", "Save", width = "100%"), style = "padding-left:0px;")
                        )
        ) %>% bsTooltip("Some features will be removed after transforming the data if Inf values are introduced.", placement = "bottom", trigger = "hover"),
        bsCollapsePanel("Merge datasets",
                        style = "primary",
                        fluidRow(selectInput("mergeFile", "Select dataset to merge with", choices = NULL, width = "100%")),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, numericInput("merge_ppm", "M/z tolerance ppm", min = 0, value = 10, width = "100%"), style = "padding-left:0px;"),
                          column(6, numericInput("merge_rt", "RT tolerance", min = 0, value = 0.1, step = 0.01, width = "100%"), style = "padding-left:0px;")
                        ),
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(6, bsButton("mergeRankings", "Edit priorities", width = "100%"), style = "padding-left:0px;"),
                          column(6, bsButton("mergeDatasets", "Run", width = "100%"), style = "padding-left:0px;")
                        )
        ) %>% bsTooltip("Merge datasets with same samples and different ion mode. The datasets must have the same number of samples.", placement = "bottom", trigger = "hover"),
        bsCollapsePanel("Remove files",
                        style = "primary",
                        fluidRow(
                          style = "margin-right: 0px;",
                          column(12, checkboxGroupInput("filesToRemove", "Select files to remove", choices = NULL, selected = NULL), style = "padding-left:0px;"),
                          column(12, bsButton("removeFiles", "Remove", width = "50%"), style = "padding-left:0px;")
                        )
        )
      )
    ),
    sidebarMenu()
  ),
  
  # Main Body
  
  dashboardBody(
    tags$head(tags$style(".modal-sm{ width:300px}
                         .modal-lg{ width:1200px}")),
    useShinyjs(),
    fluidRow(hidden(div(
      id = "buttons", style = "padding-bottom: 49px",
      column(3, bsButton("sequence",
                         label = "Sequence",
                         icon = icon("tags"),
                         style = "default",
                         block = T
      )),
      column(3, bsButton("explore",
                         label = "Explore data",
                         icon = icon("table"),
                         style = "default",
                         block = T
      )),
      column(3, bsButton("statistics_button",
                         label = "Data analysis",
                         icon = icon("clipboard"),
                         style = "default",
                         block = T
      )),
      tags$style(type = "text/css", "#plot2 {width:100%}"),
      column(3, bsButton("export",
                         label = "Export",
                         icon = icon("download"),
                         style = "default",
                         block = TRUE
      )),
    ))),
    fluidRow(
      hidden(
        div(
          id = "sequence_panel",
          column(12, box(width = NULL,
                         strong("User guide"),
                         p("Make sure the columns are labeled correctly before proceeding."),
                         p("If you see a sample column labeled '-', this usually means there are invalid characters in the column."),
                         p("Labels cannot be edited in the app to avoid crashes. Please edit the file and re-upload."))
          ),
          column(
            width = 8,
            box(
              title = textOutput("diboxtitle"), width = NULL,
              DTOutput("seq_table") %>% withSpinner(color="#0A4F8F")
            )
          ),
          column(
            width = 4,
            box(
              width = NULL, title = "Upload sequence file", status = "danger",
              fileInput("inputSequence", "Select file", accept = c("txt/csv", "text/comma-seperated-values,text/plain", ".csv"), width = "100%"),
              column(6, actionButton("updateSequence", label = "Update", width = "100%")), 
              column(6, actionButton("reuseSequence", label = "Re-use sequence", width = "100%"))
            ),
            box(
              width = NULL, title = "Edit data columns",
              actionButton("editSequence", "Edit", width = 100)
            ),
            box(
              width = NULL, title = "Group nicknames",
              actionButton("editGroups", "Edit", width = 100)
            ),
            box(
              width = NULL, title = "Download sequence file",
              downloadButton("downloadSequence", "Download", width = 100)
            )
          )
        )
      )
    ),
    fluidRow(
      hidden(
        div(
          id = "datatable_panel",
          tabsetPanel(
            tabPanel("Data table",
                     fluidRow(
                       column(12, box(width = NULL, DTOutput("dttable") %>% withSpinner(color="#0A4F8F")))
                     )
            ),
            tabPanel("Sample distribution", 
                     fluidRow(
                       column(12,
                              box(width = NULL, h4("Median across samples"), plotlyOutput("histogram")),
                       ),
                       column(12,
                              box(width = NULL, h4("Median across QCs"), uiOutput("histogram_qc"))
                       )
                     )
            ),
            tabPanel("PCA", 
                     fluidRow(
                       column(6, box(width = NULL, 
                                     selectInput("selectpca1", "", choices = NULL, width = "100%"),
                                     checkboxInput("pca1_islog", "Is data log-transformed?", value = F, width = "100%"),
                                     actionButton("run_pca1", "Run PCA", width = "50%"),
                                     plotlyOutput("plotpca1"), br(),
                                     htmlOutput("pca1Details")
                       )),
                       column(6, box(width = NULL, 
                                     selectInput("selectpca2", "", choices = NULL, width = "100%"),
                                     checkboxInput("pca2_islog", "Is data log-transformed?", value = F, width = "100%"),
                                     actionButton("run_pca2", "Run PCA", width = "50%"),
                                     plotlyOutput("plotpca2"), br(),
                                     htmlOutput("pca2Details")
                       ))
                     ),
            ),
            
            
            
            
            
            
            
            
            
            
            tabPanel("Lipid Heatmap",
                     useShinyjs(),
                     tabsetPanel(
                       
                       tabPanel("Group selection",
                                fluidRow(
                                  column(6, box(
                                    width = NULL,
                                    title = "User guide",
                                    tags$ul(
                                      tags$li("Select the data frame to use."),
                                      tags$li("Click 'Run Data Processing' to start the analysis."),
                                      tags$li("After data processing select the numerator and denominator groups."),
                                      conditionalPanel(
                                        condition = "input.run_process > 0", # Only display the following if data processing has been triggered.
                                        tableOutput("groups_table"), # Shows table of groups.
                                      )
                                    )
                                  )),
                                  column(6, box(
                                    width = NULL, 
                                    title = "Select Data Frame", 
                                    radioButtons("selected_dataset", "Select data frame:",
                                                 choices = c("Original Data" = "original", "Merged Data" = "merged"),
                                                 selected = "original", width = "50%"),
                                    
                                    column(6, 
                                           actionButton("run_process", "Run Data Processing"),
                                           uiOutput("select_group_ui_heatmap"),
                                           tableOutput("numerator_table"),
                                           tableOutput("denominator_table"),
                                           conditionalPanel(
                                             condition = "input.run_process > 0", 
                                             actionButton("show_lipid_info", "Show Lipid Summary")
                                           )
                                    )
                                  ))
                                ),
                                
                                conditionalPanel(
                                  condition = "input.run_process > 0", 
                                  column(6, box(
                                    title = "Numerator Group Table",
                                    width = NULL,
                                    DT::dataTableOutput("numerator_group_table")  
                                  )),
                                  column(6, box(
                                    title = "Denominator Group Table",
                                    width = NULL,
                                    DT::dataTableOutput("denominator_group_table")  
                                  ))
                                )
                       ),
                       
                       # HEATMAP VISUALIZATION
                       tabPanel("Lipid Visualization",
                                box(
                                  column(width = 8,
                                         uiOutput("select_lipid_ui"),
                                         tags$li("Lipids will show up both on Lipid Heatmap and Bubble plot if they are within the threshold of p-value and logFC."),
                                         tags$li("Be aware of isoforms of each lipid.")
                                  )
                                ),
                                box(
                                  column(width = 4, 
                                         uiOutput("p_value_max_ui"),
                                         uiOutput("logFC_input_ui")
                                  ),
                                  column(width = 4, 
                                         # Note: `uiOutput("logFC_input_ui")` appears twice; ensure this is intentional.
                                  )
                                ),
                                
                                div(style = "width: 100%; height: 2000px; overflow-y: scroll;",  # Adjust overflow and height
                                    uiOutput("heatmap_ui")
                                ),   
                                
                                column(width = 12,
                                       dataTableOutput("pValueTable")
                                )
                                
                       ),
                       
                       tabPanel(
                         "Bubble plot of data",
                         #checkboxInput("swap_logfc_pvalue", "Swap logFC and p-value", value = FALSE),
                         #helpText("Bubble size represents p-value ranges. Color indicates logFC."),
                         uiOutput("bubble_plot_ui"),
                         
                       )
                       
                       
                       
                       
                       
                     )
            ),
            
            
            
            
            
            
            
            
            tabPanel("Feature drift",
                     fluidRow(
                       column(3, box(width = NULL, DTOutput("dt_drift_panel"))),
                       column(9, 
                              box(
                                width = NULL,
                                fluidRow(
                                  column(4, selectizeInput("drift_select", "Select dataset to compare with", choices = NULL, width = "100%", options = list(placeholder = "Select file"))),
                                  column(2, style = "margin-top: 25px;", bsButton("drift_1", label = "Individual", block = TRUE)),
                                  column(2, style = "margin-top: 25px;", bsButton("drift_2", label = "CV variation", block = TRUE)),
                                  column(2, style = "margin-top: 25px;", bsButton("drift_3", label = "CV distribution", block = TRUE))
                                ),
                              ),
                              uiOutput("drift_ui")
                       )
                     )
            ),
            tabPanel("Feature viewer",
                     fluidRow(
                       column(3, box(width = NULL, DTOutput("dt_boxplot_panel"))),
                       column(9, box(width = NULL, 
                                     fluidRow(
                                       column(6,
                                              radioButtons(
                                                inputId = "bloxplot_log",
                                                label = "Log",
                                                choices = c("None", "ln", "log2", "log10"),
                                                selected = "None",
                                                inline = TRUE
                                              )
                                       ),
                                       column(6, radioButtons(
                                         inputId = "bloxplot_ylog",
                                         label = "Y axis log",
                                         choices = c("None", "log2", "log10"),
                                         selected = "None",
                                         inline = TRUE
                                       ))
                                     )),
                              uiOutput("boxplot_ui")
                       )
                     )
            ),
            tabPanel("Summary",
                     box(width = NULL,
                         fluidRow(
                           column(12,htmlOutput("title")),
                         ),
                         fluidRow(
                           column(6, uiOutput("info_ui")),
                           column(6, htmlOutput("cvinfo_ui"))
                         )
                     )
            )
          )
        )
      )
    ),
    fluidRow(
      hidden(
        div(
          id = "statistics_panel",
          fluidPage(
            fluidRow(
              box(width=NULL, column(5,
                                     id = "pr_c3",
                                     h4("Analysis parameters"), 
                                     fluidRow(
                                       column(12, p("Remember to log-transform and scale data before running tests.")),
                                       column(12, selectInput("testType", "Select test", choices = c("2 group comparison (unpaired)" = "GroupsUnpaired",
                                                                                                     "2 group comparison with multiple time points (paired)" = "GroupsMultipleTime",
                                                                                                     "Compare to reference group" = "CompareToReference"), selected = NULL, 
                                                              width = "100%"))
                                       # calculate fold change as the ratio between 2 group means?                                            
                                     ),
                                     conditionalPanel(
                                       condition = "input.testType == 'GroupsUnpaired'",
                                       fluidRow(
                                         column(6, selectInput("group1", "Group", choices = NULL, width = "100%")),
                                         column(6, selectInput("time1", "Time", choices = NULL, width = "100%"))
                                       ),
                                       fluidRow(
                                         column(6, selectInput("group2", "Group", choices = NULL, width = "100%")),
                                         column(6, selectInput("time2", "Time", choices = NULL, width = "100%"))
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.testType == 'CompareToReference'",
                                       fluidRow(
                                         column(12, selectInput("referenceGroup", "Select reference group", choices = NULL, width = "100%"))
                                       )
                                     ),
                                     conditionalPanel(
                                       condition = "input.testType == 'GroupsMultipleTime'",
                                       fluidRow(
                                         column(12, checkboxGroupInput("contrasts", "Select contrasts", choices = NULL, selected = NULL, inline = FALSE))
                                       )
                                     ),
                                     fluidRow(
                                       column(12, actionButton("selectTest", "Run test", width = "40%", style="float:right; margin-right: 0px;"))
                                     )    
              )),
              column(6,
                     id = "pr_c2",
                     h4("User guide"),
                     p("")
              )
            ),
            br(),
            h4("Results"),
            uiOutput("results_ui")
          )
        )
      )
    ),
    fluidRow(
      hidden(
        div(
          id = "export_panel",
          box(title = strong(".csv and .xlsx"), width = 6,
              column(12, 
                     h4(".csv"),
                     uiOutput("export_ui")
              ),
              column(12, style = "margin-top: 20px;",
                     h4(".xlsx"),
                     checkboxGroupInput("export_xml_list", "Choose sheets", choices = NULL, selected = NULL),
                     downloadButton("export_xml", "Export combined .xlsx")
              ),
              column(12, style = "margin-top: 20px;",
                     h4("Statistics results"),
                     uiOutput("export_stats")
              ),
              column(12, style = "margin-top: 20px;",
                     h4("Settings used in app"),
                     uiOutput("export_settings")
              )
          ),
          box(title = strong("Export to other apps"), width = 6,
              column(12, 
                     h4("Statistical testing"),
                     actionButton("send_polystest", "Send to PolySTest"),
                     span(textOutput("connection_polystest"), style="color:#33DD33;")
              ),
              column(12, style = "margin-top: 20px;",
                     h4("Clustering"),
                     actionButton("send_vsclust", "Send to VSClust"),
                     span(textOutput("connection_vsclust"), style="color:#33DD33;")
              ),
              column(12, style = "margin-top: 20px;",
                     h4(".csv for MetaboAnalyst"),
                     uiOutput("export_metabo")
              )
          )
        )
      )
    ),
    fluidRow(
      div(
        id = "welcome_panel",
        column(12, box(width = NULL, includeHTML("intro_text.html")))
      )
    )
  )
))