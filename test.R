library(shiny)
library(DT)

server <- function(input, output, session) {
  # Reactive for reading uploaded data immediately
  uploaded_data <- reactiveFileReader(1000, session, input$dataFile$datapath, function(filePath) {
    read.csv(filePath, stringsAsFactors = FALSE)
  })
  
  # Display uploaded data immediately
  output$uploadedTable <- renderDT({
    req(uploaded_data())
    uploaded_data()
  }, options = list(pageLength = 50))
  
  # Processed data after clicking 'Process Data'
  data_processed <- eventReactive(input$process, {
    req(input$dataFile)  # Ensure a file is uploaded
    
    # Read the data from the uploaded file
    Data_in <- read.csv(input$dataFile$datapath, stringsAsFactors = FALSE)
    # Set the number of eqQCs, MSMSs, and QCs you need
    num_eqQCs <- input$num_eqQCs # Adjust the number of eqQC rows as needed
    num_MSMSs <- input$num_MSMSs # Adjust the number of MSMS rows as needed
    num_QCs <- input$num_QCs # Adjust the number of QC rows as needed
    insert_after_samples <- input$insert_after_samples # Change this number as needed
    
    
    
    
    # Split the data into subsets: blanks, QCs, and the rest
    blanks <- Data_in[tolower(Data_in$Sample.type) == "blank", ]
    qcs <- Data_in[tolower(Data_in$Sample.type) == "qc", ]
    samples <- Data_in[tolower(Data_in$Sample.type) == "sample", ]
    non_samples <- Data_in[!(Data_in$Sample.type %in% c("Blank", "QC", "Sample")), ]
    existing_qcs <- Data_in[tolower(Data_in$Sample.type) == "qc", ]
    
    # Name to the   # Generate new eqQC names... 
    sample_type_eqQC_names <- paste0(existing_qcs$name, "_eq_QC")
    sample_type_MSMS_names <- paste0(existing_qcs$name, "_MSMS")
    new_QC_names <- paste0(existing_qcs$name, "_QC")
    new_QC_samples_names <- paste0(existing_qcs$name, "_QC")
    
    new_eqQCs <- qcs[rep(1, num_eqQCs), ]
    
    
    # Check if there are existing QCs to duplicate
    if(nrow(qcs) > 0) {
      # Create new eqQC rows by replicating the existing ones
      new_eqQCs <- qcs[rep(1, num_eqQCs), ]
      
      # Generate new eqQC names
      new_eqQC_names <- paste(sample_type_eqQC_names, sprintf("%02d", seq_len(num_eqQCs)), sep="")
      new_eqQCs$name <- new_eqQC_names
      
      # Create new MSMS rows by replicating the existing QCs
      new_MSMSs <- qcs[rep(1, num_MSMSs), ]
      
      # Generate new MSMS names
      new_MSMS_names <- paste(sample_type_MSMS_names, sprintf("%02d", seq_len(num_MSMSs)), sep="")
      new_MSMSs$name <- new_MSMS_names
      
      # Create new QC rows
      new_QCs <- qcs[rep(1, num_QCs), ]
      
      # Generate new QC names
      new_QC_names <- paste(new_QC_names, sprintf("%02d", seq_len(num_QCs)), sep="")
      new_QCs$name <- new_QC_names
      
      # Combine the new eqQCs, MSMSs, and QCs with the blanks
      qcs_combined <- rbind(new_eqQCs, new_MSMSs, new_QCs)
      
    } else {
      stop("No QC samples available to duplicate.")
    }
    
    
    # Randomize the sample rows (true randomization)
    samples_randomized <- samples[sample(nrow(samples)), ]
    
    # Initialize an empty data frame for the samples with interspersed QCs
    samples_with_qcs <- data.frame(Position = character(0),
                                   name = character(0),
                                   Sample.type = character(0),
                                   stringsAsFactors = FALSE)
    
    # Start numbering new QC entries from num_QCs + 1
    qc_counter <- num_QCs + 1
    
    # Add the samples and intersperse new QCs
    for (i in 1:nrow(samples_randomized)) {
      # Add a sample row
      samples_with_qcs <- rbind(samples_with_qcs, samples_randomized[i, ])
      
      # After every insert_after_samples, insert a new QC row
      if (i %% insert_after_samples == 0) {
        new_qc_row <- data.frame(Position = qcs$Position,
                                 name = paste(new_QC_samples_names, sprintf("%02d", qc_counter)),
                                 Sample.type = "QC",
                                 stringsAsFactors = FALSE)
        samples_with_qcs <- rbind(samples_with_qcs, new_qc_row)
        qc_counter <- qc_counter + 1
      }
    }
    
    Data_out <- rbind(blanks, qcs_combined, samples_with_qcs, new_qc_row)
    
    return(Data_out)  # Ensure you have this line to return the processed data
  })
  
  output$table <- renderDT({
    req(data_processed())  # Ensure data has been processed
    data_processed()
  }, options = list(pageLength = 50))
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("data-output-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(data_processed())  # Ensure data has been processed
      write.csv(data_processed(), file, row.names = FALSE)
    }
  )
}

ui <- fluidPage(
  titlePanel("Hygge Project seq generator"),
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Choose CSV File",
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
      numericInput("num_eqQCs", "Number of eqQCs:", 1),
      numericInput("num_MSMSs", "Number of MSMSs:", 6),
      numericInput("num_QCs", "Number of QCs before samples:", 4),
      numericInput("insert_after_samples", "Insert QC After Every N Samples:", 5),
      actionButton("process", "Process Data"),
      downloadButton("downloadData", "Download CSV")
    ),
    mainPanel(
      DTOutput("uploadedTable"),
      DTOutput("table")
    )
  )
)

shinyApp(ui = ui, server = server)
