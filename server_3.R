# Server for lipid_heatmap #3



# Lipid Heatmap
# Values is used for HTML message display before and after data process
values <- reactiveValues(runProcessClicked = FALSE)

# When bottum clicked in interface, all the following will be processed
observeEvent(input$run_process, {
  values$runProcessClicked <- TRUE
  
  # Accessing sequence and data from active files
  sequence <- rv$sequence[[rv$activeFile]]
  data <- rv$data[[rv$activeFile]]
  
  # Removes anything that are not part of the data of the samples and name. 
  data <- data[, sequence[ , 'labels'] %in% c("Name","Sample")]
  sequence <- sequence[sequence[ , 'labels'] %in% c("Name","Sample"), ]
  
  # Check if the sequence is uploaded before proceeding
  if (is.null(sequence)) {
    return(NULL)  # Stop the observeEvent if no file is uploaded
  }
  
  # Capture the number of rows before filtering
  number_of_rows_before <- nrow(data)
  
  
  # Cleaning the noise of the data, by calling functions
  # Apply the `extract_pattern` function to the first column of the dataset
  # Removes all noise from compound name, so name and length is the only left: eg. going from "CAR 14:1'CAR'[M+H]+" to "CAR 14:1"
  data[, 1] <- sapply(data[, 1], extract_pattern)
  
  # Puts the length and double bonds numbers into a (). Eg "CAR 14:1" to "CAR(14:1)"
  data[, 1] <- sapply(data[, 1], format_strings)
  
  # Function to filter rows based on the specified pattern, meaning removes any data that are not on X(C:D) format.
  data <- filter_data_by_pattern(data)
  
  
  
  
  # This will make it possible to switch between original data and merged data. OG data: using _1 _2 ... _n. Merged will sum the values of the "duplicated" data. 
  if (input$selected_dataset == "original") {
    # Call a function to process the original data
    data <- unique_compound_names(data)
  } else if (input$selected_dataset == "merged") {
    # Call a function to process the merged data
    data <- merge_duplicates(data)
  }
  
  # Some data are in NA, calculations cannot read this, therefore NA values are set low. 
  data[is.na(data)] <- 0.000001
  
  # For data counting, used in display of how many rows are removed. 
  # Capture the number of rows after filtering
  number_of_rows_after <- nrow(data)
  
  # Calculate the number of rows removed
  rows_removed <- number_of_rows_before - number_of_rows_after
  
  # Output the number of rows removed to the console
  print(paste("Rows removed:", rows_removed))
  
  # Alternatively, you can display this information in the UI using a textOutput
  output$rows_removed_text <- renderText({
    paste("Rows removed after data cleaning are:", rows_removed, ". The removal is be due to the names in the first column of the data file not being in the X(C:D) format. Keep in mind, that a merged data will also count as a removed row.")
  })
  
  # Notification text for the UI
  output$limitation_notice <- renderText({
    data <- rv$data[[rv$activeFile]]  # or however you access your raw data
    paste("Displaying first 25 observations and 5 variables out of", 
          nrow(data), "observations and", ncol(data), "variables in the dataset before data cleaning.")
  })
  
  
  # The following is used in the tab: 'Data of groups in heatmap'.
  # Call the process_data function and use the result directly within observeEvent
  processed_data <- process_data(sequence, data)
  
  output$groups_table <- renderTable({
    if (is.null(processed_data)) {
      return()
    }
    
    # Extract the 'Sample' labels and corresponding 'class' from 'sequence'
    sample_rows <- sequence[sequence$labels == "Sample", ]
    unique_groups <- unique(sample_rows$class)
    
    # Create the dataframe to be displayed as a table
    df_processed_data <- data.frame(
      Group = unique_groups,
      Samples = sapply(unique_groups, function(group) {
        sample_identifiers <- rownames(sample_rows)[sample_rows$class == group]
        paste(sample_identifiers, collapse = ", ")
      })
    )
    # Return the data frame to be rendered as a table
    df_processed_data
  })
  
  # Add this to render the raw data table, being displayed in "Data of groups in Heatmap"
  output$raw_data_table <- renderTable({
    
    # Notify the user of the display limitation and the total size of the data
    total_obs <- nrow(data)
    total_vars <- ncol(data)
    
    # limteing the displayed data
    limited_data <- head(data, 25)[, 1:5]
    
    return(limited_data)
  })
  
  
  # Table output of the table in the tab: 'Heatmap' used for testing. 
  observeEvent(input$run_process, {
    # Process your data here
    processed_results <- process_data(sequence, data)
    grouped_data_frames <- create_grouped_data_frames(sequence, data)
    
    # Add the first column of "data" to each grouped data frame
    compound_names <- data[[1]]  # Extract the first column which contains compound names
    
    # Assuming that each grouped data frame has rows in the same order as "data"
    for (i in seq_along(grouped_data_frames)) {
      grouped_data_frames[[i]] <- cbind(Compound_Name = compound_names, grouped_data_frames[[i]])
    }
    
    # Update the names of the grouped_data_frames if they're not already set
    names(grouped_data_frames) <- paste("Group", seq_along(grouped_data_frames))
    
    # Dynamically generate selectInput for group selection
    output$select_group_ui <- renderUI({
      selectInput("selected_group", "Select Group:",
                  choices = names(grouped_data_frames))  # Use group names as choices
    })
    
    # Dynamically generate table output for the selected group
    output$selected_group_table <- renderTable({
      req(input$selected_group)  # Ensure a group is selected
      grouped_data_frames <- grouped_data_frames[[input$selected_group]]
      if (is.null(grouped_data_frames)) {
        return(data.frame())  # Return an empty data frame if group data is not available
      }
      grouped_data_frames
    })
  })
  
  
  # Heatmap input selection  
  observeEvent(input$run_process, {
    
    # Process your data here
    processed_results <- process_data(sequence, data)
    grouped_data_frames <- create_grouped_data_frames(sequence, data)
    grouped_data_frames_with_means <- calculate_means_for_grouped_data(grouped_data_frames)
    
    # Add the first column of "data" to each grouped data frame
    compound_names <- data[[1]]  # Extract the first column which contains compound names
    
    # Assuming that each grouped data frame has rows in the same order as "data"
    for (i in seq_along(grouped_data_frames_with_means)) {
      grouped_data_frames_with_means[[i]] <- cbind(Compound_Name = compound_names, grouped_data_frames_with_means[[i]])
    }
    
    # Update the names of the grouped_data_frames if they're not already set
    names(grouped_data_frames_with_means) <- paste("Group", seq_along(grouped_data_frames_with_means))
    
    # Create UI for group selection
    output$select_group_ui_heatmap <- renderUI({
      tagList(
        selectInput("selected_group_for_numerator", "Select Group for numerator:",
                    choices = names(grouped_data_frames_with_means)),
        selectInput("selected_group_for_denominator", "Select Group for denominator:",
                    choices = names(grouped_data_frames_with_means))  # Use group names as choices
      )
    })
    
    # Message shown when hovering over Original data and merged data. # Remember to change this outside of the observe event, Search for addTooltip
    observe({
      addTooltip(session, "selected_dataset", 
                 "Choose 'Original Data' to work with the data as it was initially collected. Select 'Merged Data' for a combined and cleaned dataset.", 
                 placement = "bottem", 
                 trigger = "hover")
    })
    
    # Render UI for maximum p-value input
    output$p_value_max_ui <- renderUI({
      numericInput("p_value_max", 
                   "Maximum p-value:", 
                   value = 0.05, 
                   min = 0, 
                   step = 0.01)
    })
    
    # Render UI for logFC input
    output$logFC_input_ui <- renderUI({
      numericInput("logFC_input", 
                   "Enter logFC value:", 
                   value = 2)
    })
    
    # Dynamic p-values depended on interface
    reactiveFilteredData <- reactive({
      # Get the maximum p-value threshold from the input
      p_value_max <- input$p_value_max
      
      # Filter the data based on the maximum p-value
      filtered_data <- numerator_data %>%
        filter(P_Value <= p_value_max)
      
      # Now return the filtered data
      filtered_data
    })
    
    
    
    # Dynamic logFC depended on interface
    reactive_max_logFC <- reactive({
      input$logFC_input  # This will be the positive value
    })
    
    reactive_min_logFC <- reactive({
      -input$logFC_input  # This will be the negative value
    })
    
    
    
    
    # logFC calculation
    reactiveLogFC <- reactive({
      # The required data input for the data handling. 
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      req(reactive_max_logFC(), reactive_min_logFC())
      
      # Define data input, makes it more readable 
      numerator_data <- grouped_data_frames_with_means[[input$selected_group_for_numerator]]
      denominator_data <- grouped_data_frames_with_means[[input$selected_group_for_denominator]]
      
      
      # Ensure there is data to work with
      if(nrow(numerator_data) == 0 || nrow(denominator_data) == 0) {
        return(NULL)
      }
      
      # Calculate logFC
      logFC <- log2((numerator_data$Mean + 1e-6) / (denominator_data$Mean + 1e-6))
      numerator_data$logFC <- logFC
      
      
      # Filter based on the selected lipid(s), if not 'All'
      if(!"All" %in% input$selected_lipid) {
        numerator_data <- numerator_data[lipid_names$Class %in% input$selected_lipid, ]
      }
      
      # Filter based on the input logFC range
      filtered_data <- numerator_data[numerator_data$logFC >= reactive_min_logFC() & numerator_data$logFC <= reactive_max_logFC(), ]
      return(filtered_data)
    })
    
    
    
    
    # Calculation of p-value, using t-test
    reactiveP_value <- reactive({
      req(input$selected_group_for_numerator, input$selected_group_for_denominator)
      
      numerator_data <- grouped_data_frames_with_means[[input$selected_group_for_numerator]]
      denominator_data <- grouped_data_frames_with_means[[input$selected_group_for_denominator]]
      
      if(nrow(numerator_data) == 0 || nrow(denominator_data) == 0) {
        return(NULL)
      }
      
      # Initialize a vector to store the p-values
      p_values <- numeric(nrow(numerator_data))
      
      # Loop through each row to perform the t-test
      for (i in 1:nrow(numerator_data)) {
        t_test_result <- t.test(numerator_data[i, -1], denominator_data[i, -1])
        p_values[i] <- t_test_result$p.value
      }
      
      # store the data
      numerator_data$p_values <- p_values
      
      # Filtrate/store the data, so it is ready to display in table in 'Heatmap'.
      filtered_data <- numerator_data %>%
        mutate(p_value = p_values) %>%
        filter(p_value <= input$p_value_max)
      
      return(filtered_data[, c("Compound_Name", "p_values")])
    })
    
    
    # Combine both logFC and p-values into one reactive expression
    reactiveFilteredData <- reactive({
      # Retrieve the filtered datasets based on logFC and p-values
      logFCData <- reactiveLogFC()
      pValuesData <- reactiveP_value()
      
      # Ensure both datasets are not NULL before proceeding
      req(logFCData, pValuesData)
      
      # Combine the datasets to have both logFC and p-value information
      # Assuming both datasets have a 'Compound_Name' column to join on
      combinedData <- merge(logFCData, pValuesData, by = "Compound_Name")
      
      # Now return the combined dataset
      combinedData
    })
    
    
    # Interface of selectinos of lipids to display
    output$select_lipid_ui <- renderUI({
      # Extract the lipid names from first column of the file 'data'
      lipid_names <<- group_lipids_by_class(data)
      
      ###### Above this in the definition of the lipid_names I instead tried to call the function. 
      
      selectizeInput("selected_lipid", "Select lipid(s) to display:",
                     choices = c("All", unique(lipid_names$Class)),
                     multiple = TRUE,
                     options = list(placeholder = 'Choose lipids...',
                                    onInitialize = I('function() { this.setValue(""); }')))
    })
    
    # Reactive expression to track the number of selected lipids or the "All" selection
    selected_lipid_count <- reactive({
      # If "All" is selected, we could set this to a value that causes the default text size to be used
      if ("All" %in% input$selected_lipid) {
        return(Inf)  # 'Inf' is used here as a flag for "All"
      } else {
        return(length(input$selected_lipid))
      }
    })
    
    
    
    
    
    # Heatmap display
    # Use the reactive expression in renderPlotly
    output$heatmapPlot <- renderPlot({
      filtered_data <- reactiveLogFC()
      
      # Take the input from user in interface and change p-value and logFC
      filtered_data <- reactiveFilteredData()
      
      
      # Ensure the data is not NULL and has rows to plot
      req(nrow(filtered_data) > 0)
      
      # Apply any necessary filtering based on logFC
      filtered_data <- filtered_data[filtered_data$logFC >= -2 & filtered_data$logFC <= 2, ]
      
      # Get the count of selected lipids
      num_of_lipids <- selected_lipid_count()
      
      
      # Adjust text size based on the count of selected lipids
      lipid_text_size <- if (num_of_lipids < 5) {
        10  # Smaller number of lipids, larger text
      } else if (num_of_lipids == Inf) {
        4  # "All" is selected, use default text size
      } else {
        5   # More than 5 lipids, smaller text
      }        
      # Ensure compound_names are available. If compound_names were defined earlier,
      # make sure they are accessible here, either as reactive values or as global variables.
      names.mapping <- map_lipid_names(x = filtered_data$Compound_Name)
      
      heatmap_plot <- heatmap_lipidome(
        x = filtered_data[ , c("Compound_Name", "logFC")],
        names.mapping = names.mapping,
        class.facet = "wrap",
        x.names = "Compound_Name",
        fill.limits = c(-2.5, 2.5), # Set limits to include -2.5 to 2.5
        fill.midpoint = 2.5, # Set midpoint explicitly to 2.5
        melt.value.name = "logFC",
        scales = "free"
      ) +
        scale_fill_gradient2(
          low = "blue",       # Color for low values
          mid = "white",     # Color for midpoint values
          high = "red",     # Color for high values
          midpoint = 0,      # Midpoint value, adjust as needed (e.g., the neutral point in your data)
          limit = c(-2.5, 2.5),  # Limits for the scale
          space = "Lab",     # Color space in which to calculate gradient
          name = "logFC"     # Legend title
        ) +
        scale_x_continuous(
          breaks = scales::pretty_breaks(n = lipid_text_size)  # Adjust n to control label frequency
        ) +
        scale_y_continuous(
          breaks = scales::pretty_breaks(n = lipid_text_size)  # Adjust n as needed for the y-axis
        ) 
      # Return the heatmap plot
      heatmap_plot
    })
    
    
    
    #Display af p_Values and logFC values under Heatmap
    output$pValueTable <- renderDataTable({
      # Access the logFC and p_values data
      logFCData <- reactiveLogFC()
      pValuesData <- reactiveP_value()
      
      
      # Ensure both are not NULL before attempting to merge
      req(logFCData, pValuesData)
      
      # Merge the dataframes based on the common "Compound_Name" column
      combinedData <- merge(logFCData, pValuesData, by = "Compound_Name")
      
      # Select only the columns you want to display
      dataTableToShow <<- combinedData[, c("Compound_Name", "logFC", "p_values")]
      
      
      # Round 'logFC' and 'p_values' to the desired number of decimal places
      dataTableToShow$logFC <- round(dataTableToShow$logFC, 5) # 2 decimal places for logFC
      dataTableToShow$p_values <- round(dataTableToShow$p_values, 5) # 4 decimal places for p-values
      
      # Render the selected data in a DataTable
      datatable(dataTableToShow, options = list(pageLength = 10, scrollX = TRUE))
    })
  })
  
  
  # Update the UI message 
  output$table_message <- renderUI({
    if (values$runProcessClicked) {
      HTML('<p>Data processing complete, see tables below.</p>')
    }
  })
  
}) # This finishes the first 'observeEvent' when 'Run data processing' is clicked

# Outside of the observeEvent, based on whether runProcessClicked is TRUE or FALSE, the message display will be placed on this: 
output$table_message <- renderUI({
  if (!values$runProcessClicked) {
    HTML('<p>Make sure sequences file is uploaded, when uploaded: Press "Run Data Processing" to get a display of data</p>')
  }
})

# Outside of the observeEvent, so the message both are shown before and after runProcessClicked is clicked. 
observe({
  addTooltip(session, "selected_dataset", 
             "Choose 'Original Data' to work with the data as it was initially collected. Select 'Merged Data' for a combined and cleaned dataset.", 
             placement = "bottom", 
             trigger = "hover")
})

# User guide inside 'Heatmap'
observeEvent(input$show_help, {
  showModal(modalDialog(
    title = "User Guide for the Lipid Heatmap",
    tags$ul(tags$li(tags$b("This Heatmap is designed for comparative lipidomic analysis."))),
    tags$ul("The user guide will explain how to use the Heatmap visualization.",
            tags$li("Upon clicking 'Run Data Processing', the app cleans the data and prepares it for plotting. This process may take a few seconds due to extensive data handling. The 'Data of groups in Heatmap' tab allows you to view the distribution of samples across groups."),
            tags$li("Original Data vs. Merged Data: Switch between 'Original Data' and 'Merged Data' by selecting the desired option and then clicking 'Run Data Processing'. This action will automatically reset the heatmap plot and any selected scales to their default settings."),
            tags$li(tags$b("When using the Heatmap:")),
            tags$li("Select groups for comparison from the dropdown menus."),
            tags$li("The heatmap will automatically update to display the data for the chosen groups."),
            tags$li("To display specific lipids, use 'Select lipid(s) to display:' to type or click your selections. To remove selections, use 'Backspace' on your keyboard. If there are no display of the selected lipids, if may be due to the thresholds of logFC and P-value"),
            tags$li("Adjustments to logFC and P-values are available. Click on 'Enter max logFC value' or 'Maximum P-value:' to enter new thresholds. The application accepts both comma and period as decimal separators. Lipids will be displayed within the specified logFC and P-value ranges. For instance, entering a logFC of 2 will automatically consider a range from -2 to 2."),
            tags$li("For further details on lipid data, scroll down to the table beneath the heatmap."),
            tags$li(tags$b("Data Calculations:")),
            tags$li("The color scale of the heatmap represents log-fold change (logFC) values. LogFC is computed as log2((numerator_data$Mean + 1e-6) / (denominator_data$Mean + 1e-6)), where numerator_data and denominator_data correspond to the groups selected. An offset of 1e-6 is included to avoid division by zero. The logFC scale is set to span from -2.5 to 2.5."),
            tags$li("P-values are determined using the 't.test' function in R, which conducts a statistical comparison between corresponding rows of lipid data from the selected groups."),
            
            
    ),
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})  

