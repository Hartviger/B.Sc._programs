# Lipid heatmap 2, Working Shiny app. 
# Here the Interactive Heatmap is possible to run, but beware that some data is lost doing the process. 

#Locate data storage #fixing 
getwd()

#Set data location
setwd("/Users/jakobhartvig/Desktop/GitHub/B.Sc.-project") 

#Load data
lab_dataset <-read.csv("lab_data.csv")


#Check the data
#View(lab_dataset)
#View(seq)ª

#string manipulating 
#Removes all other than name and length: going from CAR 14:1'CAR'[M+H]+ to CAR 14:1
# Define a function to extract the desired pattern from a name 
extract_pattern <- function(name) {
  # Pattern to find first part consisting of letters and numbers with a colon or a letter before the numbers
  pattern <- "([A-Za-z]+\\s[0-9]+:[0-9]+)|([A-Za-z]+\\s[[:alpha:]]?-?[0-9]+:[0-9]+)"
  matches <- regmatches(name, gregexpr(pattern, name))
  
  # Returns the first match, or the hole name if no match
  if (length(matches[[1]]) > 0) {
    return(matches[[1]][1])
  } else {
    return(name)
  }
}

# Update the "Names" column using the extract_pattern function 
lab_dataset$Compound.Name <- sapply(lab_dataset$Compound.Name, extract_pattern)

#Check the data, make sure the function work
#View(lab_dataset)


#puts the length and double bonds numbers into a (), eg. CAR 14:1 to CAR(14:1)
format_strings <- function(input_strings) {
  # Use gsub with regular expression to remove all whitespace characters
  formatted_strings <- gsub("\\s+", "", input_strings)
  # Add parentheses around the numbers
  formatted_strings <- gsub("([A-Za-z]*)(\\d+):(\\d+)", "\\1(\\2:\\3)", formatted_strings)
  return(formatted_strings)
}

#Updates the table with function above. 
lab_dataset$Compound.Name <- format_strings(lab_dataset$Compound.Name)


#Check the data, make sure the function work
#View(lab_dataset)


# Use grepl and regular expressions to filter rows
#removes any data that are not on X(C:D) format
pattern <- "^.+\\(\\d+:\\d+\\)$"  # Regular expression pattern
lab_dataset <- lab_dataset[grepl(pattern, lab_dataset$Compound.Name), ]


#Merged data
#If any data is duplicated, this function will merge them
#Merge of duplicates. Example: going from 229 obs. to 180 obs.
# Identify duplicates in the "Compound.Name" column 
duplicated_compounds <- duplicated(lab_dataset$Compound.Name)

#Removes duplicates by summing them
library(dplyr)
lab_dataset_merged <- lab_dataset %>%
  group_by(Compound.Name) %>%
  summarise_all(sum)


#Original data structure 
#add _1 _2 _3 to duplicated names depending on how many duplicated names there are. 
lab_dataset$Compound.Name <- ave(lab_dataset$Compound.Name, lab_dataset$Compound.Name, FUN = function(x) if (length(x) > 1) paste0(sub("\\(.*\\)", "", x), "_", seq_along(x), sub(".*\\(", "(", x)) else x)


#Calculation of data

#mean of the hole datasheet
#Add a new column in Lab_dataset.
lab_dataset$Mean <- rowMeans(lab_dataset[, -1, -2], na.rm = TRUE) 

# Calculate the mean for class 2
class_2_columns <- grep("X\\d+_B_2", names(lab_dataset), value = TRUE)
lab_dataset$Mean_class_2 <- rowMeans(lab_dataset[class_2_columns], na.rm = TRUE)

# Calculate the mean for class 3
class_3_columns <- grep("X\\d+_B_3", names(lab_dataset), value = TRUE)
lab_dataset$Mean_class_3 <- rowMeans(lab_dataset[class_3_columns], na.rm = TRUE)

# Calculate the mean for class 4
class_4_columns <- grep("X\\d+_B_4", names(lab_dataset), value = TRUE)
lab_dataset$Mean_class_4 <- rowMeans(lab_dataset[class_4_columns], na.rm = TRUE)

# Calculate the mean for class 5
class_5_columns <- grep("X\\d+_C_1", names(lab_dataset), value = TRUE)
lab_dataset$Mean_class_5 <- rowMeans(lab_dataset[class_5_columns], na.rm = TRUE)

# Calculate the mean for class 6
class_6_columns <- grep("X\\d+_C_2", names(lab_dataset), value = TRUE)
lab_dataset$Mean_class_6 <- rowMeans(lab_dataset[class_6_columns], na.rm = TRUE)

# Calculate the mean for class 7
class_7_columns <- grep("X\\d+_C_3", names(lab_dataset), value = TRUE)
lab_dataset$Mean_class_7 <- rowMeans(lab_dataset[class_7_columns], na.rm = TRUE)


#original data
#Frame the data in the different class samples
#Data needed for p-value calculation # Utillzes the function from: Calculate the mean for class ... 
class_2 <- lab_dataset[, c("Compound.Name", class_2_columns)]
class_3 <- lab_dataset[, c("Compound.Name", class_3_columns)]
class_4 <- lab_dataset[, c("Compound.Name", class_4_columns)]
class_5 <- lab_dataset[, c("Compound.Name", class_5_columns)]
class_6 <- lab_dataset[, c("Compound.Name", class_6_columns)]
class_7 <- lab_dataset[, c("Compound.Name", class_7_columns)]

#merged data
#Frame the data in the different class samples 
#Data needed for p-value calculation # Utillzes the function from: Calculate the mean for class ... 
class_2_merged <- lab_dataset_merged[, c("Compound.Name", class_2_columns)]
class_3_merged <- lab_dataset_merged[, c("Compound.Name", class_3_columns)]
class_4_merged <- lab_dataset_merged[, c("Compound.Name", class_4_columns)]
class_5_merged <- lab_dataset_merged[, c("Compound.Name", class_5_columns)]
class_6_merged <- lab_dataset_merged[, c("Compound.Name", class_6_columns)]
class_7_merged <- lab_dataset_merged[, c("Compound.Name", class_7_columns)]



#P_value calculation local on R
# Assuming class_5 and class_6 are data frames that contain the same lipids 
# in the same order and that the first column is "Compound.Name"

# Initialize a vector to store the p-values
p_values <- numeric(nrow(class_5))

# Loop through each lipid
for (i in 1:nrow(class_5)) {
  # Perform the t-test on the ith lipid
  t_test_result <- t.test(class_5[i, -1], class_6[i, -1])  # -1 to exclude the "Compound.Name" column
  # Store the p-value
  p_values[i] <- t_test_result$p.value
}

# Add the p-values to the class_2 data frame for demonstration purposes
class_5$p_values <- p_values

print(p_values)

names.mapping <- map_lipid_names(x = lab_dataset$Compound.Name)


#shiny app

library(shiny)
library(dplyr)
library(plotly)

# Define UI for application
ui <- fluidPage(
  titlePanel("Interactive Heatmap Visualization"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("dataset", "Select data frame:",
                   choices = c("Original Data" = "lab_dataset", "Merged Data" = "lab_dataset_merged"),
                   selected = "lab_dataset"),
      selectInput("group1", "Select group for numerator:",
                  choices = c("Class 2" = "class_2", 
                              "Class 3" = "class_3",
                              "Class 4" = "class_4",
                              "Class 5" = "class_5",
                              "Class 6" = "class_6",
                              "Class 7" = "class_7")),
      
      selectInput("group2", "Select group for denominator:",
                  choices = c("Class 2" = "class_2", 
                              "Class 3" = "class_3",
                              "Class 4" = "class_4",
                              "Class 5" = "class_5",
                              "Class 6" = "class_6",
                              "Class 7" = "class_7")),
      actionButton("show_help", "Show User Guide")
    ),
    mainPanel(
      plotlyOutput("heatmapPlot", width = "100%", height = "650px"),
      dataTableOutput("pValueTable")  # Use dataTableOutput for DT::renderDataTable #shows the table
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  current_data <- reactiveVal(lab_dataset)
  
  observeEvent(input$dataset, {
    if (input$dataset == "lab_dataset") {
      current_data(lab_dataset)
    } else {
      current_data(lab_dataset_merged)
    }
  })
  
  # Define a function to calculate mean class data
  calculate_class_data <- function(lab_dataset, class_columns) {
    class_data <- lab_dataset[, c("Compound.Name", class_columns)]
    class_data$Mean <- rowMeans(class_data[, -1], na.rm = TRUE)
    return(class_data)
  }
  
  # Reactive expression to store class data frames
  reactive_class_data <- reactive({
    list(
      class_2 = calculate_class_data(current_data(), class_2_columns),
      class_3 = calculate_class_data(current_data(), class_3_columns),
      class_4 = calculate_class_data(current_data(), class_4_columns),
      class_5 = calculate_class_data(current_data(), class_5_columns),
      class_6 = calculate_class_data(current_data(), class_6_columns),
      class_7 = calculate_class_data(current_data(), class_7_columns)
    )
  })
  
  # Reactive expression for logFC
  reactive_logFC <- reactive({
    df_group1 <- reactive_class_data()[[input$group1]]
    df_group2 <- reactive_class_data()[[input$group2]]
    
    # Ensure that the data frames have the same number of rows
    min_rows <- min(nrow(df_group1), nrow(df_group2))
    
    df_group1$Mean <- rowMeans(df_group1[, -1], na.rm = TRUE)
    df_group2$Mean <- rowMeans(df_group2[, -1], na.rm = TRUE)
    
    df_group1$logFC <- log2((df_group1$Mean + 1e-6) / (df_group2$Mean + 1e-6))
    
    return(df_group1$logFC[1:min_rows])
  })
  
  # Reactive expression for p-values
  reactive_p_values <- reactive({
    group1_data <- reactive_class_data()[[input$group1]]
    group2_data <- reactive_class_data()[[input$group2]]
    
    # Initialize a vector to store the p-values
    p_values <- numeric(nrow(group1_data))
    
    # Loop through each lipid
    for (i in 1:nrow(group1_data)) {
      # Perform the t-test on the ith lipid
      t_test_result <- t.test(group1_data[i, -1], group2_data[i, -1])  # -1 to exclude "Compound.Name"
      # Store the p-value
      p_values[i] <- t_test_result$p.value
    }
    
    # Add p-values to the data frame
    df <- current_data()
    df$p_values <- p_values
    
    # Return the p-values
    return(df$p_values)
  })
  
  # Render heatmap with the updated logFC values
  output$heatmapPlot <- renderPlotly({
    df <- current_data()
    df$logFC <- reactive_logFC()
    df$p_values <- reactive_p_values()  
    df <- df[df$logFC >= -2 & df$logFC <= 2, ]
    names.mapping <- map_lipid_names(x = df$Compound.Name)
    p <- heatmap_lipidome(
      x = df[, c("Compound.Name", "logFC")],
      names.mapping = names.mapping,
      class.facet = "wrap",
      x.names = "Compound.Name",
      fill.limits = c(-2.5, 2.5),
      fill.midpoint = 0,
      melt.value.name = "logFC",
      scales = "free"
    )
  })
  
  # Render the table of p-values and logFC in the UI
  output$pValueTable <- renderDataTable({
    df <- current_data()
    df$logFC <- reactive_logFC()
    df$p_values <- reactive_p_values()
    
    # Ensure that the data frames have the same number of rows
    min_rows <- min(nrow(df), length(df$logFC))
    
    df <- df[1:min_rows, c("Compound.Name", "logFC", "p_values")]
    df
  }, options = list(pageLength = 10, scrollX = TRUE))  # Options for DataTable
}

# Observe event for help button click
observeEvent(input$show_help, {
  showModal(modalDialog(
    title = "User Guide for Heatmap",
    tags$ul("This app is designed for comparative lipidomic analysis."),
    tags$ul("This user guide will explain how to use the interactive heatmap visualization and interpret the results.",
            tags$li("Select groups for comparison using the dropdown menus."),
            tags$li("The heatmap will automatically update to reflect the selected groups."),
            tags$li("Hover over the heatmap to see detailed information about each lipid."),
            tags$li("The color scale represents the log-fold change (logFC) values."),
            tags$li("A logFC value close to zero suggests no significant change."),
            tags$li("Positive logFC values indicate higher concentrations in the numerator group, meaning it is upregulated compared to the class selected in the denominator."),
            tags$li("Negative logFC values suggest higher concentrations in the denominator group, indicating downregulation in the numerator group compared to the denominator."),
            tags$li("The p-value table below the heatmap provides insight into the statistical significance of the observed differences in lipid concentrations between the selected groups. A p-value is a measure of the probability that the observed data would occur if there were no real effect or difference (null hypothesis). In other words, it indicates how likely it is that any observed difference in lipid concentrations happened by random chance. A low p-value (typically less than 0.05) suggests that the difference is unlikely to be due to chance, and therefore, may be considered statistically significant. Conversely, a high p-value suggests that the differences observed could easily happen by random chance and are not necessarily due to any specific cause or effect. It is important to consider the p-value in the context of the study design, sample size, and other relevant factors when interpreting the results.")
    ),
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})

# Run the application
shinyApp(ui = ui, server = server)
