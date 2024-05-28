# Lipid Heatmap 1, first standalone script that works. 
# This code is my fist probably working script, that gives the intented output. 
# Run all lines before any display. 


#Locate data storage 
getwd()

#Set data location
setwd("/Users/jakobhartvig/Desktop/GitHub/B.Sc._programs") 

#Load data
lab_dataset <-read.csv("lipidomics_data.csv")

#Check the data
#View(lab_dataset)


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

#View(lab_dataset)


# Use grepl and regular expressions to filter rows
#removes any data that are not on X(C:D) format
pattern <- "^.+\\(\\d+:\\d+\\)$"  # Regular expression pattern
lab_dataset <- lab_dataset[grepl(pattern, lab_dataset$Compound.Name), ]


#### duplicated lipids #### 

library(dplyr)

# Merge of duplicates. Going from 229 obs. to 180 obs.
# Identify duplicates in the "Compound.Name" column 
duplicated_compounds <- duplicated(lab_dataset$Compound.Name)
#Removes duplicates by summing them
lab_dataset <- lab_dataset %>%
  group_by(Compound.Name) %>%
  summarise_all(sum)

# Used for future programming. 
# If testing this, do not run the lines between this line and "#### duplicated lipids ####"
# After running the line below, the data frame should be fixed to _n
# add _1 _2 _3 to duplicated names
lab_dataset$Compound.Name <- ave(lab_dataset$Compound.Name, lab_dataset$Compound.Name, FUN = function(x) if (length(x) > 1) paste0(sub("\\(.*\\)", "", x), "_", seq_along(x), sub(".*\\(", "(", x)) else x)

print(duplicated_compounds)

#mean of the data
lab_dataset$Mean <- rowMeans(lab_dataset[, -1, -2], na.rm = TRUE)  # Use -1 to omit the first column of names


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

#log fold change 
#average of a group (class 4) divided by the average of another group (class 6)
lab_dataset$logFC_4_6 <- (lab_dataset$Mean_class_4/lab_dataset$Mean_class_6) 



#Frames the data as desired
selected_columns <- lab_dataset[, c("Compound.Name", "Mean", "Mean_class_4", "Mean_class_6", "logFC_4_6")]

View(selected_columns)

#code from colleagues used at our data
library( "lipidomeR" )

# Enumerate the lipid names into values.
names.mapping <- map_lipid_names( x = selected_columns$Compound.Name )

# Create the lipidomeR heatmap of lipid concentrations.
heatmap_lipidome(
  x = selected_columns[ , c( "Compound.Name", "logFC_4_6" ) ],
  names.mapping = names.mapping,
  class.facet = "wrap",
  x.names = "Compound.Name",
  fill.limits =
    range(
      x = selected_columns$"logFC_4_6",
      na.rm = TRUE
    ),
  fill.midpoint =
    sum(
      range(
        x = selected_columns$"logFC_4_6",
        na.rm = TRUE
      )
    ) / 2,
  melt.value.name = "logFC_4_6",
  scales = "free"
)




