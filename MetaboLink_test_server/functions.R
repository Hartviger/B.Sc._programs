checkSequence <- function(sequence) {
  columnsToCheck <- c("sample", "batch", "order", "class", "time", "paired", "amount")
  missingColumns <- setdiff(columnsToCheck, colnames(sequence))
  if (length(missingColumns) > 0) {
    sequence[missingColumns] <- lapply(seq_along(missingColumns), function(x) sequence[missingColumns[x]] <- NA)
  }
  sequence <- sequence[, c("sample", "batch", "order", "class", "time", "paired", "amount")]
  return(sequence)
}

checkColumns <- function(columns, labels) {
  if(any(labels == "-")) {
    sendSweetAlert(title = "Info", text = paste("Columns ", paste(columns[labels == "-"], collapse = ", "), " labeled '-'.\n
        If this is wrong, check file for invalid characters in these columns."), type = "info")
  }
}

blankFiltration <- function(data, sequence, signalStrength, keepIs) {
  data[sequence[, 1] %in% "Blank"][is.na(data[sequence[, 1] %in% "Blank"])] <- 0
  bf <- apply(data[sequence[, 1] %in% "Blank"], 1, mean) * signalStrength < 
                    apply(data[sequence[, 1] %in% "QC"], 1, mean, na.rm = TRUE)
  if (keepIs) { #TODO
    is <- grepl("\\(IS\\)", toupper(data[sequence[, 1] %in% "Name"][, 1]))
    data <- data[bf | is, ]
  } else {
    data <- data[bf, ]
  }
  return(data)
}

normalizationIS <- function(data, sequence, is, method, qc) {
  rt <- which(sequence[, 1] == "RT")
  isname <- is
  is <- as.numeric(gsub(" .*$", "", is))
  sel <- if (qc) c("Sample", "QC") else "Sample"
  sdat <- data[sequence[, 1] %in% sel]
  sdat[sdat == 0] <- NA
  is <- is[complete.cases(sdat[is, ])]
  near <- sapply(data[, rt], function(y) {
    which.min(abs(data[is, rt] - y))
  })
  if (method == "Same lipid structure") {
    name <- data[sequence[, 1] %in% "Name"]
    istype <- gsub(" .*$", "", name[is, ])
    near <- sapply(seq(name[, 1]), function(x) {
      if (gsub(" .*$", "", name[x, 1]) %in% istype) {
        which(istype %in% gsub(" .*$", "", name[x, 1]))
      } else {
        near[x]
      }
    })
  }
  sdat <- sapply(seq(ncol(sdat)), function(j) {
    sapply(seq(nrow(sdat)), function(i) {
      sdat[i, j] <- sdat[i, j] / sdat[is, j][near[i]]
    })
  })
  isnorm <- sapply(seq(nrow(sdat)), function(x) {
    isname[near[x]]
  })
  data[sequence[, 1] %in% sel] <- sdat
  data <- cbind(data, data.frame(isnorm = isnorm))
  return(data)
}

optimizeIS <- function(data, sequence, is, method, qc) {
  iscomb <- Map(combn, list(is), seq_along(is), simplify = FALSE)
  iscomb <- lapply(rapply(iscomb, enquote, how = "unlist"), eval)
  progressSweetAlert(id = "pbis", title = "Finding best IS combination", value = 0, total = length(iscomb), striped = T, display_pct = T)
  islow <- lapply(iscomb, function(x) {
    isdat <- normalizationIS(data, sequence, x, method, qc)
    mean(apply(isdat[, sequence[, 1] %in% "QC"], 1, sd, na.rm = T) / apply(isdat[, sequence[, 1] %in% "QC"], 1, mean, na.rm = T) * 100)
    # updateProgressBar(id = "pbis", value = which(iscomb %in% iscomb[x]), total = length(iscomb))
  })
  closeSweetAlert()
  return(unlist(iscomb[which.min(unlist(islow))]))
}

findInternalStandards <- function(data) {
  isIndex <- grepl("\\(IS\\)", toupper(data[, 1]))
  if (sum(isIndex) > 0) {
    featureName <- as.vector(data[isIndex, 1])
    internalStandards <- paste(which(isIndex), " - ", featureName)
    return(internalStandards)
  } else {
    return(character(0))
  }
}

identifyLabels <- function(data) {
  labels <- sapply(names(data), function(x) {
    if (grepl("BLANK", toupper(x), fixed = TRUE) && is.numeric(data[, x])) {
      "Blank"
    } else if (grepl("QC", toupper(x), fixed = TRUE) && is.numeric(data[, x])) {
      "QC"
    } else if (grepl("NAME", toupper(x), fixed = TRUE)) {
      "Name"
    } else if (grepl("MASS|M/Z|M.Z", toupper(x)) && is.numeric(data[, x])) {
      "Mass"
    } else if (grepl("RT|TIME|RETENTION", toupper(x)) && is.numeric(data[, x])) {
      "RT"
    } else if (grepl("ADDUCT_POS", toupper(x), fixed = TRUE)) {
      "Adduct_pos"
    } else if (grepl("ADDUCT_NEG", toupper(x), fixed = TRUE)) {
      "Adduct_neg"
    } else if (grepl("ADDUCT", toupper(x)) && grepl("\\]\\+", data[, x])) {
      "Adduct_pos"
    } else if (grepl("ADDUCT", toupper(x)) && grepl("\\]\\-", data[, x])) {
      "Adduct_neg"
    } else if (grepl("[[:digit:]]", toupper(x)) && is.numeric(data[, x])) {
      "Sample"
    }  else {
      "-"
    }
  })
  if (sum(labels == "Name") > 1) {
    labels[labels == "Name" & duplicated(labels)] <- "-"
  }
  labels <- factor(labels, levels = c("Name", "Blank", "QC", "Sample", "RT", "Mass", "Adduct_pos", "Adduct_neg", "-"))
  return(labels)
}

imputation <- function(data, seq, method, minx = 1, onlyqc, remaining) {
  if (onlyqc) {
    qcData <- data[seq[, 1] %in% "QC"]
    qcSequence <- seq[seq[, 1] %in% "QC", ]
  } else {
    qcData <- data[seq[, 1] %in% c("Sample", "QC")]
    qcSequence <- seq[seq[, 1] %in% c("Sample", "QC"), ]
  }
  qcData[qcData == 0] <- NA
  qcSequence[qcSequence[, 1] %in% c("QC"), ]$class <- "QC"
  qcSequence$class[is.na(qcSequence$class)] <- "Sample"

  if (method == "KNN") {
    qcData <- as.matrix(qcData)
    knndat <- impute.knn(qcData, k = 10, rowmax = .99, colmax = .99, maxp = 15000) # TODO skal k sættes anderledes evt mindste class -1?
    impsqdat <- as.data.frame(knndat$data)
  } else if (method == "Median") {
    impsqdat <- imp_median(qcData, qcSequence)
  } else if (method == "Min/X") {
    impsqdat <- imp_minx(qcData, qcSequence, minx)
  }

  if (sum(is.na(impsqdat)) > 0) {
    if (remaining == "Min/X") {
      for (i in 1:nrow(impsqdat)) {
        impsqdat[i, is.na(impsqdat[i, ])] <- min(impsqdat[i, ], na.rm = T) / minx
      }
    } else if (remaining == "zero") {
      for (i in 1:nrow(impsqdat)) {
        impsqdat[i, is.na(impsqdat[i, ])] <- 0
      }
    } else if (remaining == "Median") {
      for (i in 1:nrow(impsqdat)) {
        impsqdat[i, is.na(impsqdat[i, ])] <- median(as.numeric(impsqdat[i, ]), na.rm = T)
      }
    }
  }
  if (onlyqc) {
    data[seq[, 1] %in% "QC"] <- impsqdat
  } else {
    data[seq[, 1] %in% c("Sample", "QC")] <- impsqdat
  }
  return(data)
}

cutoffrm <- function(data, seq, cutoff, method) {
  cutoff <- cutoff / 100
  if ("entire data" %in% method) {
    datm <- data[seq[, 1] %in% "Sample"]
    datm[datm == 0] <- NA
    keep <- rowSums(!is.na(datm)) / ncol(datm) >= cutoff
  }
  if ("in QC" %in% method) {
    datm <- data[seq[, 1] %in% "QC"]
    datm[datm == 0] <- NA
    keep <- rowSums(!is.na(datm)) / ncol(datm) >= cutoff
  }
  if ("in group" %in% method) {
    datm <- data[seq[, 1] %in% "Sample"]
    datm[datm == 0] <- NA
    classes <- factor(seq[, 4], exclude = NA)
    nseq <- seq[seq[, 1] %in% "Sample", ]
    keep_m <- matrix(FALSE, nrow(datm), ncol = length(levels(classes)))
    for (cl in 1:length(levels(classes))) {
      cl_f <- datm[, nseq[, 4] %in% levels(classes)[cl]]
      keep_m[, cl] <- rowSums(!is.na(cl_f)) / ncol(cl_f) >= cutoff
    }
    keep <- apply(keep_m, 1, function(x) any(x))
  }
  data <- data[keep, ]
  return(data)
}

imp_median <- function(data, seq) {
  datm <- data.frame(seq$class, t(data))
  datm <- datm %>%
    group_by(seq.class) %>%
    mutate_if(
      is.numeric,
      function(x) {
        ifelse(is.na(x), median(x, na.rm = T), x)
      }
    )
  datm <- data.frame(t(datm[, -1]))
  colnames(datm) <- colnames(data)
  row.names(datm) <- row.names(data)
  return(datm)
}

imp_minx <- function(data, seq, minx) {
  datm <- data.frame(seq$class, t(data))
  datm <- datm %>%
    group_by(seq.class) %>%
    mutate_if(
      is.numeric,
      function(x) {
        ifelse(is.na(x),
          min(x, na.rm = T) / minx,
          x
        )
      }
    )
  datm <- data.frame(t(datm[, -1]))
  colnames(datm) <- colnames(data)
  row.names(datm) <- row.names(data)
  datm[datm == "Inf"] <- NA
  return(datm)
}

cvmean <- function(data) {
  mean(apply(data, 1, sd, na.rm = T) / apply(data, 1, mean, na.rm = T), na.rm = T) * 100
}

cv <- function(data) {
  round(apply(data, 1, sd, na.rm = T) / apply(data, 1, mean, na.rm = T) * 100, 2)
}

pcaplot <- function(data, class, islog) {
  data[data == 0] <- NA
  data <- data[complete.cases(data), ]
  class[!is.na(class)] <- class[!is.na(class)]
  class[is.na(class)] <- "No class"
  ifelse(islog, data <- t(data), data <- log(t(data)))
  prin <- prcomp(data, rank. = 2, center = T, scale = F)
  pov <- summary(prin)[["importance"]]["Proportion of Variance", ]
  pov <- round(pov * 100, 2)
  components <- prin[["x"]]
  components <- data.frame(components)
  label <- paste0(row.names(components), ": ", class)
  col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(class)))
  pca <- plot_ly(components, x = ~PC1, y = ~PC2, type = "scatter", mode = "markers", text = label, hoverinfo = "text", color = class, colors = col)
  pca <- pca %>% layout(
    legend = list(title = list(text = "color")),
    plot_bgcolor = "#e5ecf6",
    xaxis = list(
      title = paste0("PC1 (", pov[1], "% explained var.)"),
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor = "#ffff"
    ),
    yaxis = list(
      title = paste0("PC2 (", pov[2], "% explained var.)"),
      zerolinecolor = "#ffff",
      zerolinewidth = 2,
      gridcolor = "#ffff"
    )
  )
  return(pca)
}





# Lipid Heatmap
# Grouping the data
process_data <- function(sequence, data) {
  # Initialize an empty list to store the sample names by group and other information
  results <- list()
  
  sample_identifiers <- rownames(sequence)[sequence[, "labels"] == "Sample"]
  groups <- sequence[sample_identifiers, "class"]
  names <- sample_identifiers
  
  # Extract the 'Sample' labels and corresponding 'class' from 'sequence'
  sample_rows <- sequence[sequence$labels == "Sample", ]
  results$sample_rows <- sample_rows  # Store the sample rows in the results list
  
  unique_groups <- unique(sample_rows$class)
  results$unique_groups <- unique_groups  # Store the unique groups in the results list
  
  # Initialize an empty list within results for grouped_samples
  results$grouped_samples <- list()
  
  # Iterate over each group to get the corresponding sample names
  for (group in unique_groups) {
    if (!is.na(group)) {
      # Get the sample names for the current group
      samples_in_group <- sample_rows[sample_rows$class == group, 1]
      
      # Add the sample names to the list, named by their group
      results$grouped_samples[[paste("group", group)]] <- samples_in_group
      results[[paste("Group", group)]] <- samples_in_group  # Store each group separately in the results list
    }
  }
  
  # Check if any of the components are NULL or empty and handle accordingly
  if (length(results$grouped_samples) == 0) {
    results$grouped_samples <- list(group1 = NA)  # Placeholder if no groups found
  }
  
  
  # Return the list of results that now includes grouped samples, unique groups, sample rows, and each group
  return(results)
}



# Function to create grouped data frames based on sequence and data
create_grouped_data_frames <- function(sequence, data) {
  # Initialize an empty list to store data frames for each group
  grouped_data_frames <- list()
  
  # Extract the 'Sample' labels and corresponding 'class' from 'sequence'
  sample_rows <- sequence[sequence$labels == "Sample", ]
  unique_groups <- unique(sample_rows$class)
  
  # Iterate over each unique group to create data frames
  for (group in unique_groups) {
    if (!is.na(group)) {
      # Get the sample identifiers for the current group
      sample_identifiers <- rownames(sample_rows)[sample_rows$class == group]
      
      # Find the matching column indices, excluding NA values
      matching_indices <- match(sample_identifiers, colnames(data))
      matching_indices <- matching_indices[!is.na(matching_indices)]
      
      # Check if we have any matching columns at all
      if (length(matching_indices) > 0) {
        # Select only the columns for the current group
        group_data <- data[, matching_indices, drop = FALSE]
        
        # Store the filtered data frame in the list, named by the group
        grouped_data_frames[[paste("group", group)]] <- group_data
      } else {
        warning(paste("Group", group, "contains column names that are not in the data. Skipping this group."))
      }
    }
  }
  
  return(grouped_data_frames)
}



calculate_means_for_grouped_data <- function(grouped_data_frames) {
  # Initialize a new list to store the modified data frames
  new_grouped_data_frames <- list()
  
  # Iterate over each group's data frame in the list
  for (group_name in names(grouped_data_frames)) {
    # Clone the current group's data frame to avoid modifying the original
    group_data <- grouped_data_frames[[group_name]]
    
    # Assuming the first column is not numeric and should be excluded from the mean calculation
    # Calculate the mean for each row across all other columns
    means <- rowMeans(group_data[, drop = FALSE], na.rm = TRUE)
    
    # Append the calculated means as a new column to the cloned data frame
    group_data$Mean <- means
    
    # Add the modified data frame to the new list
    new_grouped_data_frames[[group_name]] <- group_data
  }
  
  # Return the new list of grouped data frames with means calculated
  return(new_grouped_data_frames)
}


# Function to group lipids by their class prefix (e.g., "CAR", "LP", etc.)
group_lipids_by_class <- function(data) {
  # Assuming the first column of 'data' contains the lipid names like "CAR(18:1)"
  lipid_names <- data[[1]]  # Replace 1 with the actual column name or index if different
  
  # Use a regular expression to extract the class prefix from lipid names
  # This matches any consecutive alphabetic characters at the beginning of the string
  lipid_classes <- sub("\\(([0-9]+:[0-9]+)\\).*", "", lipid_names)
  
  # Create a data frame that maps lipid names to their class
  class_mapping <- data.frame(Lipid_Name = lipid_names, Class = lipid_classes, stringsAsFactors = FALSE)
  
  # Optionally, if you want to return a list that names each group by its class
  # names(grouped_data) <- unique(lipid_classes)
  
  return(class_mapping)
}


# Data cleaning

# Function to extract patterns from compound names
# Removes all noise from compound name, so name and length is the only left: eg. going from "CAR 14:1'CAR'[M+H]+" to "CAR 14:1"
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

# Function to format strings
# Puts the length and double bonds numbers into a (). Eg "CAR 14:1" to "CAR(14:1)"
format_strings <- function(input_strings) {
  # Step 1: Remove all whitespace characters
  formatted_strings <- gsub("\\s+", "", input_strings)
  
  # Step 2: Identify and move the prefix (e.g., O-, P-) to the suffix position
  formatted_strings <- gsub("([A-Za-z]+)(O-|P-)(\\d+):(\\d+)", "\\1-\\2(\\3:\\4)", formatted_strings)
  
  # Step 3: Add parentheses around the numbers for strings without the O- or P- prefix
  formatted_strings <- gsub("([A-Za-z]+)(\\d+):(\\d+)", "\\1(\\2:\\3)", formatted_strings)
  
  return(formatted_strings)
}


### pattern_column is that corret? 
# Function to filter rows based on the specified pattern, meaning removes any data that are not on X(C:D) format.
filter_data_by_pattern <- function(data) {
  # Define the regular expression pattern
  pattern <- "^.+\\(\\d+:\\d+\\)$"
  
  # Check if the first column of data matches the pattern
  filtered_data <- data[grepl(pattern, data[[1]]), ]
  
  return(filtered_data)
}


#merge duplicated names of the data
merge_duplicates <- function(data) {
  # Ensure the first column is treated as the Compound Name
  compound_name_col <- names(data)[1]
  
  # Group by the first column and then summarise all other columns by summing
  data_merged <- data %>%
    group_by(.data[[compound_name_col]]) %>%
    summarise(across(everything(), sum, na.rm = TRUE), .groups = 'drop')
  
  return(data_merged)
}

#duplicated names have add _1, _2 and _3 depending on how many duplicates. 
unique_compound_names <- function(data) {
  # Ensure that 'data' is a data frame and has at least one column
  if (!is.data.frame(data) || ncol(data) < 1) {
    stop("The input must be a data frame with at least one column.")
  }
  
  # Apply the processing to the first column of 'data'
  data[[1]] <- ave(data[[1]], data[[1]], FUN = function(x) {
    if (length(x) > 1) {
      # Extract the base name without parentheses
      base_name <- sub("\\(.*\\)", "", x)
      # Extract the part within parentheses
      suffix <- sub(".*\\(", "(", x)
      # Combine base name with sequence number and the part within parentheses
      paste0(base_name, "_", seq_along(x), suffix)
    } else {
      x
    }
  })
  
  # Return the modified data
  return(data)
}









driftcorrection <- function(dat, seq, method, ntree = 500, degree = 2, QCspan) {
  seqsq <- seq[seq[, 1] %in% c("Sample", "QC"), ]
  datsq <- dat[, seq[, 1] %in% c("Sample", "QC")]
  datsqsorted <- datsq %>% select(order(seqsq$order))
  qcid <- as.numeric(sort(seqsq[seqsq[, 1] %in% "QC", ]$order)) # order of QCs
  frame <- data.frame("qcid" = 1:ncol(datsq)) 
  dcdat <- as.matrix(datsqsorted)
  dcdat1 <- dcdat
  progressSweetAlert(id = "pbdc", title = "Correcting drift", value = 0, total = nrow(dcdat), striped = T, display_pct = T)
  if (method == "QC-RFSC (random forrest)") {
    for (i in 1:nrow(dcdat)) {
      forest <- randomForest(data.frame(qcid), as.numeric(dcdat[i, qcid]), ntree = ntree)
      pv <- predict(forest, frame)
      dcdat[i, ] <- as.numeric(dcdat[i, ]) / pv
      updateProgressBar(id = "pbdc", value = i, total = nrow(dcdat))
    }
  }
  if (method == "QC-RLSC (robust LOESS)") {
    for (i in 1:nrow(dcdat)) {
      loess <- loess(dcdat[i, qcid] ~ qcid,
        span = QCspan,
        degree = degree
      )
      pv <- predict(loess, frame)
      dcdat[i, ] <- as.numeric(dcdat[i, ]) / pv
      updateProgressBar(id = "pbdc", value = i, total = nrow(dcdat))
    }
  }
  dat[seq[, 1] %in% c("Sample", "QC")] <- dcdat[, seqsq$order]
  closeSweetAlert()
  return(dat)
}

findAdduct <- function(data, sequence) {
  names <- data[, sequence[, 1] %in% "Name"]
  reg <- gregexpr("_\\[(.*?)_", names)
  adduct <- regmatches(names, reg)
  adduct <- as.character(adduct)
  adduct[adduct > 0] <- NA
  adduct <- gsub("^.|.$", "", adduct)
  adduct <- gsub(" ", "", adduct)
  return(adduct)
}

mergeDatasets <- function(dataset1, sequence1, dataset2, sequence2, ppmTolerance, rtTolerance) {
  progressSweetAlert(id = "pb", title = "Work in progress", display_pct = T, value = 0, striped = T)
  mass <- dataset1[, sequence1[, 1] %in% "Mass"]
  if ("Adduct_pos" %in% sequence1[, 1]) {
    adduct <- dataset1[, sequence1[, 1] %in% "Adduct_pos"]
    ionmode1 <- "pos"
  } else {
    adduct <- dataset1[, sequence1[, 1] %in% "Adduct_neg"]
    ionmode1 <- "neg"
  }
  mass <- monomass(adduct, mass, ionmode1)
  rt <- dataset1[, sequence1[, 1] %in% "RT"]
  first <- data.frame(mass, rt)
  updateProgressBar(id = "pb", value = 10)
  mass <- dataset2[, sequence2[, 1] %in% "Mass"]
  if ("Adduct_pos" %in% sequence2[, 1]) {
    adduct <- dataset2[, sequence2[, 1] %in% "Adduct_pos"]
    ionmode2 <- "pos"
  } else {
    adduct <- dataset2[, sequence2[, 1] %in% "Adduct_neg"]
    ionmode2 <- "neg"
  }
  mass <- monomass(adduct, mass, ionmode2)
  rt <- dataset2[, sequence2[, 1] %in% "RT"]
  second <- data.frame(mass, rt)
  updateProgressBar(id = "pb", value = 20)
  comb <- rbind(first, second)
  distmz <- as.matrix(dist(comb$mass))
  updateProgressBar(id = "pb", value = 30)
  v <- rep(comb$mass, each = dim(distmz)[1])
  updateProgressBar(id = "pb", value = 40)
  distp <- distmz / v
  updateProgressBar(id = "pb", value = 50)
  distppm <- distp * 10^6
  updateProgressBar(id = "pb", value = 60)
  distrt <- as.matrix(dist(comb$rt))
  updateProgressBar(id = "pb", value = 70)
  adj <- distppm <= ppmTolerance & distrt <= rtTolerance
  updateProgressBar(id = "pb", value = 80)
  graph <- graph.adjacency(adj)
  updateProgressBar(id = "pb", value = 90)
  mergeid <- clusters(graph)$membership
  updateProgressBar(id = "pb", value = 100)
  colnames(dataset2) <- colnames(dataset1)
  combineddat <- as.data.frame(rbind(dataset1, dataset2))
  combineddat[, "mergeID"] <- mergeid
  ionmode1 <- rep(ionmode1, nrow(dataset1))
  ionmode2 <- rep(ionmode2, nrow(dataset2))
  combineddat[, "ionmode"] <- c(ionmode1, ionmode2)
  closeSweetAlert()
  return(combineddat)
}

monomass <- function(adduct, mz, ionmode) {
  masscorrection <- read.csv("./csvfiles/adducts.csv")
  monomass <- sapply(seq(adduct), function(i) {
    if (adduct[i] %in% masscorrection[, 1]) {
      j <- which(masscorrection[, 1] %in% adduct[i])
      cm <- masscorrection[j, 3]
      mass <- masscorrection[j, 2]
      monomass <- mz[i] * cm - mass
      return(monomass)
    } else if (ionmode == "pos") {
      return(mz[i] - 1.007276)
    } else {
      return(mz[i] + 1.007276)
    }
    return(monomass)
  })
}

finddup <- function(out_dup, rankings) {
  prio <- apply(out_dup, 1, function(x) {
    duplicaterank(x, rankings)
  })
  dup_prio <- cbind(out_dup, prio)
  rows <- NULL
  for (i in unique(dup_prio[, 2])) {
    lowp <- dup_prio[i == dup_prio[, 2], ][, 9] %in% min(dup_prio[i == dup_prio[, 2], ][, 9])
    mincv <- min(dup_prio[i == dup_prio[, 2], ][lowp, 8])
    keeprow <- rownames(dup_prio[i == dup_prio[, 2], ][lowp, ][dup_prio[i == dup_prio[, 2], ][lowp, 8] == mincv, ])
    keeprow <- keeprow[1]
    rows <- c(rows, keeprow)
  }
  numb <- which(rownames(out_dup) %in% rows)
  return(numb)
}

duplicaterank <- function(duplicate, rankings) {
  j <- sapply(1:10, function(x) {
    if (rankings[x, 1] == "") {
      FALSE
    } else {
      grepl(toupper(rankings[x, 1]), toupper(duplicate[5]))
    }
  })
  if (sum(j) > 0) {
    return(min(rankings[j, 2]))
  } else {
    return(10)
  }
}

# PolySTest
addEmptyCols <- function(data, sequence, groups, replicates) {
  processed <- data[, 1] # feature names
  rgroup <- c("") # group vector
  rtime <- c("") # time vector

  for(group in 1:length(groups)) {
    groupCols <- data[, sequence[, 4] %in% groups[group]]
    time <- sequence[sequence[, 4] %in% groups[group], 5]
    processed <- cbind(processed, groupCols)
    if(length(groupCols) < replicates) {
      missing <- t(rep(NA, replicates - length(groupCols)))
      processed <- cbind(processed, missing)
    }
    rgroup <- append(rgroup, rep(paste("g", groups[group], sep = ""), replicates))
    if(any(complete.cases(sequence[, 5])))
      rtime <- append(rtime, t(time))
  }
  if(any(complete.cases(sequence[, 5])))
    colnames(processed) <- paste(colnames(processed), rgroup, paste("t", rtime, sep=""), sep = "_")
  else
    colnames(processed) <- paste(colnames(processed), rgroup, sep = "_")
  return(processed)
}

windowselect <- function(input) {
  hide("welcome_panel")
  if (input == "pca") show("pca_panel") else hide("pca_panel")
  if (input == "drift") show("drift_panel") else hide("drift_panel")
  if (input == "feature") show("boxplot_panel") else hide("boxplot_panel")
  if (input == "sequence") show("sequence_panel") else hide("sequence_panel")
  if (input == "datatable") show("datatable_panel") else hide("datatable_panel")
  if (input == "statistics") show("statistics_panel") else hide("statistics_panel")
  if (input == "export") show("export_panel") else hide("export_panel")
  if (input == "info") show("info_panel") else hide("info_panel")
}
