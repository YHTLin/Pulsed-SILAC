### Analysis of pulsed SILAC proteomics data (FOR OZLEM)

#############################################
# Step 1 - Set working directory
#############################################
setwd("C:/Users/Tony Lin/Desktop/Wiita_lab/Projects/Proteomics_project/Pulsed_silac_Ozlem/Proteomics/analysis/")


#############################################
# Step 2 - Read in proteinGroups file
#############################################
raw.data = read.csv("proteinGroups.csv", header = TRUE, colClasses = "character")


#############################################
# Step 3 - Data clean up
#############################################
# Remove reverse proteins, contaminants, and only identified by site
require(dplyr)
data = raw.data
data = filter(data, data$Reverse != "+")
data = filter(data, data$Potential.contaminant != "+")
data = filter(data, data$Only.identified.by.site != "+")


#############################################
# Step 4 - Data pre-processing
#############################################
# Assign gene symbols according to Fasta headers
require(stringr)
first_fasta = strsplit(data$Fasta.headers, ";")   # Extract the first fasta header in each row
data$Symbol = sapply(first_fasta, function(x) 
  str_extract(x[1], "(?<=\\|.{1,10}\\|).*(?=_MOUSE)"))

# Assign gene names according to Fasta headers
data$Gene.name = sapply(first_fasta, function(x)
  str_extract(x[1], "(?<=\\|.{1,10}\\|.{1,30}_MOUSE).*(?=OS\\=)"))

# Reduce Protein.IDs and Majority.protein.IDs to UniProtID
extractID = function(df, col.name, sep = ";") {
  # df = data frame containing data of interest
  # col.name = a character of length 1 describing the name of a column in df
  # sep = pattern for splitting string up into pieces
  require(stringr)
  ids = str_split(df[, col.name], ";")
  df[, col.name] = sapply(ids, # Extracts UniProtID and stitches the results together separated by semicolon
                          function(x) paste(str_extract(x, "(?<=\\|).*(?=\\|)"), collapse = ";"))
  return(df)
}
data = extractID(data, "Protein.IDs")
data = extractID(data, "Majority.protein.IDs")

# Cast "Intensity" data as numeric
intensity.names = names(select(data, starts_with("Intensity")))
data[intensity.names] = lapply(data[intensity.names], function(x) as.numeric(x))

# Compute LOG2 on intensity
log2Int.names = sub("^Intensity", "LOG2", intensity.names)
data[log2Int.names] = log2(data[intensity.names])

# Cast "Ratio.H.L." data as numeric
ratio.names = names(select(data, starts_with("Ratio.H.L.")))
ratio.names = ratio.names[!grepl("type", ratio.names)]  # Exlude Ratio.H.L.type columns
data[ratio.names] = lapply(data[, ratio.names], function(x) as.numeric(x))

# Compute LOG2 on Ratio.H.L
ratio.names = names(select(data, matches("^Ratio\\.H\\.L\\.(DMSO|M1071|Rapa|Heavy_incorp)")))
log2Ratio.names = sub("^Ratio", "LOG2", ratio.names)
data[log2Ratio.names] = log2(data[ratio.names])


#############################################
# Step 5 - Data filtering
#############################################
# Organize P4heavy data (heavy incorporation after passage 4)
protein.info = c("Protein.IDs", "Majority.protein.IDs", "Number.of.proteins", "Peptides",
                 "Symbol", "Gene.name")
P4heavy.names = c("LOG2.L.Heavy_incorp", "LOG2.H.Heavy_incorp", 
                  "LOG2.H.L.Heavy_incorp", "Ratio.H.L.count.Heavy_incorp")
P4heavy = select(data, c(protein.info, P4heavy.names))

# Organize light data
require(gtools)   #for invoking mixedsort to order column labels
light.names = mixedsort(grep("^LOG2\\.L\\.(DMSO|M1071|Rapa)", 
                             names(data), value = TRUE))
light = select(data, c(protein.info, light.names))

# Organize heavy data
heavy.names = mixedsort(grep("^LOG2\\.H\\.(DMSO|M1071|Rapa)", 
                             names(data), value = TRUE))
heavy = select(data, c(protein.info, heavy.names))

# Organize H/L ratio data
HL.names = c(mixedsort(grep("^LOG2\\.H\\.L\\.(DMSO|M1071|Rapa)", names(data), value = TRUE)),
             mixedsort(grep("^Ratio\\.H\\.L\\.count\\.(DMSO|M1071|Rapa)", names(data), value = TRUE)))
ratio = select(data, c(protein.info, HL.names))
rm(protein.info, light.names, heavy.names, HL.names)

# Filter out missing values
filter_valids = function(df, conditions, min_count, 
                         is_infinite = TRUE, at_least_one = FALSE) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # is_infinite = Boolean indicating the nature of missing data
  #     if TRUE, counts Inf/-Inf values as missing values
  #     if FALSE, counts NaN values as missing values
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  require(dplyr)
  
  log2.names = names(select(df, starts_with("LOG2")))   # Extract LOG2 columns
  cond.names = sapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE))
  
  if (is_infinite) {
    cond.filter = sapply(1:length(cond.names), function(i) {
      df2 = df[cond.names[[i]]]   # Extract columns of interest
      df2 = as.matrix(df2)   # Cast as matrix for the following command
      sums = rowSums(!is.infinite(df2)) # count the number of valid values for each condition
      sums >= min_count[i]   # Calculates whether min_count requirement is met
    })
  } else {
    cond.filter = sapply(1:length(cond.names), function(i) {
      df2 = df[cond.names[[i]]]   # Extract columns of interest
      df2 = as.matrix(df2)   # Cast as matrix for the following command
      sums = rowSums(!is.nan(df2)) # count the number of valid values for each condition
      sums >= min_count[i]   # Calculates whether min_count requirement is met
    })
  }
  
  if (at_least_one) {
    df$KEEP = apply(cond.filter, 1, any)
  } else {
    df$KEEP = apply(cond.filter, 1, all)
  }
  
  return(df)  # No rows are omitted, filter rules are listed in the KEEP column
}
lightF = filter_valids(light, c("DMSO", "M1071", "Rapa"), c(5, 4, 4), # F for filtered
                       is_infinite = TRUE, at_least_one = FALSE) # contains infinite as invalids
heavyF = filter_valids(heavy, c("DMSO", "M1071", "Rapa"), c(5, 4, 4), # F for filtered
                       is_infinite = TRUE, at_least_one = FALSE) # contains infinite as invalids
ratioF = filter_valids(ratio, c("DMSO", "M1071", "Rapa"), c(5, 4, 4), # F for filtered
                       is_infinite = FALSE, at_least_one = FALSE) # contains infinite as invalids


#############################################
# Step 6 - Data normalization (light and heavy channels)
#############################################
# Perform global median normalization to correct for variable sample concentration/injection volume
# Adjust the light and heavy intensity by a normalization factor calculated from the combined intensity
sample_centering = function(df, df2) {
  # df = dataframe containing filtered light or heavy channel data (output filter_valids)
  # df2 = dataframe containing the summed H + L intensity data
  require(dplyr)
  filter.df = select(df2, matches("^LOG2\\.(DMSO|M1071|Rapa)"))
  norm.factor = sapply(filter.df, median)  # Include the -Infs to use all ID'd proteins (-Inf represents lowly expressed proteins)
  names(norm.factor) = sub("^LOG2\\.", "", names(norm.factor))
  for (sample in names(norm.factor)) {
    label = grep(sample, names(df), value = TRUE)
    df[label] = df[label] - norm.factor[sample]   # Subtract sample median from all observations
  }
  return(df)
}
lightFC = sample_centering(lightF, data)
heavyFC = sample_centering(heavyF, data)

# Normalize light and heavy channel data with DMSO 0hr
DMSO_norm = function(df) {
  # df = dataframe containing filtered and normalized light or heavy channel data
  log.names = grep("^LOG2.*(3hr|6hr|12hr|24hr)", names(df), value = TRUE)
  zero.name = grep("^LOG2.*0hr", names(df), value = TRUE)
  
  matr = as.matrix(df[log.names])  # Casting required for the following operation
  for (i in 1:nrow(matr)) {
    if (!is.infinite(df[i, zero.name])) {
      matr[i, ] = matr[i, ] - df[i, zero.name]
    }
  }
  
  norm.names = sub("^LOG2", "NORM", log.names)
  df[norm.names] = matr
  return(df)
}
lightFCD = DMSO_norm(lightFC)
heavyFCD = DMSO_norm(heavyFC)


#############################################
# Step 7 - Score treatment/control rate of change in intensity over time
#############################################
# Use the LOG2 or NORM columns to calculate slope of linear regression
score = function(df, sample.names, suffix) {
  require(gtools)   # For mixedsort
  # df = data frame containing LOG2 columns for computing differences in rate of change
  # sample.names = vector of regex patterns for extracting column names for each condition
  # suffix = character vector describing the name of each condition for output
  
  log.names = lapply(sample.names, function(x) grep(x, names(df), value = TRUE))
  df2 = lapply(log.names, function(x) {   # Loop through each condition
    x = mixedsort(x)  # order the columns alphanumerically
    apply(df[x], 1, function(y) data.frame(time = 1:length(y), value = y))  # Loop through each row
  })
  # Apply linear model after excluding -Inf rows
  slopes = lapply(df2, function(x) {   # Loop through each condition
    sapply(x, function(y) {
      y = y[!is.infinite(y$value), ]   # Remove -Inf rows
      if (nrow(y) > 1) {   # Only perform linear regression if two or more data points are available
        lm(value ~ time, y)$coefficients[2]   # Extract slope from linear model
      } else {   # Otherwise return NA
        NA
      }
    })
  })
  # Output columns
  df[paste0("slope_", suffix)] = as.data.frame(slopes)
  return(df)
}
# score_ is the difference in slope between treatment and DMSO control
lightFCDS = score(lightFCD, c("^LOG2.*DMSO", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
lightFCDS = mutate(lightFCDS, score_M1071 = slope_M1071 - slope_DMSO)
lightFCDS = mutate(lightFCDS, score_Rapa = slope_Rapa - slope_DMSO)
heavyFCDS = score(heavyFCD, c("^LOG2.*DMSO", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
heavyFCDS = mutate(heavyFCDS,score_M1071 = slope_M1071 - slope_DMSO)
heavyFCDS = mutate(heavyFCDS, score_Rapa = slope_Rapa - slope_DMSO)


#############################################
# Step 8 - Data visualization
#############################################
#INSTALL SUPERHEAT FOR VISUALIZATION
#devtools::install_github("rlbarter/superheat")
### Heatmaps for light and heavy data (NOT RATIO)
plot_heatmap = function(df, pattern) {
  # df = data frame containing filtered/normalized data for plotting
  # pattern = regular expression pattern for selecting column names to plot
  require(superheat)
  require(RColorBrewer)
  require(gtools)   #for invoking mixedsort to order column labels
  require(stringr)
  # Extract column names
  label = mixedsort(grep(pattern, names(df), value = TRUE))
  
  # Turn invalid values into NA
  plot.df = as.matrix(df[df$KEEP, label])  # Use only filtered rows
  plot.df[is.infinite(plot.df) | is.nan(plot.df)] = NA   # Required for superheat
  colnames(plot.df) = sub(paste0(pattern, "\\."), "", colnames(plot.df))   # name columns
  rownames(plot.df) = df[df$KEEP, "Symbol"]   # name rows
  
  # Plot figure
  superheat(plot.df,
            
            # set heatmap color map for missing values
            heat.na.col = "black",
            
            # top line plot: sample medians
            # Calculate column median based on valid values
            # Ideally includes the NAs for light and heavy data but this is a close enough approximation
            yt = apply(plot.df, 2, function(x) median(x, na.rm = TRUE)), 
            yt.plot.type = "scatterline",
            yt.point.size = 3,
            yt.line.size = 1,
            yt.axis.name = "Sample median\n(log-2 transformed)",
            yt.cluster.col = c("#a6d854", "#fc8d62", "#8da0cb"),
            
            # Order columns
            membership.cols = sub("_.*", "", colnames(plot.df)),  # Important for placing vertical lines
            order.cols = 1:ncol(plot.df),
            
            # Separate conditions by white lines
            grid.vline.col = "white",
            grid.vline.size = 2,
            
            # Hierarchial clustering for rows
            pretty.order.rows = TRUE,
            dist.method = "euclidean",
            
            # bottom labels
            bottom.label = "variable",
            bottom.label.size = 0.25,
            bottom.label.text.size = 3.25,
            bottom.label.col = "white",
            bottom.label.text.angle = 90,
            bottom.label.text.alignment = "center",
            
            # left labels
            force.left.label = FALSE,
            left.label = "none",
            left.label.size = 0.1,
            left.label.text.size = 1,
            left.label.col = "white")
}
#plot_heatmap(lightFCDS, "^NORM")
#plot_heatmap(lightFCDS, "^LOG2")
#plot_heatmap(heavyFCDS, "^NORM")
#plot_heatmap(heavyFCDS, "^LOG2")
#plot_heatmap(ratioF, "^LOG2")


### Pair plot
plot_pairs = function(df, pattern, use_keep = FALSE) {
  # df = data frame carrying data of interest
  # pattern = regex pattern to select column of interest
  # use_keep = TRUE means to construct plot on filtered values; FALSE uses all available values
  require(gpairs)
  require(scales)
  if (use_keep) {
    plot.df = df[df$KEEP, grep(pattern, names(df), value = TRUE)]
  } else {
    plot.df = df[grep(pattern, names(df), value = TRUE)]
  }
  
  gpairs(plot.df,
         upper.pars = list(scatter = "lm"),
         scatter.pars = list(pch = 20,
                             col = alpha("black", 0.3)),
         lower.pars = list(scatter = "stats"),
         stat.pars = list(verbose = FALSE, fontsize = 15))
}
## Pairs plot: light label after DMSO normalization
#plot_pairs(lightFCDS, "^NORM.*DMSO")
#plot_pairs(lightFCDS, "^NORM.*M1071")
#plot_pairs(lightFCDS, "^NORM.*Rapa")
## Pairs plot: light label before DMSO normalization
#plot_pairs(lightFCDS, "^LOG2.*DMSO")
#plot_pairs(lightFCDS, "^LOG2.*M1071")
#plot_pairs(lightFCDS, "^LOG2.*Rapa")
## Pairs plot: heavy label after DMSO normalization
#plot_pairs(heavyFCDS, "^NORM.*DMSO")
#plot_pairs(heavyFCDS, "^NORM.*M1071")
#plot_pairs(heavyFCDS, "^NORM.*Rapa")
## Pairs plot: heavy label before DMSO normalization
#plot_pairs(heavyFCDS, "^LOG2.*DMSO")
#plot_pairs(heavyFCDS, "^LOG2.*M1071")
#plot_pairs(heavyFCDS, "^LOG2.*Rapa")


### Histograms of heavy vs light in each sample
plot_hist = function(df, heavy, light, compare_time = FALSE) {
  # df = data frame containing Log2 intensity
  # heavy, light = regex patterns dictating rules for extracting light and heavy columns
  # compare_time = logical where TRUE compares time within each plot and FALSE compares SILAC labeling
  require(ggplot2)
  require(tidyr)
  require(stringr)
  
  light.names = grep(light, names(df))
  heavy.names = grep(heavy, names(df))
  # Re-organize data frame
  df = gather(df[, c(light.names, heavy.names)], "key", "value", 
              1:length(c(light.names, heavy.names)))
  df$value[is.infinite(df$value)] = NA
  # Denote SILAC label
  df$label = str_extract(df$key, "(?<=^.{5}).")  
  # Denote treatment
  df$condition = str_extract(df$key, "(?<=^.{7}).{4,5}(?=_)")   
  # Denotes time point
  df$time = str_extract(df$key, "(?<=_).*(?=_)")
  df$time = factor(df$time, levels = c("3hr", "6hr", "12hr", "24hr"))   # Maintain chronological order
  
  if (compare_time) {
    ggplot(df, aes(x = value, fill = time)) + 
      geom_histogram(alpha = 0.2, binwidth = 0.2, 
                     position = "identity", na.rm = TRUE) +
      facet_grid(condition~label, switch = "both") +
      labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
      scale_fill_discrete(name = "Timepoint")
  } else {
    ggplot(df, aes(x = value, fill = label)) + 
      geom_histogram(alpha = 0.3, binwidth = 0.2,
                     position = "identity", na.rm = TRUE) +
      facet_grid(condition~time, switch = "both") +
      labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
      scale_fill_discrete(name = "SILAC",
                          breaks = c("H", "L"),
                          labels = c("Heavy", "Light"))
  }
}
## Compare Time
#plot_hist(cbind(lightFCDS, heavyFCDS), compare_time = TRUE,
#          "^LOG2\\.L.*(3hr|6hr|12hr|24hr)", "^LOG2\\.H.*(3hr|6hr|12hr|24hr)")
#plot_hist(cbind(lightFCDS, heavyFCDS), compare_time = TRUE,
#          "^NORM\\.L.*(3hr|6hr|12hr|24hr)", "^NORM\\.H.*(3hr|6hr|12hr|24hr)")
## Compare SILAC labeling
#plot_hist(cbind(lightFCDS, heavyFCDS), compare_time = FALSE,
#          "^LOG2\\.L.*(3hr|6hr|12hr|24hr)", "^LOG2\\.H.*(3hr|6hr|12hr|24hr)")
#plot_hist(cbind(lightFCDS, heavyFCDS), compare_time = FALSE,
#          "^NORM\\.L.*(3hr|6hr|12hr|24hr)", "^NORM\\.H.*(3hr|6hr|12hr|24hr)")


### Histograms of heavy vs light in fully-incorporated light (DMSO_0hr) and heavy (Heavy_incorp) samples
plot_hist2 = function(df) {
  # df = data frame containing Log2 intensity
  # Designed specifically to create one histogram
  require(ggplot2)
  require(tidyr)
  require(stringr)
  sample.names = c("LOG2.H.DMSO_0hr_bR1", "LOG2.H.Heavy_incorp", 
                   "LOG2.L.DMSO_0hr_bR1", "LOG2.L.Heavy_incorp")
  # Re-organize data frame
  df = gather(df[, sample.names], "key", "value", 
              1:length(sample.names))
  df$value[is.infinite(df$value)] = NA
  # Denote SILAC label
  df$label = str_extract(df$key, "(?<=^.{5}).")  
  # Denote treatment (assign plot title)
  df$condition = str_extract(df$key, "(?<=^.{7}).*")
  df$condition = factor(df$condition)
  levels(df$condition) <- list(`Light incorporation` = "DMSO_0hr_bR1",  # Renaming levels of factor
                               `Heavy incorporation` = "Heavy_incorp")

  ggplot(df, aes(x = value, fill = label)) + 
    geom_histogram(alpha = 0.3, binwidth = 0.2,
                   position = "identity", na.rm = TRUE) + 
    facet_grid(.~condition) +
    labs(x = expression("log"[2]*"-transformed Intensity"), y = "Frequency") +
    scale_fill_discrete(name = "SILAC",
                        breaks = c("H", "L"),
                        labels = c("Heavy", "Light"))
}
## Perform global median normalization (centering) on P4heavy
P4heavyC = P4heavy
P4heavyC$LOG2.L.Heavy_incorp = P4heavy$LOG2.L.Heavy_incorp - median(data$LOG2.Heavy_incorp)
P4heavyC$LOG2.H.Heavy_incorp = P4heavy$LOG2.H.Heavy_incorp - median(data$LOG2.Heavy_incorp)
#plot_hist2(cbind(lightFCDS, heavyFCDS, P4heavyC))


### Table of summary statistics
summary_stat = function(df) {
  # df = data frame containing LOG2 intensity
  
  require(dplyr)
  df2 = select(df, starts_with("LOG2"))
  data.frame(
    Sample = sub("(LOG2.L.|LOG2.H.|LOG2.H.L.)", "", names(df2)),
    Raw = sapply(df2, function(x) {
      paste0(sum(!is.infinite(x) & !is.na(x)), " / ", length(x), 
             " (", round(sum(!is.infinite(x) & !is.na(x)) / length(x) * 100, 1), "%)")
    }),
    Filtered = sapply(df2[df$KEEP,], function(x) {
      paste0(sum(!is.infinite(x) & !is.na(x)), " / ", length(x), 
             " (", round(sum(!is.infinite(x) & !is.na(x)) / length(x) * 100, 1), "%)")
    })
  )
}
#summary_stat(lightFCDS)
#summary_stat(heavyFCDS)
#summary_stat(ratioF)


#############################################
# Step 9 - Save workspace as Rmd file
#############################################
save.image("pulsed_SILAC.RData")

