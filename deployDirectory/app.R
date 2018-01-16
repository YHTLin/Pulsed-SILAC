### ShinyApp accompanying pulsed_SILAC_ozlem.R
## Some functions from pulsed_SILAC_ozlem.R has been reordered for app development

#############################################
# Step 1 - Set working directory
#############################################
#setwd('C:/Users/Tony Lin/Desktop/Wiita_lab/Projects/Proteomics_project/Pulsed_silac_Ozlem/Proteomics/analysis/deployDirectory')


#############################################
# Step 2 - Read in proteinGroups file
#############################################
raw.data = read.csv("./Data/proteinGroups.csv", header = TRUE, colClasses = "character")


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
# Step 5 - Data organization
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


#############################################
# Step 6 - Data normalization (light and heavy channels)
#############################################
# Perform global median normalization to correct for variable sample concentration/injection volume
# Adjust the light and heavy intensity by a normalization factor calculated from the combined intensity
sample_centering = function(df, df2) {
  # df = dataframe containing filtered light or heavy channel data (output filter_valids)
  # df2 = dataframe containing the summed heavy and light intensity data
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
lightC = sample_centering(light, data)
heavyC = sample_centering(heavy, data)

# Normalize light and heavy channel data with DMSO 0hr
DMSO_norm = function(df) {
  # df = dataframe containing filtered and normalized light or heavy channel data
  log.names = grep("^LOG2.*(3hr|6hr|12hr|24hr)", names(df), value = TRUE)
  zero.name = grep("^LOG2.*0hr", names(df), value = TRUE)
  
  matr = as.matrix(df[log.names])  # Casting required for the following operation
  for (i in 1:nrow(matr)) {
    if (!is.infinite(df[i, zero.name]) & !is.na(df[i, zero.name])) {
      matr[i, ] = matr[i, ] - df[i, zero.name]
    }
  }
  
  norm.names = sub("^LOG2", "NORM", log.names)
  df[norm.names] = matr
  return(df)
}
lightCD = DMSO_norm(lightC)
heavyCD = DMSO_norm(heavyC)
ratioD = DMSO_norm(ratio)


#############################################
# Step 7 - Score treatment/control rate of change in intensity over time
#############################################
# Use the LOG2 or NORM columns to calculate slope of linear regression
score = function(df, sample.names, suffix) {
  require(gtools)   # For mixedsort
  # df = data frame containing LOG2 columns for computing differences in rate of change
  # sample.names = vector of regex patterns for extracting column names for each condition
  # suffix = character vector describing the name of each condition for output
  
  log.names = lapply(sample.names, function(x) grep(x, names(df), value = TRUE, perl = TRUE)) #perl = TRUE enables negative lookahead
  df2 = lapply(log.names, function(x) {   # Loop through each condition
    x = mixedsort(x)  # order the columns alphanumerically
    apply(df[x], 1, function(y) data.frame(time = 1:length(y), value = y))  # Loop through each row
  })
  # Apply linear model after excluding -Inf rows
  slopes = lapply(df2, function(x) {   # Loop through each condition
    sapply(x, function(y) {
      y = y[!is.infinite(y$value) & !is.na(y$value), ]   # Remove -Inf or NA rows
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
# IGNORE DMSO_0HR IN LINEAR REGRESSION BECAUSE HEAT MAP SUGGESTS THAT THE SAMPLE COULD BE MISLABLED
lightCDS = score(lightCD, c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
lightCDS = mutate(lightCDS, score_M1071 = slope_M1071 - slope_DMSO)
lightCDS = mutate(lightCDS, score_Rapa = slope_Rapa - slope_DMSO)
heavyCDS = score(heavyCD, c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
heavyCDS = mutate(heavyCDS, score_M1071 = slope_M1071 - slope_DMSO)
heavyCDS = mutate(heavyCDS, score_Rapa = slope_Rapa - slope_DMSO)
ratioDS = score(ratioD, c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
ratioDS = mutate(ratioDS, score_M1071 = slope_M1071 - slope_DMSO)
ratioDS = mutate(ratioDS, score_Rapa = slope_Rapa - slope_DMSO)


#############################################
# Step 8 - Data Filtering Function
#############################################
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


#############################################
# Step 8 - Data visualization functions
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


#############################################
# Step 9 - ShinyApp code
#############################################
library(shinythemes)
library(shiny)
library(DT)  # For interactive data tables

# USER INTERFACEs ------------------------------------------
ui = navbarPage(
  "Pulsed SILAC Proteomics",
  theme = shinytheme("sandstone"),
  #theme = shinytheme("simplex"),
  
  # ABOUT PAGE ----------------------------------------------
  tabPanel("About",
           h2("Interactive Proteomics"),
           hr(),
           fluidRow(
             column(8,
               h3("Motivation"),
               p(HTML("The aim of this study is to identify proteins whose biosyntheses are disrupted 
                      by mTOR inhibitors
                      rapamycin and RapaLink-1 (M1071) in vivo. To this end, we have employed mass 
                      spectrometry-based proteomics in conjunction with a technique called stable isotope 
                      labeling with amino acids in cell culture (SILAC). This webpage helps showcase
                      and visualize the proteomics data in an interactive way.")),
               br(),
               h3("Experimental Design"),
               p(HTML("To start, GTML5 mouse cells were conditioned in light SILAC media (Arg-0 and Lys-0). 
                      When fully light-labeled, the cells were split into three flasks containing 
                      heavy SILAC media (Arg-10 and Lys-8) and immediately treated with DMSO, M1071, 
                      or rapamycin. Cultures from each condition were harvested at the 3-hour, 6-hour,
                      12-hour, and 24-hour timepoint. The process was repeated for
                      a biological replicate.")),
               br(),
               h3("Sample Preparation"),
               p(HTML("A standard proteomics workflow was followed to prepare the samples for liquid
                      chromatography-tandem mass spectrometry (LC-MS/MS) analysis on a Thermo
                      Q-Exactive Plus. Cell lysis, protein digestion, and peptide clean-up were 
                      performed prior to data acquisition on a 4-hour gradient. <b>Both biological 
                      replicates were prepped but only one was analyzed by LC-MS/MS and included
                      in the current analysis.</b>")),
               br(),
               h3("Data Analysis"),
               p(HTML("The acquired data were analyzed by the proteomics software MaxQuant to
                      match the identified spectra to peptides and to derive quantitative measures
                      of protein abundance. Appropriate SILAC settings were used to detect both 
                      light- and heavy-labeled proteins across all samples. Please visit the 
                      <b>Download</b> page for the raw MaxQuant output.")),
               p(HTML("The MaxQuant output was processed in a series of steps prior to visualization.
                      First, the data were cleaned and organized into three data types: 
                      <b>Heavy</b>, <b>Light</b>, and <b>Ratio</b> representing the protein abundance
                      in the heavy channel, protein abundance in the light channel, and the
                      heavy-to-light ratio, respectively. For more information on how MaxQuant
                      determines these values, please see the "),
                 a("documentation", href = "https://www.nature.com/articles/nprot.2016.136", target="_blank"),
                 HTML("online.")),
               p(HTML("Next, the protein intensity values were log<sub>2</sub>-transformed and 
                      and corrected
                      to account for the variability in sample injection on the mass spectrometer.
                      These values are reported in the <b>LOG2</b> columns in the downloadable CSVs.
                      Additionally, in order to account for the varying baseline protein abundances, the 
                      <b>NORM</b> columns were computed by taking the difference between the 
                      log<sub>2</sub>-transformed intensity of a protein in a sample and its 
                      corresponding value in the DMSO 0-hour control.")),
               p(HTML("Finally, a filter to remove missing values can be applied independently to 
                      each of the analysis pages above. The controls are initialized at the most selective 
                      threshold and can be adjusted by the user. Take note that some graphics may 
                      throw an error at low thresholds."))
               )
             ),
           hr(),
           helpText("Please direct any questions to Tony Lin (lin.yu.hsiu@ucsf.edu).")
  ),
  
  # HEAVY PAGE ----------------------------------------------
  tabPanel("Heavy",
    h2("Protein Synthesis Analysis"),
    hr(),
    # CONSTRUCT FILTER PANEL
    wellPanel(
      h4("Filter Settings"),
      
      p("Adjust the slider below to indicate the minimum number of valid values",
        "required for retention."),
      
      fluidRow(
        # Input for missing value filter: DMSO
        column(4, sliderInput("DMSO_valid_heavy", label = "DMSO",
                              min = 0, max = 5, value = 5)
        ),
        # Input for missing value filter: M1071
        column(4, sliderInput("M1071_valid_heavy", label = "M1071", 
                              min = 0, max = 4, value = 4)
        ),
        # Input for missing value filter: Rapa
        column(4, sliderInput("Rapa_valid_heavy", label = "Rapamycin",
                              min = 0, max = 4, value = 4))
      ),
        
      # Input for EACH or AT LEAST ONE
      radioButtons("filter_type_heavy", label = "Filter Type",
                   choices = list("Keep observations for which ALL conditions above are met" = FALSE, 
                                  "Keep observations for which ANY condition above is met" = TRUE),
                   selected = FALSE),
        
      # Help text
      helpText(HTML("Note: All tabs below will automatically apply the filter. 
                    Low thresholds could interfere with heatmap construction."))
    ),
    
    # CONSTRUCT PAGE TABS
    tabsetPanel(
      tabPanel("Summary",
               fluidRow(
                 column(3, br(),  # Vertical spacing
                        p("Summary statistics are displayed for each sample in the form of ",
                          tags$b("X"), "/", tags$b("Y"), "(", tags$b("Z"), "%), where", tags$b("X"), 
                          " is the number of proteins quantified,", tags$b("Y"), "is the number identified, and",
                          tags$b("Z") , "is the expression in percentages.")
                        ),
                 column(9, br(), tableOutput("summary_table_heavy")))
               ),
      tabPanel("Pairwise Plot",
               fluidRow(
                 column(3, br(),  # Vertical spacing
                        p(HTML("Pairwise plots show three types of information.
                               Histograms display the distribution of
                               log<sub>2</sub>-transformed protein intensities along the diagonal.
                               Additionally, scatter plots between any two 
                               timepoints are shown at the upper right corner along with lines of best fit 
                               in red. The corresponding correlation coefficients are found at the lower left 
                               corner. Proteins quantified in one sample but not the other are counted as \"missing\".
                               Data are organized by treatment groups.")),
                        hr(),
                        selectInput('pair_plot_heavy_display', 'Display Type:',
                                    c("NORM", "LOG2"), selectize = FALSE),
                        helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                 to the DMSO 0-hour timepoint.")
                        ),
                 column(9, br(),
                        h4("DMSO"),
                        plotOutput("pair_plot_heavy_DMSO", height = "700px", width = "700px"),
                        hr(),
                        h4("M1071"),
                        plotOutput("pair_plot_heavy_M1071", height = "700px", width = "700px"),
                        hr(),
                        h4("Rapamycin"),
                        plotOutput("pair_plot_heavy_Rapa", height = "700px", width = "700px")
                        ))
               ),
      tabPanel("Heatmap",
               fluidRow(
                 column(3, br(),
                        p(HTML("Heatmaps represent the values of the <b>filtered</b> data table in color.
                               The samples are organized into columns, and each row shows the 
                               log<sub>2</sub>-transformed intensities for a protein. The sample medians 
                               are plotted as a line chart above the figure. See legend below for 
                               color scale. Missing values are colored in black.")),
                        selectInput('heatmap_heavy_display', 'Display Type:',
                                    c("NORM", "LOG2"), selectize = FALSE),
                        helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                 to the DMSO 0-hour timepoint. NORM is recommended.")
                        ),
                 column(9, br(),
                        plotOutput("heatmap_heavy", height = "700px", width = "450px")
                        ))
               ),
      tabPanel("Filtered Table", br(),
               p(HTML("The data table below expresses the pulsed SILAC proteomics output as slopes 
                      and scores for ease of evaluation. A linear regression model is fitted on the
                      log<sub>2</sub>-transformed intensities for each condition. Here, timepoints are used
                      as the independent variable while the intensities as the dependent. At least two
                      data points are required for modeling. No slopes are reported for those that fail 
                      to meet this criterion. Scores are computed by subtracting the slopes of the 
                      DMSO control from those of the treatment group. A negative score suggests that the rate of 
                      synthesis decreases with treatment and a positive score indicates the opposite effect. 
                      <b>All values are extremely sensitive to outliers. Please interpret with caution.</b> 
                      A complete data table including the measured intensities can be downloaded below.")),
               p(HTML("<b>NOTE: There is reason to believe that DMSO-0hr was mishandled or mislabeled (see 
                             heatmap under \"Ratio\" page). It is thus omitted before applying the linear model.</b>")),
               downloadButton("DL_filtered_heavy", "Download Filtered CSV"),
               hr(),
               DT::dataTableOutput('filtered_table_heavy')),
      tabPanel("Raw Table", br(),
               p(HTML("The data table below expresses the pulsed SILAC proteomics output as slopes 
                      and scores for ease of evaluation. A linear regression model is fitted on the
                      log<sub>2</sub>-transformed intensities for each condition. Here, timepoints are used
                      as the independent variable while the intensities as the dependent. At least two
                      data points are required for modeling. No slopes are reported for those that fail 
                      to meet this criterion. Scores are computed by subtracting the slopes of the 
                      DMSO control from those of the treatment group. A negative score suggests that the rate of 
                      synthesis decreases with treatment and a positive score indicates the opposite effect. 
                      <b>All values are extremely sensitive to outliers. Please interpret with caution.</b> 
                      A complete data table including the measured intensities can be downloaded below.")),
               p(HTML("<b>NOTE: There is reason to believe that DMSO-0hr was mishandled or mislabeled (see 
                             heatmap under \"Ratio\" page). It is thus omitted before applying the linear model.</b>")),
               downloadButton("DL_raw_heavy", "Download Raw CSV"),
               hr(),
               DT::dataTableOutput('raw_table_heavy'))
    )
  ),
  
  # LIGHT PAGE --------------------------------------------------
  tabPanel("Light",
           h2("Protein Degradation Analysis"),
           hr(),
           # CONSTRUCT FILTER PANEL
           wellPanel(
             h4("Filter Settings"),
             
             p("Adjust the slider below to indicate the minimum number of valid values",
               "required for retention."),
             
             fluidRow(
               # Input for missing value filter: DMSO
               column(4, sliderInput("DMSO_valid_light", label = "DMSO",
                                     min = 0, max = 5, value = 5)
               ),
               # Input for missing value filter: M1071
               column(4, sliderInput("M1071_valid_light", label = "M1071", 
                                     min = 0, max = 4, value = 4)
               ),
               # Input for missing value filter: Rapa
               column(4, sliderInput("Rapa_valid_light", label = "Rapamycin",
                                     min = 0, max = 4, value = 4))
             ),
             
             # Input for EACH or AT LEAST ONE
             radioButtons("filter_type_light", label = "Filter Type",
                          choices = list("Keep observations for which ALL conditions above are met" = FALSE, 
                                         "Keep observations for which ANY condition above is met" = TRUE),
                          selected = FALSE),
             
             # Help text
             helpText(HTML("Note: All tabs below will automatically apply the filter. 
                           Low thresholds could interfere with heatmap construction."))
             ),
           
           # CONSTRUCT PAGE TABS
           tabsetPanel(
             tabPanel("Summary",
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p("Summary statistics are displayed for each sample in the form of ",
                                 tags$b("X"), "/", tags$b("Y"), "(", tags$b("Z"), "%), where", tags$b("X"), 
                                 " is the number of proteins quantified,", tags$b("Y"), "is the number identified, and",
                                 tags$b("Z") , "is the expression in percentages.")
                        ),
                        column(9, br(), tableOutput("summary_table_light")))
             ),
             tabPanel("Pairwise Plot",
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p(HTML("Pairwise plots show three types of information.
                                      Histograms display the distribution of
                                      log<sub>2</sub>-transformed protein intensities along the diagonal.
                                      Additionally, scatter plots between any two 
                                      timepoints are shown at the upper right corner along with lines of best fit 
                                      in red. The corresponding correlation coefficients are found at the lower left 
                                      corner. Proteins quantified in one sample but not the other are counted as \"missing\".
                                      Data are organized by treatment groups.")),
                               hr(),
                               selectInput('pair_plot_light_display', 'Display Type:',
                                           c("NORM", "LOG2"), selectize = FALSE),
                               helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                        to the DMSO 0-hour timepoint.")
                               ),
                        column(9, br(),
                               h4("DMSO"),
                               plotOutput("pair_plot_light_DMSO", height = "700px", width = "700px"),
                               hr(),
                               h4("M1071"),
                               plotOutput("pair_plot_light_M1071", height = "700px", width = "700px"),
                               hr(),
                               h4("Rapamycin"),
                               plotOutput("pair_plot_light_Rapa", height = "700px", width = "700px")
                        ))
                               ),
             tabPanel("Heatmap",
                      fluidRow(
                        column(3, br(),
                               p(HTML("Heatmaps represent the values of the <b>filtered</b> data table in color.
                                      The samples are organized into columns, and each row shows the 
                                      log<sub>2</sub>-transformed intensities for a protein. The sample medians 
                                      are plotted as a line chart above the figure. See legend below for 
                                      color scale. Missing values are colored in black.")),
                               selectInput('heatmap_light_display', 'Display Type:',
                                           c("NORM", "LOG2"), selectize = FALSE),
                               helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                        to the DMSO 0-hour timepoint. NORM is recommended.")
                               ),
                        column(9, br(),
                               plotOutput("heatmap_light", height = "700px", width = "450px")
                        ))
                               ),
             tabPanel("Filtered Table", br(),
                      p(HTML("The data table below expresses the pulsed SILAC proteomics output as slopes 
                             and scores for ease of evaluation. A linear regression model is fitted on the
                             log<sub>2</sub>-transformed intensities for each condition. Here, timepoints are used
                             as the independent variable while the intensities as the dependent. At least two
                             data points are required for modeling. No slopes are reported for those that fail 
                             to meet this criterion. Scores are computed by subtracting the slopes of the 
                             DMSO control from those of the treatment group. A negative score suggests that the rate of 
                             degradation increases with treatment while a positive score indicates the opposite effect. 
                             <b>All values are extremely sensitive to outliers. Please interpret with caution.</b> 
                             A complete data table including the measured intensities can be downloaded below.")),
                      p(HTML("<b>NOTE: There is reason to believe that DMSO-0hr was mishandled or mislabeled (see 
                             heatmap under \"Ratio\" page). It is thus omitted before applying the linear model.</b>")),
                      downloadButton("DL_filtered_light", "Download Filtered CSV"),
                      hr(),
                      DT::dataTableOutput('filtered_table_light')),
             tabPanel("Raw Table", br(),
                      p(HTML("The data table below expresses the pulsed SILAC proteomics output as slopes 
                             and scores for ease of evaluation. A linear regression model is fitted on the
                             log<sub>2</sub>-transformed intensities for each condition. Here, timepoints are used
                             as the independent variable while the intensities as the dependent. At least two
                             data points are required for modeling. No slopes are reported for those that fail 
                             to meet this criterion. Scores are computed by subtracting the slopes of the 
                             DMSO control from those of the treatment groups. A negative score suggests that the rate of 
                             degradation increases with treatment while a positive score indicates the opposite effect. 
                             <b>All values are extremely sensitive to outliers. Please interpret with caution.</b> 
                             A complete data table including the measured intensities can be downloaded below.")),
                      p(HTML("<b>NOTE: There is reason to believe that DMSO-0hr was mishandled or mislabeled (see 
                             heatmap under \"Ratio\" page). It is thus omitted before applying the linear model.</b>")),
                      downloadButton("DL_raw_light", "Download Raw CSV"),
                      hr(),
                      DT::dataTableOutput('raw_table_light'))
                      )
             ),

  # RATIO PAGE --------------------------------------------------
  tabPanel("Ratio",
           h2("Protein Turnover Analysis"),
           hr(),
           # CONSTRUCT FILTER PANEL
           wellPanel(
             h4("Filter Settings"),
             
             p("Adjust the slider below to indicate the minimum number of valid values",
               "required for retention."),
             
             fluidRow(
               # Input for missing value filter: DMSO
               column(4, sliderInput("DMSO_valid_ratio", label = "DMSO",
                                     min = 0, max = 5, value = 5)
               ),
               # Input for missing value filter: M1071
               column(4, sliderInput("M1071_valid_ratio", label = "M1071", 
                                     min = 0, max = 4, value = 4)
               ),
               # Input for missing value filter: Rapa
               column(4, sliderInput("Rapa_valid_ratio", label = "Rapamycin",
                                     min = 0, max = 4, value = 4))
             ),
             
             # Input for EACH or AT LEAST ONE
             radioButtons("filter_type_ratio", label = "Filter Type",
                          choices = list("Keep observations for which ALL conditions above are met" = FALSE, 
                                         "Keep observations for which ANY condition above is met" = TRUE),
                          selected = FALSE),
             
             # Help text
             helpText(HTML("Note: All tabs below will automatically apply the filter. 
                           Low thresholds could interfere with heatmap construction."))
             ),
           
           # CONSTRUCT PAGE TABS
           tabsetPanel(
             tabPanel("Summary",
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p("Summary statistics are displayed for each sample in the form of ",
                                 tags$b("X"), "/", tags$b("Y"), "(", tags$b("Z"), "%), where", tags$b("X"), 
                                 " is the number of quantified H/L ratios,", tags$b("Y"), "is the total possible, and",
                                 tags$b("Z") , "is the expression in percentages.")
                        ),
                        column(9, br(), tableOutput("summary_table_ratio")))
             ),
             tabPanel("Pairwise Plot",
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p(HTML("Pairwise plots show three types of information.
                                      Histograms display the distribution of
                                      log<sub>2</sub>-transformed H/L ratios along the diagonal.
                                      Additionally, scatter plots between any two 
                                      timepoints are shown at the upper right corner along with lines of best fit 
                                      in red. The corresponding correlation coefficients are found at the lower left 
                                      corner. Ratios quantified in one sample but not the other are counted as \"missing\". 
                                      Data are organized by treatment groups.")),
                               hr(),
                               selectInput('pair_plot_ratio_display', 'Display Type:',
                                           c("NORM", "LOG2"), selectize = FALSE),
                               helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                        to the DMSO 0-hour timepoint.")
                               ),
                        column(9, br(),
                               h4("DMSO"),
                               plotOutput("pair_plot_ratio_DMSO", height = "700px", width = "700px"),
                               hr(),
                               h4("M1071"),
                               plotOutput("pair_plot_ratio_M1071", height = "700px", width = "700px"),
                               hr(),
                               h4("Rapamycin"),
                               plotOutput("pair_plot_ratio_Rapa", height = "700px", width = "700px")
                        ))
                               ),
             tabPanel("Heatmap",
                      fluidRow(
                        column(3, br(),
                               p(HTML("Heatmaps represent the values of the <b>filtered</b> data table in color.
                                      The samples are organized into columns, and each row shows the 
                                      log<sub>2</sub>-transformed H/L ratios for a protein. The sample medians 
                                      are plotted as a line chart above the figure. See legend below for 
                                      color scale. Missing values are colored in black.")),
                               selectInput('heatmap_ratio_display', 'Display Type:',
                                           c("NORM", "LOG2"), selectize = FALSE),
                               helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                        to the DMSO 0-hour timepoint. NORM is recommended.")
                               ),
                        column(9, br(),
                               plotOutput("heatmap_ratio", height = "700px", width = "450px")
                        ))
                               ),
             tabPanel("Filtered Table", br(),
                      p(HTML("The data table below expresses the pulsed SILAC proteomics output as slopes 
                             and scores for ease of evaluation. A linear regression model is fitted on the
                             log<sub>2</sub>-transformed H/L ratios for each condition. Here, timepoints are used
                             as the independent variable while the measured ratios as the dependent. At least two
                             data points are required for modeling. No slopes are reported for those that fail 
                             to meet this criterion. Scores are computed by subtracting the slopes of the 
                             DMSO control from those of the treatment groups. A negative score suggests that the rate of 
                             protein turnover is reduced with treatment while a positive score indicates the opposite effect. 
                             <b>All values are extremely sensitive to outliers. Please interpret with caution.</b> 
                             A complete data table including the ratios can be downloaded below.")),
                      p(HTML("<b>NOTE: There is reason to believe that DMSO-0hr was mishandled or mislabeled (see 
                             heatmap under \"Ratio\" page). It is thus omitted before applying the linear model.</b>")),
                      downloadButton("DL_filtered_ratio", "Download Filtered CSV"),
                      hr(),
                      DT::dataTableOutput('filtered_table_ratio')),
             tabPanel("Raw Table", br(),
                      p(HTML("The data table below expresses the pulsed SILAC proteomics output as slopes 
                             and scores for ease of evaluation. A linear regression model is fitted on the
                             log<sub>2</sub>-transformed H/L ratios for each condition. Here, timepoints are used
                             as the independent variable while the measured ratios as the dependent. At least two
                             data points are required for modeling. No slopes are reported for those that fail 
                             to meet this criterion. Scores are computed by subtracting the slopes of the 
                             DMSO control from those of the treatment groups. A negative score suggests that the rate of 
                             protein turnover is reduced with treatment while a positive score indicates the opposite effect. 
                             <b>All values are extremely sensitive to outliers. Please interpret with caution.</b> 
                             A complete data table including the ratios can be downloaded below.")),
                      p(HTML("<b>NOTE: There is reason to believe that DMSO-0hr was mishandled or mislabeled (see 
                             heatmap under \"Ratio\" page). It is thus omitted before applying the linear model.</b>")),
                      downloadButton("DL_raw_ratio", "Download Raw CSV"),
                      hr(),
                      DT::dataTableOutput('raw_table_ratio'))
                      )
           ),
  
  # LIGHT + HEAVY PAGE ------------------------------------------
  tabPanel("Light + Heavy",
           h2("SILAC Incorporation Analysis"),
           hr(),
           tabsetPanel(
             tabPanel("Full Incorporation",
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p(HTML("Histograms comparing the distribution of SILAC heavy- and 
                                      light-labeled proteins are shown here for two samples. 
                                      The <b>Heavy incorporation</b> sample was processed 
                                      after passaging GTML5 cells four times in heavy 
                                      SILAC media. <b>Light incorporation</b> refers to the DMSO 
                                      0-hr sample, which was conditioned in light SILAC media in 
                                      preparation for pulsed SILAC proteomics."))
                        ),
                        column(9, br(), plotOutput("HL_full_incorp", height = "350px", width = "800px")))
             ),
             tabPanel("Partial Incorporation",
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p(HTML("Validation of heavy label incorporation and light label reduction
                                      across timepoints is shown on the right. Histograms are organized 
                                      by conditions and timepoints.")),
                               hr(),
                               selectInput('HL_partial_label_display', 'Display Type:',
                                           c("NORM", "LOG2"), selectize = FALSE),
                               helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                 to the DMSO 0-hour timepoint.")
                        ),
                        column(9, br(),
                               h4("Compare SILAC Labeling Across Timepoints"),
                               plotOutput("HL_partial_label", height = "600px", width = "900px")
                        )
                      ),
                      hr(),
                      fluidRow(
                        column(3, br(),  # Vertical spacing
                               p(HTML("Validation that changes in protein ratio can be expressed as a  
                                      function with a greater weight on protein synthesis than degradation. 
                                      Histograms are organized by conditions and SILAC labels.")),
                               hr(),
                               selectInput('HL_partial_time_display', 'Display Type:',
                                           c("NORM", "LOG2"), selectize = FALSE),
                               helpText("Option to view the data prior to (LOG2) or following (NORM) normalization 
                                        to the DMSO 0-hour timepoint.")
                               ),
                        column(9, br(),
                               h4("Compare Timepoints Between SILAC Labels"),
                               plotOutput("HL_partial_time", height = "600px", width = "900px")
                        )
                      )
             )
          )
  ),
  
  # DOWNLOAD PAGE ------------------------------------------
  tabPanel("Download",
           h2("Download MaxQuant Output"),
           hr(),
           p(HTML("This interactive webpage is built on the <b>Protein Groups</b> output from MaxQuant.")),
           p(HTML("The data file is available below.")),
           downloadButton("DL_maxquant", "Download CSV")
  )
)



# SERVER ------------------------------------------------------------------
server = function(input, output) {
  
  # HEAVY PAGE --------------------------------------------------
  # HEAVY: Data tables after applying filter
  heavy_df = reactive({ filter_valids(heavyCDS, 
                                      c("DMSO", "M1071", "Rapa"), is_infinite = TRUE, 
                                      c(input$DMSO_valid_heavy, input$M1071_valid_heavy, input$Rapa_valid_heavy), 
                                      at_least_one = input$filter_type_heavy) 
    })
  heavy_table_raw = reactive({ 
    outTable = filter_valids(heavyCDS,
                             c("DMSO", "M1071", "Rapa"), is_infinite = TRUE,
                             c(input$DMSO_valid_heavy, input$M1071_valid_heavy, input$Rapa_valid_heavy),
                             at_least_one = input$filter_type_heavy)
    outTable = outTable[c("Gene.name", "Symbol", "slope_DMSO", "slope_M1071", "slope_Rapa", 
                          "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
    })
  heavy_table_filtered = reactive({ 
    outTable = filter_valids(heavyCDS,
                             c("DMSO", "M1071", "Rapa"), is_infinite = TRUE,
                             c(input$DMSO_valid_heavy, input$M1071_valid_heavy, input$Rapa_valid_heavy),
                             at_least_one = input$filter_type_heavy)
    outTable = outTable[outTable$KEEP, c("Gene.name", "Symbol", "slope_DMSO", "slope_M1071", "slope_Rapa",
                                         "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
    })
  
  # HEAVY: Output summary table statistics
  output$summary_table_heavy <- renderTable({ summary_stat(heavy_df()) })
  
  # HEAVY: Output pairwise plot
  output$pair_plot_heavy_DMSO <- renderPlot({ 
    if (input$pair_plot_heavy_display == "LOG2") {
      plot_pairs(heavy_df(), "^LOG2.*DMSO")
    } else {
      plot_pairs(heavy_df(), "^NORM.*DMSO")
    }
  })
  output$pair_plot_heavy_M1071 <- renderPlot({ 
    if (input$pair_plot_heavy_display == "LOG2") {
      plot_pairs(heavy_df(), "^LOG2.*M1071")
    } else {
      plot_pairs(heavy_df(), "^NORM.*M1071")
    }
  })
  output$pair_plot_heavy_Rapa <- renderPlot({ 
    if (input$pair_plot_heavy_display == "LOG2") {
      plot_pairs(heavy_df(), "^LOG2.*Rapa")
    } else {
      plot_pairs(heavy_df(), "^NORM.*Rapa")
    }
  })
  
  # HEAVY: Output heatmap
  output$heatmap_heavy <- renderPlot({ 
    if (input$heatmap_heavy_display == "LOG2") {
      plot_heatmap(heavy_df(), "^LOG2")
    } else {
      plot_heatmap(heavy_df(), "^NORM")    
    }
  })
  
  # HEAVY: Filtered Table
  output$filtered_table_heavy <- DT::renderDataTable(
    DT::datatable(
      heavy_table_filtered(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      )
    )
  )
  
  # HEAVY: Raw Table
  output$raw_table_heavy <- DT::renderDataTable(
    DT::datatable(
      heavy_table_raw(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      )
    )
  )
  
  # HEAVY: Download Filtered Table
  output$DL_filtered_heavy <- downloadHandler(
    filename = function() {"filtered_heavy.csv"},
    content = function(file) {
      write.csv(heavy_df()[heavy_df()$KEEP, -which(names(heavy_df()) == "KEEP")], file, row.names = FALSE)
    }   # Remove KEEP column
  )
  
  # HEAVY: Download Raw Table
  output$DL_raw_heavy <- downloadHandler(
    filename = function() {"raw_heavy.csv"},
    content = function(file) {
      write.csv(heavy_df()[-which(names(heavy_df()) == "KEEP")], file, row.names = FALSE)
    }   # Remove KEEP column
  )
  
  
  # LIGHT PAGE --------------------------------------------------
  # LIGHT: Data tables after applying filter
  light_df = reactive({ filter_valids(lightCDS, 
                                      c("DMSO", "M1071", "Rapa"), is_infinite = TRUE, 
                                      c(input$DMSO_valid_light, input$M1071_valid_light, input$Rapa_valid_light), 
                                      at_least_one = input$filter_type_light) 
  })
  light_table_raw = reactive({ 
    outTable = filter_valids(lightCDS,
                             c("DMSO", "M1071", "Rapa"), is_infinite = TRUE,
                             c(input$DMSO_valid_light, input$M1071_valid_light, input$Rapa_valid_light),
                             at_least_one = input$filter_type_light)
    outTable = outTable[c("Gene.name", "Symbol", "slope_DMSO", "slope_M1071", "slope_Rapa", 
                          "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  light_table_filtered = reactive({ 
    outTable = filter_valids(lightCDS,
                             c("DMSO", "M1071", "Rapa"), is_infinite = TRUE,
                             c(input$DMSO_valid_light, input$M1071_valid_light, input$Rapa_valid_light),
                             at_least_one = input$filter_type_light)
    outTable = outTable[outTable$KEEP, c("Gene.name", "Symbol", "slope_DMSO", "slope_M1071", "slope_Rapa",
                                         "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  
  # LIGHT: Output summary table statistics
  output$summary_table_light <- renderTable({ summary_stat(light_df()) })
  
  # LIGHT: Output pairwise plot
  output$pair_plot_light_DMSO <- renderPlot({ 
    if (input$pair_plot_light_display == "LOG2") {
      plot_pairs(light_df(), "^LOG2.*DMSO")
    } else {
      plot_pairs(light_df(), "^NORM.*DMSO")
    }
  })
  output$pair_plot_light_M1071 <- renderPlot({ 
    if (input$pair_plot_light_display == "LOG2") {
      plot_pairs(light_df(), "^LOG2.*M1071")
    } else {
      plot_pairs(light_df(), "^NORM.*M1071")
    }
  })
  output$pair_plot_light_Rapa <- renderPlot({ 
    if (input$pair_plot_light_display == "LOG2") {
      plot_pairs(light_df(), "^LOG2.*Rapa")
    } else {
      plot_pairs(light_df(), "^NORM.*Rapa")
    }
  })
  
  # LIGHT: Output heatmap
  output$heatmap_light <- renderPlot({ 
    if (input$heatmap_light_display == "LOG2") {
      plot_heatmap(light_df(), "^LOG2")
    } else {
      plot_heatmap(light_df(), "^NORM")    
    }
  })
  
  # LIGHT: Filtered Table
  output$filtered_table_light <- DT::renderDataTable(
    DT::datatable(
      light_table_filtered(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      )
    )
  )
  
  # LIGHT: Raw Table
  output$raw_table_light <- DT::renderDataTable(
    DT::datatable(
      light_table_raw(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      )
    )
  )
  
  # LIGHT: Download Filtered Table
  output$DL_filtered_light <- downloadHandler(
    filename = function() {"filtered_light.csv"},
    content = function(file) {
      write.csv(light_df()[light_df()$KEEP, -which(names(light_df()) == "KEEP")], file, row.names = FALSE)
    }   # Remove KEEP column
  )
  
  # LIGHT: Download Raw Table
  output$DL_raw_light <- downloadHandler(
    filename = function() {"raw_light.csv"},
    content = function(file) {
      write.csv(light_df()[-which(names(light_df()) == "KEEP")], file, row.names = FALSE)
    }   # Remove KEEP column
  )

  
  # RATIO PAGE --------------------------------------------------
  # RATIO: Data tables after applying filter
  ratio_df = reactive({ filter_valids(ratioDS, 
                                      c("DMSO", "M1071", "Rapa"), is_infinite = FALSE, 
                                      c(input$DMSO_valid_ratio, input$M1071_valid_ratio, input$Rapa_valid_ratio), 
                                      at_least_one = input$filter_type_ratio) 
  })
  ratio_table_raw = reactive({ 
    outTable = filter_valids(ratioDS,
                             c("DMSO", "M1071", "Rapa"), is_infinite = FALSE,
                             c(input$DMSO_valid_ratio, input$M1071_valid_ratio, input$Rapa_valid_ratio),
                             at_least_one = input$filter_type_ratio)
    outTable = outTable[c("Gene.name", "Symbol", "slope_DMSO", "slope_M1071", "slope_Rapa", 
                          "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  ratio_table_filtered = reactive({ 
    outTable = filter_valids(ratioDS,
                             c("DMSO", "M1071", "Rapa"), is_infinite = FALSE,
                             c(input$DMSO_valid_ratio, input$M1071_valid_ratio, input$Rapa_valid_ratio),
                             at_least_one = input$filter_type_ratio)
    outTable = outTable[outTable$KEEP, c("Gene.name", "Symbol", "slope_DMSO", "slope_M1071", "slope_Rapa",
                                         "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  
  # RATIO: Output summary table statistics
  output$summary_table_ratio <- renderTable({ summary_stat(ratio_df()) })
  
  # RATIO: Output pairwise plot
  output$pair_plot_ratio_DMSO <- renderPlot({ 
    if (input$pair_plot_ratio_display == "LOG2") {
      plot_pairs(ratio_df(), "^LOG2.*DMSO")
    } else {
      plot_pairs(ratio_df(), "^NORM.*DMSO")
    }
  })
  output$pair_plot_ratio_M1071 <- renderPlot({ 
    if (input$pair_plot_ratio_display == "LOG2") {
      plot_pairs(ratio_df(), "^LOG2.*M1071")
    } else {
      plot_pairs(ratio_df(), "^NORM.*M1071")
    }
  })
  output$pair_plot_ratio_Rapa <- renderPlot({ 
    if (input$pair_plot_ratio_display == "LOG2") {
      plot_pairs(ratio_df(), "^LOG2.*Rapa")
    } else {
      plot_pairs(ratio_df(), "^NORM.*Rapa")
    }
  })
  
  # RATIO: Output heatmap
  output$heatmap_ratio <- renderPlot({ 
    if (input$heatmap_ratio_display == "LOG2") {
      plot_heatmap(ratio_df(), "^LOG2")
    } else {
      plot_heatmap(ratio_df(), "^NORM")    
    }
  })
  
  # RATIO: Filtered Table
  output$filtered_table_ratio <- DT::renderDataTable(
    DT::datatable(
      ratio_table_filtered(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      )
    )
  )
  
  # RATIO: Raw Table
  output$raw_table_ratio <- DT::renderDataTable(
    DT::datatable(
      ratio_table_raw(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      )
    )
  )
  
  # RATIO: Download Filtered Table
  output$DL_filtered_ratio <- downloadHandler(
    filename = function() {"filtered_ratio.csv"},
    content = function(file) {
      write.csv(ratio_df()[ratio_df()$KEEP, -which(names(ratio_df()) == "KEEP")], file, row.names = FALSE)
    }   # Remove KEEP column
  )
  
  # RATIO: Download Raw Table
  output$DL_raw_ratio <- downloadHandler(
    filename = function() {"raw_ratio.csv"},
    content = function(file) {
      write.csv(ratio_df()[-which(names(ratio_df()) == "KEEP")], file, row.names = FALSE)
    }   # Remove KEEP column
  )
  
  
  # HEAVY + LIGHT PAGE --------------------------------------------------
  # HEAVY + LIGHT: Full incorporation
  output$HL_full_incorp <- renderPlot({ plot_hist2(cbind(lightCDS, heavyCDS, P4heavyC)) })
  
  # HEAVY + LIGHT: Partial incorporation (By Labeling)
  output$HL_partial_label = renderPlot({
    if (input$HL_partial_label_display == "LOG2") {
      plot_hist(cbind(lightCDS, heavyCDS), compare_time = FALSE,
                "^LOG2\\.L.*(3hr|6hr|12hr|24hr)", "^LOG2\\.H.*(3hr|6hr|12hr|24hr)")
    } else {
      plot_hist(cbind(lightCDS, heavyCDS), compare_time = FALSE,
                "^NORM\\.L.*(3hr|6hr|12hr|24hr)", "^NORM\\.H.*(3hr|6hr|12hr|24hr)")
    }
  })
  
  # HEAVY + LIGHT: Partial incorporation (By Time)
  output$HL_partial_time = renderPlot({
    if (input$HL_partial_time_display == "LOG2") {
      plot_hist(cbind(lightCDS, heavyCDS), compare_time = TRUE,
                "^LOG2\\.L.*(3hr|6hr|12hr|24hr)", "^LOG2\\.H.*(3hr|6hr|12hr|24hr)")
    } else {
      plot_hist(cbind(lightCDS, heavyCDS), compare_time = TRUE,
                "^NORM\\.L.*(3hr|6hr|12hr|24hr)", "^NORM\\.H.*(3hr|6hr|12hr|24hr)")
    }
  })
  
  
  # DOWNLOAD PAGE --------------------------------------------------
  # DOWNLOAD: Link
  output$DL_maxquant <- downloadHandler(
    filename = function() {"proteinGroups_maxquant.csv"},
    content = function(file) {
      write.csv(raw.data, file, row.names = FALSE)
    } 
  )
  
}

shinyApp(ui = ui, server = server)


#############################################
# Step 9 - Deploy App!
#############################################
# Copy and run code in console
#library(rsconnect)
#setwd('C:/Users/Tony Lin/Desktop/Wiita_lab/Projects/Proteomics_project/Pulsed_silac_Ozlem/Proteomics/analysis/deployDirectory')
#deployApp()
