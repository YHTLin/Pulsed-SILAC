### ShinyApp for pulsed SILAC data on GTML5 cells
## Adapted app.R in deployDirectory-v1 for LFQ!
## Some functions from pulsed_SILAC_ozlem.R has been reordered for app development

#############################################
# Step 1 - Set working directory
#############################################
#setwd('C:/Users/Tony Lin/Desktop/Wiita_lab/Projects/Proteomics_project/Pulsed_silac_Ozlem/Proteomics/analysis/LFQ-newDMSO')


#############################################
# Step 2 - Read in proteinGroups file
#############################################
raw.data = read.delim("./Data/proteinGroups_newDMSO.txt", header = TRUE, stringsAsFactors = FALSE)


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
data$Protein.ID = sapply(first_fasta, function(x)
  str_extract(x[1], "(?<=\\|.{1,10}\\|).*(?=_MOUSE)"))

# Assign gene names according to Fasta headers
data$Gene.name = sapply(first_fasta, function(x)
  str_extract(x[1], "(?<=\\|.{1,10}\\|.{1,30}_MOUSE).*(?=OS\\=)"))

# Compute LOG2 on intensity
LFQ.names = names(dplyr::select(data, starts_with("LFQ.intensity")))
log2LFQ.names = sub("^LFQ.intensity", "LOG2", LFQ.names)
data[log2LFQ.names] = log2(data[LFQ.names])

# Compute LOG2 on Ratio.H.L
ratio.names = grep("^Ratio\\.H\\.L\\.(DMSO|M1071|Rapa|Heavy_incorp)", names(data), value = T)
log2Ratio.names = sub("^Ratio", "LOG2", ratio.names)
data[log2Ratio.names] = log2(data[ratio.names])


#############################################
# Step 5 - Data organization
#############################################
# Organize P4heavy data (heavy incorporation after passage 4)
protein.info = c("Protein.IDs", "Gene.names", "Protein.ID", "Gene.name",
                 "Number.of.proteins", "Peptides")
P4heavy.names = c("LOG2.L.Heavy_incorp", "LOG2.H.Heavy_incorp", 
                  "LOG2.H.L.Heavy_incorp", "Ratio.H.L.count.Heavy_incorp")
P4heavy = dplyr::select(data, c(protein.info, P4heavy.names))

# Organize light data
require(gtools)   #for invoking mixedsort to order column labels
light.names = mixedsort(grep("^LOG2\\.L\\.(DMSO|M1071|Rapa)", 
                             names(data), value = TRUE))
light = dplyr::select(data, c(protein.info, light.names))

# Organize heavy data
heavy.names = mixedsort(grep("^LOG2\\.H\\.(DMSO|M1071|Rapa)", 
                             names(data), value = TRUE))
heavy = dplyr::select(data, c(protein.info, heavy.names))

# Organize H/L ratio data
HL.names = c(mixedsort(grep("^LOG2\\.H\\.L\\.(DMSO|M1071|Rapa)", names(data), value = TRUE)),
             mixedsort(grep("^Ratio\\.H\\.L\\.count\\.(DMSO|M1071|Rapa)", names(data), value = TRUE)))
ratio = dplyr::select(data, c(protein.info, HL.names))
rm(protein.info, light.names, heavy.names, HL.names)


#############################################
# Step 6 - Score treatment/control rate of change in intensity over time
#############################################
# Use the LOG2 or NORM columns to calculate slope of linear regression
score = function(df, sample.names, suffix) {
  require(gtools)   # For mixedsort
  # df = data frame containing LOG2 columns for computing differences in rate of change
  # sample.names = vector of regex patterns for extracting column names for each condition
  # suffix = character vector describing the name of each condition for output
  
  # Organize sample names by treatment
  log.names = lapply(sample.names, function(x) grep(x, names(df), value = TRUE, perl = TRUE)) #perl = TRUE enables negative lookahead
  
  # Create time vs intensity matrix for each protein in every condition
  df2 = lapply(log.names, function(x) {   # Loop through each condition
    x = mixedsort(x)  # order the columns alphanumerically
    apply(df[x], 1, function(y) data.frame(time = 1:length(y), value = y))  # Loop through each row
  })
  
  # Apply linear model after excluding -Inf or NA rows
  slopes = lapply(df2, function(x) {   # Loop through each condition
    sapply(x, function(y) {
      y = filter(y, is.finite(value))   # Remove -Inf or NA rows
      if (nrow(y) > 1) {   # Only perform linear regression if two or more data points are available
        x = lm(value ~ time, y)
        c(x$coefficients[2], summary(x)$r.squared) # Extract slope and r-squared from linear model
      } else {   # Otherwise return NA
        c(NA, NA)
      }
    })
  })
  slopes = lapply(slopes, t)
  
  # Output columns
  for (i in 1:length(suffix)) {
    df[paste0(c("slope_", "rSquared_"), suffix[i])] = slopes[[i]]
  }
  
  return(df)
}
# score_ is the difference in slope between treatment and DMSO control
# IGNORE DMSO_0HR IN LINEAR REGRESSION BECAUSE HEATMAP SUGGESTS THAT THE SAMPLE COULD BE MISLABLED
lightS = score(light, c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
lightS = mutate(lightS, score_M1071 = slope_M1071 - slope_DMSO)
lightS = mutate(lightS, score_Rapa = slope_Rapa - slope_DMSO)
heavyS = score(heavy, c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
heavyS = mutate(heavyS, score_M1071 = slope_M1071 - slope_DMSO)
heavyS = mutate(heavyS, score_Rapa = slope_Rapa - slope_DMSO)
ratioS = score(ratio, c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), c("DMSO", "M1071", "Rapa"))
ratioS = mutate(ratioS, score_M1071 = slope_M1071 - slope_DMSO)
ratioS = mutate(ratioS, score_Rapa = slope_Rapa - slope_DMSO)


#############################################
# Step 8 - Data Filtering Function
#############################################
# Filter out missing values
filter_valids = function(df, conditions, min_count, at_least_one = FALSE,
                         special_cond = NA, special_filter = NA) {
  # df = data frame containing LOG2 data for filtering and organized by data type
  # conditions = a character vector dictating the grouping
  # min_count = a numeric vector of the same length as "conditions" indicating the minimum 
  #     number of valid values for each condition for retention
  # at_least_one = TRUE means to keep the row if min_count is met for at least one condition
  #     FALSE means min_count must be met across all conditions for retention
  # special_cond = a character vector of regex patterns for isolating certain samples/conditions
  # special_filter = a logical vector of same length as special_cond indicating whether to
  #     require quantification in the conditions listed in special_cond
  require(dplyr)
  
  log2.names = grep("^LOG2", names(df), value = TRUE)   # Extract LOG2 column names
  if (!is.na(special_cond[1])) {
    # Extract names of samples requiring their own filtering
    special.names = unlist(sapply(special_cond, 
                                  function(x) grep(x, log2.names, value = TRUE)))
    log2.names = log2.names[!(log2.names %in% special.names)]
  }
  cond.names = lapply(conditions, # Group column names by conditions
                      function(x) grep(x, log2.names, value = TRUE))
  
  # Apply filter to "conditions"
  cond.filter = sapply(1:length(cond.names), function(i) {
    df2 = df[cond.names[[i]]]   # Extract columns of interest
    df2 = as.matrix(df2)   # Cast as matrix for the following command
    sums = rowSums(is.finite(df2)) # count the number of valid values for each condition
    sums >= min_count[i]   # Calculates whether min_count requirement is met
  })
  
  if (at_least_one) {
    KEEP = apply(cond.filter, 1, any)
  } else {
    KEEP = apply(cond.filter, 1, all)
  }
  
  # Apply filter to "special"
  if (!is.na(special_cond[1])) {
    special.names = special.names[special_filter]
    special_df = df[special.names]
    special_df = as.matrix(special_df)
    KEEP = cbind(KEEP, apply(special_df, 1, 
                             function(x) all(is.finite(x))))
    KEEP = apply(KEEP, 1, all)
  }
  
  df$KEEP = KEEP
    
  return(df)  # No rows are omitted, filter rules recorded in the KEEP column
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
  plot.df[!is.finite(plot.df)] = NA   # Required for superheat
  colnames(plot.df) = sub(paste0(pattern, "\\."), "", colnames(plot.df))   # name columns
  rownames(plot.df) = df[df$KEEP, "Protein.ID"]   # name rows
  
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



### Table of summary statistics
summary_stat = function(df) {
  # df = data frame containing LOG2 intensity
  
  require(dplyr)
  df2 = dplyr::select(df, starts_with("LOG2"))
  data.frame(
    Sample = sub("(LOG2.L.|LOG2.H.|LOG2.H.L.)", "", names(df2)),
    Raw = sapply(df2, function(x) {
      paste0(sum(is.finite(x)), " / ", length(x), 
             " (", round(sum(is.finite(x)) / length(x) * 100, 1), "%)")
    }),
    Filtered = sapply(df2[df$KEEP,], function(x) {
      paste0(sum(is.finite(x)), " / ", length(x), 
             " (", round(sum(is.finite(x)) / length(x) * 100, 1), "%)")
    })
  )
}


### Line Plot of LOG2 intensity vs Time
##Derived from function "score"
line_plot = function(df, ID, sample.names, label) {
  # df = data frame containing LOG2 intensity
  # ID = a numeric vector of length one indicating which row of data to use for plotting
  # sample.names = vector of regex patterns for extracting column names for each condition
  # label = character vector describing the name of each condition for output (same length as sample.names)
  
  require(ggplot2)
  log.names = lapply(sample.names, function(x) grep(x, names(df), value = TRUE, perl = TRUE)) #perl = TRUE enables negative lookahead
  
  df2 = lapply(1:length(log.names), function(i) {   # Loop through each condition
    x = mixedsort(log.names[[i]])  # order the columns alphanumerically
    data.frame(time = 1:length(df[ID, x]), value = as.numeric(df[ID, x]), lab = label[i])
  })
  df2 = Reduce(rbind, df2)
  df2 = df2[!is.infinite(df2$value) & !is.na(df2$value), ] # Remove -Inf or NA rows

  if (nrow(df2) == 0) {
    return("No valid values available for plotting!")
  } else {
    ggplot(df2, aes(x = time, y = value, color = lab, alpha = 0.3)) +
      geom_line(size = 1.5) +
      geom_point(size = 3) +
      ylab(expression(log[2]-transformed~values)) +
      scale_x_continuous("Time",
                         breaks = c(1,2,3,4),
                         labels = c("3H", "6H", "12H", "24H"),
                         limits = c(1,4)) +
      scale_color_discrete(name = "Legend") +
      guides(alpha = FALSE)
  }
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
                      labeling with amino acids in cell culture (SILAC). This webpage showcases
                      the proteomic data using interactive visualization.")),
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
               p(HTML("Next, the protein intensity values were log<sub>2</sub>-transformed.
                      These values are reported in the <b>LOG2</b> columns in the downloadable CSVs.")),
               p(HTML("Finally, a filter to remove missing values can be applied independently to 
                      each of the analysis pages above. The controls are initialized at the most selective 
                      threshold and can be adjusted by the user. Note that some graphics may 
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
                              min = 0, max = 4, value = 4)
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
      
      # Input checkbox for DMSO-0hr quantification
      checkboxInput("DMSO_quant_heavy", label = "Require quantification in DMSO-0hr", 
                    value = FALSE),
        
      # Help text
      helpText(HTML("Note: All tabs below will automatically apply the filter. 
                    Low thresholds could interfere with visualization."))
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
                               Data are organized by treatment groups."))
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
                               color scale. Missing values are colored in black."))
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
                      <b>All values are extremely sensitive to outliers. Please interpret with caution and
                      select rows to check for linearity.</b> 
                      A complete data table including the measured intensities can be downloaded below.")),
               p(HTML("<b>NOTE: DMSO-0hr is omitted when applying the linear model.</b>")),
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
                      <b>All values are extremely sensitive to outliers. Please interpret with caution and
                      select rows to check for linearity.</b>  
                      A complete data table including the measured intensities can be downloaded below.")),
               p(HTML("<b>NOTE: DMSO-0hr is omitted when applying the linear model.</b>")),
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
                                     min = 0, max = 4, value = 4)
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
             
             # Input checkbox for DMSO-0hr quantification
             checkboxInput("DMSO_quant_light", label = "Require quantification in DMSO-0hr", 
                           value = FALSE),
             
             # Help text
             helpText(HTML("Note: All tabs below will automatically apply the filter. 
                           Low thresholds could interfere with visualization."))
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
                                      Data are organized by treatment groups."))
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
                                      color scale. Missing values are colored in black."))
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
                             <b>All values are extremely sensitive to outliers. Please interpret with caution and
                             select rows to check for linearity.</b>  
                             A complete data table including the measured intensities can be downloaded below.")),
                      p(HTML("<b>NOTE: DMSO-0hr is omitted when applying the linear model.</b>")),
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
                             <b>All values are extremely sensitive to outliers. Please interpret with caution and
                             select rows to check for linearity.</b>  
                             A complete data table including the measured intensities can be downloaded below.")),
                      p(HTML("<b>NOTE: DMSO-0hr is omitted when applying the linear model.</b>")),
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
                                     min = 0, max = 4, value = 4)
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

             # Input checkbox for DMSO-0hr quantification
             checkboxInput("DMSO_quant_ratio", label = "Require quantification in DMSO-0hr", 
                           value = FALSE),
             
             # Help text
             helpText(HTML("Note: All tabs below will automatically apply the filter. 
                           Low thresholds could interfere with visualization."))
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
                                      Data are organized by treatment groups."))
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
                                      color scale. Missing values are colored in black."))
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
                             <b>All values are extremely sensitive to outliers. Please interpret with caution and
                             select rows to check for linearity.</b>  
                             A complete data table including the ratios can be downloaded below.")),
                      p(HTML("<b>NOTE: DMSO-0hr is omitted when applying the linear model.</b>")),
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
                             <b>All values are extremely sensitive to outliers. Please interpret with caution and
                             select rows to check for linearity.</b>  
                             A complete data table including the ratios can be downloaded below.")),
                      p(HTML("<b>NOTE: DMSO-0hr is omitted when applying the linear model.</b>")),
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
                                      by conditions and timepoints."))
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
                                      Histograms are organized by conditions and SILAC labels."))
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
           p(HTML("The raw data file is available below.")),
           downloadButton("DL_maxquant", "Download CSV")
  )
)



# SERVER ------------------------------------------------------------------
server = function(input, output, session) {
  
  # HEAVY PAGE --------------------------------------------------
  # HEAVY: Data tables after applying filter
  heavy_df = reactive({ filter_valids(heavyS, 
                                      c("DMSO", "M1071", "Rapa"),
                                      c(input$DMSO_valid_heavy, input$M1071_valid_heavy, input$Rapa_valid_heavy), 
                                      at_least_one = input$filter_type_heavy,
                                      special_cond = "_0hr",
                                      special_filter = input$DMSO_quant_heavy) 
    })
  heavy_table_raw = reactive({ 
    outTable = heavyS
    outTable = outTable[c("Gene.name", "Protein.ID", "slope_DMSO", "slope_M1071", "slope_Rapa", 
                          "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
    })
  heavy_table_filtered = reactive({ 
    outTable = filter_valids(heavyS,
                             c("DMSO", "M1071", "Rapa"),
                             c(input$DMSO_valid_heavy, input$M1071_valid_heavy, input$Rapa_valid_heavy),
                             at_least_one = input$filter_type_heavy,
                             special_cond = "_0hr",
                             special_filter = input$DMSO_quant_heavy)
    outTable = outTable[outTable$KEEP, c("Gene.name", "Protein.ID", "slope_DMSO", "slope_M1071", "slope_Rapa",
                                         "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
    })
  
  # HEAVY: Output summary table statistics
  output$summary_table_heavy <- renderTable({ summary_stat(heavy_df()) })
  
  # HEAVY: Output pairwise plot
  output$pair_plot_heavy_DMSO <- renderPlot({
    plot_pairs(heavy_df(), "^LOG2.*DMSO", use_keep = TRUE)
  })
  output$pair_plot_heavy_M1071 <- renderPlot({ 
    plot_pairs(heavy_df(), "^LOG2.*M1071", use_keep = TRUE)
  })
  output$pair_plot_heavy_Rapa <- renderPlot({ 
    plot_pairs(heavy_df(), "^LOG2.*Rapa", use_keep = TRUE)
  })
  
  # HEAVY: Output heatmap
  output$heatmap_heavy <- renderPlot({ 
    plot_heatmap(heavy_df(), "^LOG2")
  })
  
  # HEAVY: Filtered Table
  output$filtered_table_heavy <- DT::renderDataTable(
    DT::datatable(
      heavy_table_filtered(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # HEAVY: Raw Table
  output$raw_table_heavy <- DT::renderDataTable(
    DT::datatable(
      heavy_table_raw(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )

  # HEAVY: Filtered Table Line Plot (in modal dialog)
  ID_filtered_heavy = reactive({ rownames(heavy_table_filtered()[input$filtered_table_heavy_rows_selected, ]) })
  output$line_plot_heavy_filtered = renderPlot({
    line_plot(heavy_df(), ID_filtered_heavy(),
              c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"),
              c("DMSO", "M1071", "Rapa"))
  })
  observeEvent(input$filtered_table_heavy_rows_selected,
               {
                 showModal(modalDialog(
                   title = heavy_df()[ID_filtered_heavy(), "Gene.name"],
                   plotOutput("line_plot_heavy_filtered"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
  # HEAVY: Raw Table Line Plot (in modal dialog)
  ID_raw_heavy = reactive({ rownames(heavy_table_raw()[input$raw_table_heavy_rows_selected, ]) })
  output$line_plot_heavy_raw = renderPlot({
    line_plot(heavy_df(), ID_raw_heavy(), 
              c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), 
              c("DMSO", "M1071", "Rapa"))
  })
  observeEvent(input$raw_table_heavy_rows_selected,
               {
                 showModal(modalDialog(
                   title = heavy_df()[ID_raw_heavy(), "Gene.name"],
                   plotOutput("line_plot_heavy_raw"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
  # HEAVY: Download Filtered Table
  output$DL_filtered_heavy <- downloadHandler(
    filename = function() {"filtered_heavy.csv"},
    content = function(file) {
      write.csv(heavy_df()[heavy_df()$KEEP, -which(names(heavy_df()) == "KEEP")], 
                file, row.names = FALSE)
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
  light_df = reactive({ filter_valids(lightS, 
                                      c("DMSO", "M1071", "Rapa"),
                                      c(input$DMSO_valid_light, input$M1071_valid_light, input$Rapa_valid_light), 
                                      at_least_one = input$filter_type_light,
                                      special_cond = "_0hr",
                                      special_filter = input$DMSO_quant_light) 
  })
  light_table_raw = reactive({ 
    outTable = lightS
    outTable = outTable[c("Gene.name", "Protein.ID", "slope_DMSO", "slope_M1071", "slope_Rapa", 
                          "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  light_table_filtered = reactive({ 
    outTable = filter_valids(lightS,
                             c("DMSO", "M1071", "Rapa"),
                             c(input$DMSO_valid_light, input$M1071_valid_light, input$Rapa_valid_light),
                             at_least_one = input$filter_type_light,
                             special_cond = "_0hr",
                             special_filter = input$DMSO_quant_light)
    outTable = outTable[outTable$KEEP, c("Gene.name", "Protein.ID", "slope_DMSO", "slope_M1071", "slope_Rapa",
                                         "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  
  # LIGHT: Output summary table statistics
  output$summary_table_light <- renderTable({ summary_stat(light_df()) })
  
  # LIGHT: Output pairwise plot
  output$pair_plot_light_DMSO <- renderPlot({ 
    plot_pairs(light_df(), "^LOG2.*DMSO", use_keep = TRUE)
  })
  output$pair_plot_light_M1071 <- renderPlot({ 
    plot_pairs(light_df(), "^LOG2.*M1071", use_keep = TRUE)
  })
  output$pair_plot_light_Rapa <- renderPlot({ 
    plot_pairs(light_df(), "^LOG2.*Rapa", use_keep = TRUE)
  })
  
  # LIGHT: Output heatmap
  output$heatmap_light <- renderPlot({ 
    plot_heatmap(light_df(), "^LOG2")
  })
  
  # LIGHT: Filtered Table
  output$filtered_table_light <- DT::renderDataTable(
    DT::datatable(
      light_table_filtered(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # LIGHT: Raw Table
  output$raw_table_light <- DT::renderDataTable(
    DT::datatable(
      light_table_raw(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # LIGHT: Filtered Table Line Plot (in modal dialog)
  ID_filtered_light = reactive({ rownames(light_table_filtered()[input$filtered_table_light_rows_selected, ]) })
  output$line_plot_light_filtered = renderPlot({
    line_plot(light_df(), ID_filtered_light(),
              c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"),
              c("DMSO", "M1071", "Rapa"))
  })
  observeEvent(input$filtered_table_light_rows_selected,
               {
                 showModal(modalDialog(
                   title = light_df()[ID_filtered_light(), "Gene.name"],
                   plotOutput("line_plot_light_filtered"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
  # LIGHT: Raw Table Line Plot (in modal dialog)
  ID_raw_light = reactive({ rownames(light_table_raw()[input$raw_table_light_rows_selected, ]) })
  output$line_plot_light_raw = renderPlot({
    line_plot(light_df(), ID_raw_light(), 
              c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), 
              c("DMSO", "M1071", "Rapa"))
  })
  observeEvent(input$raw_table_light_rows_selected,
               {
                 showModal(modalDialog(
                   title = light_df()[ID_raw_light(), "Gene.name"],
                   plotOutput("line_plot_light_raw"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
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
  ratio_df = reactive({ filter_valids(ratioS, 
                                      c("DMSO", "M1071", "Rapa"),
                                      c(input$DMSO_valid_ratio, input$M1071_valid_ratio, input$Rapa_valid_ratio), 
                                      at_least_one = input$filter_type_ratio,
                                      special_cond = "_0hr",
                                      special_filter = input$DMSO_quant_ratio) 
  })
  ratio_table_raw = reactive({ 
    outTable = ratioS
    outTable = outTable[c("Gene.name", "Protein.ID", "slope_DMSO", "slope_M1071", "slope_Rapa", 
                          "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  ratio_table_filtered = reactive({ 
    outTable = filter_valids(ratioS,
                             c("DMSO", "M1071", "Rapa"),
                             c(input$DMSO_valid_ratio, input$M1071_valid_ratio, input$Rapa_valid_ratio),
                             at_least_one = input$filter_type_ratio,
                             special_cond = "_0hr",
                             special_filter = input$DMSO_quant_ratio)
    outTable = outTable[outTable$KEEP, c("Gene.name", "Protein.ID", "slope_DMSO", "slope_M1071", "slope_Rapa",
                                         "score_M1071", "score_Rapa")]
    outTable[3:7] = round(as.matrix(outTable[3:7]), 3)
    return(outTable)
  })
  
  # RATIO: Output summary table statistics
  output$summary_table_ratio <- renderTable({ summary_stat(ratio_df()) })
  
  # RATIO: Output pairwise plot
  output$pair_plot_ratio_DMSO <- renderPlot({ 
    plot_pairs(ratio_df(), "^LOG2.*DMSO", use_keep = TRUE)
  })
  output$pair_plot_ratio_M1071 <- renderPlot({ 
    plot_pairs(ratio_df(), "^LOG2.*M1071", use_keep = TRUE)
  })
  output$pair_plot_ratio_Rapa <- renderPlot({ 
    plot_pairs(ratio_df(), "^LOG2.*Rapa", use_keep = TRUE)
  })
  
  # RATIO: Output heatmap
  output$heatmap_ratio <- renderPlot({ 
    plot_heatmap(ratio_df(), "^LOG2")
  })
  
  # RATIO: Filtered Table
  output$filtered_table_ratio <- DT::renderDataTable(
    DT::datatable(
      ratio_table_filtered(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # RATIO: Raw Table
  output$raw_table_ratio <- DT::renderDataTable(
    DT::datatable(
      ratio_table_raw(), 
      options = list(
        lengthMenu = list(c(10, 20, 50, 100, -1), c("10", "20", "50", "100", "ALL")),
        pageLength = 10
      ),
      selection = "single",
      colnames = c('ID' = 1)
    )
  )
  
  # RATIO: Filtered Table Line Plot (in modal dialog)
  ID_filtered_ratio = reactive({ rownames(ratio_table_filtered()[input$filtered_table_ratio_rows_selected, ]) })
  output$line_plot_ratio_filtered = renderPlot({
    line_plot(ratio_df(), ID_filtered_ratio(),
              c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"),
              c("DMSO", "M1071", "Rapa"))
  })
  observeEvent(input$filtered_table_ratio_rows_selected,
               {
                 showModal(modalDialog(
                   title = ratio_df()[ID_filtered_ratio(), "Gene.name"],
                   plotOutput("line_plot_ratio_filtered"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
  # RATIO: Raw Table Line Plot (in modal dialog)
  ID_raw_ratio = reactive({ rownames(ratio_table_raw()[input$raw_table_ratio_rows_selected, ]) })
  output$line_plot_ratio_raw = renderPlot({
    line_plot(ratio_df(), ID_raw_ratio(), 
              c("^LOG2.*DMSO(?!_0hr)", "^LOG2.*M1071", "^LOG2.*Rapa"), 
              c("DMSO", "M1071", "Rapa"))
  })
  observeEvent(input$raw_table_ratio_rows_selected,
               {
                 showModal(modalDialog(
                   title = ratio_df()[ID_raw_ratio(), "Gene.name"],
                   plotOutput("line_plot_ratio_raw"),
                   easyClose = TRUE,
                   footer = NULL
                 ))
               })
  
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
  output$HL_full_incorp <- renderPlot({ plot_hist2(cbind(lightS, heavyS, P4heavy)) })
  
  # HEAVY + LIGHT: Partial incorporation (By Labeling)
  output$HL_partial_label = renderPlot({
    plot_hist(cbind(lightS, heavyS), compare_time = FALSE,
                "^LOG2\\.L.*(3hr|6hr|12hr|24hr)", "^LOG2\\.H.*(3hr|6hr|12hr|24hr)")
  })
  
  # HEAVY + LIGHT: Partial incorporation (By Time)
  output$HL_partial_time = renderPlot({
    plot_hist(cbind(lightS, heavyS), compare_time = TRUE,
                "^LOG2\\.L.*(3hr|6hr|12hr|24hr)", "^LOG2\\.H.*(3hr|6hr|12hr|24hr)")
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
#setwd('C:/Users/Tony Lin/Desktop/Wiita_lab/Projects/Proteomics_project/Pulsed_silac_Ozlem/Proteomics/analysis/LFQ-newDMSO')
#deployApp()
