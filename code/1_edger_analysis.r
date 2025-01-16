# This code performs EdgeR analysis for pre-processed data.
## EdgeR is a form of differential expression analysis which uses a generalized linear model (GLM). 
## The GLM uses a negative binomial distribution to model lipid ion count data (or mRNA count data) accounting for the technical and biological variability. 

# References
citation("edgeR")

## Pre-processing is required to ensure there is no duplicate lipid entries in the input data and remove NaN/NA values. 

# This code generates output files that look like this: KO_vs_WT.csv
## Columns in this output file include:
### lipid: lipid name in this general format: Headgroup(sn1/sn2/sn3)
### type: lipid type (abbreviations can be found in the supplement materials: ________), 
### mean1 & mean2: mean ion intensity for the two groups (Group1: Pla2g2f cKO, Group2: WT)
### logFC: Log2 transformed fold change of Group1 (Pla2g2f cKO) over Group2 (WT)
### logCPM: Log2 transformed average expression of the lipid species across all samples
### LR: likelihood ratio test statistic
### PValue: nominal p-value derived from qualsi-likelihood ratio test without multiple testing correction
### FDR: false discovery rate calculated using the Benjamini–Hochberg method
### Z1_sum, Z2_sum...K3_sum: ion intensity values detected in each sample (Z1, Z2... K3)

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(readxl)
library(edgeR)
library(tidyverse)

# Function to get package version
get_package_version <- function(package_name) {
  package_version <- as.character(packageVersion(package_name))
  return(package_version)
}
packages <- c("dplyr", "readxl", "edgeR", "tidyverse")
package_versions <- sapply(packages, get_package_version)
print(package_versions)

# set WD
getwd()
setwd("")

# read in pre-processed data
file_path = "./data/" # make sure processed data are in this folder
file_list = list.files(path=file_path, pattern=NULL, all.files=FALSE, full.names=FALSE)
file_list

read_lipids <- function(file_path, filename) {
  data <- read_csv(paste0(file_path, filename), show_col_types = FALSE) %>%
    select(-transition) %>%  # Exclude the 'transition' column
    dplyr::mutate_if(is.numeric, ~replace(., is.nan(.), 0))  # Replace NaN values with 0
  return(data)
}

data_pmg <- read_lipids(file_path, "1_PMG_ion_intensity.csv")
data_mix <- read_lipids(file_path, "2_MIX_ion_intensity.csv")

dim(data_pmg) # 1497 lipids, 16 samples (N=4 for each group: Impact D1, Sham D1, Impact D7, Sham D7)
dim(data_mix) # 1497 lipids, 24 samples (N=5 for 1-Pre-Impact, N=5 for 2-Impact 4H, N=5 for 3-Impact D1, N=5 for 4-Impact D7, N=4 for 5-Impact D7+phAB)


# set experiments

expr1 <- "I1 vs S1"
expr2 <- "I7 vs S7"
expr3 <- "I7 vs I1"
expr_list_pmg <- c(expr1, expr2, expr3)

expr4 <- "Pre vs 4H Post-Impact"
expr5 <- "Pre vs 3D Post-Impact"
expr6 <- "Pre vs 7D Post-Impact"
expr7 <- "7D Post-Impact vs +pHAB"
expr_list_mix <- c(expr4, expr5, expr6, expr7)

expr_list <- c(expr_list_pmg, expr_list_mix)
expr_list

# create a list of data for looping through analysis
data_list <- list(data_pmg, data_pmg, data_pmg, data_mix, data_mix, data_mix, data_mix, data_mix)

# set EdgeR groups
## PMG experimental groups
blank_name = "Blank" # column name for blank

gr1_name = "I1" # group1: Impact Day 1 PMG
gr1_length = 4
gr1 = rep(gr1_name, gr1_length)
gr2_name = "S1" # group2: Sham Day 1 PMG
gr2_length = 4 
gr2 = rep(gr2_name, gr2_length)
gr3_name = "I7" # group3: Impact Day 7 PMG
gr3_length = 4
gr3 = rep(gr3_name, gr3_length)
gr4_name = "S7" # group4: Sham Day 7 PMG
gr4_length = 4
gr4 = rep(gr4_name, gr4_length)

## MIX experimental groups
blank_name = "SolventBlank"

gr5_name = "Pre_I" # group5: Pre-Impact MIX
gr5_length = 5
gr5 = rep(gr5_name, gr5_length)
gr6_name = "I_4H" # group6: Impact 4H MIX
gr6_length = 5
gr6 = rep(gr6_name, gr6_length)
gr7_name = "I_D3" # group7: Impact D3 MIX
gr7_length = 5
gr7 = rep(gr7_name, gr7_length)
gr8_name = "I_D7" # group8: Impact D7 MIX
gr8_length = 5
gr8 = rep(gr8_name, gr8_length)
gr9_name = "I_D7_pHAB" # group9: Impact D7+pHAB MIX
gr9_length = 4
gr9 = rep(gr9_name, gr9_length)
gr10_name = "I_D7_" # group10: Impact D7-pHAB MIX corresponding to the +pHAB group
gr10_length = 4
gr10 = rep(gr10_name, gr10_length)

## for experiment design 1, comparing Impact D1 vs Sham D1 PMG
## since N number is the same, this expr design can be applied to expr 2, 3, and 7
groups_expr1 = c(gr1, gr2, blank_name) %>% # create factor from groups
              factor(levels = c(blank_name, gr1_name, gr2_name)) 
design_expr1 = model.matrix(~groups_expr1) # create design matrix with experimental groups as column names
colnames(design_expr1) = c("Intercept", gr1_name, gr2_name) 
contrast_expr1 = makeContrasts( # create contrast
    H = I1 - S1, 
    levels = groups_expr1
)

## for experiment design 2, comparing Impact 4H vs Pre-Impact MIX
## since N number is the same, this expr design can be applied to expr 4-6
groups_expr2 = c(gr6, gr5, blank_name) %>% # create factor from groups
              factor(levels = c(blank_name, gr6_name, gr5_name)) 
design_expr2 = model.matrix(~groups_expr2) # create design matrix with experimental groups as column names
colnames(design_expr2) = c("Intercept", gr6_name, gr5_name) 
contrast_expr2 = makeContrasts( # create contrast
    H = I_4H - Pre_I, 
    levels = groups_expr2
) 

# functions folder
functions_filepath = "./code/functions/"
source(paste0(functions_filepath, "edger.r")) # load functions for EdgeR analysis

# perform EdgeR analysis

colnames(data_pmg)
lipid_cols = c(1,2) # where lipid names & types are in the raw data frame
counts_expr1 = c(3:6, 11:14, 19) # where ion intensities are in the raw data frame for exper 1
counts_expr2 = c(7:10, 15:18, 19) # where ion intensities are in the raw data frame for exper 2
counts_expr3 = c(7:10, 3:6, 19) # where ion intensities are in the raw data frame for exper 3

colnames(data_mix)
counts_expr4 = c(27, grep("2", colnames(data_mix)), grep("1", colnames(data_mix))) # Impact 4H vs Pre-Impact MIX
counts_expr5 = c(27, grep("3", colnames(data_mix)), grep("1", colnames(data_mix))) # Impact D3 vs Pre-Impact MIX
counts_expr6 = c(27, grep("4", colnames(data_mix)), grep("1", colnames(data_mix))) # Impact D7 vs Pre-Impact MIX
counts_expr7 = c(27, grep("5", colnames(data_mix)), grep("4", colnames(data_mix))[c(1, 3:5)]) # Impact D7+pHAb vs Impact D7

# using a loop since there are many comparisons
counts_expr_list <- list(counts_expr1, counts_expr2, counts_expr3, counts_expr4, counts_expr5, counts_expr6, counts_expr7) 
groups_expr_list <- list(groups_expr1, groups_expr1, groups_expr1, groups_expr2, groups_expr2, groups_expr2, groups_expr1) 
design_expr_list <- list(design_expr1, design_expr1, design_expr1, design_expr2, design_expr2, design_expr2, design_expr1) 
contrast_expr_list <- list(contrast_expr1, contrast_expr1, contrast_expr1, contrast_expr2, contrast_expr2, contrast_expr2, contrast_expr1)

sapply(counts_expr_list, length) == sapply(groups_expr_list, length) # check the length of count columns is equal to the length of groups

# Initialize an empty list to store the results
results_list <- list()

# Loop through each combination of counts_expr, groups_expr, and design_expr
for (i in 1:length(expr_list)) {
  
  # Extract the expressions for the current iteration
  data <- data_list[[i]]
  counts_expr <- counts_expr_list[[i]]  # Use [[i]] to extract the vector from the list
  groups_expr <- groups_expr_list[[i]]
  design_expr <- design_expr_list[[i]]
  contrast_expr <- contrast_expr_list[[i]]  

  # Perform the analysis and store the result in the list
  expr_results <- data %>% # use all lipids
    perform_raw_analysis(counts_expr, lipid_cols, groups_expr, design_expr) %>%
    calculate_significance(design_expr, contrast_expr) %>%
    topTags(5000) %>%
    as.data.frame() %>%
    remove_rownames()
      
  # Add the results to the list with a meaningful name
  results_list[[paste0("result_expr", i)]] <- expr_results
}

head(results_list[[1]]) # preview results
dim(results_list[[1]])
length(unique(results_list[[1]]$lipid))

# merge results with ion intensity data
merged_list <- list()

for (i in 1:length(expr_list)) {
  counts_expr <- counts_expr_list[[i]]
  result_df <- results_list[[i]] 
  data <- data_list[[i]]

  # Merge the result with the lipid name & types
  merged_df <- merge(result_df, data[,c(1, counts_expr)], by = "lipid") 
  merged_list[[paste0("result_expr", i)]] <- merged_df
}

head(merged_list[[1]])
dim(merged_list[[1]])
length(unique(merged_list[[1]]$lipid))

# output file path

out_filepath = "./edger_results/"
dir.create(path = out_filepath, F)

# add group mean to merged results 

## this function calculate means for 2 groups and add them to the dataframe
add_means <- function(data, g1, g2){ 
  g1_loc <- grep(g1, colnames(data), ignore.case = FALSE) 
  g2_loc <- grep(g2, colnames(data), ignore.case = FALSE)

  data <- data %>%
    mutate(mean1 = rowMeans(select(., all_of(g1_loc))),
           mean2 = rowMeans(select(., all_of(g2_loc)))) # calculate row means for each group
}

# output full result tables

for (i in 1:length(expr_list)) {
  df_name <- paste0("result_expr", i)
  data <- merged_list[[i]]

  if (i == 1) {result_df <- add_means(data, "I1", "S1")}
  else if (i == 2) {result_df <- add_means(data, "I7", "S7")}
  else if(i == 3) {result_df <- add_means(data, "I7", "I1")}
  else if(i == 4) {result_df <- add_means(data, "1", "2")}
  else if(i == 5) {result_df <- add_means(data, "1", "3")}
  else if(i == 6) {result_df <- add_means(data, "1", "4")}
  else if(i == 7) {result_df <- add_means(data, "4", "5")}

  assign(df_name, result_df)
  write.csv(result_df, file = paste0(out_filepath, i, " ", expr_list[[i]], '.csv'), row.names=FALSE)
}

# save FDR<0.1 lipids list

out_filepath = "./edger_results/FDR0.1 lipids/"
dir.create(path = out_filepath, F)

for (i in 1:length(expr_list)) {
  data <- merged_list[[i]]
  expr <- expr_list[[i]]

  FDR_lipids <- data %>%
    filter(FDR<0.1) 
  results_df <- table(FDR_lipids$type) %>% as.data.frame()
  colnames(results_df) <- c("Lipid Class", "Number of Species")

  write.csv(results_df, file = paste0(out_filepath, i, " ", expr_list[[i]], " FDR0.1 lipids", '.csv'), row.names = FALSE)
}

print("Mission Accomplished!")
