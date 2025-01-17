# This code filters analyzed data for significant lipids (FDR<0.1) and compare to LIPEA's lipid list
# The output of this code can be uploaded to LIPEA (https://hyperlipea.org/analyze) for pathway enrichment analysis

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(tidyverse)
library(readxl)

# Get package versions
get_package_version <- function(package_name) {
  package_version <- as.character(packageVersion(package_name))
  return(package_version)
}
packages <- c("dplyr", "tidyverse", "readxl")
package_versions <- sapply(packages, get_package_version)
print(package_versions)

# set WD
getwd()
setwd("")

# read in analyzed data
analyzed_filepath = "./edgeR_results/"
analyzed_file_list = list.files(path=analyzed_filepath, pattern=".csv", all.files=FALSE, full.names=FALSE)
expr_list <- str_sub(analyzed_file_list, 3, -5)
expr_list

data_list <- list()
for (i in 1:length(expr_list)) {
   data <- read_csv(paste0(analyzed_filepath, analyzed_file_list[[i]]), show_col_types = FALSE)
   data <- filter(data, FDR < 0.1)
   data_list[[expr_list[[i]]]] <- data
}
sapply(data_list, dim) # preview the number of significant lipids from each experiment

# read in LIPEA lipid class list
# this file was obtained from LIPEA website (https://hyperlipea.org/analyze) and placed in the data folder
lipea <- read_xls("data/LIPEA - Lipid classes.xls")
dim(lipea)

# functions folder
functions_filepath = "./code/functions/"
source(paste0(functions_filepath, "lipea naming.r")) # load function

# filter for lipid classes in LIPEA lipid class list
lipea_data <- list()
for (i in seq_along(expr_list)) {
    data <- data_list[[i]]
    lipea_data[[expr_list[[i]]]] <- lipea_change(data)
}
sapply(lipea_data, head) # preview data

# export LIPEA lipids list
out_filepath = "./lipea_analysis/"
dir.create(path = out_filepath, F)

for (i in seq_along(expr_list)) {
  lipea_matches <- lipea_data[[i]]
  expr <- expr_list[[i]]
  results_df <- as.data.frame(lipea_matches)
  print(results_df)
  write.csv(results_df, file = paste0(out_filepath, i, "_", expr_list[[i]], " FDR0.1 lipids in LIPEA classification.csv"), row.names = FALSE)
}
