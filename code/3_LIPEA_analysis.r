# This code filters analyzed data for significant lipids (FDR<0.1) and compare to LIPEA's lipid list
# The output of this code can be uploaded to LIPEA (https://hyperlipea.org/analyze) for pathway enrichment analysis

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(tidyverse)
library(readxl)
library(ggplot2)


# Get package versions
get_package_version <- function(package_name) {
  package_version <- as.character(packageVersion(package_name))
  return(package_version)
}
packages <- c("dplyr", "tidyverse", "readxl", "ggplot2")
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

# read in LIPEA analysis results
# analyses were performed against pre-defined background of mouse (mus musculus)
lipea_filepath <- "./lipea_analysis/"
lipea_file_list = list.files(path=lipea_filepath, pattern=".xls", all.files=FALSE, full.names=FALSE)
lipea_file_list

lipea_results <- list()
for (i in seq_along(expr_list)) {
   data <- read_xls(paste0(lipea_filepath, lipea_file_list[[i]])) %>%
    rename(pathway = "Pathway name", 
           lipids = "Pathway lipids", 
           p = "p-value",
           converted_number = "Converted lipids (number)", 
           converted_percent = "Converted lipids (percentage)", 
           Benjamini_corr = "Benjamini correction", 
           Bonferroni_corr = "Bonferroni correction")
   lipea_results[[expr_list[[i]]]] <- data
}
sapply(lipea_results, dim) # preview
head(lipea_results[[1]])
colnames(lipea_results[[1]])

# visualize lipea pathways
# Create the bubble plot using ggplot2
out_filepath = "./plots/lipea/"
dir.create(path = out_filepath, F)

for (i in seq_along(expr_list)){
  data <- lipea_results[[i]] %>% filter(p < 0.05)
  title <- expr_list[[i]]
  
  data$pathway <- reorder(data$pathway, data$converted_percent) # arrange by decreasing converted lipid %
  data %>%
    ggplot(aes(x = converted_percent, y = pathway, size = -log10(p))) +
    geom_point(alpha = 0.6, fill = "light blue") +
    scale_size_continuous(range = c(3, 15)) +
    theme_classic(base_size = 20) +
    labs(title = "Top Lipidomic Pathways (P<0.05)", 
         x = "% Converted Lipids",
         size = "-Log10(P-value)") +
    theme(aspect.ratio = 2, legend.key.size = unit(0.01, "cm"))  
  ggsave(filename = paste0(out_filepath, i, " ", title, " top lipid pathways.pdf"), width = 10, height = 10)
}
