# This code visualizes lipidomics data in the following ways
## 1. Principal component analysis of variation between samples
## 2. Ridge Plots
## 3. Scatter plot
## 4. Scree plot
## 5. QC Bar plot
## 6. Bar plot per class
## 7. Heatmaps

# clean environment
rm(list=ls())

# import packages
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(pheatmap)
library(reshape2)

# Get package versions
get_package_version <- function(package_name) {
  package_version <- as.character(packageVersion(package_name))
  return(package_version)
}
packages <- c("dplyr", "tidyverse", "ggplot2", "ggrepel", "ggridges", "pheatmap", "reshape2")
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
   data$lipid <- make.unique(as.character(data$lipid))
   data_list[[paste0("result_expr", i)]] <- data
}
anyDuplicated(data_list[[1]]$lipid) # check for duplicated lipids
head(data_list[[1]])

# output folder
out_filepath = "./plots/"
dir.create(out_filepath, F)

# functions folder
functions_filepath = "./code/functions/"

# lipid colors
## color palette was generated using https://medialab.github.io/iwanthue/
lipid_classes <- unique(data_list[[1]]$type)
lipid_classes
length(lipid_classes) # need 11 colors
lipid_colors <- c("AC" = "#a6cee3", 
                  "Cer" = "#1f78b4", 
                  "CE" = "#b2df8a", 
                  "FFA" = "#33a02c", 
                  "SM" = "#fb9a99", 
                  "PC" = "#e31a1c", 
                  "PE" = "#fdbf6f", 
                  "PG" = "#ff7f00", 
                  "PI" = "#808080", 
                  "PS" = "#cab2d6", 
                  "TAG" = "#6a3d9a")
length(lipid_colors) == length(lipid_classes)

# bar plots for QC
## These bar plots show total ion intensity from each class in a stack for each sampmle
## If sample was preweighed/measured before analysis, expect similar total ion intensity across samples
## If samples are similar in nature (e.g. same tissue type), expect similar distribution of lipid classes

out_filepath = "./plots/QC/"
dir.create(out_filepath, F)

source(paste0(functions_filepath, "QC plot functions.r")) # load functions for QC plots

for (i in 1:length(expr_list)) {
   data <- data_list[[i]] %>%
      select(-"logFC",-"logCPM",-"LR",-"PValue",-"FDR", -"mean1", -"mean2")
   expr_title <- expr_list[[i]]
   
   # Save the plot with specified size
   QC_plot(data)
   ggsave(filename = paste0(out_filepath, "QC bar plot for ", expr_title, " .pdf"), width = 10, height = 10)
}

# dotplots for QC
## for each lipid class, examine if the ion intensity is 10 or more folds higher than blank

for (i in 1:length(expr_list)) {
   data <- data_list[[i]] %>%
      select(-"logFC",-"logCPM",-"LR",-"PValue",-"FDR", -"mean1", -"mean2")
   expr_title <- expr_list[[i]]
      for (j in lipid_classes) {
      QC_plot2(data, j)
      ggsave(filename = paste0(out_filepath, "QC dotplot for ", expr_title, " ", j, ".pdf"), width = 10, height = 10)
   }
}

# FDR<0.1 lipid breakdown
out_filepath = "./plots/breakdown/"
dir.create(out_filepath, F)

source(paste0(functions_filepath, "lipid breakdown.r")) # load functions

for (i in 1:length(expr_list)) {
   data <- data_list[[i]] 
   expr_title <- expr_list[[i]]
   breakdown_lipids2
      for (j in lipid_classes) {
      breakdown_lipids2(data, expr_title)
      ggsave(filename = paste0(out_filepath, "FDR breakdown for ", expr_title, ".pdf"), width = 10, height = 10)
   }
}

# PCA
source(paste0(functions_filepath, "PCA plot functions.r")) # load functions for PCA analysis & visualization

## where to save the pca plots
out_filepath = "./plots/pca/"
dir.create(out_filepath, F)

## calculate principal component
colnames(data_list[[1]])
colnames(data_list[[4]])
colnames(data_list[[7]])

counts_expr1 = 8:15 # where ion intensity values are in the dataframe, excluding blank
counts_expr2 = 8:17
counts_list <- list(counts_expr1, counts_expr1, counts_expr1, counts_expr2, counts_expr2, counts_expr2, counts_expr1)
length(counts_list) == length(expr_list)

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- data_list[[i]]
   counts <- counts_list[[i]]

   data %>% # Scree plot with all lipids
    calculate_pc(counts) %>%
    make_scree_plot(paste0("Scree plot for ", title))    
    ggsave(filename = paste0(out_filepath, i, " Scree plot for ", title, ".pdf"), width = 10, height = 10)

   filter(data, FDR<0.1) %>% # Scree plot with FDR<0.1 cut off
    calculate_pc(counts) %>%
    make_scree_plot(paste0("Scree plot for ", title, " (FDR<0.1)"))    
    ggsave(filename = paste0(out_filepath, i, " Scree plot for ", title, " FDR_0.1.pdf"), width = 10, height = 10)
}

## make PCA plot
## Calculate the pca results
pca_results_list <- list()

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- data_list[[i]]
   counts <- counts_list[[i]]
   # calculate principal component for each experiment
   pca_results <- calculate_pc(data, counts)
   # Add the results to the list with a meaningful name
   pca_results_list[[paste0("result_expr", i)]] <- pca_results
}

## visualize by genotype first
## groups for visualization
expr_list
grp_expr1 = c(rep("Impact D1", 4), rep("Sham D1", 4))
grp_expr2 = c(rep("Impact D7", 4), rep("Sham D7", 4))
grp_expr3 = c(rep("Impact D7", 4), rep("Impact D1", 4))
grp_expr4 = c(rep("Impact 4H", 5), rep("Pre-Impact", 5))
grp_expr5 = c(rep("Impact D3", 5), rep("Pre-Impact", 5))
grp_expr6 = c(rep("Impact D7", 5), rep("Pre-Impact", 5))
grp_expr7 = c(rep("Impact D7 +pHAb", 4), rep("Impact D7 -pHAb", 4))
grp_list <- list(grp_expr1, grp_expr2, grp_expr3, grp_expr4, grp_expr5, grp_expr6, grp_expr7)

## PCA visualization
for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   pca_results <- pca_results_list[[i]]
   groups <- grp_list[[i]]
   make_pca_plot(pca_results, groups, title)
   ggsave(filename = paste0(out_filepath, i, " PCA for ", title, ".pdf"), width = 10, height = 10)
}

## repeat for lipids with FDR<0.1
pca_results_list_FDR <- list()

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- data_list[[i]] %>%
    filter(FDR<0.1)
   counts <- counts_list[[i]]
   # calculate principal component for each experiment
   pca_results <- calculate_pc(data, counts)
   # Add the results to the list with a meaningful name
   pca_results_list_FDR[[paste0("result_expr", i)]] <- pca_results
}

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   pca_results <- pca_results_list_FDR[[i]]
   groups <- grp_list[[i]]
   make_pca_plot(pca_results, groups, title)
   ggsave(filename = paste0(out_filepath, i, " PCA for ", title, " FDR_0.1.pdf"), width = 10, height = 10)
}

# Ridgeplots  
## Ridgeplot shows fold change as a density distribution of lipids for each class
## x-axis is Log2(Fold change KO vs WT)
## y-axis is number of lipids (not to scale as different lipid classes has different total number of lipids)
## Density distribution cannot be generated for lipid classes with <3 lipids

source(paste0(functions_filepath, "ridgeplot functions.r")) # load functions for ridge plots

out_filepath = "./plots/ridgeplots/"
dir.create(out_filepath, F)

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- data_list[[i]]

   ridge_plot(data, title) # Ridgeplot for all lipids in each experiment
   ggsave(paste0(out_filepath, i, " Ridge plot for ", title, ".pdf"), width = 10, height = 10)
   filter(data, FDR<0.1) %>% ridge_plot(title) # Ridgeplot for lipids with FDR<0.1 for each experiment
   ggsave(paste0(out_filepath, i, " Ridge plot for ", title, " FDR_0.1.pdf"), width = 10, height = 10)
}


# Heatmaps
source(paste0(functions_filepath, "heatmaps functions.r")) # load functions for heatmaps

out_filepath = "./plots/heatmaps/"
dir.create(out_filepath, F)

## calculate z scores for heatmaps
heatmap_results_list <- list()
heatmap_results_list_FDR <- list()

for (i in 1:length(expr_list)) {
   data <- data_list[[i]]
   counts <- counts_list[[i]]
   # calculate heatmap statistics for each experiment
   heatmap_results <- calculate_z_scores(data, cols = counts)
   heatmap_results_FDR <- calculate_z_scores(data, FDR_cutoff = 0.1, cols = counts)
   # Add the results to the list with a meaningful name
   heatmap_results_list[[paste0("result_expr", i)]] <- heatmap_results
   heatmap_results_list_FDR[[paste0("results_expr", i)]] <- heatmap_results_FDR
}

dim(heatmap_results_list[[1]]) # preview :P
head(heatmap_results_list[[1]])
dim(heatmap_results_list_FDR[[1]]) # preview :P
head(heatmap_results_list_FDR[[1]])

## make heatmaps
expr_list
grp_expr1 <- c("Impact D1" = "#A7535A", "Sham D1" = "#144A74")
grp_expr2 <- c("Impact D7" = "#A7535A", "Sham D7" = "#144A74")
grp_expr3 <- c("Impact D7" = "#A7535A", "Impact D1" = "#DFC243")
grp_expr4 <- c("Impact 4H" = "#49442B", "Pre-Impact" = "#144A74")
grp_expr5 <- c("Impact D3" = "#DFC243", "Pre-Impact" = "#144A74")
grp_expr6 <- c("Impact D7" = "#A7535A", "Pre-Impact" = "#144A74")
grp_expr7 <- c("Impact D7 +pHAb" = "#49442B", "Impact D7 -pHAb" = "#A7535A")
grp_col_list <- list(grp_expr1, grp_expr2, grp_expr3, grp_expr4, grp_expr5, grp_expr6, grp_expr7)

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- heatmap_results_list[[i]]
   data_fdr <- heatmap_results_list_FDR[[i]]  
   grp <- grp_list[[i]]
   grp_colors <- grp_col_list[[i]]

   rownames(data) <- data$lipid
   rownames(data_fdr) <- data_fdr$lipid
   
   sample_ann <- data.frame(Group = grp, row.names = colnames(data)[-c(1,2)])
   lipid_ann <- data.frame(Class = data$type, row.names = data$lipid)
   lipid_ann_fdr <- data.frame(Class = data_fdr$type, row.names = data_fdr$lipid)
   ann_colors <- list(LipidClass = lipid_colors, Group = grp_colors)
   make_heatmap(data, title, lipid_ann, sample_ann, ann_colors, path = paste0(out_filepath, i, " "))
   make_heatmap(data_fdr, paste0(title, "_FDR_0.1"), lipid_ann_fdr, sample_ann, ann_colors, path = paste0(out_filepath, i, " "))
}

dev.off()
graphics.off()

# scatter plot
## scatterplot shows the average intensity of lipids detected in WT vs. Pla2g2fKO samples
## x-axis is Log10 transformed fold-change of WT/Blank
## y-axis is Log10 transformed fold-change of Pla2g2fKO/Blank
## colored dots are lipids with FDR < 0.1, and they are colored based on fold-change of Pla2g2fKO/WT

source(paste0(functions_filepath, "scatterplot functions.r")) # load functions for scatterplots

out_filepath = "./plots/scatter/"
dir.create(out_filepath, F)

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- data_list[[i]]
   make_scatter(data, title, 14)
   ggsave(filename = paste0(out_filepath, i, " Scatter plot for ", title, ".pdf"), width = 10, height = 10)
}

# bar plots for each class
## bar plots show normalized average intensity of each lipid
## After average is calculated for each group, lipids are normalized by dividing over the larger group
out_filepath = "./plots/bar/"
dir.create(out_filepath, F)

for (i in 1:length(expr_list)) {
   title <- expr_list[i]
   data <- data_list[[i]]

   bar_data <- data[c("mean1","mean2")] %>% 
    apply(1, function(x)(x/max(x))) %>% # range normalization by dividing over the maximum
    t() %>%
    as_tibble() %>%
    mutate(lipid = data$lipid,
            type = data$type) %>%
    as.data.frame() %>%
    melt(id = c("lipid", "type"))

   for (lipid_class in lipid_classes) {
   bar_data %>%
        filter(type == lipid_class) %>%
        ggplot(aes(x = lipid, y = value, fill = variable)) +
        geom_col(position = "dodge") +
        xlab("") +
        ylab("Normalized Intensity") +
        ggtitle(paste0("Normalized intensity for ", lipid_class)) +
        theme_classic(base_size = 20) +
        theme(axis.text.x = element_text(angle = 90)) +
        scale_fill_manual(values = c("#4E4E4E", "#53B3CB"))
    # readline(prompt="Press [enter] to proceed")
    ggsave(filename = paste0(out_filepath, i, " Normalized intensity for ", lipid_class, " in ", title, ".pdf"), width = 10, height = 10)
   }
}

# some lipid classes have too many lipids, so their names get squished on the x-axis
# these can be fixed by changing the font size of the text in a figure editor (e.g. Affinity Designer)