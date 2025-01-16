# Heatmaps

## calculate z scores for heatmaps
calculate_z_scores <- function(data, FDR_cutoff=0, cols) {
    if (FDR_cutoff==0) {
       data <- data 
    } else {
       data <- filter(data, FDR<FDR_cutoff) # filter for FDR < cutoff value
    }    
    result <- data %>%
       column_to_rownames('lipid') %>%
       select(cols-1) %>% 
       apply(1, scale) %>% # scale function calculates z-scores using (x - mean(x)) / sd(x)
       t() %>%
       as.data.frame() %>%
       rownames_to_column(var = "lipid")
    # cleaning up the result
    colnames(result) <- c("lipid", colnames(data[counts]))
    merged_result <- data %>%
        select(c(1,2)) %>%
        merge(result, by = "lipid") %>%
        arrange(type)
    return(merged_result)
}

## make heatmaps
make_heatmap <- function(data, title, lipid_ann, sample_ann, ann_colors, path){
    pheatmap(data[-c(1,2)], 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE, 
    show_rownames = FALSE,
    border_color = NA, 
    fontsize = 10,
    legend_breaks = -2:2,
    annotation_row = lipid_ann, 
    annotation_col = sample_ann, 
    annotation_colors = ann_colors,
    main = paste0("Heatmap for ", title), 
    filename = paste0(path, "Heatmap for ", title, ".pdf"))
}

make_heatmap2 <- function(data, title, lipid_ann, sample_ann, ann_colors, path){
    pheatmap(data[-c(1,2)], 
    cluster_rows = FALSE, 
    cluster_cols = FALSE, 
    show_colnames = FALSE, 
    show_rownames = TRUE, # show lipid names here
    border_color = NA, 
    fontsize = 10,
    legend_breaks = -2:2,
    annotation_row = lipid_ann, 
    annotation_col = sample_ann, 
    annotation_colors = ann_colors,
    main = paste0("Heatmap for ", title), 
    filename = paste0(path, "Heatmap for ", title, ".pdf"))
}