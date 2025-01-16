# EdgeR functions

perform_raw_analysis <- function(data, counts_col, lipid_col, group, design) {
    y <- DGEList(counts = data[,counts_col],
                 genes = data[,lipid_col], 
                 group = group) # create DGEList object from data
    y <- calcNormFactors(y, method = "TMM") # calculate scaling factor for raw data
    y <- estimateCommonDisp(y, design = design) # estimate a common dispersion value across all genes/lipids
    y
}

calculate_significance <- function(data, design, contrast){
    data %>% 
    glmFit(design = design) %>% # fit NB model to the counts for each gene/lipid
    glmLRT(contrast = contrast) # conduct likelihood ratio tests for the given contrasts (to be tested equal to 0)
}