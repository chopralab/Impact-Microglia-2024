# code that change our lipid naming classification to LIPEA naming system

lipea_change <- function(df) {
    lipea_class <- c(unique(df$type)) # carry over existing lipid classes

    # list of lipids to scan
    lipids <- unlist(df$lipid)
    # Vector of key words
    keywords <- c("Lysp_PC", "PCo", "PCp", 
                  "Lyso_PG", 
                  "Lyso_PE", "PEo", 
                  "Lyso_PI",
                  "Lyso_PS", 
                  "steryl ester")
    # LIPEA values corresponding to each keyword
    lipea_vals <- c("Lysp_PC" = "LPC", "PCo" = "PC O-", "PCp" = "PC P-", 
                    "Lyso_PG" = "LPG", 
                    "Lyso_PE" = "LPE", "PEo" = "PE O-", 
                    "Lyso_PI" = "LPI",
                    "Lyso_PS" = "LPS", 
                    "steryl ester" = "SE")

    # Loop through the vector and check for matches
    for (lp in lipids) {
    # For each keyword, check if it matches the current text
    for (i in seq_along(keywords)) {
        if (str_detect(lp, keywords[i])) {
        # Append the corresponding custom value to result_vector
        lipea_class <- c(lipea_class, lipea_vals[i])
        }}}
        
    lipea_class <- unique(lipea_class)
    lipea_class[lipea_class == "FFA"] = "FA" # change FA naming
    return(lipea_class)
}
