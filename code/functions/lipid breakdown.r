
breakdown_lipids <- function(df, title) {
    levels <- table(df$ClassKey) %>%
        as.data.frame() %>%
        arrange(Freq) %>%
        select(Var1) %>%
        unlist()
        
    df %>%
    ggplot(aes(x = ClassKey, fill = TotalGrade)) +
    geom_bar() +
    theme_classic(base_size = 20) +
    labs(x = "Lipid", y = "Frequency") +
    theme_classic(base_size = 20) +
    ggtitle(paste0("Lipid Breakdown of ", title)) +
    scale_x_discrete(limits = levels) +
    theme(aspect.ratio = 2) + 
    coord_flip() +
    scale_fill_manual(values = c("#00a8e9", "#ff4268", "#a9d291"))
}

breakdown_lipids2 <- function(df, title) {
    levels <- table(df$type) %>%
        as.data.frame() %>%
        arrange(Freq) %>%
        select(Var1) %>%
        unlist()

    df$FDR_0.1 <- df$FDR < 0.1
        
    df %>%
    ggplot(aes(x = type, fill = FDR_0.1)) +
    geom_bar(position = "stack") +
    theme_classic(base_size = 20) +
    labs(x = "Lipid", y = "Frequency") +
    theme_classic(base_size = 20) +
    ggtitle(paste0("Lipid Breakdown of ", title)) +
    scale_x_discrete(limits = levels) +
    scale_fill_manual(values = c("#bca3ea", "#6de5cc")) +
    theme(aspect.ratio = 1) + 
    coord_flip()
}

