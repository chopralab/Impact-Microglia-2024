# QC stacked bar plots

QC_plot <- function(data) {
  data %>%
    arrange(type) %>%
    select(-lipid) %>%
    melt(id.vars = "type", variable.name = "variable") %>%
    ggplot(aes(x = variable, y = value, fill = type)) +
    geom_col() +
    theme_classic(base_size = 20) +
    # labs(color = "group") +
    ggtitle("Breakdown of lipid classes across samples") +
    # theme(axis.text = element_blank()) +
    xlab("Sample") +
    ylab("Total Ion Intensity") +
    theme(aspect.ratio = 1) +
    scale_fill_manual(values = lipid_colors)
}

QC_plot2 <- function(data, lipidclass) {
  blk_loc <- grep("blank", colnames(data), ignore.case = TRUE) 
  data <- data %>%
    filter(type == lipidclass) 
  blk <- data[blk_loc] %>% unlist() %>% log()
  data %>%
    mutate_if(is.numeric, log) %>%
    mutate_if(is.numeric, ~.-blk) %>%
    # mutate_if(is.factor, as.character) %>%
    select(-type, -blk_loc) %>%
    melt(id.vars = "lipid", variable.name = "Variable", value.name = "Value") %>%
    ggplot(aes(x = lipid, y = Value, color = Variable)) +
    geom_point(size = 1) +  # Add dots
    labs(x = "Lipid", y = "Value", title = "Dotplot of Lipid vs. Variables") +
    theme_classic(base_size = 20) +
    ggtitle("Ion Intensity of ", lipidclass) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    xlab("Lipid") +
    ylab("Log2(FC over blank)") +
    theme(aspect.ratio = 1)
}
