# Lipidomics of impacted microglia and mixed histiotypical cells

This project was generated with the manuscript "Impact Induces Phagocytic Defect in Reactive Microglia" [insert publication link].

Biological samples are generated by Ruilin Yu and Edmond A. Rogers.
Lipidomics data is collected by Caitlin E. Randolph and Connor H. Beveridge using an Agilent 6410 Triple Quadrupole MS and pre-processed into csv files. 
Code is written by Ruilin Yu and Connor H. Beveridge.

## Usage

Download everything in a folder and run the code (.r files) in code folder sequentially:

1. 1_edger_analysis.r
2. 2_data_visualization.r

1_edger_analysis.r will perform EdgeR analysis with generalized linear model on the ion intensity data in data folder.
2_data_visualization.r will visualize ion intensity data & EdgeR outputs in several formats, e.g. PCA, dot plot, bar plots, etc.

## Acknowledgments

The authors would like to thank the funding agency and any individual who assisted in the completion of this project. 
