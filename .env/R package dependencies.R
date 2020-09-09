# System
options("repos" = c(CRAN = "https://cran.revolutionanalytics.com"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref = "7ea9440", upgrade = "never")

# Microbiome analysis
install_bioc("3.10/microbiome", upgrade = "never")
install("phyloseq", update = FALSE)
install_bioc("3.8/DESeq2", upgrade = "never")
install_version("phylogram", version = "2.1.0", upgrade = "never")      # prune function
install_version("usedist",   version = "0.4.0", upgrade = "never")
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", ref = "6e09713", upgrade = "never")
install_github("ropensci/taxa",   ref = "49969dc", upgrade = "never")
install_github("jbisanz/qiime2R", ref = "08e3b30", upgrade = "never")

# Pathway analysis
install_bioc("3.9/ALDEx2", upgrade = "never")
install_github("ggloor/CoDaSeq/CoDaSeq", ref = "73f77ce", upgrade = "never")
