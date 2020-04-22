# System
options("repos" = c(CRAN = "https://cran.revolutionanalytics.com"))
library("devtools")
library("BiocManager")

# # Generic
# install_version("plyr",         version="1.8.4",   upgrade="never")
# install_version("DescTools",    version="0.99.28", upgrade="never")
# install_version("rcompanion",   version="2.0.0",   upgrade="never")
# install_version("tidyr",        version="0.8.1",   upgrade="never")
# install_version("reshape2",     version="1.4.3",   upgrade="never")
# install_version("multcompView", version="0.1-7",   upgrade="never")

# # Microbiome analysis
# install("microbiome", update=FALSE)
# install("phyloseq",   update=FALSE)
# install("picante",    update=FALSE)          # PD function
# install_version("phylogram",  version="2.1.0",  upgrade="never")      # prune function
# install_version("dendextend", version="1.13.4", upgrade="never")
# install_version("pvclust",    version="2.2-0",  upgrade="never")
# install_version("vegan",      version="2.5-2",  upgrade="never")
# install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis", ref="6e09713", upgrade="never")
# install_github("ropensci/taxa", ref="49969dc", upgrade="never")
# install_github("jbisanz/qiime2R", ref="08e3b30", upgrade="never")

# # Pathway analysis
# install("ALDEx2", update=FALSE)
# install_version("pbkrtest", version="0.4-7", upgrade="never")
# install_github("ggloor/CoDaSeq/CoDaSeq", ref="73f77ce", upgrade="never")

# Graph
install_version("gridBase",   version="0.4-7",   upgrade="never")
install_version("plotrix",    version="3.7-5",   upgrade="never")
install_version("ggplot2",    version="3.1.0",   upgrade="never")
install_version("ggpubr",     version="0.2.5",   upgrade="never")
install_version("ggdendro",   version="0.1-20",  upgrade="never")
install_version("colorspace", version="1.3-2",   upgrade="never")
install_version("cowplot",    version="1.0.0",   upgrade="never")
install_version("grid",       version="3.5.0",   upgrade="never")
install_version("igraph",     version="1.2.4.1", upgrade="never")
install_version("lemon",      version="0.4.3",   upgrade="never")
