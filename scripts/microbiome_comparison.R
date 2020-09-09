#!/usr/bin/env Rscript
# Title: microbiome_comparison.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-09-02



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    # System
    library("magrittr")
    library("tidyr")

    # Microbiome
    library("qiime2R")
    library("phyloseq")
})



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Working directory
setwd(file.path(getwd(), "scripts"))

# Suppress warning messages
options(warn=-1)

# Qiime files
asv.f  <- "../results/1-qiime/table.qza"
taxa.f <- "../results/1-qiime/rep-seqs_taxa.qza"
tree.f <- "../results/1-qiime/rooted-tree.qza"

# Metadata file
md.f   <- "../data/sample-metadata.tsv"

# Comparison file
comp.f <- "../data/whole_snails_main_bacteria.tsv"

# Population order and color
pop.data <- matrix(c("Ba",      "#ffac40",
                     "Bg26",    "#bf4d4d",
                     "Bg36",    "#a469d6",
                     "Bg121",   "#ff61e8",
                     "BgBRE",   "#6be3a7",
                     "BgBS90",  "#87d687",
                     "BgNMRI",  "#ab7a78",
                     "Water",   "#7798e0"), byrow=TRUE, ncol=2)



#=================#
# Data processing #
#=================#

# Loading ASV data
asv    <- read_qza(asv.f)$data
asv.tb <- otu_table(asv, taxa_are_rows=TRUE)

# Loading taxa data and reorganizing table
taxa <- read_qza(taxa.f)$data

## Split info
cln.nm <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa   <- separate(taxa, 2, cln.nm, sep=";")

## Cleaning names
taxa[, (1:length(cln.nm))+1] <- apply(taxa[,(1:length(cln.nm))+1], 2, function(x) gsub(".*__","", x))

## Polishing data
rownames(taxa) <- taxa[,1]
taxa <- taxa[,-c(1,ncol(taxa))]

## Assign "Unassigned" to all level of unassigned ASV
taxa[ taxa[,1] == "Unassigned", ] <- "Unassigned"

## Convert table
taxa.tb             <- tax_table(as.matrix(taxa))
taxa_names(taxa.tb) <- rownames(taxa)


# Loading tree
tree <- read_qza(tree.f)$data

# Loading metadata and reorder factors
md <- read.delim(md.f, row.names=1)
md[,"Population"] <- factor(md[,"Population"], pop.data[,1])

# Building phyloseq object
mybiom <- phyloseq(asv.tb, taxa.tb, tree, sample_data(md))


# Load comparison file
comp <- read.delim(comp.f, header=TRUE)

# Generate names of genus or family from the table
comp.ls        <- apply(comp, 1, function(x) strsplit(x[1], "\ ") %>% unlist() %>% .[1] %>% grep(., tax_table(mybiom)[,x[2]]))
names(comp.ls) <- strsplit(as.character(comp[,1]), "\ ") %>% lapply(., function(x) x[1]) %>% unlist()

# Identify populations carrying the genus or family
comp.pop <- lapply(comp.ls, function(x) if (length(x) > 0) { (otu_table(mybiom)[x,] %>% colSums() > 0) %>% sample_data(mybiom)[., "Population"] %>% as.matrix() %>% as.vector() %>% unique()})
comp.pop <- comp.pop[names(comp.pop) %>% unique()]

# Print the results on screen
print(comp.pop)
