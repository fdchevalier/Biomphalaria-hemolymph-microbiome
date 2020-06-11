#!/usr/bin/env Rscript
# Title: microbiome_density.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-03-25



#==========#
# Packages #
#==========#

suppressMessages({
    library("qiime2R")
    library("magrittr")
    library("ggpubr")

    library("rcompanion")
    library("plyr")
    library("multcompView")
})



#===========#
# Functions #
#===========#

# Compute geom_boxplot upper whisker
## source: https://stackoverflow.com/a/53795076
q <- function(x) {
    y   <- sort(x)
    iqr <- quantile(y, 0.75) - quantile(y, 0.25)
    myq <- tail(y[which(y <= quantile(y, 0.75) + 1.5*iqr)],1)
    return(myq)
}



#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Qiime file
asv.f  <- "../results/1-qiime/table.qza"

# qPCR file
qpcr.f <- "../data/qPCR/16S_v4_quantification.csv"

# Picrust file
picrust.f <- "../results/2-picrust2/marker_predicted_and_nsti.tsv.gz"

# Population order and color
pop.data <- matrix(c("Ba",      "#ffac40",
                     "Bg26",    "#bf4d4d",
                     "Bg36",    "#a469d6",
                     "Bg121",   "#ff61e8",
                     "BgBRE",   "#6be3a7",
                     "BgBS90",  "#87d687",
                     "BgNMRI",  "#ab7a78",
                     "Water",   "#7798e0"), byrow=TRUE, ncol=2)

# Graph folder
graph.d <- "../graphs/"

# GGplot options
theme_set(theme_classic())

# Volumes
elu.vol <- 50
hml.vol <- 40



#=================#
# Data processing #
#=================#

# Load data
mydata    <- read.csv(qpcr.f, header=FALSE, stringsAsFactors=FALSE)
taxa      <- read_qza(asv.f)$data
mrkr.pred <- read.delim(picrust.f)

# Note: no data cleaning because contaminant are also quantified with qPCR
asv <- apply(taxa, 2, function (x) sum(x > 0) ) %>% data.frame(asv=.)

# Reset 16S count prediction for ASV with nsis > 2 because of incertainty of annotation (see PiCRUST documentation)
mrkr.pred[mrkr.pred[, 3] > 2, 2] <- 1

# Adjust marker count
mrkr.count <- apply(mrkr.pred[,1:2], 1, function(x) taxa[ x[1], ] * as.numeric(x[2]))
mrkr.sum   <- as.data.frame(rowSums(mrkr.count))

# Bacterial index
mrkr.idx <- colSums(taxa) / mrkr.sum

# Merge data
mydata <- merge(mydata, mrkr.idx, by.x=1, by.y="row.names")


# Correct and normlize data
## Dilution factor to reflect sample (hemolymph or water) concentration (instead of elution concentration)
dil.fct <- elu.vol/hml.vol

## Correct qPCR values with dilution factor
mydata[,2] <- mydata[,2] * dil.fct

## Normalize qPCR count
mydata[,ncol(mydata)+1] <- mydata[,2] * mydata[,3]

## Merge ASV number
mydata <- merge(mydata, asv, by.x=1, by.y="row.names")


# Population labels
mydata[,ncol(mydata)+1] <- lapply(strsplit(mydata[,1], ".", fixed=TRUE), function(x) x[1]) %>% unlist()
mydata[ grepl("Water", mydata[,1]), ncol(mydata)] <- "Water"
colnames(mydata)[ncol(mydata)] <- "Population"
mydata[,"Population"] <- factor(mydata[,"Population"], pop.data[,1])
mydata <- mydata[order(match(mydata$Population,pop.data[,1])),] # Reorder

# Sample labels
mydata[,ncol(mydata)+1] <- lapply(strsplit(mydata[,1], ".[12]$"), function(x) x[1]) %>% unlist()
colnames(mydata)[ncol(mydata)] <- "ReplicateID"

# Replicate
mydata[,ncol(mydata)+1] <- lapply(strsplit(mydata[,1],".",fixed=TRUE), function(x) rev(x)[1]) %>% unlist()
colnames(mydata)[ncol(mydata)] <- "Replicate"


# Raw 16S copy number between populations
cat("Raw 16S copy number\n\n")
mycln <- 2

## Mean and range
cat("16S copy number:\n\t- mean: ", mean(mydata[,mycln]), "\n\t- range (min max): ", range(mydata[,mycln]), "\n\n")

## Test for normal distribution without replicates
not.norm <- any(sapply(pop.data[,1], function(x) shapiro.test(mydata[mydata[,"Population"] == x, mycln])$p.value) < 0.05)
if (not.norm) cat("Data per population is not following normal distribution\n")

## Formating data
qpcr.d <- mydata[, c(7,mycln,6)]

## Kruskal-Wallis test
qpcr.d.k <- kruskal.test(qpcr.d[,2], qpcr.d[,3])
print(qpcr.d.k)

## Statistical comparison
qpcr.d.stat    <- pairwise.wilcox.test(qpcr.d[,2], qpcr.d[,3], p.adjust="none")
qpcr.d.tb      <- fullPTable(qpcr.d.stat$p.value)
qpcr.d.letters <- multcompLetters(qpcr.d.tb, compare="<=", threshold=0.05, Letters=letters)$Letters

## Print results
cat("Pairwise comparisons of 16S copy number between population (mean of the replicates): ", qpcr.d.letters, "\n\n\n")


# Estimated bacteria between populations
cat("Estimateed bacteria\n\n")
mycln <- 4

## Mean and range
cat("Bacteria number:\n\t- mean: ", mean(mydata[,mycln]), "\n\t- range (min max): ", range(mydata[,mycln]), "\n\n")

## Test for normal distribution without replicates
not.norm <- any(sapply(pop.data[,1], function(x) shapiro.test(mydata[mydata[,"Population"] == x, mycln])$p.value) < 0.05)
if (not.norm) cat("Data per population is not following normal distribution\n")

## Formating data
asv.d <- mydata[, c(7,mycln,6)]

## Kruskal-Wallis test
asv.d.k <- kruskal.test(asv.d[,2], asv.d[,3])
print(asv.d.k)

## Statistical comparison
asv.d.stat    <- pairwise.wilcox.test(asv.d[,2], asv.d[,3])
asv.d.tb      <- fullPTable(asv.d.stat$p.value)
asv.d.letters <- multcompLetters(asv.d.tb, compare="<=", threshold=0.05, Letters=letters)$Letters

## Print results
cat("Pairwise comparisons of estimated bacteria between population (mean of the replicates): ", asv.d.letters, "\n\n")

## Mean bacteria in snails
snail.b <- asv.d[asv.d[,"Population"] != "Water",2]
cat("Mean bacteria number in snails:\n\t- mean: ", mean(snail.b), "\n\t- standard error: ", sd(snail.b)/sqrt(length(snail.b)), "\n\t- range (min max): ", range(snail.b), "\n\n")

## Mean bacteria in water
water.b <- asv.d[asv.d[,"Population"] == "Water",2]
cat("Mean bacteria number in waters:\n\t- mean: ", mean(water.b), "\n\t- standard error: ", sd(water.b)/sqrt(length(water.b)), "\n\t- range (min max): ", range(water.b), "\n\n")



#=========#
# Figures #
#=========#

# Check graph folder
if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

#--------------------#
# Microbiome density #
#--------------------#

p.ls <- vector("list", 3)

ylim1 <- boxplot(mydata$V2 ~ mydata$Population, outline=FALSE, plot=FALSE)$stats[1:4,] %>% max() %>% signif(., 1)
ylim2 <- boxplot(mydata$V4 ~ mydata$Population, outline=FALSE, plot=FALSE)$stats[1:4,] %>% max() %>% signif(., 1)
ylim <- max(c(ylim1, ylim2))
# Raw data
p.ls[[1]] <- ggplot(mydata, aes(x = Population, y = V2, fill = Population)) +
        geom_boxplot(outlier.shape=NA) + xlab("") + ylab("16S copy number per µL") +
        scale_y_continuous(limits = c(0, ylim)) +
        scale_fill_manual(values = pop.data[,2]) +
        stat_summary(geom = 'text', label = qpcr.d.letters, fun.y = q, vjust = -1) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Normalized data
p.ls[[2]] <- ggplot(mydata, aes(x = Population, y = V4, fill = Population)) +
        geom_boxplot(outlier.shape=NA) + xlab("") + ylab("Estimated number of bacteria per µL") + 
        scale_y_continuous(limits = c(0, ylim)) +
        scale_fill_manual(values = pop.data[,2]) +
        stat_summary(geom = 'text', label = asv.d.letters, fun.y = q, vjust = -1) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none")

# Correlation
p.ls[[3]] <- ggplot(mydata, aes(x = V4, y = asv)) +
        geom_point(aes(color = Population)) + 
        stat_smooth(method = lm, se = FALSE) +
        xlab("Bacterial number") + ylab("Number of observed ASVs") +
        scale_color_manual(values = pop.data[,2]) +
        theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
        stat_cor(method = "kendall", cor.coef.name = "tau", label.x.npc = "right", label.y.npc = "top", hjust = 1)

p <- ggarrange(plotlist=p.ls, ncol=1, labels=LETTERS[1:length(p.ls)])
ggsave(paste0(graph.d,"Fig. 6 - qPCR-diversity.pdf"), p, height=15, width=5)


# Clean tmp file
if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }
