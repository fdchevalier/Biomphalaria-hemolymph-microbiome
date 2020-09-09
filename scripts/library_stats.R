#!/usr/bin/env Rscript
# Title: library_stats.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-03-25


#==========#
# Packages #
#==========#

library("magrittr")



#===========#
# Variables #
#===========#

# Set working directory
setwd(file.path(getwd(), "scripts"))

# Results and graphs directories
res.d <- "../results/1-qiime/"
graph.d <- "../graphs/"

# Sequencing stat file
myqzv.f <- paste0(res.d, "denoising-stats.qzv")



#=================#
# Data processing #
#=================#

# Check for folder presence
if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

# Loading data
## Identification of the metadata file within the archive file
myfile <- unzip(myqzv.f, list=TRUE)[1] %>% unlist() %>% grep("metadata.tsv", ., value=T)
## Reading the file from the archive file
mydata <- read.delim(unz(myqzv.f, myfile), stringsAsFactors=FALSE)

# Inditifying numeric columns
num <- which(mydata[1,] == "numeric")

# Clean data frame
mydata <- mydata[-1,]

# Convert values to numeric
mydata[ ,num] <- sapply(mydata[,num], function(x) as.numeric(x))

# Export table
mytb <- mydata[, c("Population", "Species", "id", "input", "filtered", "denoised", "merged", "non.chimeric", "Replicate") ]
colnames(mytb) <- c("Population", "Species", "Sample ID", "Number of input sequences", "Number of sequences after filtering", "Number of sequences after denoising", "Number of sequencing after merging", "Number of non chimeric sequences", "Replicate")
mytb <- mytb[order(mytb[,3]),]
cat("Output: ", paste0(res.d, "lib-seq-stat (supp. table 1).tsv"), "\n")
write.table(mytb, paste0(res.d, "lib-seq-stat (supp. table 1).tsv"), sep="\t", quote=FALSE, row.names=FALSE)


for (i in 1:2) {
    cat("Output: ", paste0(graph.d,"seq.ct_",i,".pdf"), "\n")
    pdf(paste0(graph.d,"seq.ct_",i,".pdf"), width=10)
    boxplot(mydata[mydata$Replicate == i,"input"] ~ mydata[mydata$Replicate == i,"Population"], ylab="Number of sequences")
    dev.off()
    
    cat("Output: ", paste0(graph.d,"seq.final_",i,".pdf"), "\n")
    pdf(paste0(graph.d,"seq.final_",i,".pdf"), width=10)
    boxplot(mydata[mydata$Replicate == i,"non.chimeric"] ~ mydata[mydata$Replicate == i,"Population"], ylab="Number of sequences")
    dev.off()
}


#--------------#
# Data summary #
#--------------#

# Per library
mylb.mn <- mean(mydata[,"input"])
mylb.sd <- sd(mydata[,"input"])

mycln <- c("input", "filtered", "denoised", "merged", "non.chimeric")
mylb.tb   <- apply(mydata[,mycln], 1, function(x) x/x[1])
mylb.tb.s <- apply(mylb.tb, 1, function(x) c(mean(x), sd(x), max(x), min(x)) )
rownames(mylb.tb.s) <- c("Mean", "SD", "Max", "Min")


# Per population and replicate (mean and proportion)
mytb.s.mn <- split(mydata, paste(mydata[,"Population"], mydata[,"Replicate"])) %>% lapply(., function(x) sapply(x[,mycln], mean) ) %>% as.data.frame()
mytb.s.sd <- split(mydata, paste(mydata[,"Population"], mydata[,"Replicate"])) %>% lapply(., function(x) sapply(x[,mycln], sd) )   %>% as.data.frame()
mytb.s.p  <- apply(mytb.s.mn, 2, function(x) x/x[1])

mytb.s <- sapply(1:ncol(mytb.s.mn), function(i) paste(round(mytb.s.mn[,i],1), round(mytb.s.sd[,i],1), sep=" ± "))
colnames(mytb.s) <- colnames(mytb.s.mn)
rownames(mytb.s) <- rownames(mytb.s.mn)

cat("Output: ", paste0(res.d, "lib-stat.tsv"), "\n")
write.table(mytb.s, paste0(res.d, "lib-stat.tsv"), sep="\t", quote=FALSE)
