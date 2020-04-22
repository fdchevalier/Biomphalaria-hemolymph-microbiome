#!/usr/bin/env Rscript
# Title: pathway_analysis.R
# Version: 1.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-03-25



#==========#
# Packages #
#==========#

library("ALDEx2")
library("CoDaSeq")
library("zCompositions")
library("igraph")


#===========#
# Variables #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

# Metadata and pathway files
md.f <- "../data/sample-metadata.tsv"
pw.f <- c("../results/2-picrust2/rep1/pathways_out/path_abun_unstrat_descrip.tsv.gz", "../results/2-picrust2/rep2/pathways_out/path_abun_unstrat_descrip.tsv.gz")

# Population order and color
pop.data <- matrix(c("Ba",      "#ffac40",
                     "Bg26",    "#bf4d4d",
                     "Bg36",    "#a469d6",
                     "Bg121",   "#ff61e8",
                     "BgBRE",   "#6be3a7",
                     "BgBS90",  "#87d687",
                     "BgNMRI",  "#ab7a78",
                     "Water",   "#7798e0"), byrow=TRUE, ncol=2)

# Result folder
res.d <- "../results/3-pathways/"



#=================#
# Data processing #
#=================#

if (!dir.exists(res.d)) { dir.create(res.d, recursive=TRUE) }

# Loading metadata
md.t <- read.delim(md.f, header=TRUE, stringsAsFactors=FALSE)

nb.rep <- unique(md.t[,"Replicate"])


mycomb  <- cbind(c("snail", "Water"), combn(pop.data[,1][-8],2)) # Remove water
mydenom <- c("zero", rep("all", ncol(mycomb)-1))

pw.ls <- vector("list", 3)

for (j in 1:ncol(mycomb)) {

    mycomb.tmp <- mycomb[,j]

    if (j == 1) {
        md.tmp <- md.t
    } else {
        md.tmp <- md.t[ grepl(paste(mycomb.tmp, collapse="|"), md.t[,"Population"]), ]
    }

    for (i in nb.rep) {

        # Select first replicate
        md <- md.tmp[ md.tmp[,"Replicate"] == i, ]

        # Loading pathway information
        pw <- read.delim(pw.f[i], header=TRUE, stringsAsFactors=FALSE, check.names=FALSE, row.names=1)
        pw.r <- pw[,-1]

        pw.r <- pw[, md[,1]]

        # Rounding counts
        pw.r <- round(pw.r)

        # Count threshold
        count <- 10

        # Input matrix
        aldex.in <- data.frame(pw.r[rowMeans(pw.r) > count,], check.names=FALSE)

        # Condition vector
        conds <- sapply(colnames(aldex.in), function(x) md[md[,1] == x, 2])

        if (j == 1 ) {
            ## Snail vs water
            conds[ conds != "Water" ] <- "snail"
        } else {
            conds <- sapply(colnames(aldex.in), function(x) md[md[,1] == x, "Population"])
        }

        # Running ALDEx
        ## Denominator must be "zero" because "all" or "iqrl" generate a binomial not 0 centered. This is likely due to very specific pathway between water and snails.
        x.all <- aldex(aldex.in, conds, mc.samples=500, test="t", effect=TRUE, verbose=TRUE, denom=mydenom[j])

        pw.ls[[1]] <- x.all 

        # FDR seting
        myfdr <- 0.05

        # Pathway in snails
        myvec <- rownames(x.all)[which(x.all$we.eBH < myfdr & x.all$effect < -1)]
        pw.ls[[2]] <- cbind(pw[myvec,1, drop=FALSE], x.all[myvec, c("effect", "we.eBH")])

        # Pathway in water
        myvec <- rownames(x.all)[which(x.all$we.eBH < myfdr & x.all$effect > 1)]
        pw.ls[[3]] <- cbind(pw[myvec,1, drop=FALSE], x.all[myvec, c("effect", "we.eBH")])

        assign(paste(paste0(mycomb.tmp, collapse="_"), "rep", i, sep="."), pw.ls)
    }
}


# Generate tables
for (j in 1:ncol(mycomb)) {
    mycomb.tmp <- mycomb[,j]
    myobj      <- paste(paste0(mycomb.tmp, collapse="_"), "rep.", sep=".")
    
    # First condition
    ## if length(nb.rep ) > 2, run a loop to merge successively each replicate (but merge the two first replicates first)
    mytb      <- merge(get(paste0(myobj,1))[[2]], get(paste0(myobj,2))[[2]], by="row.names", all=TRUE)
    mydsp.cln <- grepl("description", colnames(mytb))
    mydsp     <- apply(mytb[,mydsp.cln], 1, function(x) unique(na.omit(x)))
    mytb      <- cbind(mytb[,1], mydsp, mytb[,!mydsp.cln][,-1])

    write.table(mytb, paste0(res.d, myobj, mycomb.tmp[1], ".tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    
    # Second condition
    mytb <- merge(get(paste0(myobj,1))[[3]], get(paste0(myobj,2))[[3]], by="row.names", all=TRUE)
    mydsp.cln <- grepl("description", colnames(mytb))
    mydsp     <- apply(mytb[,mydsp.cln], 1, function(x) unique(na.omit(x)))
    mytb      <- cbind(mytb[,1], mydsp, mytb[,!mydsp.cln][,-1])

    write.table(mytb, paste0(res.d, myobj, mycomb.tmp[2], ".tsv"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
}




