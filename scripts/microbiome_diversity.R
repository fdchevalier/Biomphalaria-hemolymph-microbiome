#!/usr/bin/env Rscript
# Title: microbiome_diversity.R
# Version: 1.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-08-25



#==========#
# Packages #
#==========#

cat("Loading packages...\n")

suppressMessages({
    # System
    library("magrittr")
    library("plyr")             # function count
    library("tidyr")

    # Microbiome
    library("DESeq2")
    library("microbiome")       # function core
    library("qiime2R")
    library("phyloseq")
    library("picante")          # PD function
    library("pairwiseAdonis")   # see github
    library("phylogram")        # prune function
    library("pvclust")
    library("dendextend")       # pvclust_edges function
    library("usedist")          # dist_multi_centroids function

    # Graph
    library("colorspace")       # function rainbow_hcl
    library("cowplot")
    library("ggpubr")
    library("ggdendro")
    library("grid")             # function grid.rect
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

cat("Loading functions...\n")
source("microbiome_diversity_func.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Suppress warning messages
options(warn=-1)

# Qiime files
asv.f  <- "../results/1-qiime/table.qza"
taxa.f <- "../results/1-qiime/rep-seqs_taxa.qza"
tree.f <- "../results/1-qiime/rooted-tree.qza"

# Metadata file
md.f   <- "../data/sample-metadata.tsv"

# Contaminant file
blast.f <- "../results/1-qiime/rep-seqs_unassigned/unassigned.blastn.tsv"

# Population order and color
pop.data <- matrix(c("Ba",      "#ffac40",
                     "Bg26",    "#bf4d4d",
                     "Bg36",    "#a469d6",
                     "Bg121",   "#ff61e8",
                     "BgBRE",   "#6be3a7",
                     "BgBS90",  "#87d687",
                     "BgNMRI",  "#ab7a78",
                     "Water",   "#7798e0"), byrow=TRUE, ncol=2)

# Results and Graph output directory
res.d <- "../results/1-qiime/"
graph.d <- "../graphs/"

# GGplot options
theme_set(theme_classic())

# Seed for reproducibility
myseed <- 26017
set.seed(myseed)

# Alpha threshold for DESeq analysis
myalpha.lim <- 0.05



#=================#
# Data processing #
#=================#

if (! dir.exists(graph.d)) { dir.create(graph.d, recursive=TRUE) }

#----------------#
# Preparing data #
#----------------#

# Loading ASV data
asv    <- read_qza(asv.f)$data
asv.tb <- otu_table(asv, taxa_are_rows=TRUE)

# Loading contaminants
blast <- read.delim(blast.f, header=FALSE)

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

# Populations 
mypop <- cbind(rownames(md), md[,2:4])

# Number of replicates
rep    <- unclass(sample_data(mybiom)[,"Replicate"])[[1]]
nb.rep <- unique(rep)


# Cleaning data from contaminants
nb.asv   <- nrow(tax_table(mybiom))
nb.unass <- (taxa[,1] == "Unassigned") %>% sum()

## Mitochondria and chloroplast
my.mt <- apply(tax_table(mybiom), 1, function(x) any(grepl("mitochondria", x, ignore.case=TRUE)))
my.cp <- apply(tax_table(mybiom), 1, function(x) any(grepl("chloroplast", x, ignore.case=TRUE)))

## Eukaryote contaminants among the unassigned
my.euk <- blast[ blast[,14] == "Eukaryota", 1]
my.euk <- rownames(tax_table(mybiom)) %in% my.euk

## Exporting contaminent list
write.table(rownames(tax_table(mybiom))[(my.mt | my.cp | my.euk)], "../data/contaminants", row.names=FALSE, col.names=FALSE, quote=FALSE)

## Remove contaminants
mybiom <- subset_taxa(mybiom, ! (my.mt | my.cp | my.euk))

## Final numbers
cat("\nNumber of ASVs removed:", 
    "\n\t- mitochondria and chloroplasts:", sum(my.mt | my.cp), "/", signif(sum(my.mt | my.cp) / nb.asv * 100, 2), "% of the total ASVs",
    "\n\t- eukaryotes:", sum(my.euk), "/", signif(sum(my.euk) / nb.unass * 100, 2), "% of the unassigned ASVs", "/", signif(sum(my.euk) / nb.asv * 100, 2), "% of the total ASVs",
    "\n\t- total:", sum(my.mt | my.cp | my.euk), "/", signif(sum(my.mt | my.cp | my.euk)/ nb.asv * 100, 2), "% of the total ASVs", "\n\n")


#-------------------#
# Rarefaction curve #
#-------------------#

p <- ggrare(mybiom, step = 1000, color = "Population", se = FALSE, plot=FALSE) +
        scale_color_manual(values = pop.data[,2]) +
        theme(legend.position="none") +
        facet_wrap(~Population)

# # Adjust legend
# ## Comment legend line above if legend needed
# pdf(NULL)
# p <- p + guides(col = guide_legend(ncol=2))
# p <- shift_legend(p)
# dev.off()

ggsave(paste0(graph.d,"Supp. Fig. 2 - rarefaction.pdf"), p, height=7, width=7, useDingbats=FALSE)


#-----------------#
# Alpha Diversity #
#-----------------#

myalpha <- microbiome::alpha(mybiom)
mypd    <- pd(t(asv), tree)
myalpha <- merge(myalpha, mypd, by="row.names")
myalpha <- merge(myalpha, mypop, by.x=1, by.y="row.names")

myidx <- matrix(c("observed",           "Observed richness",
                  "chao1",              "Chao1 richness",
                  "diversity_shannon",  "Shannon diversity",
                  "PD",                 "Faith's Phylogenetic diversity",
                  "evenness_simpson",   "Simpson evenness",
                  "dominance_simpson",  "Simpson dominance"), ncol=2, byrow=TRUE)


# Snail and water diversity
snail.mean <- mean(myalpha[ myalpha[,"Population"] != "Water", 2])
snail.se   <- se(myalpha[ myalpha[,"Population"] != "Water", 2]) / sqrt(length(myalpha[ myalpha[,"Population"] != "Water", 2]))

water.mean <- mean(myalpha[ myalpha[,"Population"] == "Water", 2])
water.se   <- se(myalpha[ myalpha[,"Population"] == "Water", 2]) / sqrt(length(myalpha[ myalpha[,"Population"] == "Water", 2]))

cat("\nObserved number of ASVs in snails:\n\t- mean:", snail.mean, "\n\t- standard error:", snail.se, "\n")
cat("\nObserved number of ASVs in water:\n\t- mean:",  water.mean, "\n\t- standard error:", water.se, "\n\n")


# Statistical tests
mya.test    <- vector("list", length(nb.rep))
mya.letters <- vector("list", length(nb.rep))

for (r in nb.rep) { 
    mya.test[[r]]   <- vector("list", length(myidx[,1]))
    for (i in myidx[,1]) {
        mya.test[[r]][[ match(i,myidx[,1]) ]] <- pairwise.wilcox.test(myalpha[myalpha["Replicate"]==r,i], myalpha[myalpha["Replicate"]==r,"Population"], p.adjust.method="none")
    }

    mya.letters[[r]] <- vector("list", length(myidx[,1]))
    for (i in 1:length(mya.test[[r]])) {
        d <- as.vector(mya.test[[r]][[i]]$p.value)
        names(d) <- sapply(colnames(mya.test[[r]][[i]]$p.value), function(x) paste0(x, "-", rownames(mya.test[[r]][[i]]$p.value))) %>% as.vector()
        d <- d[!is.na(d)]
        mya.letters[[r]][[i]] <- generate.label.df(d)
    }
}


#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

p <- vector("list", length(nb.rep))
for (r in nb.rep) {
    p[[r]] <- vector("list", length(myidx[,1]))
    for (i in myidx[,1]) {
        df.tmp <- data.frame(i=myalpha[myalpha["Replicate"]==r,i], Population=myalpha[myalpha["Replicate"]==r,"Population"])

        p[[r]][[match(i,myidx)]] <- ggplot(df.tmp, aes(x=Population, y=i, fill=Population)) +
                geom_blank(data=df.tmp, aes(x=Population, y=i*1.1)) +
                geom_boxplot() + xlab("") + ylab("") + ggtitle(myidx[myidx[,1] == i,2]) +
                stat_summary(geom = 'text', label = mya.letters[[r]][[match(i,myidx[,1])]][,1], fun.y = max, vjust = -1) +
                scale_fill_manual(values=pop.data[,2]) +
                theme(plot.title = element_text(hjust =0.5), legend.position="none")

    }
}

p.ls <- vector("list", 2)
p.ls[[1]] <- ggarrange(plotlist=p[[1]], labels=LETTERS[1:length(p[[r]])], ncol=2, nrow=ceiling(length(p[[r]])/2))

p.ls[[2]] <- ggarrange(ggarrange(plotlist=p[[1]], labels=LETTERS[1], ncol=2, nrow=ceiling(length(p[[r]])/2)),
                       as_ggplot(grid.lines(c(0.5,0.5), c(0,1), gp=gpar(lwd=1.5, col="black"))),
                       ggarrange(plotlist=p[[2]], labels=LETTERS[2], ncol=2, nrow=ceiling(length(p[[r]])/2)),
                       ncol=length(nb.rep)+1, nrow=1, widths=c(1,0.01,1))

ggsave(paste0(graph.d,"Fig. 2 - alpha-div.pdf"), p.ls[[1]], width=10, height=10, useDingbats=FALSE)
ggsave(paste0(graph.d,"alpha-div_all.pdf"), p.ls[[2]], width=10*length(nb.rep), height=10, useDingbats=FALSE)

# Correlation between replicates
p <- vector("list", length(myidx[,1]))
for (i in myidx[,1]) {
    df.tmp <- spread(myalpha[,c(i, "ReplicateID", "Replicate", "Population")] , "Replicate", i)
    colnames(df.tmp)[3:4] <- c("V1", "V2")

    p[[match(i,myidx)]] <- ggplot(df.tmp, aes(x=V1, y=V2)) +
                geom_point(aes(color=Population)) + ggtitle(myidx[myidx[,1] == i,2]) +
                stat_smooth(method=lm, se=FALSE) +
                xlab("Replicate 1") + ylab("Replicate 2") +
                scale_color_manual(values=pop.data[,2]) +
                theme(plot.title = element_text(hjust = 0.5), legend.position="none") +
                stat_cor(method="pearson", label.x.npc="left", label.y.npc="top", hjust=0)
}

p <- ggarrange(plotlist=p, ncol=2, nrow=ceiling(length(p)/2), common.legend=TRUE, legend="bottom", align="hv")
ggsave(paste0(graph.d,"Supp. Fig. 3 - alpha-div_cor.pdf"), p, width=10, height=10, useDingbats=FALSE)


#----------------#
# Beta diversity #
#----------------#

# Distance indices
mydist <- c("jaccard", "Unweighted UniFrac", "bray", "Weighted UniFrac")

# Lists for saving statistics results
mybiom.r.cdis   <- vector("list", length(nb.rep))
mybiom.r.pcoa   <- vector("list", length(nb.rep))
mybiom.r.stat   <- vector("list", length(nb.rep))
mybiom.r.stat.p <- vector("list", length(nb.rep))
mybiom.r.stat.b <- vector("list", length(nb.rep))
mybiom.r.stat.T <- vector("list", length(nb.rep))
mybiom.r.stat.m <- vector("list", length(nb.rep))

# Statistical tests for each replicate
for (r in nb.rep) {
    # Data
    mybiom.tmp <- subset_samples(mybiom, Replicate == r)
    mybiom.r   <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose=FALSE)

    # Lists for saving results
    mybiom.r.cdis[[r]]   <- vector("list", length(mydist))
    mybiom.r.pcoa[[r]]   <- vector("list", length(mydist))
    mybiom.r.stat[[r]]   <- vector("list", length(mydist))
    mybiom.r.stat.p[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.b[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.T[[r]] <- vector("list", length(mydist))
    mybiom.r.stat.m[[r]] <- vector("list", length(mydist))

    # Test for each distance index
    for (i in mydist) {
        myls.idx <- match(i, mydist)

        if (i == "jaccard") {
            mydist.tmp <- distance(mybiom.r, i, binary=TRUE)
        } else {
            mydist.tmp <- distance(mybiom.r, i)
        }

        # Centroid distance between popuulations
        mybiom.r.cdis[[r]][[myls.idx]] <- dist_multi_centroids(mydist.tmp, sample_data(mybiom.r)$Population)
        
        # PCoA
        mybiom.r.pcoa[[r]][[myls.idx]] <- ordinate(mybiom.r, "PCoA", mydist.tmp)

        # Permanova stat
        ## Reorder table and select only the current samples
        mypop.tmp <- mypop[ order(mypop[,"Population"]), ]
        mypop.tmp <- mypop.tmp[ mypop.tmp[,1] %in% labels(mydist.tmp),]

        ## Reorder the distance matrix
        mydist.tmp2 <- as.matrix(mydist.tmp)[ as.vector(mypop.tmp[,1]), as.vector(mypop.tmp[,1])] %>% as.dist()

        ## Run test
        mybiom.r.stat[[r]][[myls.idx]]   <- adonis(mydist.tmp ~ Population, data = get_variable(mybiom.r), permutations = 1000)
        set.seed(myseed)
        mybiom.r.stat.p[[r]][[myls.idx]] <- pairwise.adonis(mydist.tmp2, mypop.tmp[,"Population"], perm = 1000)

        # Dispersion stat
        mybiom.r.stat.b[[r]][[myls.idx]] <- betadisper(mydist.tmp, get_variable(mybiom.r)[,"Population"], type = "centroid")
        mybiom.r.stat.T[[r]][[myls.idx]] <- TukeyHSD(mybiom.r.stat.b[[r]][[myls.idx]])

        # P-values matrix
        mymat <- matrix(NA, ncol=nrow(pop.data), nrow=nrow(pop.data))
        colnames(mymat) <- rownames(mymat) <- pop.data[,1]
        mymat[lower.tri(mymat)] <- mybiom.r.stat.p[[r]][[myls.idx]][,6]
        mymat <- t(mymat)
        mymat[lower.tri(mymat)] <- mybiom.r.stat.T[[r]][[myls.idx]]$group[,4]
        mybiom.r.stat.m[[r]][[myls.idx]] <- mymat
    }

    # Naming slots
    names(mybiom.r.cdis[[r]])   <- mydist
    names(mybiom.r.pcoa[[r]])   <- mydist
    names(mybiom.r.stat[[r]])   <- mydist
    names(mybiom.r.stat.p[[r]]) <- mydist
    names(mybiom.r.stat.b[[r]]) <- mydist
    names(mybiom.r.stat.T[[r]]) <- mydist
    names(mybiom.r.stat.m[[r]]) <- mydist
}


# Centroid distances between species and populations

## List species and remove water
sp.vec <- sample_data(mybiom)$Species %>% as.character() %>% unique()
sp.vec <- sp.vec[! sp.vec %in% "Water"]

## List samples per species
spl.ls        <- sapply(sp.vec, function(x) sample_data(mybiom)[ sample_data(mybiom)$Species == x, "Population"] %>% unique())
names(spl.ls) <- sp.vec

for (r in nb.rep) {
    cat("\nCentroid distances for replicate ", r, ":\n", sep="")

    for(i in mydist) {
        mymt <- mybiom.r.cdis[[r]][[i]] %>% as.matrix()

        # Add NA
        mymt[upper.tri(mymt, diag = TRUE)] <- NA

        # Remove water
        myspl <- ! rownames(mymt) %in% "Water"
        mymt  <- mymt[myspl, myspl]

        # Test selection
        norm <- lapply(spl.ls, function(x) shapiro.test(as.vector(mymt[, x]))$p.value > 0.05) %>% unlist() %>% all()

        # Test
        if (norm) {
            p.val <- t.test(mymt[, spl.ls[[1]]], mymt[, spl.ls[[2]]])$p.value %>% signif(., 2)
            test  <- "T test"
        } else {
            p.val <- wilcox.test(mymt[, spl.ls[[1]]], mymt[, spl.ls[[2]]])$p.value %>% signif(., 2)
            test  <- "Wilcoxon test"
        }

        # Average distance
        average.d <- lapply(spl.ls, function(x) mean(mymt[, x], na.rm=TRUE) %>% signif(., 2))

        cat("\t- ", tools::toTitleCase(i), " distance: the average distance between Ba and Bg is ", average.d[[1]], " and the average distance within Bg is ", average.d[[2]], " (", test, "; p-value :", p.val, ").\n", sep="")

    }
}

#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Axes to plot
myax <- list(c(1,2), c(1,3), c(2,3))

# List to save plots
p.ls.r <- vector("list", length(nb.rep))

for (r in nb.rep) {

    # Data
    mybiom.tmp <- subset_samples(mybiom, Replicate == r)
    mybiom.r <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose=FALSE)

    p.ls.r[[r]] <- vector("list", length(mydist) * length(myax))

    k <- 0
    for (i in 1:length(mydist)) {
        for (j in 1:length(myax)) {
            mytitle <- ""
            if (j == 1) { mytitle <- tools::toTitleCase(mydist[i]) }

            k <- k + 1
            p.ls.r[[r]][[k]] <- plot_ordination(mybiom.r, mybiom.r.pcoa[[r]][[i]], color = "Population", axes=myax[[j]]) +
                           scale_color_manual(values=pop.data[,2]) +
                           stat_ellipse(type = "norm") +
                           ggtitle(mytitle) +
                           theme(axis.text = element_text(size = rel(0.45))) +
                           theme(axis.title = element_text(size = rel(0.6))) +
                           theme(legend.position="none")
        }
    }
}

myletters <- lapply(1:length(mydist), function(x) c(LETTERS[x], rep("",length(myax)-1))) %>% unlist()

mylg <- length(p.ls.r[[1]])
p <- ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom", labels=myletters)

ggsave(paste0(graph.d,"b-div.pdf"), p, width=3.5*3, height=3.5*length(p.ls.r[[1]])/3, useDingbats=FALSE)

# Panel with replicates
p <- ggarrange(ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
               as_ggplot(grid.lines(c(0.5,0.5), c(0,1), gp=gpar(lwd=1.5, col="black"))),
               ggarrange(plotlist=p.ls.r[[2]], nrow=length(p.ls.r[[2]])/3, ncol=3, common.legend=TRUE, legend="bottom"),
               ncol=length(nb.rep)+1, nrow=1, widths=c(1,0.01,1), labels=c(LETTERS[1], "", LETTERS[2]))

ggsave(paste0(graph.d,"Supp. Fig. 4 - b-div_all.pdf"), p, width=3.5*3*2, height=3.5*length(p.ls.r[[1]])/3, useDingbats=FALSE)


# List to save plots
p.ls.r <- vector("list", length(nb.rep))
p.ls.t <- vector("list", length(nb.rep))

for (r in nb.rep) {

    # Data
    mybiom.tmp <- subset_samples(mybiom, Replicate == r)
    mybiom.r   <- rarefy_even_depth(mybiom.tmp, rngseed=myseed, verbose=FALSE)

    p.ls.r[[r]] <- vector("list", length(mydist))
    p.ls.t[[r]] <- vector("list", length(mydist)) 

    for (i in 1:length(mydist)) {
        mytitle <- tools::toTitleCase(mydist[i])

        p.ls.r[[r]][[i]] <- plot_ordination(mybiom.r, mybiom.r.pcoa[[r]][[i]], color = "Population", axes=myax[[1]]) +
                       scale_color_manual(values=pop.data[,2]) +
                       stat_ellipse(type = "norm") +
                       ggtitle(mytitle) +
                       theme(legend.position="none")

        mymat <- mybiom.r.stat.m[[r]][[i]]
        mybreaks <- c(10^(seq(log10(min(mymat, na.rm=T)), log10(0.05), length.out=4)), 1)

        p.ls.t[[r]][[i]] <- plotPvalues(mymat, border.col="black", pval.col="black", font.size=3) +
                    xlab(expression(paste(beta,"-dispersion test"))) + ylab("Permanova test") +
                    scale_fill_gradient2(low = "red", mid = "yellow", high = "white", limits = c(NA, 0.05), guide = "colorbar", na.value = "transparent", trans = "log10", breaks = mybreaks, label = format(mybreaks, digit = 2, scientific = TRUE)) +
                    guides(fill = guide_colourbar(barwidth = 9, title.vjust = 0.75)) +
                    theme(axis.line = element_blank(), axis.ticks = element_blank(), legend.position = "none")
    }
}

mylg <- length(p.ls.r[[1]])
p <- ggarrange(ggarrange(plotlist=p.ls.r[[1]], nrow=length(p.ls.r[[1]]), ncol=1, common.legend=TRUE, legend="bottom", labels=LETTERS[1:mylg]),
               ggarrange(plotlist=p.ls.t[[1]], nrow=length(p.ls.t[[1]]), ncol=1, common.legend=TRUE, legend="bottom"),
               ncol=2, nrow=1)

ggsave(paste0(graph.d,"Fig. 3 - b-div.pdf"), p, width=5*2, height=3.5*length(p.ls.r[[1]]), useDingbats=FALSE)


#---------------------#
# Taxonomic diversity #
#---------------------#

# Phyla only
mybiom.P <- mytax_glom(mybiom, "Phylum")

# Dominant OTUs
otu.dom <- vector("list", length(nb.rep))
for (r in nb.rep) {
    mybiom.tmp <- subset_samples(mybiom.P, Replicate == r)

    otu.dom[[r]] <- vector("list", nrow(pop.data))
    names(otu.dom[[r]]) <- pop.data[,1]

    for (p in pop.data[,1]) {
        myind <- rownames(sample_data(mybiom.tmp))[sample_data(mybiom.tmp)$Population == p]
        otu.dom[[r]][[p]] <- tax_table(mybiom.P)[,2][otu_table(mybiom.tmp)[,myind] %>% rowSums() %>% order() %>% tail(., 1)] %>% as.vector()

        # Message
        cat(p, " population from replicate ", r, " is dominated by ", otu.dom[[r]][[p]], ".\n", sep="")
    }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Differentially represented phyla #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Comparison matrices
mycomp.pop <- sample_data(mybiom)$Population %>% as.character() %>% unique() %>% combn(., 2)
mycomp.sp  <- sample_data(mybiom)$Species %>% as.character() %>% unique() %>% combn(., 2)

# Lists to store statistical results
mybiom.tax.pop <- vector("list", length(nb.rep))
mybiom.tax.sp  <- vector("list", length(nb.rep))

# Statistical tests for each replicate
for (r in nb.rep) {

    # Subset by replicate
    mybiom.tmp <- subset_samples(mybiom.P, Replicate == r)

    ### Populations ###

    # Data conversion
    pop.dds <- phyloseq_to_deseq2(mybiom.tmp, ~ Population)

    # Counts
    pop.dds <- DESeq(pop.dds, test="Wald", fitType="parametric", parallel=TRUE)

    # Lists to store results
    res        <- vector("list", ncol(mycomp.pop))
    names(res) <- apply(mycomp.pop, 2, function(x) paste(x[1], "vs", x[2]))
    sig.tb     <- res

    for (i in 1:ncol(mycomp.pop)) {
        res[[i]] <- results(pop.dds, contrast=c("Population", mycomp.pop[1,i], mycomp.pop[2,i]), cooksCutoff=FALSE, parallel=TRUE)

        sig.tb.tmp <- res[[i]][which(res[[i]]$padj < myalpha.lim), ]

        if (nrow(sig.tb.tmp) > 0) {
            sig.tb[[i]] <- cbind(tax_table(mybiom.tmp)[rownames(sig.tb.tmp),], sig.tb.tmp)
        } else {
            sig.tb[[i]] <- as.matrix(sig.tb.tmp)
        }
    }

    # Save results
    mybiom.tax.pop[[r]] <- list(res, sig.tb)

    ### Species ###

    # Data conversion
    pop.dds <- phyloseq_to_deseq2(mybiom.tmp, ~ Species)

    # Counts
    pop.dds <- DESeq(pop.dds, test="Wald", fitType="parametric", parallel=TRUE)

    # Lists to store results
    res        <- vector("list", ncol(mycomp.sp))
    names(res) <- apply(mycomp.sp, 2, function(x) paste(x[1], "vs", x[2]))
    sig.tb     <- res

    for (i in 1:ncol(mycomp.sp)) {
        res[[i]] <- results(pop.dds, contrast=c("Species", mycomp.sp[1,i], mycomp.sp[2,i]), cooksCutoff=FALSE, parallel=TRUE)

        sig.tb.tmp <- res[[i]][which(res[[i]]$padj < myalpha.lim), ]

        if (nrow(sig.tb.tmp) > 0) {
            sig.tb[[i]] <- cbind(tax_table(mybiom.tmp)[rownames(sig.tb.tmp),], sig.tb.tmp)
        } else {
            sig.tb[[i]] <- as.matrix(sig.tb.tmp)
        }
    }

    # Save results
    mybiom.tax.sp[[r]] <- list(res, sig.tb)
}

# Export results
exp.d.m <- "differentially represented phyla"
cat("\nExporting differentially represent phyla in ", exp.d.m, "...\n\n", sep="")
## Populations
for (i in 1:length(mybiom.tax.pop)) {
    exp.d <- paste0(res.d, exp.d.m, "/Populations/rep_", i, "/")
    if (! dir.exists(exp.d)) {dir.create(exp.d, recursive=TRUE)}

    mynames <- names(mybiom.tax.pop[[i]][[2]])
    for (j in 1:length(mybiom.tax.pop[[i]][[2]])) {
        write.table(data.frame(mybiom.tax.pop[[i]][[2]][[j]]), file=paste0(exp.d, mynames[j], ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
}

## Species
for (i in 1:length(mybiom.tax.sp)) {
    exp.d <- paste0(res.d, exp.d.m, "/Species/rep_", i, "/")
    if (! dir.exists(exp.d)) {dir.create(exp.d, recursive=TRUE)}

    mynames <- names(mybiom.tax.sp[[i]][[2]])
    for (j in 1:length(mybiom.tax.sp[[i]][[2]])) {
        write.table(data.frame(mybiom.tax.sp[[i]][[2]][[j]]), file=paste0(exp.d, mynames[j], ".tsv"), sep="\t", quote=FALSE, row.names=FALSE)
    }
}

#~~~~~~~~~#
# Cluster #
#~~~~~~~~~#

# Cluster and bootstrap
mybiom.d <- apply(otu_table(mybiom), 1, function(x){ x <- x+1; log(x) - mean(log(x))})
mypvc    <- pvclust(t(mybiom.d), method.hclust="ward.D2", method.dist="euclidian", parallel=TRUE, iseed=myseed)

# Rotate water cluster
myordr <- c(grep("Water", label(mypvc$hclust)[mypvc$hclust$order], invert=TRUE, value=TRUE), grep("Water", label(mypvc$hclust)[mypvc$hclust$order], value=TRUE))
mypvc$hclust <- dendextend::rotate(mypvc$hclust, myordr)

# Bootstrap table
mybs <- pvclust_edges(mypvc)


#~~~~~~~~~#
# Figures #
#~~~~~~~~~#

# Data storage for graphs
p.ls   <- vector("list", length(nb.rep))
pop.wd <- vector("list", length(nb.rep))

# Color
set.seed(myseed)
mybiom.P <- mytax_glom(mybiom, "Phylum") %>% core(., detection=0, prevalence=0)
myclr.P <- tax_table(mybiom.P) %>% nrow() %>% distinctColorPalette(., seed=myseed - 223) 
names(myclr.P) <- tax_table(mybiom.P)[,2] %>% sort()

# Figure for each replicate
for (r in nb.rep) {

    if (r ==1 ) {
        # Remove replicates
        mybiom.tmp <- subset_samples(mybiom, Replicate == r)
        mypvc.tmp  <- mypvc
        mypvc.tmp$hclust <- dendextend::prune(mypvc$hclust, grep(".2$", mypvc$hclust$labels, value=TRUE)) ## Remove phylogram

        # Update bootstrap table
        mybs.tmp <- mybs[ ! (grepl("\\.2$", mybs[,1]) | grepl("\\.2$", mybs[,2])), ]
        mybs.tmp <- mybs.tmp[ -1, ]
        mybs.tmp <- cbind(mybs.tmp, pvclust:::hc2axes(mypvc.tmp$hclust))

        # Update pvclust object
        mypvc.tmp$edges <- mypvc.tmp$edges[rownames(mybs.tmp),]
        
        # Dendrogram conversion
        mydd.tmp <- as.dendrogram(mypvc.tmp$hclust)

        myhj <- 1
    } else {
        mybiom.tmp <- mybiom
        mybs.tmp   <- mybs
        mybs.tmp   <- cbind(mybs.tmp, pvclust:::hc2axes(mypvc$hclust))
        mydd.tmp   <- as.dendrogram(mypvc$hclust)
        myhj       <- 0.5
    }

    # Remove bootstrap of the first node
    mybs.tmp <- mybs.tmp[ -nrow(mybs.tmp), ]

    # Remove non present taxa
    mybiom.P <- mytax_glom(mybiom.tmp, "Phylum")
    mybiom.P <- core(mybiom.P, detection=0, prevalence=0)

    # Sample biome
    mybiom.P.s <- transform(mybiom.P, transform="compositional")

    # Population biome
    mybiom.P.p <- mysample_glom(mybiom.P, "Population")
    mybiom.P.p <- transform(mybiom.P.p, transform="compositional")

    # Color vector
    myclr <- myclr.P[ tax_table(mybiom.P)[,2] ]

    # List to store graphs
    p <- vector("list", 3)

    # Tree
    mylabels <- as.vector(dendro_data(mydd.tmp)$labels$label)
    p[[1]] <- ggdendrogram(mydd.tmp) +
        scale_x_continuous(expand = c(1/(42*2),0), labels = mylabels, breaks = 1:length(mylabels)) +
        annotate(geom = "text", x = mybs.tmp[,"x.axis"], y = mybs.tmp[,"y.axis"], label = signif(mybs.tmp[,"bp"],2)*100, hjust = 0.5, vjust = 1.2, size = 2, col = "grey70") +
        theme(axis.text.y = element_blank(), plot.margin = margin(b = -10))

    # Barplot at sample level
    mybiom.P.m <- psmelt(mybiom.P.s)
    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]
    p[[2]] <- ggplot(mybiom.P.m, aes(x = Sample, y = Abundance, fill = Phylum)) + 
        geom_bar(stat = "identity", color = "black", size = 0.05) +
        scale_fill_manual(values = myclr) +
        scale_x_discrete(limits = mylabels) +
        ## Update legend look (because of the bar borders
        guides(fill = guide_legend(override.aes = list(colour = NA, size = 0.5))) +
        ## Remove x axis title and legend
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), legend.position = "none") + 
        ylab("Relative Abundance\n")


    # Barplot at population level
    mybiom.P.m <- psmelt(mybiom.P.p)
    mybiom.P.m <- mybiom.P.m[ order(mybiom.P.m$Abundance, decreasing=TRUE), ]

    ## Number of samples per population
    spl.count   <- plyr::count(sample_data(mybiom.P)[,2])
    pop.ls      <- as.vector(mybiom.P.m$Population) %>% sapply(., function(x) spl.count[ spl.count[,1] == x, 2])
    pop.wd[[r]] <- pop.ls

    ## Correction factor to make sure to have odd number under bars
    if (r == 1) myftr <- 1
    if (r == 2) myftr <- 2

    ## Break vector
    mylabels.b <- rep("", nsamples(mybiom.tmp) / myftr)

    ## Filling appropriate breaks with labels
    mylabels.u  <- sapply(as.vector(mylabels), function(x) sample_data(mybiom.tmp)[x, "Population"]) %>% unlist() %>% as.vector() %>% unique()
    j <- 0
    for (l in mylabels.u) {
        i <- spl.count[ spl.count[,1] == l,2]/2 / myftr
        mylabels.b[ j+ceiling(i) ] <- l
        j <- j+i*2
    }

    ## Width of the population bar
    pop.wd[[r]] <- pop.wd[[r]]/myftr

    p[[3]] <- ggplot(mybiom.P.m, aes(x = Population, y = Abundance, fill = Phylum, width = pop.wd[[!!r]]*0.98)) +   # Use of quasiquotation !!r otherwise lastest slot of the list taken for plotting
        geom_bar(stat = "identity", color = "black", size =0.05) +
        scale_x_discrete(limits=mylabels.b, labels=mylabels.b) +
        scale_fill_manual(values = myclr) +
        ## Update legend look (because of the bar borders
        guides(fill = guide_legend(override.aes = list(colour = NA, size = 0.5))) +
        ## Remove x axis title and legend
        theme(axis.title.x = element_blank(), legend.position = "none", axis.ticks = element_line(linetype = c("blank", rep("solid", length(mylabels.u))))) +     # axis.ticks needed to remove the first tick
        ylab("Relative Abundance\n")

    # Saving all plots
    p.ls[[r]] <- ggarrange(plotlist = p, labels = LETTERS[1:length(p[[1]])], ncol = 1, common.legend = TRUE, legend = "bottom", align = "v")
}

ggsave(paste0(graph.d,"Fig. 4 - taxo-div.pdf"), p.ls[[1]], width=10, height=9, useDingbats=FALSE)
ggsave(paste0(graph.d,"Supp. Fig. 6 - taxo-div_all.pdf"), p.ls[[2]], width=15, height=9, useDingbats=FALSE)


#-----------------#
# Unassigned taxa #
#-----------------#

# Proportion of unassigned
unass.t <- tax_table(mybiom)[,1] == "Unassigned"
unass.p <- apply(otu_table(mybiom), 2, function(x) sum(x[unass.t]) / sum(x))

snail.u <- unass.p[ ! grepl("Water", names(unass.p))]
water.u <- unass.p[ grepl("Water", names(unass.p))]

cat("\nUnassigned bacteria proportion in snails:\n\t- mean:", mean(snail.u), "\n\t- range (min max):", range(snail.u), "\n")
cat("\nUnassigned bacteria proportion in water:\n\t- mean:",  mean(water.u), "\n\t- range (min max):", range(water.u), "\n\n")


#---------#
# Heatmap #
#---------#

p.ls <- vector("list", length(nb.rep))
for (r in nb.rep) {

    if (r == 1 ) {
        mybiom.tmp <- subset_samples(mybiom, Replicate == r)
    } else {
        mybiom.tmp <- mybiom
    }

    myotu.tb <- otu_table(mybiom.tmp)

    # Identify sample in each populations
    mygrp    <- split(get_variable(mybiom.tmp), get_variable(mybiom.tmp, "Population")) %>% lapply(., rownames)
    mygrp.ls <- mygrp %>% lapply(., function(x) match(x, colnames(myotu.tb)))
    mygrp.ls[[length(mygrp.ls) + 1 ]] <- do.call(c, mygrp.ls[1:7]) %>% unname()

    # List specific ASVs
    p.ls.otu <- vector("list", length(mygrp))
    for (i in 1:length(mygrp.ls)) {
        mycln  <- mygrp.ls[[i]]
        row.w  <- rowMeans(myotu.tb[, mycln]) > 50
        row.wo <- rowSums(myotu.tb[, -mycln]) == 0
        row.sl <- rowSums(cbind(row.w, row.wo)) == 2

        p.ls.otu[[i]] <- myotu.tb[row.sl, ]
    }

    # Convert list to data.frame
    p.ls.otu <- do.call(rbind, p.ls.otu)

    # Reduce to unique ASVs
    mygrp.otu.uniq <- rownames(p.ls.otu) %>% unique()

    # Prune table
    mybiom.o  <- prune_taxa(mygrp.otu.uniq, mybiom.tmp)

    # Identify interval to add blank bars
    spl.count <- plyr::count(sample_data(mybiom.tmp)[,2])
    myint     <- cumsum(spl.count[-nrow(spl.count),2])

    # Plot
    p.ls[[r]] <- plot_heatmap(mybiom.o, method=NULL, distance="bray", sample.order=unlist(mygrp), taxa.order=mygrp.otu.uniq) +
                    geom_vline(xintercept = myint+0.5, size = 2, color="white") +
                    labs(y="ASV") +
                    theme(axis.ticks.y=element_blank(), axis.line=element_blank(), axis.text.y=element_blank())
}

ggsave(paste0(graph.d,"taxo-div_hm.pdf"), p.ls[[1]], width=10, height=20, useDingbats=FALSE)
ggsave(paste0(graph.d,"Supp. Fig. 5 - taxo-div_hm_all.pdf"), p.ls[[2]], width=10, height=20, useDingbats=FALSE)

rare.taxa.count <- p.ls[[2]]$data[,1] %>% unique() %>% length()
rare.taxa.range <- p.ls[[2]]$data[,3] %>% .[. > 0] %>% range()
rare.taxa.prop  <- (sapply(p.ls[[2]]$data[,1] %>% unique(), function(x) sum(p.ls[[2]]$data[p.ls[[2]]$data[, 1] == x, 3] > 0) / length(p.ls[[2]]$data[p.ls[[2]]$data[, 1] == x, 3]) ) %>% unlist() %>% range() * 100) %>% round(., digit=3) %>% paste0(., "%")

cat("\nRare bacteria in snails and water:\n\t- total number of taxa:", rare.taxa.count, "\n\t- Count range (min max):", rare.taxa.range, "\n\t- Prevalence range (min max):", rare.taxa.prop, "\n")

# Clean tmp file
if (file.exists("Rplots.pdf")) { unlink("Rplots.pdf") }
