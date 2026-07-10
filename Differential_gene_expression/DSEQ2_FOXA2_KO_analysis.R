library("pasilla")
library("tidyverse")
library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("apeglm")
library("DEGreport")
library("tximport")
library("readr") 


files <- list.files("/Users/j.roberts12/Desktop/FOXA2_KO_paper", pattern = "genes.results$", full.names = TRUE)

names(files) <- basename(files)

# https://support.bioconductor.org/p/92763/ 
# Reads the rsem expected count column and the effective length column
txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)

# editing 0 lengths to be 1 
txi.rsem$length[txi$length == 0] <- 1 #https://support.bioconductor.org/p/92763/

# Define the Genotype based on the file namesv - this sorts the levels for the condition 
# factor in the sampleTable
genotype <- c("Control", "Control", "FOX_KO", "FOX_KO", "FOX_KO", 
              "Control", "Control", "FOX_KO", "FOX_KO", "FOX_KO")

# Also defining the levels for the factor of the stage 
stage <- c(rep("EP", 5), rep("PP2", 5))

# Create the sample information table and using row.names = names(files) 
sampleTable <- data.frame(
  condition = factor(genotype, levels = c("Control", "FOX_KO")),
  stage = factor(stage),
  row.names = names(files) # This ensures names match your txi object
)

#-------It is critical that the columns of the count matrix and the 
# rows of the column data (sample information) are in the same order ---------

rownames(sampleTable)

colnames(txi.rsem$counts)

# first check if the row names in sampleTable are the same as the column names of txi.rsem$counts data
all(rownames(sampleTable) %in% colnames(txi.rsem$counts))

# check if order is the same
all(rownames(sampleTable) == colnames(txi.rsem$counts))

# constructing DESeqDataSet
dds <- DESeqDataSetFromTximport(txi = txi.rsem,
                              colData = sampleTable,
                              design = ~ stage + condition)
# Pre-filtering - this is not necessary but removing rows in which there are very few reads, 
# we reduce the memory size of the dds data object, and we increase the speed of the 
# transformation and testing functions within DESeq2

# Note that more strict filtering to increase power is automatically applied via independent 
# filtering on the mean of normalized counts within the results function

keep <- rowSums(counts(dds)) >=10

dds <- dds[keep,]

# By default R will choose a reference level for factors based on alphabetical order 
# If you dont tell DESeq2 functions which level you want to compare against, the ocmaprisons will be based on the alphabetical order of the levels 

#Here we set the factor levels explicitly - https://introtogenomics.readthedocs.io/en/latest/2021.11.11.DeseqTutorial.html

dds$condition <- factor(dds$condition, levels = c("Control", "FOX_KO"))

dds$stage <- factor(dds$stage, levels = c("EP", "PP2"))

# Running the function - we can look at the differential expressed functions using the 
# Results() function. 

#By default DESeq2 perform pair-wise comparison of the first and the last variable in the experimental design variables and provide a result table

dds <- DESeq(dds)

res <- results(dds)

res

# Note that we could have specified the coefficient or contrast we want to build a results table for, 
# using either of the following equivalent commands:

res_stage <- results(dds, name="stage_PP2_vs_EP")
res_condition_FOXKOvscontrol <- results(dds, contrast = c("condition", "FOX_KO", "Control"))

resultsNames(dds)

resLFC <- lfcShrink(dds, coef="condition_FOX_KO_vs_Control", type="apeglm")

resLFC

resOrdered <- res[order(res$pvalue),]

summary(res)

res05 <- results(dds, alpha=0.05)

summary(res05)

sum(res05$padj < 0.05, na.rm=TRUE)

plotMA(res, ylim=c(-2,2))

plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

d <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", 
                returnData=TRUE)

library("ggplot2")
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

mcols(res)$description

vsd <- vst(dds)
head(assay(vsd), 3)


