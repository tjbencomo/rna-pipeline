#!/usr/bin/env Rscript


# Author: Tomas Bencomo
# R script to analyze gene expression data produced by HTSeq-count.
# Final step of the rna-pipeline
# Meant to be run after all HTSeq-count jobs have finished for a given collection of samples
# Script performs the following:
# 1) Checks if DESeq2 installed - if not, installs it
# 2) Parse command line to get the directory of the counts data, the list of sample conditions, and the directory to write output
# 3) Analyze RNA expression and write 3 output files:
# 	-Unfiltered output list with all data
#	-Filtered down-regulated file, filtered by pvalue <= .05
#	-Filtered up-regulated file, filted by pvalue <= .05

if(!require(DESeq2)) {
	print("DESeq2 not installed, installing now")
	source("https://bioconductor.org/biocLite.R")
	y #no idea if this will work - needed to answer prompts when running source()
	y
	biocLite("DESeq2")
	library(DESeq2)
} else {
	print("DESeq2 installed, continuing")
}
if(!require(AnnotationDbi)) {
	print("AnnotationDbi not installed, installing now")
        source("https://bioconductor.org/biocLite.R")
        #idk if this is needed - no idea if this will work - needed to answer prompts when running source()
        biocLite("AnnotationDbi")
        library(AnnotationDbi)
} else {
	print("AnnotationDbi installed, continuing")
}
if(!require(org.Hs.eg.db)) {
	print("org.Hs.eg.db not installed, installing now")
        source("https://bioconductor.org/biocLite.R")
        #idk if this is needed - no idea if this will work - needed to answer prompts when running source()
        biocLite("org.Hs.eg.db")
        library(org.Hs.eg.db)
} else {
	print("org.Hs.eg.db installed, continuing")
}
if(!require(dplyr)) {
	print("dplyr not installed, installing now")
        #idk if this is needed - no idea if this will work - needed to answer prompts when running source()
        packages.install("dplyr")
        library(dplyr)
} else {
	print("dplyr installed, continuing")
}





#Args we are expecting
# --countsDir directory with htseq counts files
# --sampleTypes list of sample conditions separated spaces
# --resultsDir directory to output results

args = commandArgs(trailingOnly=TRUE)

countsDir <- args[match("--countsDir", args) + 1]
directory <- countsDir

sampleFiles <- grep("*counts*", list.files(directory), value=TRUE)
print(sampleFiles)


sampleTypesIndex <- match("--sampleTypes", args) + 1
sampleCondition <- args[sampleTypesIndex : (sampleTypesIndex + length(sampleFiles) - 1)]

resultsDir = args[match("--resultsDir", args) + 1]
setwd(resultsDir)

sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = sampleCondition)


#Build DESeq2 data object and run analysis
#create data object
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~condition)

#display info about data object
ddsHTSeq

#ensure factor levels are correct
ddsHTSeq$condition <- factor(ddsHTSeq$condition, levels = c("Normal","Tumor"))

#call DESeq analysis object and store results in res
dds <- DESeq(ddsHTSeq)
res <- results(dds)

gene_ids <- rownames(res)
gene_id_no_isoforms <- sapply(gene_ids, substr, 0, which(strsplit(gene_ids, "")[[1]]==".")-1)

res$ENSEMBL_CODE <- gene_id_no_isoforms
res$symbol <- mapIds(org.Hs.eg.db, keys=res$ENSEMBL_CODE, column="SYMBOL", keytype = "ENSEMBL", multiVals = "first")


sum(res$pvalue < .05, na.rm = TRUE)
sum(res$padj < .05, na.rm = TRUE)

print("Minimum pvalues, first for pvalue and then for padj")

min(res$pvalue, na.rm = TRUE)
min(res$padj, na.rm = TRUE)


#Analysis finished. Filter and save results
#First, write to csv file an unfiltered output list
resOrdered <- res[order(res$pvalue), ]
resOrdered <- data.frame(resOrdered)
resOrdered <- resOrdered[, c(8, 7, 1, 2, 3, 4, 5, 6)]
write.csv(as.data.frame(resOrdered), file='unfiltered_normal_tumor_results.csv', row.names = FALSE)

#Now write upregulated file with padj <= .05
upregulated <- filter(resOrdered, log2FoldChange > 0 & padj <= .05)
head(upregulated)
dim(upregulated)
write.csv(upregulated, file='upregulated_genes.csv', row.names = FALSE)


#Next write downregulated file with padj <= .05
downregulated <- filter(resOrdered, log2FoldChange < 0 & padj <= .05)
head(downregulated)
dim(downregulated)
write.csv(downregulated, file = 'downregulated_genes.csv', row.names = FALSE)




