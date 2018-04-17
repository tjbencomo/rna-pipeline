#!/usr/bin/env Rscript

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
res

sum(res$pvalue < .05, na.rm = TRUE)
sum(res$padj < .05, na.rm = TRUE)
sum(res$padj < .8, na.rm = TRUE)

print("Minimum pvalues, first for pvalue and then for padj")

min(res$pvalue, na.rm = TRUE)
min(res$padj, na.rm = TRUE)


#Analysis finished. Filter and save results
#First, write to csv file an unfiltered output list
resOrdered <- res[order(res$pvalue), ]
write.csv(as.data.frame(resOrdered), file='unfiltered_normal_tumor_results.csv')


