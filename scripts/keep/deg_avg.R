# Based on https://learn.gencore.bio.nyu.edu/rna.seq.analysis/deseq.2/
# Install the latest version of DEseq2

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DESeq2", version = "3.8")
#BiocManager::install("stringr")

qvalue=0.25
# load the library
library(DESeq2)
library(stringr)

# Read in the raw read counts
rawCounts <- read.table("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/dedup_featurescounts.tsv", row.names=1, header=TRUE)
head(rawCounts)

gene_info <- rawCounts[1:5]

# ======= Convert names to epiRIL

# delete part (_S and after) of column names
for (col in 6:ncol(rawCounts)){
        colname = colnames(rawCounts)[col]
        print(colname)
        colnames(rawCounts)[col] <- sub("_S.*", "", colname)
}

old <- c("At.test.1", "At.test.2", "At.test.3", "At.test.4", "At.test.5", "At.test.6", "At.test.7", "At.test.8", "At.test.9", "At.test.10","At.test.11", "At.test.12", "At.test.13", "At.test.14", "At.test.15")

new <- c("epiRIL53_Rep1", "epiRIL53_Rep2", "epiRIL53_Rep3", "epiRIL71_Rep1", "epiRIL71_Rep2", "epiRIL71_Rep3", "epiRIL93_Rep1", "epiRIL93_Rep2", "epiRIL93_Rep3", "epiRIL101_Rep1", "epiRIL101_Rep2", "epiRIL101_Rep3", "epiRIL159_Rep1", "epiRIL159_Rep2", "epiRIL159_Rep3")

varnames <- as.data.frame(cbind(old, new))
mm <- match(names(rawCounts), varnames$old)
names(rawCounts)[!is.na(mm)] <- as.character(varnames$new[na.omit(mm)])

sampleData = data.frame()
row = 1
for (i in 6:length(names(rawCounts))){
	epiRIL = strsplit(names(rawCounts)[i], "_")[[1]][1]
	sampleData[row, 1] = names(rawCounts)[i]
	sampleData[row, 2] = epiRIL
	row = row + 1
}
rownames(sampleData) <- sampleData[,1]
sampleData[,1] <- NULL
names(sampleData) <- "group"
sampleData$group <- as.factor(sampleData$group)

# Convert count data to a matrix of appropriate form that DEseq2 can read
counts <- subset(rawCounts, select = -c(Chr, Start, End, Strand, Length))
counts <- counts[row.names(sampleData)] # order data frame
print("Number of genes in original file")
print(nrow(counts))

# ======== not filtered ==============
# remove genes not expressed in all samples
counts <- counts[apply(counts[,-1], 1, function(x) !all(x==0)),]
print("Number of genes after removing all zeros")
print(nrow(counts))

# sort by columns
counts <- counts[ , order(names(counts))]

# calculate means
n <- 3
i <- seq(1, length(counts), n)

counts2 <- data.frame(nut = rownames(counts))
counts2[, paste0('NewCol', seq_along(i))] <- lapply(i, function (j) rowMeans(counts[, j:min(j+2, length(counts))]))
rownames(counts2) <- counts2$nut
counts2$nut <- NULL

# add colnames
newnames <- c()
for (j in 1:length(i)){
	name = str_split(names(counts)[i[j]], pattern="_")[[1]][1]
	newnames[j] = name
}
names(counts2) <- newnames

# rename counts2 to counts
counts <- counts2

cds_not_filtered <- DESeqDataSetFromMatrix(countData=counts, colData=sampleData, design= ~group)
cds_not_filtered = estimateSizeFactors(cds_not_filtered)
normalized_counts_not_filter <- counts(cds_not_filtered, normalized=TRUE)

# Merge counts with gene informaiton
g_normalized_counts_not_filter <- merge(x=gene_info, y=normalized_counts_not_filter, by=0, all.y=TRUE)
# Save
write.table(format(normalized_counts_not_filter, digits=2, nsmall=2), file="not_filtered_normalized_counts.tsv", sep="\t", quote=F, col.names=NA)

# ------ Calculate average and save --------
# sort

# ------ Calculate average and save --------
cds_not_filtered = estimateDispersions(cds_not_filtered)
pdf("not_filtered_dispersion.pdf")
plotDispEsts(cds_not_filtered)
dev.off()
cds_not_filtered <- DESeq(cds_not_filtered)
res_not_filtered <- results(cds_not_filtered)
print(head(res_not_filtered))
print("Number of differentially expressed genes.")
print(sum(res_not_filtered$padj < qvalue, na.rm=T)) 

# MAplot
print("Generate MAplot")
pdf("not_filtered_MAplot.pdf")
plotMA(res_not_filtered, ylim=c(-5,5))
dev.off()

## Significant genes
resSigind = res_not_filtered[ which(res_not_filtered$padj < qvalue & res_not_filtered$log2FoldChange > 0), ] # genes that are induced - up regulated
resSigrep = res_not_filtered[ which(res_not_filtered$padj < qvalue & res_not_filtered$log2FoldChange < 0), ] # genes that are repressed - down regulated
#resSig = rbind(resSigind, resSigrep)
#head(resSig)

pvalues <- format(as.data.frame(resSig), digits=2, nsmall=2)
g_pvalues <- merge(x=gene_info, y=pvalues, by=0, all.y=TRUE)
write.table(pvalues, file="not_filtered_differential_genes.tsv", sep="\t", quote=F, col.names=NA)
print("Significant genes")
print(rownames(resSigind))

# Counts of signifined genes
counts.resSigind <- as.data.frame(normalized_counts_not_filter[rownames(normalized_counts_not_filter) %in% rownames(resSigind), ])
counts.resSigind$status <- "up"
counts.resSigrep <- as.data.frame(normalized_counts_not_filter[rownames(normalized_counts_not_filter) %in% rownames(resSigrep), ])
counts.resSigrep$status <- "down"
print("Counts of significant genes.")
countsSig <- rbind(counts.resSigind, counts.resSigrep)
count_sig <- format(countsSig, digits=2, nsmall=2)
g_count_sig <- merge(x=gene_info, y=count_sig, by=0, all.y=TRUE)
write.table(g_count_sig, file="not_filtered_norm_counts_sig_genes.tsv", sep="\t", quote=F, col.names=NA)
# ======== not iltered ==============

# ======== Filtered ==============

# Filtering counts. We will remove all genes if neither of the groups epiRILs have a median count of 10 and call the new dataframe
medianCountByGroup <- t(apply(counts, 1, tapply, sampleData, median)) 
maxMedian <- apply(medianCountByGroup, 1, max) 

# call filtered data frame
counts_filtered <- counts[maxMedian >= 10,]

## Differentially expressed genes
# Create the DEseq2DataSet object

## *** Starting here there is error since sample data is from before average with 15 rows, but counts only have 5 colnames, therefore colData=sampleData is not true
cds <- DESeqDataSetFromMatrix(countData=counts_filtered, colData=sampleData, design= ~group)

# normalize by library size
cds = estimateSizeFactors(cds)
normalized_counts <- counts(cds, normalized=TRUE)

# save normalized data
nc <- format(normalized_counts, digits=2, nsmall=2)
# join with gene information
g_nc <- merge(x=gene_info, y=nc, by=0, all.y=TRUE)
write.table(g_nc, file="normlized_counts.tsv", sep="\t", quote=F, col.names=NA)

# estimate the dispersion ( or variation ) of the dat
cds = estimateDispersions( cds )

# plot the dispersion
pdf("dispersion.pdf")
plotDispEsts( cds )

# what to compare
cds <- DESeq(cds)
res <- results(cds)
head(res)

# count the number of genes that have a an adjusted p-value less than 0.05.
write.table(sum(res$padj < qvalue, na.rm=T), "filtered_numboer_of_differential.txt" )

# MAplot
pdf("MAplot.pdf")
plotMA(res, ylim=c(-5,5))
dev.off()

## Significant genes
resSigind = res[ which(res$padj < qvalue & res$log2FoldChange > 0), ]
resSigrep = res[ which(res$padj < qvalue & res$log2FoldChange < 0), ]
resSig = rbind(resSigind, resSigrep)
head(resSig)
f_resSig <- format(as.data.frame(resSig), digits=2, nsmall=2)
g_f_resSig <- merge(x=gene_info, y=f_resSig, by=0, all.y=TRUE)
write.table(f_resSig, file="differential_genes.tsv", sep="\t", quote=F, col.names=NA)

rownames(resSigind)

# Counts of signifined genes
counts.resSigind <- as.data.frame(normalized_counts_not_filter[rownames(normalized_counts_not_filter) %in% rownames(resSigind), ])
counts.resSigind$status <- "up"
counts.resSigrep <- as.data.frame(normalized_counts_not_filter[rownames(normalized_counts_not_filter) %in% rownames(resSigrep), ])
counts.resSigrep$status <- "down"

print("Counts of significant genes.")
countsSig <- rbind(counts.resSigind, counts.resSigrep)
#countsSig <- normalized_counts[rownames(normalized_counts) %in% rownames(resSig), ]
f_countsSig <- format(countsSig, digits=2, nsmall=2)
g_f_countsSig <- merge(x=gene_info, y=f_countsSig, by=0, all.y=TRUE)
write.table(g_f_countsSig, file="norm_counts_sig_genes.tsv", sep="\t", quote=F, col.names=NA)
