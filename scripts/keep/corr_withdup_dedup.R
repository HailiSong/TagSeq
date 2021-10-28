library(ggplot2)
library(gridExtra)
library(DESeq2)
library(plyr)

withdup <- read.table("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/withdup_featurescounts.tsv", header=TRUE, sep="\t", row.names=1)

dedup <- read.table("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/dedup_featurescounts.tsv", header=TRUE, sep="\t", row.names=1)

sample_names <- names(dedup)[grepl("test", names(dedup))]

# ====================================
## Normalize with DESeq2
sampleData <- data.frame(read.csv("/home/haili/data/raw/epiRIL_tagseq_Jun/expdesign.txt", row.names=1, header=TRUE))
row.names(sampleData) <- gsub("-", "\\.", row.names(sampleData))
sampleData2 <- sampleData

# ==================================
# Get counts only
counts_withdup <- subset(withdup, select = -c(Chr, Start, End, Strand, Length))
counts_dedup <- subset(dedup, select = -c(Chr, Start, End, Strand, Length))

# =================================
# order data frame
counts_withdup <- counts_withdup[row.names(sampleData)]
counts_dedup <- counts_dedup[row.names(sampleData)]

# =================================
# Create the DEseq2DataSet object
cds_withdup <- DESeqDataSetFromMatrix(countData=counts_withdup, colData=sampleData, design= ~group)
cds_dedup <- DESeqDataSetFromMatrix(countData=counts_dedup, colData=sampleData, design= ~group)

# ===============================
# Normalize by libaray size
cds_withdup <- estimateSizeFactors(cds_withdup)
cds_dedup <- estimateSizeFactors(cds_dedup)

norm_withdup =counts(cds_withdup, normalized=TRUE) 
norm_dedup =counts(cds_dedup, normalized=TRUE)

print("==================normalized==============")

col_i_withdup = norm_withdup[, "At.test.7_S55"]
col_i_dedup = norm_dedup[, "At.test.7_S55"]
mylist <- list( one=as.data.frame(col_i_withdup), two=as.data.frame(col_i_dedup) )
for(i in 1:length(mylist)){
	colnames(mylist[[i]]) <- paste0( names(mylist)[i], "_", colnames(mylist[[i]]) )
	 mylist[[i]]$GENEID  <- rownames(mylist[[i]])
}
out <- join_all( mylist, by="GENEID", type="right" )
rownames(out) <- out$GENEID
out$GENEID <- NULL
names(out) = c("withdup", "dedup")
pdf(paste0("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/norm_scatterplot_", sample_names[i], ".pdf"))
p = ggplot(out, aes(x=withdup, y=dedup)) + geom_point() + geom_smooth(method=lm)
dev.off()


for (i in 1:length(sample_names))
{
	print(sample_names[i])
	col_i_withdup = norm_withdup[, i]
	col_i_dedup = norm_dedup[, i]
	mylist <- list( one=as.data.frame(col_i_withdup), two=as.data.frame(col_i_dedup) )
	for(i in 1:length(mylist)){
		colnames(mylist[[i]]) <- paste0( names(mylist)[i], "_", colnames(mylist[[i]]) )
		mylist[[i]]$GENEID  <- rownames(mylist[[i]])
        }
        out <- join_all( mylist, by="GENEID", type="right" )
        rownames(out) <- out$GENEID
        out$GENEID <- NULL
        names(out) = c("withdup", "dedup")
        print(head(out))
	pdf(paste0("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/norm_scatterplot_", sample_names[i], ".pdf"))
        p = ggplot(out, aes(x=withdup, y=dedup)) + geom_point()
	p + geom_smooth(method=lm)
	dev.off()
#	ggsave(plot=p, filename=paste0("norm_scatterplot_", sample_names[i], ".pdf"))
}

## ===============================
# Filter genes with low expression

maxMedian_withdup <- apply(medianCountByGroup_withdup, 1, max)

# call filtered data frame
counts_filtered_withdup <- counts_withdup[maxMedian_withdup >= 10,]
# ---------------------
medianCountByGroup_dedup <- t(apply(counts_dedup, 1, tapply, sampleData, median))
maxMedian_dedup <- apply(medianCountByGroup_dedup, 1, max)

# call filtered data frame
counts_filtered_dedup <- counts_dedup[maxMedian_dedup >= 10,]

## ===============================
# Create the DEseq2DataSet object
cds_filtered_withdup <- DESeqDataSetFromMatrix(countData=counts_filtered_withdup, colData=sampleData, design= ~group)

# normalize by library size
cds_filtered_withdup = estimateSizeFactors(cds_filtered_withdup)
normalized_counts_filtered_withdup <- counts(cds_filtered_withdup, normalized=TRUE)
# ------------------

cds_filtered_dedup <- DESeqDataSetFromMatrix(countData=counts_filtered_dedup, colData=sampleData, design= ~group)

# normalize by library size
cds_filtered_dedup = estimateSizeFactors(cds_filtered_dedup)
normalized_counts_filtered_dedup <- counts(cds_filtered_dedup, normalized=TRUE)

print("==================filterred normalized==============")

# generate scatter plots
for (i in 1:length(sample_names))
{
        col_i_withdup = normalized_counts_filtered_withdup[, i]
        col_i_dedup =normalized_counts_filtered_dedup[, i]
	mylist <- list( one=as.data.frame(col_i_withdup), two=as.data.frame(col_i_dedup) )
	for(i in 1:length(mylist)){
		colnames(mylist[[i]]) <- paste0( names(mylist)[i], "_", colnames(mylist[[i]]) )
  		mylist[[i]]$GENEID  <- rownames(mylist[[i]])
	}
	out <- join_all( mylist, by="GENEID", type="right" )
	rownames(out) <- out$GENEID 
	out$GENEID <- NULL
#        df = data.frame(merge(col_i_withdup, col_i_dedup, all.y=TRUE))
#	print(head(df))
#	df2 = df[,-1]
#	rownames(df2) = df[,1]
	names(out) = c("withdup", "dedup")
	print(head(out))
        pdf(paste0("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/filter_norm_scatterplot_", sample_names[i], ".pdf"))
	p = ggplot(out, aes(x=withdup, y=dedup)) + geom_point()
	dev.off()
#        ggsave(plot=p, filename=paste0("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/filter_norm_scatterplot_", sample_names[i], ".pdf"))
}

