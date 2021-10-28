library(ggplot2)

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

counts <- subset(rawCounts, select = -c(Chr, Start, End, Strand, Length))

# order counts by epiRILs
counts <- counts[names(counts)[order(names(counts))]]

# Calcualte PCA
counts_pca <- prcomp(counts)

# plot
pdf("PCA.pdf")
plot(counts_pca$x[,1], counts_pca$x[,2])
dev.off()

