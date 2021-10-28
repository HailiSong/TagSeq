# based on: https://www.biostars.org/p/282685/
# based on: http://huboqiang.cn/2016/03/03/RscatterPlotPCA
library(ggplot2)

normCounts <- read.table("/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/qvalue0.25/norm_counts_sig_genes.tsv", row.names=1, header=TRUE)
head(normCounts)

counts <- subset(normCounts, select = -c(Chr, Start, End, Strand, Length, status))

# order counts by epiRILs
counts <- counts[names(counts)[order(names(counts))]]

# change row names
rownames(counts) <- counts$Row.names
counts$Row.names <- NULL

counts <- t(counts)

# add group
group <- sapply( strsplit(as.character(row.names(counts)), "_"), "[[", 1 )
group <- factor(group, levels = c("epiRIL101", "epiRIL159", "epiRIL53", "epiRIL71", "epiRIL93"))
colType <- c("green1", "yellow1", "blue1", "grey1", "red1")[group]
pchType <- c(18, 16, 15, 17, 19)[group]

library(ggplot2)
library(grid)
library(gridExtra)

counts_pca <- prcomp(counts)

pdf("PCA_deg_norm.pdf")
plot(
	counts_pca$x,
	col = colType,
	pch = pchType,
	cex = 3.0)
legend(
	"bottomright",
	bty = "n",
	c("epiRIL101", "epiRIL159", "epiRIL53", "epiRIL71", "epiRIL93"),
	fill = c("green1", "yellow1", "blue1", "grey1", "red1"),
	cex = 2.0)
dev.off()

#pdf("PCA_deg_norm.pdf")
#ggplot(counts_out,aes(x=PC1,y=PC2,color=group ))
#geom_point()
#dev.off()
