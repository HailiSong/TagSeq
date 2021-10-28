#!/usr/bin/env Rscript
rm(list=ls())

# See https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

setwd("/Users/haili/Documents/WorkUT/epiRIL_all_20210524/anova")

bad_epis = c(
  "epi8.1", "epi8.2", "epi8.3", "epi11.1", "epi11.2", "epi11.3", 
  "epi36.1", "epi36.2", "epi36.3", "epi60.1", "epi60.2", "epi60.3", 
  "epi94.1", "epi94.2", "epi94.3", "epi144.1", "epi144.2", "epi144.3", 
  "epi150.1", "epi150.2", "epi150.3", "epi166.1", "epi166.2", "epi166.3", 
  "epi172.1", "epi172.2", "epi172.3", "epi297.1", "epi297.2", "epi297.3", 
  "epi361.1", "epi361.2", "epi361.3", "epi362.1", "epi362.2", "epi362.3", 
  "epi363.1", "epi363.2", "epi363.3", "epi366.1", "epi366.2", "epi366.3", 
  "epi368.1", "epi368.2", "epi368.3", "epi371.1", "epi371.2", "epi371.3"
)

# Read in count file
library(plyr)
library(dplyr)

rawcounts <- read.table("/Users/haili/Documents/WorkUT/epiRIL_all_20210524/rawcounts.txt", 
                        header = T, sep="\t", row.names = 1)
rawcounts <- select(rawcounts, -c("Chr", "Start", "End", "Strand", "Length"))
genes <- select(rawcounts, 1:5)

# Keep the first part of sample names
names(rawcounts) <- sapply(strsplit(names(rawcounts),"_"), `[`, 1)

# remove duplicated columns
rawcounts <- rawcounts[, !duplicated(colnames(rawcounts))]

# remove bad samples
rawcounts <- rawcounts[,!colnames(rawcounts) %in% bad_epis]

# sort by column names
rawcounts <- rawcounts[,order(colnames(rawcounts))]


#batches <- read.table("/Users/haili/Documents/WorkUT/epiRIL_all_20210524/metadata.txt", header = F)
#colnames(batches) <- c("epi", "batch")
#batches <- batches[!batches$epi %in% bad_epis,]

rils <- as.data.frame(colnames(rawcounts))
group <- substr(rils[,1],1,nchar(rils[,1])-2)
rils$group <- group
rownames(rils) <- rils[,1]
rils[,1] <- NULL

dds <- DESeqDataSetFromMatrix(countData=rawcounts, colData=rils, design= ~group)

# ANOVA/Likelihood ratio test
# See: https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html
dds_lrt <- DESeq(dds, test="LRT", full=~group, reduced = ~ 1)
# -----------------------------------
# Extract results
res_LRT <- results(dds_lrt)

library(tidyverse)

# Subset the LRT results to return genes with padj < 0.05
padj.cutoff = 0.05
sig_res_LRT <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < padj.cutoff)

# save DE genes
write.table(sig_res_LRT, "likelyhood_ratio_test_DE_genes.tsv", sep="\t", quote = F)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)

length(sigLRT_genes)

save.image(file = "DESeq_LRT.RData")

