args <- commandArgs(trailingOnly = TRUE)

#in_dir = "/home/haili/workspace/epiRIL_tagseq_Jun/count"
#out_dir = "/home/haili/workspace/epiRIL_tagseq_Jun/results"

DF <- do.call(cbind,
         lapply( list.files(path=args[1], pattern=".*genes.results", full.names=T),
                     FUN=function(x) {
#            print(x)
            aColumn = read.table(x,header=T, sep="\t")[,c("gene_id", "expected_count")];
            colnames(aColumn)[2] = x;
            aColumn;
             }
            )
        )
DF <- DF[,!duplicated(colnames(DF))]
DF2 <- DF[,-1]
rownames(DF2) <- DF[,1]

#Save orignal data to file


write.table(DF2, file=file.path(args[2], "epiRIL_test_15.csv"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

# change names
old <- c("At-test-1", "At-test-2", "At-test-3", "At-test-4", "At-test-5", "At-test-6", "At-test-7", "At-test-8", "At-test-9", "At-test-10","At-test-11", "At-test-12", "At-test-13", "At-test-14", "At-test-15")

new <- c("epiRIL53_Rep1", "epiRIL53_Rep2", "epiRIL53_Rep3", "epiRIL71_Rep1", "epiRIL71_Rep2", "epiRIL71_Rep3", "epiRIL93_Rep1", "epiRIL93_Rep2", "epiRIL93_Rep3", "epiRIL101_Rep1", "epiRIL101_Rep2", "epiRIL101_Rep3", "epiRIL159_Rep1", "epiRIL159_Rep2", "epiRIL159_Rep3")

# delete part (_S and after) of column names
for (col in 1:ncol(DF2)){
	colname = colnames(DF2)[col]
#       print(colname)
        filename = substr(colname, nchar(args[1])+2, nchar(colname)+1) # remove directory path from colname
	colnames(DF2)[col] <- sub("_S.*", "", filename)	
}

# replace colnames
varnames <- as.data.frame(cbind(old, new))
mm <- match(names(DF2), varnames$old)
names(DF2)[!is.na(mm)] <- as.character(varnames$new[na.omit(mm)])

# order by column names
DF2<- DF2[ , order(names(DF2))]
print("Raw count:")
print(head(DF2))

# Normalized by edgeR
library("edgeR")
group = sub("_.*", "", colnames(DF2))
dge <- DGEList(counts=DF2, group=group)
nc <- data.frame(cpm(dge, normalized.lib.sizes=FALSE))
print("First few lines of normalized count:")
print(head(nc))
format_nc <- format(round(nc, 2), nsmall = 2)
# Save normalized data
write.table(format_nc, file=file.path(args[2], "epiRIL_test_15_norm.csv"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

# Compute average of replicates
# average <- data.frame(apply(array(as.matrix(DF2[,-1]), c(nrow(DF2),3, ncol(DF2)/3)),3, rowMeans)) # this one has no gene name and sample name

DF_avg <- data.frame(matrix(nrow = nrow(nc)))

print("Declare matrix of average count")
print(head(DF_avg))


print(names(nc))
for (i in unique(gsub("_.*", "", names(nc)))){
	DF_avg[,i]=apply(nc[,grepl(gsub("_.*", "", i), names(nc))],1, mean)
}
# DF_avg[,1] <- DF2[,1]

row.names(DF_avg) <- row.names(nc)
DF_avg[1] <- NULL

print(head(DF_avg))

# Save the average value
write.table(DF_avg, file=file.path(args[2], "epiRIL_test_15_norm_avg.csv"), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

# ========= Intervial ===========


