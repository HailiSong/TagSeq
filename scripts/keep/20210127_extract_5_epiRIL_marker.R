# source("/home/haili/code/epiRIL_tagseq_Jun/scripts/20210127_extract_5_epiRIL_marker.R")
epiRILs = c("epiRIL53", "epiRIL71", "epiRIL93", "epiRIL101", "epiRIL159")

marker = read.table(file="/home/haili/genomes/tair10/123epiRILs_matrix.csv", header=TRUE, row.names=1, sep=";")

wd = "/home/haili/workspace/epiRIL_tagseq_Jun_umi2/results/qvalue0.25/"

counts_file = paste0(wd, "norm_counts_sig_genes.tsv")

exp = read.table(file=counts_file, header=TRUE, row.names=1, sep="\t")
print(colnames(exp))


# scale
exp_gene <- exp[c("Row.names", "Chr", "Start", "End", "Strand", "Length", "status")]
exp_counts <- exp[grepl("epiRIL", names(exp))]
exp_counts <- exp_counts[,sort(names(exp_counts))]
# ===== average ====
counts <- exp_counts
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

# ===== average ====

# ===== scale   ====
exp_counts <- counts2
exp_scale <- as.data.frame(scale(exp_counts, scale=TRUE))
exp_scale$Row.names <- exp$Row.names
exp <- merge(exp_gene , exp_scale, by="Row.names")
# ===== scale   ====


# Function for merging two data frames
fastmerge <- function(d1, d2) 
{
	d1.names <- names(d1)
	d2.names <- names(d2)

	# columns in d1 but not in d2
	d2.add <- setdiff(d1.names, d2.names)

	# columns in d2 but not in d1
	d1.add <- setdiff(d2.names, d1.names)

	# add blank columns to d2
	if(length(d2.add) > 0) 
	{
		for(i in 1:length(d2.add)) 
		{
			d2[d2.add[i]] <- NA
		}
  	}

	# add blank columns to d1
	if(length(d1.add) > 0) 
	{
		for(i in 1:length(d1.add)) 
		{
			d1[d1.add[i]] <- NA
		}
	}

	return(rbind(d1, d2))
}

for (epiRIL in epiRILs){
	print(epiRIL)

	# slice and save marker data
	v_epiRIL = marker[,epiRIL]
        g_sliced_methyl = marker[c("chr", "start_bp", "stop_bp")]
#	g_sliced_methyl[4:6] = "--"
	g_sliced_methyl$epiRIL = v_epiRIL
	names(g_sliced_methyl) <- c("Chr", "Start", "End", paste0("methyl_",epiRIL))	

	#sliced_cols = colnames(exp)[grepl(epiRIL, colnames(exp))] # This one is used for non-average

	# the column to select
	sliced_cols = epiRIL
	print(sliced_cols)

	# slice data
	sliced_data = exp[c("Row.names", "Chr", "Start", "End", sliced_cols, "status")]
	rownames(sliced_data) <- sliced_data$Row.names

	# rename rownames
	sliced_data$Row.names <- NULL

	# remove "Chr" form the Chr column
	sliced_data$Chr <- sub("Chr", "", sliced_data$Chr)
	# add empty column
	#sliced_data$methyl <- "--"
	print(head(sliced_data))
	
	merged <- fastmerge(sliced_data, g_sliced_methyl)
	
	ordered <- merged[with(merged, order(Chr, Start, End)),]
	# change rownames to first column
	feature <- row.names(ordered)
	ordered <- cbind(feature, ordered)
	rownames(ordered) <- NULL
	write.table(format(ordered, digits=2, nsmall=2, scientific=FALSE), paste0("scale_avg_exp_with_marker_", epiRIL, ".tsv"), sep="\t", quote=FALSE, row.names=FALSE, na="")
}



# write.table(samples5, file="/home/haili/workspace/results/test5epiRILs_matrix.csv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
