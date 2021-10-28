
samples = c("epiRIL53", "epiRIL71", "epiRIL93", "epiRIL101", "epiRIL159")

marker = read.table(file="/home/haili/genomes/tair10/123epiRILs_matrix.csv", header=TRUE, row.names=1, sep=";")


exp = read.table(file="/home/haili/workspace/results/epiRIL_test_15_norm_loc_srtd.csv", header=TRUE, row.names=1, sep="\t")


exp_avg = read.table(file="/home/haili/workspace/results/epiRIL_test_15_norm_avg_loci_srtd.csv", header=TRUE, row.names=1, sep="\t")

for (i in 1:length(samples)){
	column_data = exp_avg[c("chr", "start", "end", samples[i])]
	column_data[length(column_data)+1] = "--"
	colnames(column_data) = NA
	write.table(column_data, paste("/home/haili/workspace/results/combine/exp_avg_", samples[i], ".txt", sep=""), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

	print(i)
	print(samples[i])
	sliced_cols = colnames(exp)[grepl(paste(samples[i], "_", sep=""), colnames(exp))]
	sliced_data = exp[c("chr", "start", "end", sliced_cols)]
	print(head(sliced_data))
	#sliced_data[length(sliced_data)+1] = "--"
	#colnames(sliced_cols) = NA
	write.table(sliced_data, paste("/home/haili/workspace/results/combine/exp_", samples[i], ".txt", sep=""), row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
	
	# slice and save marker data
	sliced_methyl = marker[c("chr", "start_bp", "stop_bp", samples[i])]
	sliced_methyl[length(sliced_methyl) + 1] = sliced_methyl[length(sliced_methyl)]
	sliced_methyl[length(sliced_methyl)-1] = "--"
	colnames(sliced_methyl) = NA #c("chr", "start_bp", "stop_bp", "NAcol", samples[i])
	write.table(sliced_methyl, paste("/home/haili/workspace/results/combine/marker_", samples[i], ".txt", sep=""),  row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

#	rbind(column_data, sliced_methyl)
#	rbind(sliced_data, sliced_methyl)
	
}



# write.table(samples5, file="/home/haili/workspace/results/test5epiRILs_matrix.csv", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
