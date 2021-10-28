import sys
import os
import csv
import glob
import pandas as pd
from pathlib import Path

## after some google https://mail.python.org/pipermail/tutor/2004-November/033475.html
## The idea is to keep the count column into a list.
files_dir = Path(sys.argv[1])
pattern = sys.argv[2]
output_file = sys.argv[3]

print(files_dir)
files = glob.glob(os.path.join(files_dir, "*"+pattern))
print(files)

final = pd.DataFrame()
print(final.head())


n = 1
header = ""
for file in files:
	print(n)
	data = pd.read_csv(file, sep="\t", comment="#", index_col=0)
	header= data.iloc[0]
#	print(data.iloc[:,5].head())
#	print(data.head())
	if n<=1:
		final=data
#		print(data.head())
	else:
		final = pd.concat([final, data.iloc[:,5]], axis=1)

	n = n + 1

#final.columns = final.columns.str.lstrip(files_dir).str.strip('_genome.dedup.bam') 

print(final.head())
final.to_csv(os.path.join(files_dir, output_file), sep="\t")

'''
# ======================================
files = glob.glob(os.path.join(files_dir, "*withDup.txt"))
print(files_withDup)
n = 1
header = ""
for file in files:
	print(n)
	data = pd.read_csv(file, sep="\t", comment="#", index_col=0)
	header= data.iloc[0]
#       print(data.iloc[:,5].head())
#       print(data.head())
	if n<=1:
		final=data
#               print(data.head())
	else:
		final = pd.concat([final, data.iloc[:,5]], axis=1)

	n = n + 1

final.to_csv(os.path.join(files_dir, "final_all_samples_featurecounts_withDup.tsv"), sep="\t")
'''
