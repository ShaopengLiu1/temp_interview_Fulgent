import argparse
import numpy as np
import pandas as pd

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description="Get chr11 coverage",
	                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-f', '--input_file', type=str, help="Path to file containing reads' bed file")


	### read parameters
	args = parser.parse_args()
	input_bed = args.input_file

	out_chr11_cov = np.zeros(136000000)
	df_read = pd.read_csv(input_bed, sep="\t", header=None)
	# loop rows
	for index in range(df_read.shape[0]):
	  out_chr11_cov[df_read.iloc[index, 1] : df_read.iloc[index, 2]] += 1
	  
	# store results
	np.save(file="chr11_base_level_coverage", arr=out_chr11_cov)

