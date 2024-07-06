# Converts a bpm or wpm scantools output file into a bed

# Usage e.g.
	# python3 /Users/sian_bray/Dropbox/Scripts/scantools_to_bed.py -st /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/UK_dips_Vs_tets_15_BPM.txt -b /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/UK_dips_Vs_tets_15_BPM.bed

# bed file format: chrom	start	end
	# The first base in a chromosome is numbered 0.
	# The chromEnd base is not included in the display of the feature.
		# For example, the first 100 bases of chromosome 1 are defined as
		# chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 
		# 0-99 in our software (not 0-100), but will represent the position
		# notation chr1:1-100

# So end pos is unchanged and start pos is -1

import argparse

parser = argparse.ArgumentParser(description="Converts a bpm or wpm scantools output file into a bed")
parser.add_argument('-st', type=str,  metavar='input_scantools_file', required=True, help='Path to the input scantools file')
parser.add_argument('-b', type=str, metavar='output_bed_file', required=True, help='Path to the output bed file.')
args = parser.parse_args()

in_file = open(args.st, 'r')
out_file = open(args.b, 'w')

for line in in_file:
	if 'started' in globals():
		line = line.replace('.0', '')
		line = line.split('\t')
		start = int(line[2]) - 1
		if 'Genome' not in line[1]:
			out_file.write(f'{line[1]}\t{start}\t{line[3]}\n')
	else:
		started = True

in_file.close()
out_file.close()
