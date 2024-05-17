# Takes a ScanTools bpm file, pull out genes for each window, then re-write each line preceded by the genes (once per gene such that a window may occur more than once)

# Previously used bedtools command:
	# bedtools intersect -wb -a /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/UK_dips_Vs_tets_15_BPM.bed -b /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/C_excelsa_V5_braker2_wRseq_ANNOTATION/C_excelsa_V5_braker2_wRseq.gff3 | grep gene | awk '{print $12}' | sort | uniq > /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/bpm_genes.txt

# python3 /Users/sian_bray/Dropbox/Scripts/bpm_windows_to_genes.py -i /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/UK_dips_Vs_tets_15_BPM.txt -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM.txt -g /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/C_excelsa_V5_braker2_wRseq_ANNOTATION/C_excelsa_V5_braker2_wRseq.gff3

import argparse
import subprocess

parser = argparse.ArgumentParser(description="Converts a bpm or wpm scantools output file into a bed")
parser.add_argument('-i', type=str,  metavar='input_scantools_file', required=True, help='Path to the input scantools file.')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='Path to the output file to be created.')
parser.add_argument('-g', type=str, metavar='gff3_file', required=True, help='Path to the corresponding gff3 file.')
args = parser.parse_args()

in_file = open(args.i, 'r')
out_file = open(args.o, 'w')
gff3 = args.g

bed_zones = []

for full_line in in_file:

	if 'started' in globals():
		line = full_line.replace('.0', '')
		line = line.split('\t')
		start = int(line[2]) - 1

		if 'Genome' not in line[1]:
			zone=f'{line[1]}\t{start}\t{line[3]}\n'
			bed_file = open('temp_bed.bed', 'w')
			bed_file.write(zone)
			bed_file.close()
			awkward_bit = '{print $12}'
			cmd = f"bedtools intersect -wb -a temp_bed.bed -b {gff3} | grep gene | awk '{awkward_bit}'"
			result = subprocess.check_output(cmd, shell=True).splitlines()
			
			for gene in result:
				gene=str(gene)
				gene=gene.replace("b'ID=", "")
				gene=gene.replace(";'", "")
				out_file.write(str(gene)+'\t'+full_line)

	else:
		started = True
		out_file.write('gene\t'+full_line)

in_file.close()
out_file.close()