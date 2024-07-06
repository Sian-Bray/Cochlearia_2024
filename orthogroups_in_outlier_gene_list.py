# Outputs a list of orthogroups in a list of outlier genes
# Takes an orthogroup file and a list of outlier genes
# For Venn input in R
# Written by Sian Bray on 30th July 2023

# python3 /Users/sian_bray/Dropbox/Scripts/compare_gene_lists_and_orthogroups_V2.py -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/06_Orthofinder/OrthoFinder/Results_May19/Orthogroups/Orthogroups.tsv -p /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top_1perc_minus_inversion.txt -t /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top1perc_Fst_arenosa_lyrataIDs.txt -a /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top1perc_Fst_amara_amaraIDs.txt

import argparse

parser=argparse.ArgumentParser(description="Outputs a list of orthogroups in a list of outlier genes.")
parser.add_argument('-io', type=str, metavar='input_orthogroups', required=True, help='Path to the orthogroups file. Output of orthofinder; Orthogroups.tsv')
parser.add_argument('-ig', type=str, metavar='input_genes', required=True, help='Text file containing input genes, exactly as they are in the orthogroups file, one gene per line.')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='The output file that will be created.')
args=parser.parse_args()

# Define a function that breaks an orthofinder line into a clean list
def ortho_split(ortho_line):
	ortho_line = ortho_line.replace(', ', '\t')
	ortho_line = ortho_line.replace('\n', '')
	ortho_line = ortho_line.split('\t')
	return ortho_line

# open files
in_ortho = open(args.io, 'r')
in_genes = open(args.ig, 'r')
out_file = open(args.o, 'w')

# Convert input genes into list of genes
gene_list = []
for line in in_genes:
	line = line.replace('\n', '')
	gene_list.append(line)
# print(gene_list)
# print(len(gene_list))

# Compare each line of the orthogroup file to the genes
for ortho_line in in_ortho:
	OG = ortho_split(ortho_line)
	for gene in gene_list:
		if gene in OG:
			out_file.write(OG[0]+'\n')

# Close files
in_ortho.close()
in_genes.close()
out_file.close()



