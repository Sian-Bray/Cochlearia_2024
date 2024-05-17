# Adds annotation to your candidate gene list (quite crude)
# Annotation comes from the output of add_gene_descriptions_to_rbh_and_orthogroups.py

# Example command:
# python3 /Users/sian_bray/Dropbox/Scripts/annotate_gene_list.py -i /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/02_Visualisation_Folder/01_UK_Dips_Tets/temp_files/temp4.txt -a /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Annotated_Orthogroups.tsv -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/1perc_fst_outliers_annotated.tsv

import argparse

parser=argparse.ArgumentParser(description="Adds annotation to your candidate gene list (quite crude). Annotation comes from the output of add_gene_descriptions_to_rbh_and_orthogroups.py")
parser.add_argument('-i', type=str, metavar='input_gene_list', required=True, help='Path to the input candidate gene list, one gene per line. This could be the temp4.txt file from my annotation pipe.')
parser.add_argument('-a', type=str, metavar='input_annotation', required=True, help='Path to the gene annotation file, i.e. the output of add_gene_descriptions_to_rbh_and_orthogroups.py.')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='Path to the output file.')
args=parser.parse_args()

gene_file = open(args.i, 'r')
annotation = open(args.a, 'r')
out_file = open(args.o, 'w')
gene_list = []

for line in gene_file:
	gene = line.replace('\n', '')
	gene = gene + '.'
	gene_list.append(gene)

for line2 in annotation:
	for gene2 in gene_list:
		if gene2 in line2:
			out_file.write(gene2+'\t'+line2)

