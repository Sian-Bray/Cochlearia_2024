#Add gene descriptions to recipricol-best-hits or ortholouge groups output
#Sian Bray 2nd July 2019

# Example command:
# python3 /Users/sian_bray/Dropbox/Scripts/add_gene_descriptions_to_rbh_and_orthogroups.py -b /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Orthogroups.tsv -d /Users/sian_bray/Dropbox/Yant/Cochlearia/02_Data/References/Thaliana_proteins/Araport11_functional_descriptions_20181231.txt -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Annotated_Orthogroups.tsv -bc 1

import argparse

parser=argparse.ArgumentParser(description="Add gene descriptions to recipricol-best-hits or ortholouge groups. Defaults are set to use the output orthogroups.tsv from OrthoFinder.")
parser.add_argument('-b', type=str, metavar='input_blast', required=True, help='Path to the input reciprocol-best-hits or orthologue groups file.')
parser.add_argument('-d', type=str, metavar='descriptions_file', required=True, help='Path to the gene descriptions file.')
parser.add_argument('-o', type=str, metavar='output_file', required=True, help='Name and/or path for the output file.')
parser.add_argument('-bc', type=int, metavar='blast_column', required=False, default=2, help='The column that contains the gene names of the annotated species (i.e. A thaliana) in the input reciprocol-best-hits or orthologue groups file. Column is zero-based i.e. the first column is 0, the second is 1, the third 2, etc.')
args=parser.parse_args()

input_blast=open(args.b, 'r')
output_file=open(args.o, 'w+')

for count0, line in enumerate(input_blast):
	split_line=line.replace('\n', '')
	split_line=split_line.split('\t')
	genes_to_find=split_line[args.bc].split(', ')
	try:
		genes_to_find.remove('')
	except ValueError:
		pass
	output_file.write(line.replace('\n', ''))
	for count1, gene in enumerate(genes_to_find):
		input_descriptions=open(args.d, 'r')
		for count2, line2 in enumerate(input_descriptions):
			if gene in line2:
				line2=line2.replace('\t', ', ')
				output_file.write('\t'+line2.replace('\n', ''))
		input_descriptions.close()
	output_file.write('\n')

input_blast.close()
output_file.close()
