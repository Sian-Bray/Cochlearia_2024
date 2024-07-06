# Script to count how many MAV SNPs are in genes of interest
# Then add a MAVs column to the end of the script
# Originally written on 9th May 2023
# Modified on 13th May 2023

# Testing commands:
	# cd ~/
	# python3 /Users/sian_bray/Dropbox/Scripts/MAV_count_to_file.py -e test_excel.tsv -g test.gff -m test_MAV.txt -o test_out.txt

# Example command:
	# python3 /Users/sian_bray/Dropbox/Scripts/MAV_count_to_file.py -e /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_1perc.txt -g /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/C_excelsa_V5_braker2_wRseq_ANNOTATION/genes_only.gtf -m /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/03_Pirita_Mav/british_1percent_outliers_3.5_scores.txt -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_1perc_plus_MAVs.txt

import argparse
parser = argparse.ArgumentParser(description='Script to count how many MAV SNPs are in genes of interest.')
parser.add_argument('-e', type=str, required=True, help='Candadate genes and BPM metrics list. Headder is: gene	outname	scaff	start	end	win_size	num_sites	num_snps	Rho	FstWC	dxy	AFD	FixedDiff	FstH	FstN. Any file that has a headder, is tab delimited and has genes in the first column will work.')
parser.add_argument('-g', type=str, required=True, help='The gff file. Example line: Cexcelsa_scaf_4	AUGUSTUS	gene	1	9453	1	+	.	ID=g1')
parser.add_argument('-m', type=str, required=True, help='MAV list. Example line: Cexcelsa_scaf_3 	4285395 	p.Phe168Cys 	1 	169 	0.0168890124048 	0.827730497961 	205 	3.46224754299 	169.684752082 	166.222504539')
parser.add_argument('-o', type=str, required=True, help='Output file.')
args=parser.parse_args()

# create list of candidate genes and ...
# ... create dictionary of excel file gene lines (gene as key, full line as stuff)
# gene	outname	scaff	start	end	win_size	num_sites	num_snps	Rho	FstWC	dxy	AFD	FixedDiff	FstH	FstN

in_excel = open(args.e, 'r')
gene_list = []
line_dict = {}

for line in in_excel:
	if 'gene\t' in line:
		continue
	gene = line.split('\t')
	gene = gene[0]
	gene_list.append(gene)
	line = line.replace('\n', '')
	line_dict[gene] = line

in_excel.close()

# read through the gff file until you find genes in the list
# Cexcelsa_scaf_4	AUGUSTUS	gene	1	9453	1	+	.	ID=g1

gff_dict = {}
gff_file = open(args.g, 'r')

for line in gff_file:
	line = line.split('\t')
	gff_gene = line[8]
	gff_gene = gff_gene.replace('\n', '')
	gff_gene = gff_gene.replace('ID=', '')
	if gff_gene in gene_list:
		gff_dict[gff_gene] = [line[3], line[4]]

gff_file.close()

# read through the MAVs checking to see if they fall within the gene
# Cexcelsa_scaf_3 	4285395 	p.Phe168Cys 	1 	169 	0.0168890124048 	0.827730497961 	205 	3.46224754299 	169.684752082 	166.222504539

MAV_count_dict = {}

for gene in gene_list:
	MAV_file = open(args.m, 'r')
	MAV_count_dict[gene] = 0
	for line in MAV_file:
		line = line.split('\t')
		location = int(line[1])
		gene_range = gff_dict[gene]
		gene_range[0] = int(gene_range[0])
		gene_range[1] = int(gene_range[1])
		if location >= gene_range[0] and location <= gene_range[1]:
			MAV_count_dict[gene] += 1
	MAV_file.close()

# Write out lines

out_file = open(args.o, 'w')
in_excel = open(args.e, 'r')

for line in in_excel:
	line = line.replace('\n', '')
	if 'gene\t' in line:
		out_file.write(line+'\tMAV_count\n')
	if 'gene\t' not in line:
		final_gene = line.split('\t')
		final_gene = final_gene[0]
		final_MAV_count = MAV_count_dict[final_gene]
		out_file.write(line+f'\t{final_MAV_count}\n')






