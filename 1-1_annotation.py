# Assigns 1-2-1 homologue names for genes in a non-model organism to genes in a model organism
# Written by Sian Bray on 20/21st December 2022

# Example run command:
# python3 1-1_annotation.py -f input_fasta -i input_ortho -o output_file -r input_rbh_list -a blast_file
# cd /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/1-2-1_annotation
# python3 /Users/sian_bray/Dropbox/Scripts/1-1_annotation.py -f test_cochlearia_genes.fasta -i test_orthogroups.tsv -o test_output.txt -r test_rbh.txt -a test_all_vs_all.tsc
# python3 /Users/sian_bray/Dropbox/Scripts/1-1_annotation.py -f /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta -i /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Orthogroups.tsv -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_1-2-1_annotation.tsv -r /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/1-2-1_annotation/cochlearia_vs_thaliana_rbh -a /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/1-2-1_annotation/cochlearia_thaliana_all_Vs_all.tsv

# Requires pre-runs of: BLAST, blast_rbh.py, orthofinder

import argparse

# User input
parser = argparse.ArgumentParser(description="Figures out the best student-project combos.")
parser.add_argument('-f', type=str,  metavar='input_fasta', required=True, help='Path to the input fasta file. Name of gene is exactly the same as in the orthogroup.')
parser.add_argument('-i', type=str,  metavar='input_ortho', required=True, help='Path to the input Orthofinder file. This is the Orthogroups.tsv file produced. Headder is: "Orthogroup\tspecies_1\tspecies_2\tetc."')
parser.add_argument('-o', type=str,  metavar='output_file', required=True, help='Path to the output file.')
parser.add_argument('-r', type=str,  metavar='input_rbh_list', required=True, help='Path to the input rbh file. This is the .txt file produced by blast_rbh.py. Tab delimited with the headder: "#A_id	B_id	A_length	B_length	A_qcovhsp	B_qcovhsp	length	pident	bitscore". Make sure that A is the target species (non-model) and B is the model species.')
parser.add_argument('-a', type=str,  metavar='blast_file', required=True, help='Path to the all Vs all BLAST file. It has no headder, e.g. line: "1	1	100.000	5156	0	0	1	5156	1	5156	0.0	9522", where column 0 = target gene name, column 1 = model gene name, column 10 = E-value')
args = parser.parse_args()


def clean_line(line):
	line = line.replace('\n', '')
	line = line.replace('>', '')
	return line

def find_rbh_match(in_rbh, gene_name):
	current_hit = ''
	rbh = open(in_rbh, 'r')
	for hit in rbh:
		hit = clean_line(hit)
		hit = hit.split('\t')
		if gene_name in hit[0]:
			rbh.close()
			return hit[1]
	rbh.close()

def find_blast(in_blast, gene_name):
	blast_hits = open(in_blast, 'r')
	with_gene = []
	for line in blast_hits:
		line = clean_line(line)
		line = line.split('\t')
		if gene_name in line[0]:
			with_gene.append(line)
	blast_hits.close()
	return with_gene

def lowest_e_value(blast_hits, must_contain=False): # model is a list of model genes, blast_hits is a list of blast lines split by columns
	if must_contain != False:
		for b_hit in blast_hits:
			if b_hit[1] not in must_contain:
				blast_hits.remove(b_hit)

	e_vals = []
	bit_vals = []
	for b_hit in blast_hits:
		e_vals.append(float(b_hit[10]))
		bit_vals.append(float(b_hit[11]))
	try:
		min_e = min(e_vals)
	except ValueError:
		return ['', 'NA']
	if e_vals.count(min_e) == 1:
		for b_hit in blast_hits:
			if float(b_hit[10]) == min_e:
				return b_hit
	if e_vals.count(min_e) > 1:
		best_bits = []
		for count, item in enumerate(e_vals):
			if item == min_e:
				best_bits.append(bit_vals[count])
		best_bit = max(best_bits)
		for b_hit in blast_hits:
			if (float(b_hit[10]) == min_e) and (float(b_hit[11]) == best_bit):
				return b_hit

def find_homo(gene_name, ortho_file=args.i, model_column=1, target_column=2, in_rbh=args.r, in_blast=args.a): # columns are zero bases, model = model species, target = your non-model species
	in_orto = open(ortho_file, 'r')

	# First if the genes’ OG contained only one A. thaliana gene ID, that gene ID was used.
	for ortho_group in in_orto:
		ortho_group = clean_line(ortho_group)
		ortho_group = ortho_group.split('\t')
		model = ortho_group[model_column]
		target = ortho_group[target_column]
		if gene_name in target:
			model = model.split(', ')
			if (len(model) == 1) and (model != ['']):
				in_orto.close()
				return model[0], '1-2-1 orthogroup'

			# If the OG contained more than one A. thaliana gene ID then the RBH was taken.
			if len(model) > 1:
				hit = find_rbh_match(in_rbh, gene_name)
				if hit != None:
					for OG_gene in model:
						if hit in OG_gene: # if the hit gene is in the OG
							in_orto.close()
							return hit, 'within orthogroup RBH'

			# If there was no RBH then the OG gene with the lowest E-value in a BLAST versus the TAIR10 database was taken.
			if len(model) > 1:
				blast_hits = find_blast(in_blast, gene_name)
				best_blast = lowest_e_value(blast_hits, must_contain=model)
				if best_blast != None:
					in_orto.close()
					return best_blast[1], 'within orthogroup BLAST'

	# If no OG contained the Cochlearia gene then the RBH was taken.
	hit = find_rbh_match(in_rbh, gene_name)
	if hit != None:
		in_orto.close()
		return hit, 'no orthogroup RBH'

	# Finally, if there was no OG or RBH then the gene with the lowest E-value in a BLAST versus the TAIR10 database was taken.
	blast_hits = find_blast(in_blast, gene_name)
	best_blast = lowest_e_value(blast_hits)
	in_orto.close()
	if best_blast[1] != 'NA':
		return best_blast[1], 'no orthogroup BLAST'
	elif best_blast[1] == 'NA':
		return best_blast[1], 'no hit of any kind!'

	in_orto.close()


# Run commands
in_fasta = open(args.f, 'r')
out_file = open(args.o, 'w')

# for each gene (keep gene name and sequence from fasta)
for line in in_fasta:
	if '>' in line:
		gene_name = clean_line(line)
		print('\n'+gene_name)
		homo, reason = find_homo(gene_name)
		out_file.write (f'{gene_name}\t{homo}\t{reason}\n')


