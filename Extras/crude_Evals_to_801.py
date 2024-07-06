# Add BLAST E-values to 801 list
# python3 /Users/sian_bray/Dropbox/Scripts/crude_Evals_to_801.py

gene_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/1perc_fst_outliers_annotated.tsv', 'r')
all_vs_all = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/1-2-1_annotation/cochlearia_thaliana_all_Vs_all.tsv', 'r')
out_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/801_BLAST_Evals.tsv', 'w')

genes = {}

# Make a dir of the genes
# g22736.t1	AT1G48580.2	no orthogroup BLAST	Note=None: conf_class=5: computational_description=nuclear localized protein: conf_rating=***: symbol=None: curator_summary=: 
for line in gene_file:
	split_line = line.split('\t')
	if 'BLAST' in split_line[2]:
		genes[split_line[0]] = [split_line[1], 'XXX', line]
	if 'BLAST' not in split_line[2]:
		out_file.write(line.replace('\n', ''))
		out_file.write('\tNA\n')

for line in all_vs_all:
	split_line = line.split('\t')
	for gene in genes:
		if (gene == split_line[0]) and (genes[gene][0] == split_line[1]):
			if 'XXX' == genes[gene][1]:
				genes[gene][1] = split_line[10]
			if float(genes[gene][1]) > float(split_line[10]):
				genes[gene][1] = split_line[10]

for gene in genes:
	line = genes[gene][2]
	line = line.replace('\n', '')
	out_file.write(f'{line}\t{genes[gene][1]}\n')








# Pull out the best E-val for the gene




