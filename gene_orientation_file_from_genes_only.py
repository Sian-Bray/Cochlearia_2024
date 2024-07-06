# Generates a gene_orientation_file.txt file from the genes_only.gtf file

# e.g. in:
	# Cexcelsa_scaf_4	AUGUSTUS	gene	1	9453	1	+	.	ID=g1

# e.g. target:
	# g1	scaffold_33	1	2616	+

# cd /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/02_Visualisation_Folder/00_Blank_Folder
# python3 /Users/sian_bray/Dropbox/Scripts/gene_orientation_file_from_genes_only.py

in_file = open('genes_only.gtf', 'r')
out_file = open('gene_orientation_file.txt', 'w')

for line in in_file:
	line=line.replace('\n', '')
	line = line.split('\t')
	gene=line[8]
	gene=gene.replace('ID=', '')
	new_line=f'{gene}\t{line[0]}\t{line[3]}\t{line[4]}\t{line[6]}\n'
	out_file.write(new_line)




