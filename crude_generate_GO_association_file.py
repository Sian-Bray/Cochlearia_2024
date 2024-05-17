# python3 /Users/sian_bray/Dropbox/Scripts/crude_generate_GO_association_file.py

GO_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/gene_association.tair', 'r')
assoc_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Physics_GWAS/Data/thaliana_GO_associations.txt', 'w')
gff_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Araport11_GFF3_genes_transposons.2023-01-02.gff', 'r')
all_genes_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Physics_GWAS/Data/all_thaliana_genes.txt', 'w')

gene_dict = {}

for line in GO_file:
	if '!' not in line:
		gene = line.replace('\n', '')
		gene = line.split('\t')
		GO_term = gene[4]
		gene = gene[1]
		if gene not in gene_dict:
			gene_dict[gene] = GO_term
		if gene in gene_dict:
			GO_terms = f'{gene_dict[gene]};{GO_term}'
			gene_dict[gene] = GO_terms

for key in gene_dict:
	value=gene_dict[key]
	assoc_file.write(f'{key}\t{value}\n')

all_genes=[]
for line in gff_file:
	if 'ID=' in line:
		line = line.split('\t')
		line = line[8]
		line = line.replace(':', ';')
		line = line.split(';')
		gene = line[0]
		gene = gene.replace('ID=', '')
		gene = gene.split('.')
		gene = gene[0]
		if 'AT' in gene:
			if 'ATM' not in gene:
				if 'ATC' not in gene:
					if 'TE' not in gene:
						all_genes.append(gene)

all_genes = list(dict.fromkeys(all_genes))

for gene in all_genes:
	all_genes_file.write(gene+'\n')


GO_file.close()
assoc_file.close()
gff_file.close()
all_genes_file.close()

