# cd /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data
# python3 crude_MAVs_to_genes.py

levi_genes = open('new_top_1perc_nr.txt', 'r')
sian_MAVs = open('genes_UK_dips_Vs_tets_15_BPM_1per_plus_MAVs.txt', 'r')
out_file = open('genes_1perc_MAVs.txt', 'w')

MAV_dict = {}

# gene	outname	scaff	start	end	win_size	num_sites	num_snps	Rho	FstWC	dxy	AFD	FixedDiff	FstH	FstN	MAVcount
for line in sian_MAVs:
	line = line.replace('\n', '')
	line = line.split('\t')
	MAV_dict[line[0]] = line[15]

gene_list = []

for gene in levi_genes:
	if 'CochID' in gene:
		continue
	gene = gene.replace('\n', '')
	#gene = gene.split('.')
	#gene = gene[0]
	gene_list.append(gene)


for gene in gene_list:
	try:
		out_file.write(f'{gene}\t{MAV_dict[gene]}\n')
	except KeyError:
		pass



