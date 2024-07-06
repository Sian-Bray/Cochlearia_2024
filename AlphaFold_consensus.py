# Generates batch files to create consensus sequences for AlphaFold.
# These sequences need to be manually checked against non-biallelic sites...
# ... and manually translated (reverse complimented where nessasary).

# Written by Sian Bray on Firday 19th May 2023

# Testing command (local):
	# python3 /Users/sian_bray/Dropbox/Scripts/AlphaFold_consensus.py -v my_vcf.vcf -r /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/C_excelsa_V5.fasta -g /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/C_excelsa_V5_braker2_wRseq_ANNOTATION/C_excelsa_V5_braker2_wRseq.gff3 -i /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/new_top_1perc_nr.txt -d /Users/sian_bray -p tetraploids -s ALO_001 ALO_002 ALO_003

# Real HPC command:
	# python3 /gpfs01/home/sbzsmb/Scripts/AlphaFold_consensus.py -v /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/10.filtered.depth/reheadered.F4_133.ANNOTATED.vcf -r /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ionops/ref.gen/C_excelsa_V5.fasta -g /gpfs01/home/sbzsmb/Cochlearia/C_excelsa_V5_braker2_wRseq.gff3 -i /gpfs01/home/sbzsmb/Cochlearia/new_top_1perc_nr.txt -d /gpfs01/home/sbzsmb/Cochlearia/AlphaFold/Consensus -p tetraploids -s ROT_004 ROT_006 ROT_007 ROT_013 SKN_001 SKN_002 SKN_005 SKN_008 ALO_006 ALO_007 ALO_013 ALO_017 ELI_001 ELI_002 ELI_003 ELI_004 ERS_1 ERS_2 ERS_3 ERS_4 FTW_1 FTW_2 FTW_3 FTW_5 LAL_1 LAL_2 LAL_3 LAL_4 LNL_001 LNL_002 LNL_003 LNL_008 LOS_1 LOS_2 LOS_6 LOS_7 NEI_1 NEI_3 NEI_8 NEI_9 SCU_14 SCU_15 SCU_16 SCU_19
	# python3 /gpfs01/home/sbzsmb/Scripts/AlphaFold_consensus.py -v /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/10.filtered.depth/reheadered.F4_133.ANNOTATED.vcf -r /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ionops/ref.gen/C_excelsa_V5.fasta -g /gpfs01/home/sbzsmb/Cochlearia/C_excelsa_V5_braker2_wRseq.gff3 -i /gpfs01/home/sbzsmb/Cochlearia/new_top_1perc_nr.txt -d /gpfs01/home/sbzsmb/Cochlearia/AlphaFold/Consensus -p diploids -s  BNK_21 CHA_1 CHA_2 JOR_1 JOR_12 JOR_13 JOR_3 LAB_004 LAB_1 LAB_2 LAB_300 LAB_4 LAB_400 LAB_5 LAB_500 NEN_001 NEN_003 NEN_200 NEN_300 NEN_4 NEN_5 NEN_6 ODN_10 ODN_2 ODN_4 ODN_5 ODN_6 ODN_7 ODN_9

# Take input files
import argparse
parser = argparse.ArgumentParser(description="Align by ID.")
parser.add_argument('-v', type=str, metavar='vcf_file', required=True, help='VCF file that contains all populations.')
parser.add_argument('-r', type=str, metavar='ref_fasta', required=True, help='Refernce fasta file.')
parser.add_argument('-g', type=str, metavar='gff_file', required=True, help='A gff3 file that contains exon locations.')
parser.add_argument('-i', type=str, metavar='gene_list', required=True, help='The input gene list. One gene per line.')
parser.add_argument('-d', type=str, metavar='batch_dir', required=True, help='Directory to write output files to.')
parser.add_argument('-p', type=str, metavar='pop_name', required=True, help='Population name, e.g. tetraploids.')
parser.add_argument('-s', type=str, metavar='samples', nargs='+', default=[], required=True, help='Sample (individual) names in the population, e.g. ALO_001, ALO_002, ALO_003.')
args = parser.parse_args()

# Make a dictionary of genes using the gene list and the gff
# Dictionary contains gene name as the key and a list of exon positions...
# ...in gatk -L flag format (scaffold_6:23043016-23043109). A gff line:
	# Cexcelsa_scaf_1	AUGUSTUS	exon	1512310	1513122	.	+	.	ID=g10041.t1.exon1;Parent=g10041.t1;
# Also keep a list with the gene orientation - write this to the file name later for the manual steps

# Make a list of genes
gene_file = open(args.i, 'r')
gene_list=[]

for line in gene_file:
	gene = line.replace('\n', '')
	gene_list.append(gene)

gene_file.close()

# Make the dictionary
gff = open(args.g, 'r')
position_dict = {}

for gene in gene_list:
	position_dict[gene] = []

for feature in gff:
	feature = feature.replace('\n', '')
	feature = feature.split('\t')
	info = feature[8]
	info = info.replace('ID=', '')
	info = info.replace('Parent=', '')
	info = info.replace('.', ';')
	info = info.split(';')
	info = info[0]

	if (feature[2] == 'exon') and (info in gene_list):
		pos_list = position_dict[info]
		pos_list.append(f'{feature[0]}:{feature[3]}-{feature[4]}')
		position_dict[info] = pos_list

gff.close()

# Make a new batch file for each gene.

individuals = args.s
pop_name = args.p

out_dir = args.d
if out_dir[-1] != '/':
	out_dir = out_dir+'/'

for gene in gene_list:

	gene_exons = position_dict[gene]

	batch = open(f'{out_dir}{gene}_{pop_name}_consensus_generation.sh', 'w')

	# Write the batch headder
	batch.write('#!/bin/bash\n')
	batch.write(f'#SBATCH --job-name={gene}_{pop_name}_consensus\n')
	batch.write('#SBATCH --partition=defq\n')
	batch.write('#SBATCH --nodes=1\n')
	batch.write('#SBATCH --ntasks-per-node=1\n')
	batch.write('#SBATCH --mem=16g\n')
	batch.write('#SBATCH --time=8:00:00\n')
	batch.write('#SBATCH --output=/gpfs01/home/sbzsmb/OandE/%x.out\n')
	batch.write('#SBATCH --error=/gpfs01/home/sbzsmb/OandE/%x.err\n')
	batch.write('source $HOME/.bash_profile\n')
	batch.write('conda activate ngs_pipe_env_gatk4\n')
	batch.write('ml bcftools-uoneasy/1.17-GCC-12.2.0\n')
	# use gatk to create consensus sequence for that region, filter for biallelic and AF > 0.50
	batch.write(f'gatk SelectVariants -R {args.r} -V {args.v} --restrict-alleles-to BIALLELIC --max-nocall-fraction 0.2 -select "AF > 0.5"')
	for ind in individuals:
		batch.write(f' -sn {ind}')
	for exon in gene_exons:
		batch.write(f' -L {exon}')
	batch.write(f' -O {out_dir}{gene}_{pop_name}_biallelic.vcf\n')

	# keep a vcf of the sites that are filtered out for a manual check
	batch.write(f'gatk SelectVariants -R {args.r} -V {args.v} --restrict-alleles-to MULTIALLELIC --max-nocall-fraction 0.2 --exclude-non-variants')
	for ind in individuals:
		batch.write(f' -sn {ind}')
	for exon in gene_exons:
		batch.write(f' -L {exon}')
	batch.write(f' -O {out_dir}{gene}_{pop_name}_multiallelic.vcf\n')

	# Convert vcfs to bgzipped files and create an index
	batch.write(f'bgzip {out_dir}{gene}_{pop_name}_biallelic.vcf\n')
	batch.write(f'tabix -p vcf {out_dir}{gene}_{pop_name}_biallelic.vcf.gz\n')

	# Generate consensus with bcftools
	batch.write(f'samtools faidx {args.r}')
	for exon in gene_exons:
		batch.write(f' {exon}')

	batch.write(f' | bcftools consensus -s - {out_dir}{gene}_{pop_name}_biallelic.vcf.gz -o {out_dir}{gene}_{pop_name}_biallelic.fasta\n')

	# Concatonate the exons
	# Forward
	batch.write(f'grep -v ">" {out_dir}{gene}_{pop_name}_biallelic.fasta > {out_dir}{gene}_{pop_name}_biallelic_concat.fasta\n')












