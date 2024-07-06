# Code to perform permutation test on GO overlaps from different species
# Example run command:
	# python3 /gpfs01/home/sbzsmb/Scripts/GO_iteration_test_Cochlearia_2024_HPC.py -l1 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/all_genes_transcript_numbers.tsv -l2 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/all_genes.csv -l3 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/all_genes.csv -pl1 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/genes_in_scan_clean.txt -n1 752 -n2 451 -n3 227 -go1 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/Cochlearia_Thaliana_GO_universe_restrictive.tsv -go2 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/gene_GO.tsv -go3 /gpfs01/home/sbzsmb/GO_Permutation_Test/Input_Materials/gene_GO.tsv -p 10 -x_1_2 17 -x_1_3 1 -x_2_3 5 -x_all 0 -o /gpfs01/home/sbzsmb/GO_Permutation_Test/Test_1/

# import modules
import random, argparse, datetime, os, time, subprocess


# Function to convert input file (list of gene names, one per line) into a Python list
def in_file_to_list(in_file):
	open_file = open(in_file, 'r')
	new_list = []
	for line in open_file:
		line = line.replace('\n', '')
		new_list.append(line)
	open_file.close()
	return new_list


# Has a possible gene list been provided? If not use the full gene list.
def possible_or_full(possible_list, full_list):
	if possible_list != None:
		return possible_list
	else:
		return full_list


# Random number generation
# Working test example
	# population = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
	# random.sample(population, 3)

# Function to generate a psudo-random subset of files
# population is a list of genes
# subset_size is an integer
def random_subset(population, subset_size):
	subset = random.sample(population, subset_size)
	return subset


# Save a python list as a file, one item per line
def write_py_list(listy, file_name):
	written_file = open(file_name, 'w')
	for item in listy:
		written_file.write(item+'\n')
	written_file.close()


# Paste together a file name and directory, both are strings
# The code checks that the directory end s with a /
def create_file_name(file, directory):
	if directory[-1:] == '/':
		file_name = directory+file
	else:
		file_name = directory+'/'+file
	return file_name


# Write the file containing the R code required to generate GO enrichment tables for the randomised gene sets
# Example R GO file generated for this script is called: R_for_permutation_tests.R
	# GO_universe = e.g. Cochlearia_Thaliana_GO_universe_restrictive.tsv
		# A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.
	# scan_genes = e.g. genes_in_scan_clean.txt
		# A txt file containing all possible genes in the scan (e.g. all the genes that appear in the BPM metrics file and so have enough coverage and SNPs) with a single gene per line. No headder.
	# candidate_genes = e.g. top_1perc_minus_inversion.txt
		# A txt file containing all the candidate genes detected in the selection scan. One gene per line, no headder.
# Note: gene names in the three files above must be identical, i.e. all must contain g100.t1 not g100.
	# threshold = float, p-value cutoff, default is 0.05
	# iteration = string, independent identifier for each round of GOs 
def GO_enrichment_R(R_file, GO_universe, scan_genes, candidate_genes, out_name, out_dir, threshold=0.05):
	V1 = "'V1'"
	out_R_file = open(R_file, 'w')
	out_R_file.write(f'# load topGO\n'+
	f'library("topGO")\n'+
  	f'\n'+
  	f'# Open my files\n'+
  	f'\n'+
  	f'# This function converts the file to the gene2GO format that topGO requires\n'+
  	f'go <- readMappings(file="{GO_universe}", \n'+
  	f'                   sep = "\t", IDsep = ",")\n'+
  	f'\n'+
  	f'# All the genes in my scan (i.e. all the genes in regions from the bpm file)\n'+
  	f'genesInScan <- read.csv("{scan_genes}", \n'+
  	f'                        header = FALSE)\n'+
  	f'\n'+
  	f'# My genes of interest (i.e. Fst top 1%)\n'+
  	f'top1perc <- read.csv("{candidate_genes}", \n'+
  	f'                     header = FALSE)\n'+
  	f'\n'+
  	f'# Create a factor where the indicies are the gene names and the factors are 1/0 for "candidate gene"/"scan gene but not candidate gene"\n'+
  	f'\n'+
  	f'# First convert top1perc from a dataframe to a vector\n'+
  	f'top1perc_vector <- as.vector(top1perc[,{V1}])\n'+
  	f'\n'+
  	f'# Then convert genesInScan to a vector\n'+
  	f'genesInScan_vector <- as.vector(genesInScan[,{V1}])\n'+
  	f'\n'+
  	f'# Then make the factor\n'+
  	f'allGenes_factor <- factor(as.integer(genesInScan_vector %in% top1perc_vector))\n'+
  	f'names(allGenes_factor) = genesInScan_vector\n'+
  	f'\n'+
  	f'# Website used for help: https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/\n'+
  	f'# Make the topGOdata object\n'+
  	f'go_data_BP <- new("topGOdata", \n'+
  	f'                  ontology="BP", \n'+
  	f'                  allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan\n'+
  	f'                  annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"\n'+
  	f'                  gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms\n'+
  	f'                  nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene\n'+
  	f'\n'+
  	f'go_data_CC <- new("topGOdata", \n'+
  	f'                  ontology="CC", \n'+
  	f'                  allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan\n'+
  	f'                  annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"\n'+
  	f'                  gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms\n'+
  	f'                  nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene\n'+
  	f'\n'+
  	f'go_data_MF <- new("topGOdata", \n'+
  	f'                  ontology="MF", \n'+
  	f'                  allGenes=allGenes_factor, # This needs to be a factor, where 1 is genes in my candidate list and 0 is genes in my scan\n'+
  	f'                  annot=annFUN.gene2GO, # "function which maps genes identifiers to GO terms"\n'+
  	f'                  gene2GO=go, # This is my list of named vectors where the gene name is the name and the vector is all the GO terms\n'+
  	f'                  nodeSize=5) # an integer larger or equal to 1. This parameter is used to prune the GO hierarchy from the terms which have less than nodeSize annotated gene\n'+
  	f'\n'+
  	f'\n'+
  	f'# Run the GO analysis classic fisher\n'+
  	f'res_BP_classicFisher <- runTest(go_data_BP, statistic = "fisher")\n'+
  	f'res_CC_classicFisher <- runTest(go_data_CC, statistic = "fisher")\n'+
  	f'res_MF_classicFisher <- runTest(go_data_MF, statistic = "fisher")\n'+
  	f'\n'+
  	f'# Run the GO analysis conserative elim fisher\n'+
  	f'resBP_elimFisher <- runTest(go_data_BP, algorithm = "elim", statistic = "fisher")\n'+
  	f'resCC_elimFisher <- runTest(go_data_CC, algorithm = "elim", statistic = "fisher")\n'+
  	f'resMF_elimFisher <- runTest(go_data_MF, algorithm = "elim", statistic = "fisher")\n'+
  	f'\n'+
  	f'###########################\n'+
  	f'\n'+
  	f'# Extract p-values from the results\n'+
  	f'#p_values_BP_classicFisher <- score(res_BP_classicFisher)\n'+
  	f'#p_values_CC_classicFisher <- score(res_CC_classicFisher)\n'+
  	f'#p_values_MF_classicFisher <- score(res_MF_classicFisher)\n'+
  	f'\n'+
  	f'p_values_BP_elimFisher <- score(resBP_elimFisher)\n'+
  	f'p_values_CC_elimFisher <- score(resCC_elimFisher)\n'+
  	f'p_values_MF_elimFisher <- score(resMF_elimFisher)\n'+
  	f'\n'+
  	f'# Define the significance threshold\n'+
  	f'threshold <- {threshold}\n'+
  	f'\n'+
  	f'# Count the number of significant genes\n'+
  	f'#significant_genes_BP_classicFisher <- sum(p_values_BP_classicFisher < threshold)\n'+
  	f'#significant_genes_CC_classicFisher <- sum(p_values_CC_classicFisher < threshold)\n'+
  	f'#significant_genes_MF_classicFisher <- sum(p_values_MF_classicFisher < threshold)\n'+
  	f'\n'+
  	f'significant_genes_BP_elimFisher <- sum(p_values_BP_elimFisher < threshold)\n'+
  	f'significant_genes_CC_elimFisher <- sum(p_values_CC_elimFisher < threshold)\n'+
  	f'significant_genes_MF_elimFisher <- sum(p_values_MF_elimFisher < threshold)\n'+
  	f'\n'+
  	f'###########################\n'+
  	f'\n'+
  	f'# Find the number of significant nodes at {threshold} cutoff\n'+
  	f'####\n'+
  	f'\n'+
  	f'# Show combined table for {threshold} sig cutoff\n'+
	f'if (significant_genes_BP_elimFisher == 0) {"{"}\n'
	f'  allResBP <- data.frame()\n'
	f'{"}"} else if (significant_genes_BP_elimFisher == 1) {"{"}\n'
	f'  allResBP <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,\n'
	f'                       elimFisher = resBP_elimFisher,\n'
	f'                       orderBy = "elimFisher",\n'
	f'                       topNodes = 2,\n'
	f'                       numChar = 250)\n'
	f'  allResBP <- allResBP[-nrow(allResBP), ]\n'
	f'{"}"} else if (significant_genes_BP_elimFisher > 1) {"{"}\n'
	f'  allResBP <- GenTable(go_data_BP, classicFisher = res_BP_classicFisher,\n'
	f'                       elimFisher = resBP_elimFisher,\n'
	f'                       orderBy = "elimFisher",\n'
	f'                       topNodes = significant_genes_BP_elimFisher,\n'
	f'                       numChar = 250)\n'
	f'{"}"}\n'
	f'if (significant_genes_CC_elimFisher == 0) {"{"}\n'
	f'  allResCC <- data.frame()\n'
	f'{"}"} else if (significant_genes_CC_elimFisher == 1) {"{"}\n'
	f'  allResCC <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,\n'
	f'                       elimFisher = resCC_elimFisher,\n'
	f'                       orderBy = "elimFisher",\n'
	f'                       topNodes = 2,\n'
	f'                       numChar = 250)\n'
	f'  allResCC <- allResCC[-nrow(allResCC), ]\n'
	f'{"}"} else if (significant_genes_CC_elimFisher > 1) {"{"}\n'
	f'  allResCC <- GenTable(go_data_CC, classicFisher = res_CC_classicFisher,\n'
	f'                       elimFisher = resCC_elimFisher,\n'
	f'                       orderBy = "elimFisher",\n'
	f'                       topNodes = significant_genes_CC_elimFisher,\n'
	f'                       numChar = 250)\n'
	f'{"}"}\n'
	f'if (significant_genes_MF_elimFisher == 0) {"{"}\n'
	f'  allResMF <- data.frame()\n'
	f'{"}"} else if (significant_genes_MF_elimFisher == 1) {"{"}\n'
	f'  allResMF <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,\n'
	f'                       elimFisher = resMF_elimFisher,\n'
	f'                       orderBy = "elimFisher",\n'
	f'                       topNodes = 2,\n'
	f'                       numChar = 250)\n'
	f'  allResMF <- allResMF[-nrow(allResMF), ]\n'
	f'{"}"} else if (significant_genes_MF_elimFisher > 1) {"{"}\n'
	f'  allResMF <- GenTable(go_data_MF, classicFisher = res_MF_classicFisher,\n'
	f'                       elimFisher = resMF_elimFisher,\n'
	f'                       orderBy = "elimFisher",\n'
	f'                       topNodes = significant_genes_MF_elimFisher,\n'
	f'                       numChar = 250)\n'
	f'{"}"}\n'
  	f'\n'+
  	f'# since the elim gives me better results and it was the method used in ...\n'+
  	f'# ... Madjas paper I will go ahead with this method and not classic\n'+
  	f'\n'+
  	f'# I think one of the previous problems might be the MTC, the top GO manual ...\n'+
  	f'# ... states that "in many cases adjusted p-values might be misleading".\n'+
  	f'\n'+
  	f'# Write the combined tables to files for {threshold} cutoff\n'+
  	f'write.csv(allResBP, "{out_dir}allResBP_{out_name}")\n'+
  	f'write.csv(allResCC, "{out_dir}allResCC_{out_name}")\n'+
  	f'write.csv(allResMF, "{out_dir}allResMF_{out_name}")\n')
	out_R_file.close()


# Write batch files containing R scripts
	# R_scripts is a lits of absolute file paths
	# sbatch_name is the name of the sbatch file created
# mem and time etc are designed for 10 R scripts per batch file
def create_sbatch(R_script, sbatch_name):
	file_written = open(sbatch_name, 'w')
	file_written.write(f'#!/bin/bash\n'+
	f'#SBATCH --job-name={sbatch_name[:-3]}\n'+
	f'#SBATCH --partition=shortq\n'+
	f'#SBATCH --nodes=1\n'+
	f'#SBATCH --ntasks-per-node=1\n'+
	f'#SBATCH --mem=4g\n'+
	f'#SBATCH --time=00:05:00\n'+
	f'#SBATCH --output=/gpfs01/home/sbzsmb/OandE/%x.out\n'+
	f'#SBATCH --error=/gpfs01/home/sbzsmb/OandE/%x.err\n'+
	f'	\n'+
	f'source $HOME/.bash_profile\n'+
	f'conda activate R_GO_env\n'+
	f'\n'+
	f'# Run these R files\n'+
	f'Rscript {R_script}\n')
	file_written.close()


# Required for run and monitor
def terminal_pipe(cmd): 
    return subprocess.Popen(f'{cmd}', shell=True, stdout=subprocess.PIPE).communicate()[0].decode("utf-8").strip(' \n')


# Matts run and monitor!
def run_and_monitor(sbatch_directory, max_jobs=50):
    # get list of all sbatch files in directory
    batch_jobs=os.listdir(sbatch_directory)
    # get rid of anything that is not a shell file
    for jobs in batch_jobs:
        if '.sh' not in jobs:
            batch_jobs.remove(jobs)
    user = terminal_pipe(f'echo $USER')
    while len(batch_jobs) > 0:
        running=terminal_pipe(f'squeue -u {user} -h | wc -l')
        if int(running) < max_jobs:
            # run another script
            os.system(f'sbatch {sbatch_directory}{batch_jobs[0]}')
            # remove script from batch jobs list
            batch_jobs.remove(batch_jobs[0])
        # when there are 100 scripts running wait for 10 mins
        if int(running) >= max_jobs:
            time.sleep(60)


# Convert GO results files into a single list of GO terms
# files = a list containing the names of the files
def GO_file_to_py_list(files):
	GO_terms = []
	for file in files:
		current_file = open(file, 'r')
	for count, line in enumerate(current_file):
		if (line != '') and (line != '\n') and (line != '""') and (count != 0):
			line = line.replace('"', '')
			line = line.split(',')
			GO_terms.append(line[1])
	return GO_terms


# Compare two GO_enrichment lists
def GO_two_way(GO_terms_1, GO_terms_2):
	in_both = 0
	for item in GO_terms_1:
		if item != []:
			if item in GO_terms_2:
				in_both += 1
	return in_both


# Compare three GO enrichment lists
def GO_three_way(GO_terms_1, GO_terms_2, GO_terms_3):
	in_all = 0
	for item in GO_terms_1:
		if item != []:
			if (item in GO_terms_2) and (item in GO_terms_3):
				in_all += 1
	return in_all


# Compare an original number of overlaps with the permutation test output
def permutation_comparison(original_overlap, permutation_out_list):
	greater_or_equal = 0
	for overlaps in permutation_out_list:
		if overlaps >= original_overlap:
			greater_or_equal += 1
	p_value = greater_or_equal / len(permutation_out_list)
	return greater_or_equal, p_value


# Run the code!
if __name__ == "__main__":

	# Save the start time
	start_time = datetime.datetime.now()

	# Take command line input
	parser = argparse.ArgumentParser(description="Permutation test for GO functional catagory overlaps between three selection scans. For example the code was written for the selection scans in response to WGD in Cochlearia, arenosa and Cardamine.")
	
	# Three total gene lists (Cochlearia, Cardamine and arenosa)
	parser.add_argument('-l1', type=str, metavar='list_1', required=True, help='Path to the first list of total genes. This should be a list of all the genes in the genome, one gene per line.')
	parser.add_argument('-l2', type=str, metavar='list_2', required=True, help='Path to the second list of total genes. This should be a list of all the genes in the genome, one gene per line.')
	parser.add_argument('-l3', type=str, metavar='list_3', required=True, help='Path to the third list of total genes. This should be a list of all the genes in the genome, one gene per line.')

	# Lists of genes that are possible to find in the selection scan
		# Can I make the full list the default if not included?
	parser.add_argument('-pl1', type=str, metavar='possible_list_1', required=False, help='Path to the first list of total possible genes. That is all the genes that could possibly come up as candidate gene in the scan (e.g. enough coverage and SNPs). Csv file, one gene per line. If not provided the code defaults to the total gene list (-1).')
	parser.add_argument('-pl2', type=str, metavar='possible_list_2', required=False, help='Path to the second list of total possible genes. That is all the genes that could possibly come up as candidate gene in the scan (e.g. enough coverage and SNPs). Csv file, one gene per line. If not provided the code defaults to the total gene list (-1).')
	parser.add_argument('-pl3', type=str, metavar='possible_list_3', required=False, help='Path to the third list of total possible genes. That is all the genes that could possibly come up as candidate gene in the scan (e.g. enough coverage and SNPs). Csv file, one gene per line. If not provided the code defaults to the total gene list (-1).')

	# Lengths of the initial condidate lists for all three (Cochlearia, Cardamine and arenosa)
	parser.add_argument('-n1', type=int, metavar='number_1', required=True, help='Lenght of the candidate gene lists from the first selection scan. This should be an integer.')
	parser.add_argument('-n2', type=int, metavar='number_2', required=True, help='Lenght of the candidate gene lists from the second selection scan. This should be an integer.')
	parser.add_argument('-n3', type=int, metavar='number_3', required=True, help='Lenght of the candidate gene lists from the third selection scan. This should be an integer.')

	# GO universe used for all three (Cochlearia, Cardamine and arenosa)
	parser.add_argument('-go1', type=str, metavar='GO_universe_1', required=True, help='Full absolute path to the GO universe for the first genome/selection scan. A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.')
	parser.add_argument('-go2', type=str, metavar='GO_universe_2', required=True, help='Full absolute path to the GO universe for the second genome/selection scan. A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.')
	parser.add_argument('-go3', type=str, metavar='GO_universe_3', required=True, help='Full absolute path to the GO universe for the third genome/selection scan. A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.')

	# How many permutations to run
	parser.add_argument('-p', type=int, metavar='number_of_permutations', required=False, default=1000, help='The number of permutations to run. This should be an integer. Default is 1000.')

	# P-value threshold
	parser.add_argument('-t', type=float, metavar='threshold', required=False, default=0.05, help='The significance cutoff/P-value to be used in the GO enrichment analysis.')

	# Number of GO overlaps in the initial contrasts
	parser.add_argument('-x_1_2', type=int, metavar='overlap', required=True, help='The amount of overlap (cross over) between GO terms in the original contrast between species 1 and 2. An integer.')
	parser.add_argument('-x_1_3', type=int, metavar='overlap', required=True, help='The amount of overlap (cross over) between GO terms in the original contrast between species 1 and 3. An integer.')
	parser.add_argument('-x_2_3', type=int, metavar='overlap', required=True, help='The amount of overlap (cross over) between GO terms in the original contrast between species 3 and 2. An integer.')
	
	parser.add_argument('-x_all', type=int, metavar='overlap', required=True, help='The amount of overlap (cross over) between all three GO terms in the original contrast. An integer.')

	# Output directory
	parser.add_argument('-o', type=str, metavar='output_directory', required=True, help='Path the output directory. I recomend an empty directory to contain the results of each permutation.')

	# Add arguments to run the three stages seperatly?
	# 1) Generate random lists and R scripts?
	# 2) Run R scripts with run and monitor?
	# 3) Calculate metrics from complete runs

	args = parser.parse_args()

	# Create variables to use in the following loop

	all_genes_1 = in_file_to_list(args.l1)
	all_genes_2 = in_file_to_list(args.l2)
	all_genes_3 = in_file_to_list(args.l3)

	candidate_list_1_len = args.n1
	candidate_list_2_len = args.n2
	candidate_list_3_len = args.n3

	GO_universe_1 = args.go1
	GO_universe_2 = args.go2
	GO_universe_3 = args.go3

	scan_genes_1 = possible_or_full(args.pl1, args.l1)
	scan_genes_2 = possible_or_full(args.pl2, args.l2)
	scan_genes_3 = possible_or_full(args.pl3, args.l3)

	permutations = args.p

	threshold = args.t

	out_dir = create_file_name('', args.o)

	two_way_1_2_all = []
	two_way_1_3_all = []
	two_way_2_3_all = []

	three_way_all = []

	R_scripts = []

	# For each permutaion
	while permutations > 0:

		# Generate random 'candidate gene' lists
		random_candidates_1 = random_subset(all_genes_1, args.n1)
		random_candidates_2 = random_subset(all_genes_2, args.n2)
		random_candidates_3 = random_subset(all_genes_3, args.n3)

		# Save random 'candidate gene' lists to file
		candidates_file_1 = out_dir+'random_genes_species1_round'+str(permutations)+'.txt'
		candidates_file_2 = out_dir+'random_genes_species2_round'+str(permutations)+'.txt'
		candidates_file_3 = out_dir+'random_genes_species3_round'+str(permutations)+'.txt'

		write_py_list(random_candidates_1, candidates_file_1)
		write_py_list(random_candidates_2, candidates_file_2)
		write_py_list(random_candidates_3, candidates_file_3)

		# Create names for the files you will create next
		prefix_1 = 'species_1_'+str(permutations)
		prefix_2 = 'species_2_'+str(permutations)
		prefix_3 = 'species_3_'+str(permutations)

		R_file_1 = create_file_name(prefix_1+'.R', out_dir)
		R_file_2 = create_file_name(prefix_2+'.R', out_dir)
		R_file_3 = create_file_name(prefix_3+'.R', out_dir)

		out_name_1 = prefix_1+'.csv'
		out_name_2 = prefix_2+'.csv'
		out_name_3 = prefix_3+'.csv'

		# Get GOs for each 'candidate gene' list
		GO_enrichment_R(R_file_1, GO_universe_1, scan_genes_1, candidates_file_1, out_name_1, out_dir, threshold=0.05)
		GO_enrichment_R(R_file_2, GO_universe_2, scan_genes_2, candidates_file_2, out_name_2, out_dir, threshold=0.05)
		GO_enrichment_R(R_file_3, GO_universe_3, scan_genes_3, candidates_file_3, out_name_3, out_dir, threshold=0.05)

		### Replace below with run and monitor for a HPC version of this code. ###

		# Run the R files
		# os.system('Rscript '+R_file_1)
		# os.system('Rscript '+R_file_2)
		# os.system('Rscript '+R_file_3)

		create_sbatch(R_file_1, R_file_1[:-2]+'.sh')
		create_sbatch(R_file_2, R_file_2[:-2]+'.sh')
		create_sbatch(R_file_3, R_file_3[:-2]+'.sh')

		### Break this loop in two here so that it can 1) generate R files (above), 
		### 2) run them with 'run and monitor' (next), then 3) (below in the second
		### loop) extract the info.

		# Next permutation!
		permutations -= 1

	# Run and monitor untill all the R scripts have been run	
	run_and_monitor(out_dir)

	# Wait until all the jobs are finished
	user = terminal_pipe(f'echo $USER')
	running=terminal_pipe(f'squeue -u {user} -h | wc -l')
	
	# i.e. until this is the only batch job running
	while int(running) > 1:
		time.sleep(60)
		running=terminal_pipe(f'squeue -u {user} -h | wc -l')

	# Check that all the scripts have run and report any that are missing?

	# Reset the permutations
	permutations = args.p

	# New loop
	while permutations > 0:

		# Things I need to re-generate in the new loop

		# re-create names for the files you will use next
		prefix_1 = 'species_1_'+str(permutations)
		prefix_2 = 'species_2_'+str(permutations)
		prefix_3 = 'species_3_'+str(permutations)

		out_name_1 = prefix_1+'.csv'
		out_name_2 = prefix_2+'.csv'
		out_name_3 = prefix_3+'.csv'

		# Convert output files to python lists
		GO_files_1 = [out_dir+'allResBP_'+out_name_1, out_dir+'allResCC_'+out_name_1, out_dir+'allResMF_'+out_name_1]
		GO_files_2 = [out_dir+'allResBP_'+out_name_2, out_dir+'allResCC_'+out_name_2, out_dir+'allResMF_'+out_name_2]
		GO_files_3 = [out_dir+'allResBP_'+out_name_3, out_dir+'allResCC_'+out_name_3, out_dir+'allResMF_'+out_name_3]

		GO_list_1 = GO_file_to_py_list(GO_files_1)
		GO_list_2 = GO_file_to_py_list(GO_files_2)
		GO_list_3 = GO_file_to_py_list(GO_files_3)

		# Get GO overlaps for this iteration
		two_way_1_2 = GO_two_way(GO_list_1, GO_list_2)
		two_way_1_3 = GO_two_way(GO_list_1, GO_list_3)
		two_way_2_3 = GO_two_way(GO_list_2, GO_list_3)

		three_way = GO_three_way(GO_list_1, GO_list_2, GO_list_3)

		# Save to overlaps lists
		two_way_1_2_all.append(two_way_1_2)
		two_way_1_3_all.append(two_way_1_3)
		two_way_2_3_all.append(two_way_2_3)

		three_way_all.append(three_way)

		# Next permutation!
		permutations -= 1

	# Save the end time
	end_time = datetime.datetime.now()

	# How frequent is the experimental level of overlap Vs the permutations?

	master_out_file = open(out_dir+'000_master_output_file.txt', 'w')

	master_out_file.write(f'Sumary statistics for the GO overlap permutation test run started on {start_time} and finished on {end_time}.'+
		f'\n\nOutput files can be found in {out_dir}.\n'+
		
		f'\n\tSpecies 1 input files were:'+
		f'\n\t\tTotal gene list: {args.l1}'+
		f'\n\t\tPossible candidate genes list: {args.pl1}'+
		f'\n\t\tGO universe: {args.go1}'+
		f'\n\t\tNumber of genes in the initial candidate list: {args.n1}\n'+

		f'\n\tSpecies 2 input files were:'+
		f'\n\t\tTotal gene list: {args.l2}'+
		f'\n\t\tPossible candidate genes list: {args.pl2}'+
		f'\n\t\tGO universe: {args.go2}'+
		f'\n\t\tNumber of genes in the initial candidate list: {args.n2}\n'+

		f'\n\tSpecies 1 input files were:'+
		f'\n\t\tTotal gene list: {args.l2}'+
		f'\n\t\tPossible candidate genes list: {args.pl2}'+
		f'\n\t\tGO universe: {args.go2}'+
		f'\n\t\tNumber of genes in the initial candidate list: {args.n2}\n\n')


	# How many permutations had == or > overlaps?
	greater_or_equal_1_2, pval_1_2 = permutation_comparison(args.x_1_2, two_way_1_2_all)

	master_out_file.write(f'In the original contrast between species 1 and 2 there were {args.x_1_2} overlaps.\n'+
		f'In this permutation test {greater_or_equal_1_2} out of {len(two_way_1_2_all)} had a greater or equal number of GO overlaps than the original scan.\n'+
		f'This means that there is a ~ {pval_1_2} chance of getting the original result by chance alone.')


	# How many permutations had == or > overlaps?
	greater_or_equal_1_3, pval_1_3 = permutation_comparison(args.x_1_3, two_way_1_3_all)

	master_out_file.write(f'\n\nIn the original contrast between species 1 and 3 there were {args.x_1_3} overlaps.\n'+
		f'In this permutation test {greater_or_equal_1_3} out of {len(two_way_1_3_all)} had a greater or equal number of GO overlaps than the original scan.\n'+
		f'This means that there is a ~ {pval_1_3} chance of getting the original result by chance alone.')


	# How many permutations had == or > overlaps?
	greater_or_equal_2_3, pval_2_3 = permutation_comparison(args.x_2_3, two_way_2_3_all)

	master_out_file.write(f'\n\nIn the original contrast between species 2 and 3 there were {args.x_2_3} overlaps.\n'+
		f'In this permutation test {greater_or_equal_2_3} out of {len(two_way_2_3_all)} had a greater or equal number of GO overlaps than the original scan.\n'+
		f'This means that there is a ~ {pval_2_3} chance of getting the original result by chance alone.')


	# How many permutations had == or > overlaps?
	greater_or_equal_all, pval_all = permutation_comparison(args.x_all, three_way_all)

	master_out_file.write(f'\n\nIn the original contrast between all three species there were {args.x_all} overlaps.\n'+
		f'In this permutation test {greater_or_equal_all} out of {len(three_way_all)} had a greater or equal number of GO overlaps than the original scan.\n'+
		f'This means that there is a ~ {pval_all} chance of getting the original result by chance alone.\n')



