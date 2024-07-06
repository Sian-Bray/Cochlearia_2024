# This script compares the final Cochlearia genome scan with the mini
# inland tetraploid reviewers comment list to sumarise what genes are
# in the new list Vs the old one.

# Sian Bray 04May2024

# Define functions
def file_to_dict(in_file):
	new_dict = {}
	for line in in_file:
		gene = line.split('\t')
		gene = gene[0]
		if gene in new_dict:
			new_dict[gene].append(line)
		else:
			new_dict[gene] = []
			new_dict[gene].append(line)
	return new_dict

def list_write(in_list):
	stringy = ''
	for line in in_list:
		stringy = stringy + line
	return stringy

def highest_FstH(list_of_lists):
	fstH = -2
	to_write = ''
	headder_skip = 0
	for in_list in list_of_lists:
		try:
			current_FstH = in_list.split('\t')
			current_FstH = float(current_FstH[13])
			if current_FstH > fstH:
				to_write = in_list
				fstH = current_FstH
		except(ValueError):
			return in_list
	return to_write


# Read input files
full_1perc = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_1perc.txt', 'r')
no_salt_1perc = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/genes_UK_EU_no-salt_15_BPM_top1perc.txt', 'r')
full_all = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/genes_UK_dips_Vs_tets_15_BPM_all.txt', 'r')
no_salt_all = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/genes_UK_EU_no-salt_15_BPM.txt', 'r')

# Convert to dictionaries and lists
full_1perc_dict = file_to_dict(full_1perc)
no_salt_1perc_dict = file_to_dict(no_salt_1perc)
full_all_dict = file_to_dict(full_all)
no_salt_all_dict = file_to_dict(no_salt_all)

# Close files
full_1perc.close()
no_salt_1perc.close()
full_all.close()
no_salt_all.close()

# Ask questions

# How many and what genes are in the old and new 1%?
out_old = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/in_old_and_new_OldValues.txt', 'w')
out_new = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/in_old_and_new_NewValues.txt', 'w')
for key in no_salt_1perc_dict:
	if key in full_1perc_dict:
		out_old.write(list_write(full_1perc_dict[key]))
		out_new.write(list_write(no_salt_1perc_dict[key]))

out_old.close()
out_new.close()

# What genes from the old 1% are not in the new 1%
old_missing = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/in_old_missing_from_new.txt', 'w')

for key in full_1perc_dict:
	if key not in no_salt_1perc_dict:
		old_missing.write(list_write(full_1perc_dict[key]))

old_missing.close()

# What genes are in the new 1% but not the old 1%
new_missing = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/in_new_missing_from_old.txt', 'w')

for key in no_salt_1perc_dict:
	if key not in full_1perc_dict:
		new_missing.write(list_write(no_salt_1perc_dict[key]))

new_missing.close()

# Full metrics in the original scan for all of the genes in the new 1%
old_metrics = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/new_scan_1perc_old_metrics.txt', 'w')

for key in no_salt_1perc_dict:
	old_metrics.write(list_write(full_all_dict[key]))

old_metrics.close()

# Full metrics in the new scan for all of the genes in the old 1%
new_metrics = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/old_scan_1perc_new_metrics.txt', 'w')
# What genes are missing from the old scan?
missing_genes = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/missing_genes.txt', 'w')

for key in full_1perc_dict:
	try:
		new_metrics.write(list_write(no_salt_all_dict[key]))
	except(KeyError):
		missing_genes.write(key + '\n')

new_metrics.close()
missing_genes.close()

# What is the top FstH window for all the top 1% genes in the old scan in the new scan?
top_windows_old = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/top_windows_for_old_1perc_in_new_scan.txt', 'w')

for key in full_1perc_dict:
	try:
		current_list = no_salt_all_dict[key]
		to_write = highest_FstH(current_list)
		top_windows_old.write(to_write)
	except(KeyError):
		pass

top_windows_old.close()

# What is the top FstH window for all the top 1% genes in the new scan in the old scan?

top_windows_new = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/Cochlearia_Final_2023/Reviewers_Comments/No_Salt_Tets/top_windows_for_new_1perc_in_old_scan.txt', 'w')

for key in no_salt_1perc_dict:
	try:
		current_list = full_all_dict[key]
		to_write = highest_FstH(current_list)
		top_windows_new.write(to_write)
	except(KeyError):
		pass

top_windows_new.close()

















































