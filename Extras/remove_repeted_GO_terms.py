# Removes duplicated OG's from a id2gos file
# Written by Sian Bray on 28th June 2023

# Example run command:
	# python3 /Users/sian_bray/Dropbox/Scripts/remove_repeted_GO_terms.py

in_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_GO_universe_restrictive_id2gos_no_obsolete.tsv', 'r')
out_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_GO_universe_restrictive_id2gos_no_obsolete_nonredundant.tsv', 'w')

def GO_write(GO, written_GOs, out_line, seperator):
	if GO not in written_GOs:
		out_line = out_line + (f'{GO}{seperator}')
		written_GOs.append(GO)
	return GO, written_GOs, out_line


for line in in_file:
	out_line = ''
	gene = line.split('\t')
	GOs = gene[1]
	gene = gene[0]
	GOs = GOs.replace('\n', '')
	GOs	= GOs.split(';')
	out_line = out_line + (f'{gene}\t')
	written_GOs = []
	Go_no = len(GOs) - 1
	for count, GO in enumerate(GOs): # count is 0-based, len is 1-based
		if count == Go_no:
			GO, written_GOs, out_line = GO_write(GO, written_GOs, out_line, '\n')
		if count != Go_no:
			GO, written_GOs, out_line = GO_write(GO, written_GOs, out_line, ';')
	if out_line[-1:] == '\n':
		out_file.write(out_line)
	if out_line[-1:] != '\n':	
		out_file.write(f'{out_line[:-1]}\n')

