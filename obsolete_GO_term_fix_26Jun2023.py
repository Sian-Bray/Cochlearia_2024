# replaces obsolete GO terms with updated ones
# Written on 26th June 2023 by Sian Bray (sian.bray@nottingham.ac.uk)

# example usage command:
	# python3 Users/sian_bray/Dropbox/Scripts/obsolete_GO_term_fix_26Jun2023.py

##############################################

# input file - my home made GO universe
i_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_GO_universe_restrictive_id2gos.tsv', 'r')

# output file
o_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/Cochlearia_Thaliana_GO_universe_restrictive_id2gos_no_obsolete.tsv', 'w')

# Substitution dictionary; old as key, new as value
# If a term needs to be removed it's value is an 'X' (capital)
s_dict = {'GO:0000469':'GO:0006364', 'GO:0000478':'GO:0006364', 'GO:0000737':'GO:0006308', 'GO:0000738':'GO:0006308', 'GO:0004652':'GO:1990817', 'GO:0006073':'GO:0044042', 'GO:0006165':'GO:0009142', 'GO:0006379':'X', 'GO:0006876':'GO:0098849', 'GO:0007568':'X', 'GO:0010216':'GO:0044027', 'GO:0010252':'GO:0140964', 'GO:0010529':'GO:0010526', 'GO:0010767':['GO:0006357', 'GO:0034644'], 'GO:0010941':'X', 'GO:0010942':'X', 'GO:0016508':'GO:0018812', 'GO:0018024':['GO:0140938','GO:0140939','GO:0140940'], 'GO:0030104':['GO:0003091','GO:0009992','GO:0050891'], 'GO:0031058':['GO:0006338','GO:0040029'], 'GO:0033692':'GO:0000271', 'GO:0034414':'GO:0042780', 'GO:0034637':'GO:0016051', 'GO:0034645':'GO:0009059', 'GO:0035174':['GO:0140995','GO:0140996','GO:0140997','GO:0140998'], 'GO:0035551':'X', 'GO:0042779':'GO:0042780', 'GO:0043631':'X', 'GO:0044260':'GO:0043170', 'GO:0044262':'GO:0005975', 'GO:0044264':'GO:0005976', 'GO:0044265':'GO:0009057', 'GO:0046838':'X', 'GO:0046855':'GO:0043647', 'GO:0046939':'GO:0016310', 'GO:0048478':'GO:0031297', 'GO:0050072':['GO:0140932','GO:0140933'], 'GO:0050213':'GO:0047751', 'GO:0052746':'X', 'GO:0055072':['GO:0006879','GO:0060586'], 'GO:0055073':'GO:0071585', 'GO:0060548':'X', 'GO:0070191':'GO:0033745', 'GO:0070816':'GO:0008353', 'GO:0070940':'GO:0008420', 'GO:0072684':'GO:1990180', 'GO:0072755':'X', 'GO:0080023':'GO:0018812', 'GO:0090502':'X', 'GO:0097659':'GO:0006366', 'GO:0098789':'GO:0031124', 'GO:0102158':'GO:0018812', 'GO:0102339':'GO:0141040', 'GO:0102340':'GO:0141040', 'GO:0102341':'GO:0141040', 'GO:0102342':'GO:0141040', 'GO:0102343':'GO:0018812', 'GO:0102344':'GO:0018812', 'GO:0102345':'GO:0018812', 'GO:0102419':'GO:0090447', 'GO:0102756':'GO:0009922', 'GO:1900363':'X', 'GO:1901407':'GO:0006357', 'GO:1901485':'GO:0032436', 'GO:1901972':'X', 'GO:1903506':'GO:0006357', 'GO:1903507':'GO:0000122', 'GO:1903508':'GO:0045944', 'GO:2000653':'GO:0071514', 'GO:2001020':'GO:0006974', 'GO:2001022':'GO:0006974', 'GO:2001253':'GO:0140673'}

# Below are genes that have no 'replace by' on EMBL-EBI QuickGO
# GO:0000478 is "endonucleolytic cleavage involved in rRNA processing" replaced with GO:0006364 "rRNA processing"
# GO:0006165 is "nucleoside diphosphate phosphorylation" replaced with GO:0009142 "The chemical reactions and pathways resulting in the formation of a nucleoside triphosphate, a compound consisting of a nucleobase linked to a deoxyribose or ribose sugar esterified with triphosphate on the sugar."
# GO:0006379 is a biological process term that is a molecular function - remove
# GO:0007568 the reason for obsoletion is that this represents a phenotype - remove
# GO:0010941 and GO:0010942 the reason for obsoletion is that this term represent an assay and not a GO process -remove
# GO:0035551 the reason for obsoletion is that this term represents a molecular function - could not find a good alternative, remove
# GO:0043631 the reason for obsoletion is that this term represents a molecular function - remove
# GO:0046838 this term was obsoleted because it represents a molecular function - remove
# GO:0046939 is "nucleotide phosphorylation" replace with GO:0016310 "phosphorylation"
# GO:0052746 the reason for obsoletion is this is a single-step process - remove
# GO:0060548 the reason for obsoletion is that this term represent an assay and not a GO process - remove
# GO:0072755 this term is out of scope for GO - remove
# GO:0090502 this term was obsoleted because it represents a molecular function - remove
# GO:0098789 "pre-mRNA cleavage required for polyadenylation" replaced with GO:0031124 "mRNA 3'-end processing"
# GO:1900363 the reason for obsoletion is that this term represents a molecular function - remove
# GO:1901972 this term was obsoleted because it represents a molecular function - remove

##############################################

# Check to make sure I am not repeating GO terms with my replacements
# Note that some GO terms have more than one replacement (key extracts list)
# Note that some GO terms need to simply be removed

for line in i_file:
	
	new_line = ''
	written = False
	
	for key in s_dict:

		if (key in line) and (written == False):
			
			line2 = line.replace('\n', '')
			line2 = line2.split('\t')
			gene = line2[0]
			new_line = new_line + (f'{gene}\t')
			terms = line2[1]
			terms = terms.split(';')
			
			for term in terms:
				
				if term in s_dict:

					if (isinstance(s_dict[term], str)) and (s_dict[term] != 'X'):
					
						new_line = new_line + f'{s_dict[term]};'

					if isinstance(s_dict[term], list):

						for OG in s_dict[term]:

							new_line = new_line + f'{OG};'

					if s_dict[term] == 'X':

						pass

				if term not in s_dict:

					new_line = new_line + f'{term};'

			written = True
			new_line = new_line[:-1]+'\n'
			o_file.write(new_line)

	if new_line == '':

		o_file.write(line)


















