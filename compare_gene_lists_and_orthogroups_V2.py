# Convergent Evolution Checker #
# Output will be produced in your current directory
# Modified to work on 3 by Sian Bray (who wrote it) on Sunday 30th July2023

# Command example
# cd /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/06_Orthofinder
# python3 /Users/sian_bray/Dropbox/Scripts/compare_gene_lists_and_orthogroups_V2.py -o /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/06_Orthofinder/OrthoFinder/Results_May19/Orthogroups/Orthogroups.tsv -p /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top_1perc_minus_inversion.txt -t /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top1perc_Fst_arenosa_lyrataIDs.txt -a /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/01_Data/top1perc_Fst_amara_amaraIDs.txt

# Gene names must all be exactly as in the orthogroup file and must be unique to each species
import argparse
parser = argparse.ArgumentParser(description="Convergent Evolution Checker")
parser.add_argument('-o', type=str, metavar='orthogroups_file', required=True, help='Location of the orthogroup file.')
parser.add_argument('-t', type=str, metavar='input_thalianaa', required=True, help='List of input file locations.')
##parser.add_argument('-l', type=str, metavar='input_lyrata', required=True, help='List of input file locations.')
parser.add_argument('-p', type=str, metavar='input_pyrenica', required=True, help='List of input file locations.')
parser.add_argument('-a', type=str, metavar='input_amara', required=True, help='List of input file locations.')
args = parser.parse_args()

# Break each file into list of genes
def file_process(file_in):
	in_file=open(file_in, 'r')
	variable_out=[]
	for count0, line0 in enumerate(in_file):
		line0=line0.replace('\n', '')
		variable_out.append(line0)
	variable_out=set(variable_out)
	in_file.close()
	return variable_out

# Import gene files as sets
thaliana=file_process(args.t)
##lyrata=file_process(args.l)
amara=file_process(args.a)
pyrenica=file_process(args.p)

# Break orthogroups into list of lists
ortho_file=open(args.o, 'r')
orthogroups=[]
for count1, line1 in enumerate(ortho_file):
	line1=line1.replace(',', '\t')
	line1=line1.replace(' ', '')
	line1=line1.replace('\n', '')
	line1=line1.split('\t')
	orthogroups.append(line1)
ortho_file.close()

# Make output files
##t_by_l=open('thaliana_Vs_lyrata.txt', 'w+')
t_by_a=open('thaliana_Vs_amara.txt', 'w+')
t_by_p=open('thaliana_Vs_pyrenica.txt', 'w+')
##l_by_a=open('lyrata_Vs_amara.txt', 'w+')
##l_by_p=open('lyrata_Vs_pyrenica.txt', 'w+')
a_by_p=open('amara_Vs_pyrenica.txt', 'w+')
##in_all=open('all_Vs_all.txt', 'w+')
in_three=open('three_out_of_three.txt', 'w+')

def set_to_str(match):
	match=str(match)
	match=match.replace("{'", "\t")
	match=match.replace("'}", "\t")
	match=match.replace("', '", "\t")
	match=match.replace("set()", "")
	return match

def compare_2(og1, og, g1, g2, out_file):
	if bool(og & g1) and bool(og & g2) == True:
		match1=set_to_str(og & g1)
		match2=set_to_str(og & g2)
		out_file.write(og1+match1+match2+'\n')

# Make comparisons
for count2, item2 in enumerate(orthogroups):
	group=set(item2)
	##compare_2(item2[0], group, thaliana, lyrata, t_by_l)
	compare_2(item2[0], group, thaliana, amara, t_by_a)
	compare_2(item2[0], group, thaliana, pyrenica, t_by_p)
	##compare_2(item2[0], group, lyrata, amara, l_by_a)
	##compare_2(item2[0], group, lyrata, pyrenica, l_by_p)
	compare_2(item2[0], group, amara, pyrenica, a_by_p)
	### Four way comparison
	##if bool(group & thaliana) and bool(group & lyrata) and bool(group & amara) and bool(group & pyrenica) == True:
	##	match3=set_to_str(group & thaliana)
	##	match4=set_to_str(group & lyrata)
	##	match5=set_to_str(group & amara)
	##	match6=set_to_str(group & pyrenica)
	##	in_all.write(item2[0]+match3+match4+match5+match6+'\n')
	# Three out of three
	threes=0
	if bool(group & thaliana) == True:
		threes+=1
	##if bool(group & lyrata) == True:
	##	threes+=1
	if bool(group & amara) == True:
		threes+=1
	if bool(group & pyrenica) == True:
		threes+=1
	if threes == 3:
		match3=set_to_str(group & thaliana)
		##match4=set_to_str(group & lyrata)
		match5=set_to_str(group & amara)
		match6=set_to_str(group & pyrenica)
		in_three.write(item2[0]+match3+match5+match6+'\n')

# Close output files
##t_by_l.close()
t_by_a.close()
t_by_p.close()
##l_by_a.close()
##l_by_p.close()
a_by_p.close()
##in_all.close()
in_three.close()
