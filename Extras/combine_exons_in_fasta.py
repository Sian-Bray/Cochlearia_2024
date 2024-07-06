# Concatonate fastas composed of exons
# Written by Sian Bray on 24th May 2023

# Example commands:
	# python3 /Users/sian_bray/Dropbox/Scripts/combine_exons_in_fasta.py -f /Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/07_AlphaFold/g7445_test/g7445_diploids_biallelic.fasta -o +
	# python3 /Users/sian_bray/Dropbox/Scripts/combine_exons_in_fasta.py

import argparse
parser = argparse.ArgumentParser(description="Align by ID.")
parser.add_argument('-f', type=str, metavar='fasta_file', required=True, help='Fasta file that contains several operons.')
parser.add_argument('-o', type=str, metavar='vcf_file', required=True, help='The gene orientation, + or - , + for forward, - for reverse.')
args = parser.parse_args()

in_fasta = open(args.f, 'r')
out_fasta = open(args.f+'.combined.fasta', 'w')

if args.o == '+':
	for line in in_fasta:
		if '>' not in line:
			out_fasta.write(line)

if args.o == '-':
	print('Not written yet - waiting on a manual test!')
