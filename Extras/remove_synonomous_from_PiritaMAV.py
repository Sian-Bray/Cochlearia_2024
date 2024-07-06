import sys

input_file=open(sys.argv[1], 'r')
output_file=open(sys.argv[1]+'non-synonomous', 'w+')

for count1, line1 in enumerate(input_file):
	split_line=line1.split('\t')
	residues=split_line[2]
	residue1=residues[2:5]
	residue2=residues[-4:-1]
	if residue1 != residue2:
		output_file.write(line1)
