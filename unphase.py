# removed phasing for scantools recode

# python3 /Users/sian_bray/Dropbox/Scripts/unphase.py -i 
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', type=str, metavar='input_table_file', required=True, help='path the the input file')
args = parser.parse_args()

in_file = open(args.i, 'r')
out_name = args.i+'.unphased'
out_file = open(out_name, 'w')

for line in in_file:
	line = line.replace('|', '/')
	out_file.write(line)

in_file.close()
out_file.close()


