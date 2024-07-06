#remove lines without a protein change from PiritaMAV input

input_file=open('british_diploid_tetraploid_contrast.txt', 'r')
output_file=open('british_diploid_tetraploid_contrast_changes-only_noTer.txt', 'w+')
for count, line in enumerate(input_file):
    if '#CHROM' in line:
        output_file.write(line)
    if 'p.' in line:
        if 'Ter' not in line:
            if '???' not in line:
                output_file.write(line)