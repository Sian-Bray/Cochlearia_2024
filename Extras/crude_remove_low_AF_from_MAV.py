# removes SNPs with AF < 0.25 from british_diploid_tetraploid_contrast_changes-only_noTer.txt
# this is part of the MAV pipeline

in_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/03_Pirita_Mav/british_diploid_tetraploid_contrast_changes-only_noTer.txt', 'r')
out_file = open('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/03_Pirita_Mav/british_diploid_tetraploid_contrast_changes-only_noTer_AF_0.25.txt', 'w')

for line in in_file:
	if 'CHROM' in line:
		out_file.write(line)
	if 'CHROM' not in line:
		calc = line.split(' ')
		AF1 = int(calc[5])/int(calc[6]) # col 5 and 6, AC AN
		AF2 = int(calc[7])/int(calc[8]) # col 7 and 8, AC AN
		AFD = abs(AF1 - AF2)
		if AFD >= 0.25:
			out_file.write(line)

in_file.close()
out_file.close()
