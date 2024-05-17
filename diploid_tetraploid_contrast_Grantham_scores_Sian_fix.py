#### requires python-2.7.10 (and no longer biopython-1.69 FFS)

#!/usr/bin/env python
import re
import string
import sys
# import gzip
# import Bio
 
# from Bio.SubsMat import MatrixInfo   
# blosum = MatrixInfo.blosum100
# grant = MatrixInfo.grant
# blosum['D','N']
grant = {('G', 'G'): 215, ('K', 'G'): 88, ('S', 'E'): 135, ('Y', 'E'): 93, ('W', 'R'): 114, ('V', 'M'): 194, ('N', 'R'): 129, ('W', 'Q'): 85, ('L', 'Q'): 102, ('V', 'N'): 82, ('F', 'K'): 113, ('G', 'E'): 117, ('S', 'L'): 70, ('P', 'R'): 112, ('E', 'D'): 170, ('Y', 'G'): 68, ('W', 'P'): 68, ('F', 'D'): 38, ('G', 'D'): 121, ('K', 'D'): 114, ('T', 'N'): 150, ('W', 'W'): 215, ('L', 'D'): 43, ('S', 'S'): 215, ('K', 'C'): 13, ('S', 'A'): 116, ('Y', 'I'): 182, ('V', 'I'): 186, ('Q', 'C'): 61, ('T', 'G'): 156, ('T', 'L'): 123, ('F', 'G'): 62, ('V', 'T'): 146, ('S', 'H'): 126, ('I', 'Q'): 106, ('Y', 'K'): 130, ('W', 'T'): 87, ('P', 'D'): 107, ('I', 'C'): 17, ('K', 'R'): 189, ('T', 'E'): 150, ('Q', 'R'): 172, ('K', 'Q'): 162, ('Y', 'M'): 179, ('V', 'E'): 94, ('Y', 'D'): 55, ('V', 'W'): 127, ('T', 'C'): 66, ('T', 'H'): 168, ('F', 'Q'): 99, ('L', 'I'): 210, ('M', 'Q'): 114, ('R', 'A'): 103, ('C', 'D'): 61, ('V', 'F'): 165, ('F', 'C'): 10, ('C', 'R'): 35, ('D', 'D'): 215, ('V', 'P'): 147, ('S', 'D'): 150, ('P', 'C'): 46, ('F', 'R'): 118, ('C', 'C'): 215, ('I', 'G'): 80, ('W', 'K'): 105, ('I', 'N'): 66, ('T', 'A'): 157, ('K', 'L'): 108, ('L', 'G'): 77, ('F', 'A'): 102, ('S', 'K'): 94, ('K', 'K'): 215, ('E', 'N'): 173, ('Y', 'Q'): 116, ('V', 'A'): 151, ('W', 'I'): 154, ('V', 'S'): 91, ('T', 'T'): 215, ('F', 'M'): 187, ('L', 'E'): 77, ('M', 'M'): 215, ('W', 'H'): 100, ('S', 'R'): 105, ('P', 'Q'): 139, ('P', 'N'): 124, ('H', 'A'): 129, ('P', 'G'): 173, ('F', 'N'): 57, ('H', 'N'): 147, ('P', 'K'): 112, ('T', 'M'): 134, ('K', 'H'): 183, ('T', 'R'): 144, ('L', 'C'): 17, ('W', 'N'): 41, ('E', 'Q'): 186, ('S', 'G'): 159, ('Y', 'S'): 71, ('G', 'R'): 90, ('W', 'M'): 148, ('Q', 'A'): 124, ('T', 'K'): 137, ('C', 'N'): 76, ('T', 'P'): 177, ('V', 'L'): 183, ('F', 'I'): 194, ('G', 'Q'): 128, ('L', 'A'): 119, ('M', 'I'): 205, ('W', 'L'): 154, ('S', 'N'): 169, ('I', 'R'): 118, ('H', 'E'): 175, ('Y', 'W'): 178, ('I', 'D'): 47, ('W', 'C'): 0, ('N', 'A'): 104, ('T', 'I'): 126, ('Q', 'N'): 169, ('M', 'K'): 120, ('K', 'E'): 159, ('S', 'C'): 103, ('Y', 'Y'): 215, ('V', 'Y'): 160, ('W', 'A'): 67, ('Y', 'F'): 193, ('M', 'R'): 124, ('V', 'H'): 131, ('F', 'E'): 75, ('M', 'E'): 89, ('H', 'R'): 186, ('P', 'P'): 215, ('P', 'I'): 120, ('Q', 'Q'): 215, ('P', 'F'): 101, ('I', 'A'): 121, ('F', 'F'): 215, ('I', 'H'): 121, ('W', 'G'): 31, ('Y', 'H'): 132, ('M', 'L'): 200, ('M', 'G'): 88, ('S', 'Q'): 147, ('W', 'F'): 175, ('D', 'A'): 89, ('K', 'A'): 109, ('N', 'N'): 215, ('V', 'K'): 118, ('W', 'E'): 63, ('L', 'R'): 113, ('T', 'S'): 157, ('M', 'N'): 73, ('V', 'D'): 63, ('Q', 'D'): 154, ('M', 'A'): 131, ('V', 'V'): 215, ('W', 'D'): 34, ('S', 'F'): 60, ('D', 'N'): 192, ('P', 'M'): 128, ('H', 'D'): 134, ('I', 'E'): 81, ('R', 'R'): 215, ('K', 'N'): 121, ('Y', 'L'): 179, ('T', 'Q'): 173, ('P', 'L'): 117, ('M', 'H'): 128, ('M', 'C'): 19, ('S', 'M'): 80, ('E', 'R'): 161, ('E', 'E'): 215, ('V', 'G'): 106, ('G', 'N'): 135, ('A', 'A'): 215, ('V', 'Q'): 119, ('L', 'N'): 62, ('Y', 'N'): 72, ('V', 'R'): 119, ('P', 'H'): 138, ('H', 'C'): 41, ('P', 'A'): 188, ('F', 'L'): 193, ('H', 'H'): 215, ('C', 'A'): 20, ('I', 'I'): 215, ('T', 'F'): 112, ('L', 'L'): 215, ('Y', 'P'): 105, ('D', 'R'): 119, ('M', 'D'): 55, ('G', 'C'): 56, ('S', 'I'): 73, ('Y', 'A'): 103, ('E', 'A'): 108, ('K', 'I'): 113, ('V', 'C'): 23, ('T', 'D'): 130, ('Y', 'R'): 138, ('G', 'A'): 155, ('S', 'P'): 141, ('H', 'Q'): 191, ('Y', 'C'): 21, ('E', 'C'): 45, ('H', 'G'): 117, ('P', 'E'): 122, ('F', 'H'): 115, ('W', 'S'): 38, ('L', 'H'): 116, ('Y', 'T'): 123}
 
Amino_acid_code={'Gly':'G','Ala':'A','Val':'V','Leu':'L','Ile':'I','Phe':'F','Trp':'W','Tyr':'Y','Asp':'D','Asn':'N','Glu':'E','Lys':'K','Gln':'Q','Met':'M','Ser':'S','Thr':'T','Cys':'C','Pro':'P','His':'H','Arg':'R'}
Amino_acid_code_grant={'Gly':'Y','Ala':'L','Val':'V','Leu':'U','Ile':'W','Phe':'F','Trp':'T','Tyr':'O','Asp':'A','Asn':'N','Glu':'G','Lys':'I','Gln':'Q','Met':'M','Ser':'S','Thr':'E','Cys':'C','Pro':'P','His':'H','Arg':'R'}
 
 
exponent=float((sys.argv[1]))
scaffold=(sys.argv[2]) ## and is the  infile to be used is from the prepared file from step 1, called diploid_tetraploid_contrast.txt 
 
with open ('/Users/sian_bray/Dropbox/Bray/000_Research/Cochlearia_Massive/03_Pirita_Mav/british_diploid_tetraploid_contrast_changes-only_noTer.txt'.format(scaffold),'r') as  Pop1:  

    print "Chrom", "loci", "AA", "AF/2x", "AF/4X",  "DAP_DAF1", "DAP_DAF2", "grant_matrix_score", "grant1", "grant2", "difference in scores" # Print the header line
    for aline in Pop1:
        fields=string.split(string.strip(aline))
        if fields[8][0]!="p":   ##check that you are having real amino acid changes
            continue
        else:
            derived=float(float(fields[4])+float(fields[6]))
            if float(fields[5])>8 and float(fields[7])>48 and derived!=0:  ## Filter for half the depth, and for the fact that there are differences
                dap=float((float(fields[4])/derived)**exponent+float(float(fields[6])/derived)**exponent) #Calculates the allele purity, here exponent present
                daf_1=float(float(fields[4])/float(fields[5])) # derived allele frequency for diploids
                daf_2=float(float(fields[6])/float(fields[7])) # deribed allele frequence for tetraploids
                amino_acid_change=re.search(r'.{2}([\D]{3})\d*([\D]{3})',fields[8]) ## Fine the change in amino acids
                if amino_acid_change:
                    #pair_grant = (Amino_acid_code_grant[str(amino_acid_change.group(1))],Amino_acid_code_grant[str(amino_acid_change.group(2))]) ##turn these to one letter coder
                    pair_grant = (Amino_acid_code[str(amino_acid_change.group(1))],Amino_acid_code[str(amino_acid_change.group(2))]) #Sian edit. In my (updated?) version of Bio.SubsMat the Grantham Matrix uses standard a.a codes
                    if pair_grant not in grant:
                        grant_matrix_score= 215-grant[(tuple(reversed(pair_grant)))] ## read the code from the grantham matrix info
                    else:
                        grant_matrix_score = 215-grant[pair_grant]   ## and it depends on the order, as it upper triangular matrix
                else: continue                 
                dap_daf1=float(dap*daf_1)
                dap_daf2=float(dap*daf_2)
                print fields[0],"\t",fields[1],"\t", fields[8],"\t", fields[4],"\t", fields[6],"\t", dap_daf1,"\t", dap_daf2, "\t",grant_matrix_score,"\t",float(dap_daf1*grant_matrix_score),"\t", float(dap_daf2*grant_matrix_score), "\t",abs(float(dap_daf1*grant_matrix_score)-float(dap_daf2*grant_matrix_score))
 
##Calculates scores and prints the results line by line
 
