# Kinetochore and ionomic adaptation to whole genome duplication in autotetraploid Cochlearia reveals evolutionary convergence in three autopolyploids

## Introduction

This page contains the code for the paper:

Sian M. Bray, Tuomas Hämälä, Min Zhou, Silvia Busoms, Sina Fischer, Stuart D. Desjardins, Terezie Mandáková, Chris Moore, Thomas C. Mathers, Laura Cowan, Patrick Monnahan, Jordan Koch, Eva M. Wolf1, Martin A. Lysak, Filip Kolar, James D. Higgins, Marcus A. Koch, and Levi Yant. _Kinetochore and ionomic adaptation to whole genome duplication in autotetraploid Cochlearia reveals evolutionary convergence in three autopolyploids._ 2024. **Cell Reports.**

## Environments and Inputs

All code was run using Unix. This was done either locally on an Intel chip MacBook Pro or on a Unix HPC.

### Conda Environments

Conda environments were used for much of the below code. Conda is available at: https://conda.io/projects/conda/en/latest/user-guide/getting-started.html. Conda environments for all tools used below can be found in these Conda yml files:

	ngs_pipe_env_gatk4_23May2024.yml
	pirita_mav_12Apr2023.yml
	snpsift_env_23May2024.yml
	orthofinder_env_05Jul2024.yml

These yml files can be used to replicate the environments and contain the version numbers of the code languages and programs used.

### Input files

Because I often find example files far more useful than file descriptions I have included truncated examples of all the files used below. These have the same names as in the code below with the prefix 'subsample_'. So where I have used the file `C_excelsa_V5.fasta below`, `subsample_C_excelsa_V5.fasta` is the truncated example input file made from `C_excelsa_V5.fasta`. Most of these truncations will contain the first and last few lines of the file, though others (i.e. particularly long fasta files containg reference genomes) will also contain some trunctaed lines.

## Window-based scan for selective sweep signatures

Activate the `ngs_pipe_env_gatk4` Conda environment.

	conda activate ngs_pipe_env_gatk4

**1\)** The VCF containing all individuals must be split by phenotype, in this case ploidy. Both VCFs must be generated from the same 'master' VCF because the SNPs in each VCF have to be identical for step 6 to work. If you need to filter your data (e.g. for allele frequency or allele count) do it before you split the VCF by phenotype.

To generate the tetraploid VCF:

	gatk SelectVariants -R C_excelsa_V5.fasta -V reheadered.F4_133.ann.vcf.gz -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn LOS_1 -sn LOS_2 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 -O UK_tets.vcf

To generate the diploid VCF:

	gatk SelectVariants -R C_excelsa_V5.fasta -V reheadered.F4_133.ann.vcf.gz -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O UK_dips.vcf

**2\)** The VCF is converted into a table.

For the tetraploids:

	gatk VariantsToTable -V UK_tets.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O UK_tets.table

For the diploids:

	gatk VariantsToTable -V UK_dips.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O UK_dips.table

**3\)** The tables are run through a custom python script `unphase.py` which simply converts pipes ('|') that indicate phased data into forward slashed ('/') that indicate unphased data. Phased data is not used for this part of the analysis and the pipes mess up step 4. 

For the tetraploids:

	python3 unphase.py -i UK_tets.table

For the diploids:

	python3 unphase.py -i UK_dips.table

**4\)** The data is recoded into a format that can be accepted by the ScanTools pipeline. The ScanTools pipeline was originally written by Patrick Monahan (https://github.com/pmonnahan/ScanTools) and has since been developed, modified and troubleshot by Magdalena Bohutinska (https://github.com/mbohutinska/ScanTools_ProtEvol), Sonia Celestini and myself (Sian Bray). For clarity I include code for the version I used in this analysis here.

The custom format is a tsv file with no headder and the columns: population ID, ploidy, scaffold, position, allele number, depth and then one genotype columns per individual where 0 is reference-like homozygous, 1 is a single alternate (non-reference) allele, 2 is two alternate alleles, etc (-9 is missing data).

For the tetraploids:

	python3 recode012.py -i UK_tets.table.unphased -o ~/Cochlearia/ScanTools/ -pop UK_TETS

For the diploids:

	python3 recode012.py -i UK_dips.table.unphased -o ~/Cochlearia/ScanTools/ -pop UK_DIPS

**5\)** The two tables are now merged. I use the custom code Sian_sort_for_ScanTools.py here because of a bug I discovered in the original ScanTools code that didn't properly sort contigs/scaffolds if their names contained letters and more than 9 contigs/scaffolds.

	python3 Sian_sort_for_ScanTools.py 'UK_tets.table.unphased.recode.txt UK_dips.table.unphased.recode.txt ' ~/Cochlearia/ScanTools/Temp/ UK_dips_Vs_tets

**6\)** Calculate the Between Population Metrics (BPM) using modified ScanTools code.

	python3 sian_bpm.py -i UK_dips_Vs_tets.concat.txt -o ~/Cochlearia/ScanTools/UK_bpm/ -prefix UK_dips_Vs_tets_15 -ws 1000 -ms 15 -np 2

**7\)** Match windows to genes with custom code.

	python3 bpm_windows_to_genes.py -i UK_dips_Vs_tets_15_BPM.txt -o genes_UK_dips_Vs_tets_15_BPM.txt -g C_excelsa_V5_braker2_wRseq.gff3

Finally deactivate the Conda environment.

	conda deactivate

### Controlling for Correlated Environmental Factors

In response to insightful comments from the reviewers the window based selection was repeated using only tetraploids from freshwater environments. The commands for this are below:

**0\)** Single 'master' VCF created.

	gatk SelectVariants -V reheadered.F4_133.ann.vcf.gz -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn KVA_002 -sn KVA_003 -sn KVA_009 -sn KVA_010 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -sn TRO_001 -sn TRO_003 -sn TRO_005 -sn TRO_009 -sn WOL_002 -sn WOL_006 -sn WOL_009 -sn WOL_010 -O EU_UK_dips_no-salt-tets.vcf.gz

**1\)** Split by phenotype.

	gatk SelectVariants -V EU_UK_dips_no-salt-tets.vcf.gz -sn KVA_002 -sn KVA_003 -sn KVA_009 -sn KVA_010 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn TRO_001 -sn TRO_003 -sn TRO_005 -sn TRO_009 -O EU_UK_no-salt-tets.vcf.gz

	gatk SelectVariants -V EU_UK_dips_no-salt-tets.vcf.gz -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -sn WOL_002 -sn WOL_006 -sn WOL_009 -sn WOL_010 -O /EU_UK_dips.vcf.gz

**2\)** Convert into tables.

	gatk VariantsToTable -V EU_UK_no-salt-tets.vcf.gz -F CHROM -F POS -F REF -F AN -F DP -GF GT -O EU_UK_no-salt-tets.table

	gatk VariantsToTable -V /EU_UK_dips.vcf.gz -F CHROM -F POS -F REF -F AN -F DP -GF GT -O EU_UK_dips.table

**3\)** Remove phasing.

	python3 unphase.py -i EU_UK_no-salt-tets.table

	python3 unphase.py -i EU_UK_dips.table

**4\)** Convert to custom format.

	python3 recode012.py -i EU_UK_no-salt-tets.table.unphased -o ~/Cochlearia/Reviews/ -pop EU_UK_TETS

	python3 recode012.py -i EU_UK_dips.table.unphased -o ~/Cochlearia/Reviews/ -pop EU_UK_DIPS

**5\)** Merge tables.

	python3 Sian_sort_for_ScanTools.py 'EU_UK_no-salt-tets.table.unphased.recode.txt EU_UK_dips.table.unphased.recode.txt ' ~/Cochlearia/Reviews/Temp/ UK_EU_no-salt

**6\)** Calculate BPM.

	python3 sian_bpm.py -i UK_EU_no-salt.concat.txt -o ~/Cochlearia/Reviews/ -prefix UK_EU_no-salt_15 -ws 1000 -ms 15 -np 2

**7\)** Compare the two scans:

	python3 new_vs_old_Cochlearia_04May2024.py

## MAV Analysis.

A fineMAV like analysis was performed (https://doi.org/10.1186/s13059-017-1380-2). The code for this was originally written by Pirita Paajanen (https://github.com/paajanen/meiosis_protein_evolution) but was modified by myself to remove the need to install a defunct module.

The analysis was performed in the `ngs_pipe_env_gatk4`, snpsift_env and pirita_mav Conda environments. The order of activation and deactivation of these environments is indicated below.

### Calculate the MAV values.

**1\)** Get British diploid-only and British tetraploid-only VCFs filtered for allele frequenct of 0.25 or greater.

Activate conda.

	conda activate ngs_pipe_env_gatk4

Generate 'master' VCF.

	gatk SelectVariants -R C_excelsa_V5.fasta -V reheadered.F4_133.ANNOTATED.vcf --max-nocall-fraction 0.2 -select "DP > 7" -select "AF > 0.25" -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn LOS_1 -sn LOS_2 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O UK_mav_dips_tets.vcf

Generate tetraploid VCF.

	gatk SelectVariants -R C_excelsa_V5.fasta -V UK_mav_dips_tets.vcf -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn LOS_1 -sn LOS_2 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 -O UK_mav_tets.vcf

Generate diploid VCF.

	gatk SelectVariants -R C_excelsa_V5.fasta -V UK_mav_dips_tets.vcf -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O UK_mav_dips.vcf

Deactivate the ngs_pipe_env_gatk4 Conda environment.

	conda deactivate

**2\)** Convert the VCFs to a table.

Activate the snpsift_env Conda environment.

	conda activate snpsift_env

Generate the tetraploid table.

	SnpSift extractFields UK_mav_tets.vcf CHROM POS REF ALT AC AN "ANN[*].HGVS_P" > UK_mav_tets.table

Generate the diploid table.

	SnpSift extractFields UK_mav_dips.vcf CHROM POS REF ALT AC AN "ANN[*].HGVS_P" > UK_mav_dips.table

Deactivate the snpsift_env Conda environment.

	conda deactivate snpsift_env

Paste the tables together.

	paste UK_mav_dips.table UK_mav_tets.table > british_diploid_tetraploid_contrast.table

**3\)** Remove superfluous columns.

	awk '{print " "$1" "$2" "$3" "$4" "$5" "$6" "$12" "$13" "$14}' british_diploid_tetraploid_contrast.table > british_diploid_tetraploid_contrast.txt

**4\)** Remove any lines that don't have a protein change. 

	python3 remove_no-change_and_Ter.py

**5\)** Calculate the MAV scores.

Activate the pirita_mav Conda environment.

	conda activate pirita_mav

Calculate the MAV scores.

	python2 diploid_tetraploid_contrast_Grantham_scores_Sian_fix.py 3.5 british_diploid_tetraploid_contrast_changes-only_noTer.txt  > british_dip_tet_3.5.scores

Deactivate Conda environment.

	conda deactivate

**6\)** Get the top 1% outliers.

Sort by score.

	sort -k11 -n -r british_dip_tet_3.5.scores > sorted_british_dip_tet_3.5.scores

Manually create british_1percent_outliers_3.5_scores.txt in excel (1% of 266689 is 2667).

Manually create `british_top_1%\_FstH.txt` in excel (ScanTools output with headder). This file is then copied and renamed `UK_dips_Vs_tets_15_BPM_Apr2023_FstH_1percent.tsv`.

Make bed files.

	awk '{ print $1 "\t" $2 "\t" $2+1 }' british_1percent_outliers_3.5_scores.txt > piritaMAV_1percent_outliers.bed

	awk '{ print $2 "\t" $3 "\t" $4+1 }' british_top_1%_FstH.txt > FstH_1percent_outliers.bed

Manulally remove headder from `FstH_1percent_outliers.bed`

Get a list of genes that have top 1% MAVs in them.

	bedtools intersect -wb -a piritaMAV_1percent_outliers.bed -b C_excelsa_V5_braker2_wRseq.gff3 > british_piritaMAV_1percent_genes.gff3

get a list of MAV outlier genes in 1% Fst outlier windows.

	bedtools intersect -wa -a british_piritaMAV_1percent_genes.gff3 -b FstH_1percent_outliers.bed > MAV_FstH_1percent_overlap_genes.gff

Make a plain text list, one gene per line.

	awk '{ print $NF}' MAV_FstH_1percent_overlap_genes.gff > MAV_FstH_1percent_overlap_genes.txt

Manually clean up this file, then remove duplicates.

	sort -u MAV_FstH_1percent_overlap_genes.txt > MAV_FstH_1percent_overlap_genes_unique.txt

### Generate the MAV plots.

**1\)** Generate British diploid-only and British tetraploid-only VCFs containing identical SNPs and filtered for an allele frequency of 0.25 or greater.

Activate Conda environment for GATK.

	conda activate ngs_pipe_env_gatk4

Generate 'master' VCF.

	gatk SelectVariants -R C_excelsa_V5.fasta -V reheadered.F4_133.ANNOTATED.vcf -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn LOS_1 -sn LOS_2 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O UK_mav_dips_tets_for_plots.vcf

Generate tetraploid VCF.

	gatk SelectVariants -R C_excelsa_V5.fasta -V UK_mav_dips_tets.vcf -sn ROT_004 -sn ROT_006 -sn ROT_007 -sn ROT_013 -sn SKN_001 -sn SKN_002 -sn SKN_005 -sn SKN_008 -sn ALO_006 -sn ALO_007 -sn ALO_013 -sn ALO_017 -sn ELI_001 -sn ELI_002 -sn ELI_003 -sn ELI_004 -sn ERS_1 -sn ERS_2 -sn ERS_3 -sn ERS_4 -sn FTW_1 -sn FTW_2 -sn FTW_3 -sn FTW_5 -sn LAL_1 -sn LAL_2 -sn LAL_3 -sn LAL_4 -sn LNL_001 -sn LNL_002 -sn LNL_003 -sn LNL_008 -sn LOS_1 -sn LOS_2 -sn LOS_6 -sn LOS_7 -sn NEI_1 -sn NEI_3 -sn NEI_8 -sn NEI_9 -sn SCU_14 -sn SCU_15 -sn SCU_16 -sn SCU_19 -O UK_mav_tets_for_plots.vcf
			
Generate diploid VCF.

	gatk SelectVariants -R C_excelsa_V5.fasta -V UK_mav_dips_tets.vcf -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O /UK_mav_dips_for_plots.vcf

Deactivate conda.

	conda deactivate

**2\)** Convert to tables with the required columns (chromosome, position, allele count and allele number).

Activate conda.

	conda activate snpsift_env

For the tetraploids.

	SnpSift extractFields UK_mav_tets.vcf CHROM POS AC AN > UK_mav_plots_tets_for_plots.table

For the diploids.

	SnpSift extractFields ScanTools/UK_mav_dips.vcf CHROM POS AC AN > UK_mav_plots_dips_for_plots.table

Deactivate conda.

	conda deactivate
			
**3\)** Paste the tables together.

	paste UK_mav_plots_tets_for_plots.table UK_mav_plots_dips_for_plots.table > UK_dips_and_tets_MAV.table

**4\)** Manually fix headder of `UK_dips_and_tets_MAV.table`. The headder is:

	CHROM	POS	AC	AN	CHROM	POS	AC	AN

The headder should be:

	scaff	pos	tet_AC	tet_AN	CHROM	POS	dip_AC	dip_AN

**5\)** Generate `gene_orientation_file.txt` with `gene_orientation_file_from_genes_only.py`. This file should have no headder and the columns: gene name, scaffold, start position, end position and gene orientation (+ for the forward strand, - for the reverse strand). For example a line in this file should look like this: `g1	scaffold_33	1	2616	+`.

	python3 gene_orientation_file_from_genes_only.py

**6\)** Generate the MAV plots.

	python3 Sian_MAV_plots.py

The following files are hard coded into this Python script and their names/locations will need to be changed if you want to repeat this analysis on your own data:

line 11	`UK_dips_and_tets_MAV.table` (generated above)
line 16	`gene_orientation_file.txt` (generated above)
line 17	`FstH_1Percent.txt` (A list of genes, 1 gene per line)
line 18	`british_1percent_outliers_3.5_scores.txt` (generated above)

## Count MAV SNPs

Count the number of MAVs in each window of the BPM metrics file.

	python3 MAV_count_to_file.py -e genes_UK_dips_Vs_tets_15_BPM_1perc.txt -g genes_only.gtf -m british_1percent_outliers_3.5_scores.txt -o genes_UK_dips_Vs_tets_15_BPM_1perc_plus_MAVs.txt

Associate the window based counts directly to gene IDs.

	python3 crude_MAVs_to_genes_1perc.py

## Annotation of Cochlearia genes to a single _Arabidopsis thaliana_ homologue

The matching of each Cochlearia gene to a single most-similar homologue in the model organism _Arabidopsis thaliana_ used the script `1-1_annotation.py` which requires output files from runs of blastp, blast_rbh.py and Orthofinder (https://github.com/davidemms/OrthoFinder). The orthofinder results were later used to compare gene lists from selection scans in Cochlearia, _Cardamine amara_ and _Arabidopsis arenosa_ (aligned to the _Arabidopsis lyrata_ genome), hence the inclusion of these four species.

**1\)** Protein BLAST: `blastp`

Install BLAST with Homebrew (https://brew.sh/):

	brew install blast

Make a database out of the Arabidopsis thaliana protein fasta downloaded from the TAIR database (https://www.arabidopsis.org/) on the 21st December 2022.

	makeblastdb -in Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa -dbtype prot -out Athaliana_447_Araport11_protein_primaryTranscriptOnly_db

Run a BLAST search of the Arabidopsis thaliana database against the Cochlearia excelsa reference genome.

	blastp -evalue 0.1 -db Athaliana_447_Araport11_protein_primaryTranscriptOnly_db -query C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta -out cochlearia_thaliana_all_Vs_all_eval_0.1.tsv -outfmt 6

**2\)** Recipricol Best BLAST Hit: `blast_rbh.py`

Downloaded the code from https://github.com/peterjc/galaxy_blast on 19th December 2022. Copied the three python files to bin: 

	sudo mv ~/Downloads/galaxy_blast-master/tools/blast_rbh /usr/local/bin/

I then gave execute permissions:

	chmod +x /usr/local/bin/blast_rbh/*.py

Added an alias to my `~/.zprofile` file:

	alias blast_rbh='python3 /usr/local/bin/blast_rbh/blast_rbh.py'

Then ran a recipricol best BLAST search of:

	blast_rbh -a prot -t blastp -o cochlearia_vs_thaliana_rbh.txt C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta Athaliana_447_Araport11.protein_primaryTranscriptOnly.fa

**3\)** Orthofinder.

Install Orthofinder with Conda.

	conda create --name orthofinder_env -c bioconda orthofinder

Run Orthofinder on the `06_Orthofinder/` directory:

	/Users/sian_bray/miniconda3/envs/orthofinder_env/bin/orthofinder -f ~/06_Orthofinder/ -S blast

The log file notes the species used, which are the fasta files found in the `06_Orthofinder/` directory. These fasta files contian all of the protein sequences in the respective reference genomes:

	Species used: 
		0: Alyrata_384_v2.1.protein.fasta
		1: C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta
		2: Camara_protein.fasta
		3: TAIR10_pep_20110103_representative_gene_model.fasta

**4\)** Match each Cochlearia gene with the closest single _Arabidopsis thaliana_ homologue.

	python3 1-1_annotation.py -f C_excelsa_V5_braker2_wRseq.aa.LTPG.fasta -i Orthogroups.tsv -o Cochlearia_Thaliana_1-2-1_annotation.tsv -r cochlearia_vs_thaliana_rbh.txt -a cochlearia_thaliana_all_Vs_all_eval_0.1.tsv

## Annotate gene lists with gene descriptions

This script uses a gff/gtf file from a well studied genome (in this case _Arabidopsis thaliana_) to add gene names and descriptions to the one-on-one matches between the model (_A. thaliana_) and the non-model (Cochlearia) genes generated above (see "Annotation of Cochlearia genes to a single _Arabidopsis thaliana_ homologue"). The file `Araport11_GFF3_genes_transposons.2023-01-02.gff` is the _A. thaliana_ gff downloaded from the TAIR database (https://www.arabidopsis.org/) on 17th January 2023.

	python3 add_gene_descriptions_to_1-2-1_hits.py -d Araport11_GFF3_genes_transposons.2023-01-02.gff -o 1-2-1_hits_all_gene_descriptions.tsv -i Cochlearia_Thaliana_1-2-1_annotation.tsv

## GO Enrichment Analysis

Make the GO Universe using custom code. The file `gene_association.tair` is the gene ontology annotation for Arabidopsis thaliana downloaded from the TAIR database on 13th January 2023.

	python3 custom_GO_universe_restrictive.py -i Cochlearia_Thaliana_1-2-1_annotation.tsv -o Cochlearia_Thaliana_GO_universe_restrictive.tsv -og gene_association.tair

Run the GO enrichment analysis. The file `genes_in_scan_clean.txt` contains all the genes in the scan i.e. all genes that overlap windows with good enough coverage to generate BPM metrics (see "window-based scan for selective sweep signatures" above). This is slightly fewer than all the genes in the genome. The file `top_1perc_minus_inversion.txt` contains the genes in the top 1% outlier list from our selection scan (excluding genes that were present in a large inversion, becuase the inversion blocks recombination in that region).

	Rscript Cochlearia_GO_analysis.R

## Gene, orthologue and GO Overlaps between species

### Overlap Venn Diagrams and Statistics

**1\)** First gather the input files for the R script.

Top 1% outlier gene lists from the three selection scans, 1 gene per line. For Arabidopsis arenosa and Cardamine amara the top 1% outlier lists come from https://doi.org/10.1093/molbev/msab096.

	thalianaIDs_top_1perc_minus_inversion.txt
	top1perc_Fst_arenosa.txt
	top1perc_Fst_amara.txt

These are the output files from the above "GO Enrichment Analysis" section:

	GO_terms_cochlearia_0.05.txt
	GO_terms_arenosa_0.05.txt
	GO_terms_amara_0.05.txt

Match orthogroups to the genes in each of the three 1% outlier lists:

	python3 orthogroups_in_outlier_gene_list.py -io Orthogroups.tsv -ig top_1perc_minus_inversion.txt -o orthogroups_in_top1perc_cochlearia.txt

	python3 orthogroups_in_outlier_gene_list.py -io Orthogroups.tsv -ig top1perc_Fst_arenosa_lyrataIDs.txt -o orthogroups_in_top1perc_arenosa.txt

	python3 orthogroups_in_outlier_gene_list.py -io Orthogroups.tsv -ig top1perc_Fst_amara_amaraIDs.txt -o orthogroups_in_top1perc_amara.txt

Remove duplicate orthogroups from each of the files (NR stands for Non-Redundant).

	sort -u orthogroups_in_top1perc_cochlearia.txt > /orthogroups_in_top1perc_cochlearia_NR.txt

	sort -u orthogroups_in_top1perc_arenosa.txt > orthogroups_in_top1perc_arenosa_NR.txt

	sort -u orthogroups_in_top1perc_amara.txt > orthogroups_in_top1perc_amara_NR.txt

**2\)** Run the R script.

	Rscript Cochlearia_Venn.R

### Permutation test for GO overlaps

To test the statistical significance of the GO overlaps permutation tests were performed with custom code `GO_iteration_test_Cochlearia_2024_HPC.py`. This code ran GO analyses on random lists the same size as the candidate gene lists from the 3 species 10'000 times, then asked if the same or more overlaps were observed compared to the real data. It uses the following input:

`-l1`	Required. Path to the first list of total genes. This should be a list of all the genes in the genome, one gene per line.

`-l2`	Required. Path to the second list of total genes. This should be a list of all the genes in the genome, one gene per line.

`-l3`	Required. Path to the third list of total genes. This should be a list of all the genes in the genome, one gene per line.

`-pl1`	Optional. Path to the first list of total possible genes. That is all the genes that could possibly come up as candidate gene in the scan (e.g. enough coverage and SNPs). Csv file, one gene per line. If not provided the code defaults to the total gene list (-1).

`-pl2`	Optional. Path to the second list of total possible genes. That is all the genes that could possibly come up as candidate gene in the scan (e.g. enough coverage and SNPs). Csv file, one gene per line. If not provided the code defaults to the total gene list (-1).

`-pl3`	Optional. Path to the third list of total possible genes. That is all the genes that could possibly come up as candidate gene in the scan (e.g. enough coverage and SNPs). Csv file, one gene per line. If not provided the code defaults to the total gene list (-1).

`-n1`	Required. Lenght of the candidate gene lists from the first selection scan. This should be an integer.

`-n2`	Required. Lenght of the candidate gene lists from the second selection scan. This should be an integer.

`-n3`	Required. Lenght of the candidate gene lists from the third selection scan. This should be an integer.

`-go1`	Required. Full absolute path to the GO universe for the first genome/selection scan. A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.

`-go2`	Required. Full absolute path to the GO universe for the second genome/selection scan. A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.

`-go3`	Required. Full absolute path to the GO universe for the third genome/selection scan. A tsv file with two columns, one is the gene name, the other is all the GO terms associated with that gene. GO terms in the second column are seperated by commas. No headder.

`-p`	Default of 1'000. The number of permutations to run. This should be an integer. Default is 1000.

`-t`	Default of 0.05. The significance cutoff/P-value to be used in the GO enrichment analysis.

`-x_1_2`	Required. The amount of overlap (cross over) between GO terms in the original contrast between species 1 and 2. An integer.

`-x_1_3`	Required. The amount of overlap (cross over) between GO terms in the original contrast between species 1 and 3. An integer.

`-x_2_3`	Required. The amount of overlap (cross over) between GO terms in the original contrast between species 3 and 2. An integer.

`-x_all`	Required. The amount of overlap (cross over) between all three GO terms in the original contrast. An integer.

`-o`	Required. Path the output directory. I recomend an empty directory to contain the results of each permutation.`

Example run command:

	python3 GO_iteration_test_Cochlearia_2024_HPC.py -l1 all_genes_transcript_numbers.tsv -l2 all_genes.csv -l3 all_genes.csv -pl1 genes_in_scan_clean.txt -n1 752 -n2 451 -n3 227 -go1 Cochlearia_Thaliana_GO_universe_restrictive.tsv -go2 gene_GO.tsv -go3 gene_GO.tsv -p 10 -x_1_2 17 -x_1_3 1 -x_2_3 5 -x_all 0 -o ~/GO_Permutation_Test/

An example of the R code generated is found in `species_1_1.R`.

## Generation of Consensus Sequences and AlphaFold Structures

The Python script `AlphaFold_consensus.py` generates batch files to create consensus sequences for AlphaFold. These sequences need to be manually checked against non-biallelic sites and are manually translated using the Expasy Translate tool (https://web.expasy.org/translate/).

For the tetraploids:

	python3 AlphaFold_consensus.py -v reheadered.F4_133.ANNOTATED.vcf -r C_excelsa_V5.fasta -g C_excelsa_V5_braker2_wRseq.gff3 -i new_top_1perc_nr.txt -d ~/AlphaFold/Consensus -p tetraploids -s ROT_004 ROT_006 ROT_007 ROT_013 SKN_001 SKN_002 SKN_005 SKN_008 ALO_006 ALO_007 ALO_013 ALO_017 ELI_001 ELI_002 ELI_003 ELI_004 ERS_1 ERS_2 ERS_3 ERS_4 FTW_1 FTW_2 FTW_3 FTW_5 LAL_1 LAL_2 LAL_3 LAL_4 LNL_001 LNL_002 LNL_003 LNL_008 LOS_1 LOS_2 LOS_6 LOS_7 NEI_1 NEI_3 NEI_8 NEI_9 SCU_14 SCU_15 SCU_16 SCU_19

For the diploids:

	python3 AlphaFold_consensus.py -v reheadered.F4_133.ANNOTATED.vcf -r C_excelsa_V5.fasta -g C_excelsa_V5_braker2_wRseq.gff3 -i new_top_1perc_nr.txt -d ~/AlphaFold/Consensus -p diploids -s  BNK_21 CHA_1 CHA_2 JOR_1 JOR_12 JOR_13 JOR_3 LAB_004 LAB_1 LAB_2 LAB_300 LAB_4 LAB_400 LAB_5 LAB_500 NEN_001 NEN_003 NEN_200 NEN_300 NEN_4 NEN_5 NEN_6 ODN_10 ODN_2 ODN_4 ODN_5 ODN_6 ODN_7 ODN_9

Examples of the batch files created are:

	g7445_tetraploids_consensus_generation.sh
	g7445_diploids_consensus_generation.sh

There is a record of all manual changes in `consensus_sequence_record_of_manual_changes.txt`.

## Protein modeling

AlphaFold structures were generated on the Czech national HPC Metacentrum on 10th August 2023. The batch files used to do so are in the folder `AlphaFold`.

## Extras

This project has taken many years to reach its final state, passing through several 'final' versions on the way. In the folder `Extras/` there are some scripts that were written for this project, but not used in the final (final_actually_final_really_final_FINAL) version of the analysis that appears in the published manuscript. They are included here in case they are ever useful to anyone.

**A\)** An older messier method of generating AFD (Allele Frequency Difference) plots without the MAVs that the above MAV-plot code was based on. The R parts of this code were originally written by Christian Sailer (https://github.com/SailerChristian) and I added the Python wrapper and the readme (`00_instructions.sh`) to streamline it.

File in `00_Blank_Folder/`:

	Find_BPM_outliers.py
	BPM_to_outlier_dot_plots.py
	Annotated_Hits.txt
	temp_files
	genes_only.gtf
	gene_orientation_file.txt
	00_instructions.sh
	R.gene_window_plot.py

**B\)** An alternative method to compare gene lists for overlaps that was not used in the end.

	compare_gene_lists_and_orthogroups_V2.py

**C\)** A script that can add the annotation (gene name and description) to a list of genes using a reciprocol-best-blast-hits or orthologue groups file

	python3 add_gene_descriptions_to_rbh_and_orthogroups.py -b Orthogroups.tsv -d Araport11_functional_descriptions_20181231.txt -o Annotated_Orthogroups.tsv -bc 1

**D\)** A script that will annotate you candidate gene list using the output from `add_gene_descriptions_to_rbh_and_orthogroups.py`.

	python3 annotate_gene_list.py -i UK_0.1perc.txt -a Annotated_Orthogroups.tsv -o 0.1perc_fst_outliers_annotated.tsv

**E\)** If you forgot to filter by AF before generating your table of MAV SNPs you can do this afterwards with:

	python3 crude_remove_low_AF_from_MAV.py

**F\)** Generates one VCF file from another with sites no closer than a set value, to use as 'unlinked' data for fastStructure.

	python fastStructure_vcf_generation_diversity.py -i test.vcf -o structure_test.vcf -r structure_test.fai -g 2 -f 0.015

**G\)** Create a file in a fastSTRUCTURE format from a VCF file. Original code written by Patrick Monnahan and Jordan Koch then modified for the Cochlearia data.

	python Cochlearia_create_structure_file.py -v ~/Fast_Structure_100K/ -o 10alt_4missing -s true

**H\)** Remove duplicate orthogroups from an id2gos file.

	python3 remove_repeted_GO_terms.py

**I\)** Remove obsolete GO terms in an id2gos file.

	python3 obsolete_GO_term_fix_26Jun2023.py

**J\)** Run all pairwise between population metrics (BPMs) in the Cochlearia data:

	python3 pairwise_bpm_all_V2.py

**K\)** Concatonate fastas composed of several exons

	python3 combine_exons_in_fasta.py

**L\)** Add BLAST E-values to a list of genes:

	python3 crude_Evals_to_801.py

**M\)** Convert a ScanTools output file to a bed file:

	python3 scantools_to_bed.py -st UK_dips_Vs_tets_15_BPM.txt -b UK_dips_Vs_tets_15_BPM.bed

**N\)** Remove synonomous SNPs from MAV output:

	python3 remove_synonomous_from_PiritaMAV.py

## Spelling

I am dyslexic and the code provided will be riddled with spelling errors (no spell checker in Sublime Text haha). For information and context you can visit the British Dyslexia Association (https://www.bdadyslexia.org.uk/) or Google "dyslexic scientists".

