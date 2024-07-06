# runs all pairwise bpm for the Cochlearia samples (Dec 2022)

individuals = ['AAH_1','AAH_2','AAH_3','AAH_4','ALO_006','ALO_007','ALO_013','ALO_017','BEA_002','BEA_004','BEA_010','BNK_21','BRE_1','BRI_002','BRI_005','BRI_006','BRI_009','CHA_1','CHA_2','CUM_1','DAR_1','DAR_3','ELI_001','ELI_002','ELI_003','ELI_004','ERS_1','ERS_2','ERS_3','ERS_4','FOR_1','FRE_013','FTW_1','FTW_2','FTW_3','FTW_5','GEO_2','GEO_6','JON_001','JOR_1','JOR_12','JOR_13','JOR_3','KVA_002','KVA_003','KVA_009','KVA_010','LAB_004','LAB_1','LAB_2','LAB_300','LAB_4','LAB_400','LAB_5','LAB_500','LAL_1','LAL_2','LAL_3','LAL_4','LNL_001','LNL_002','LNL_003','LNL_008','LOS_1','LOS_2','LOS_6','LOS_7','MEL_001','MEL_002','NEI_1','NEI_3','NEI_8','NEI_9','NEN_001','NEN_003','NEN_200','NEN_300','NEN_4','NEN_5','NEN_6','ODN_10','ODN_2','ODN_4','ODN_5','ODN_6','ODN_7','ODN_9','ROT_004','ROT_006','ROT_007','ROT_013','RUZ_001','RYE_1','SAL_003','SAL_004','SAL_006','SCO_1','SCU_14','SCU_15','SCU_16','SCU_19','SKF_002','SKF_003','SKF_005','SKF_009','SKI_004','SKI_005','SKN_001','SKN_002','SKN_005','SKN_008','SPU_006','SPU_008','SPU_009','SPU_010','TET_002','TET_004','TET_006','TET_008','TRO_001','TRO_003','TRO_005','TRO_009','VAG_003','VAG_004','VAG_007','VAG_102','VEG_003','VEG_004','WOL_002','WOL_006','WOL_009','WOL_010']
pops = ['AAH','ALO','BEA','BNK','BRE','BRI','CHA','CUM','DAR','ELI','ERS','FOR','FRE','FTW','GEO','JON','JOR','KVA','LAB','LAL','LNL','LOS','MEL','NEI','NEN','ODN','ROT','RUZ','RYE','SAL','SCO','SCU','SKF','SKI','SKN','SPU','TET','TRO','VAG','VEG','WOL']
in_vcf = '/gpfs01/home/mbzly/2022.Cochlearia.MS1/VCFs/Cochlearia_112_dip_tet_4dg.purged.ann.vcf' # full path to vcf
shell_location = '/gpfs01/home/sbzsmb/Cochlearia/all_vs_all_bpm/' # full path, make sure you include the last slash
#shell_location = '/Users/sian_bray/' # full path, make sure you include the last slash

# Make all combinations
pop_combos = []

# Creates a list of lists, each sub list contains two pops
for pop1 in pops:
	for pop2 in pops:
		if pop1 != pop2:
			if [pop1, pop2] not in pop_combos:
				if [pop2, pop1] not in pop_combos:
					pop_combos.append([pop1, pop2])

# write a shell file for each combo!
# Give the shell 2 hours and 16g
for combo in pop_combos:
	combo_name=combo[0]+'_Vs_'+combo[1]+'_bpm'
	shell = open(shell_location+combo_name+'.sh', 'w')
	shell.write(f'#!/bin/bash\n#SBATCH --job-name={combo_name}\n#SBATCH --partition=defq\n#SBATCH --nodes=1\n#SBATCH --ntasks-per-node=1\n#SBATCH --mem=16g\n#SBATCH --time=04:00:00\n#SBATCH --output=/gpfs01/home/sbzsmb/OandE/%x.out\n#SBATCH --error=/gpfs01/home/sbzsmb/OandE/%x.err\nsource $HOME/.bash_profile\nconda activate ngs_pipe_env_gatk4\n')

	# Make a vcf and table for each pop
	# Pop 1 gatk (51 mins 7'576'832 K) e.g. command:
	# gatk SelectVariants -R /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ionops/ref.gen/C_excelsa_V5.fasta -V /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/10.filtered.depth/reheadered.F4_133.ann.vcf.gz -sn BNK_21 -sn CHA_1 -sn CHA_2 -sn JOR_1 -sn JOR_12 -sn JOR_13 -sn JOR_3 -sn LAB_004 -sn LAB_1 -sn LAB_2 -sn LAB_300 -sn LAB_4 -sn LAB_400 -sn LAB_5 -sn LAB_500 -sn NEN_001 -sn NEN_003 -sn NEN_200 -sn NEN_300 -sn NEN_4 -sn NEN_5 -sn NEN_6 -sn ODN_10 -sn ODN_2 -sn ODN_4 -sn ODN_5 -sn ODN_6 -sn ODN_7 -sn ODN_9 -O /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_dips.vcf
	shell.write(f'gatk SelectVariants -R /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ionops/ref.gen/C_excelsa_V5.fasta -V {in_vcf} --max-nocall-fraction 0.2')
	for pop in combo:
		for ind in individuals:
			if pop in ind:
				shell.write(f' -sn {ind}')
	shell.write(f' -O {shell_location}{combo[0]}_{combo[1]}.vcf\n')

	for pop in combo:
		shell.write(f'gatk SelectVariants -R /gpfs01/home/mbzly/2022.Cochlearia.MS1/ngs_pipe/ionops/ref.gen/C_excelsa_V5.fasta -V {shell_location}{combo[0]}_{combo[1]}.vcf')
		for ind in individuals:
			if pop in ind:
				shell.write(f' -sn {ind}')
		shell.write(f' -O {shell_location}{pop}.vcf\n')
		# pop 1 recode (26 mins 10'746'568 K)
		# gatk VariantsToTable -V /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_tets.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_tets.table
		shell.write(f'gatk VariantsToTable -V {shell_location}{pop}.vcf -F CHROM -F POS -F REF -F AN -F DP -GF GT -O {shell_location}{pop}.table\n')
		# python3 /gpfs01/home/sbzsmb/Scripts/recode012.py -i /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_tets.table.unphased -o /gpfs01/home/sbzsmb/Cochlearia/ScanTools/ -pop UK_TETS
		shell.write(f'python3 /gpfs01/home/sbzsmb/Scripts/recode012.py -i {shell_location}{pop}.table -o {shell_location} -pop {pop}\n')

	# merge (2 mins 3968K)
	# python3 /gpfs01/home/sbzsmb/Scripts/Sian_sort_for_ScanTools.py '/gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_tets.table.unphased.recode.txt /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_dips.table.unphased.recode.txt ' /gpfs01/home/sbzsmb/Cochlearia/ScanTools/Temp/ UK_dips_Vs_tets
	shell.write(f"python3 /gpfs01/home/sbzsmb/Scripts/Sian_sort_for_ScanTools.py '{shell_location}{combo[0]}.table.recode.txt {shell_location}{combo[1]}.table.recode.txt ' {shell_location} {combo_name}\n")

	# bpm (3 h 35 mins 7564K)
	# python3 /gpfs01/home/sbzsmb/Scripts/sian_bpm.py -i /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_dips_Vs_tets.concat.txt -o /gpfs01/home/sbzsmb/Cochlearia/ScanTools/UK_bpm/ -prefix UK_dips_Vs_tets_15 -ws 1000 -ms 15 -np 2
	shell.write(f'python3 /gpfs01/home/sbzsmb/Scripts/sian_bpm.py -i {shell_location}{combo_name}.concat.txt -o {shell_location} -prefix {combo_name} -ws 10000 -ms 1 -np 2')

	shell.close()











