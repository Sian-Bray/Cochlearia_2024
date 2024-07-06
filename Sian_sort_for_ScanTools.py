# Example command:
#   python3 Sian_sort_for_ScanTools.py 'British_diploids.table.recode.txt British_tetraploids.table.recode.txt ' VCF_Cochlearia_76_Filtered_DP8.M0.2/ British_dips_Vs_tets
import sys, os

file_string=sys.argv[1]
tmpdir=sys.argv[2]
output_name=sys.argv[3]

scaffold_list=[]
removal_list=[]

#create scaffold list from population files
split_file_string=file_string.split(' ')
try:
    split_file_string.remove('')
except ValueError:
    pass
for count1, file1 in enumerate(split_file_string):
    open_file=open(file1, 'r')
    for count2, line in enumerate(open_file):
        if line != '\n':
            line=line.split('\t')
            scaffold=line[2]
            if scaffold not in scaffold_list:
                scaffold_list.append(scaffold)

#create files containing one scaffold for each file with grep
for count3, scaff1 in enumerate(scaffold_list):
    for count4, file2 in enumerate(split_file_string):
        os.system('grep -w '+scaff1+' '+file2+' > '+file2+'.'+scaff1+'.ThIsNeeDSRem0vinG')
        removal_list.append(file2+'.'+scaff1+'.ThIsNeeDSRem0vinG')

#Sort and merge each scaffold file with sort
for count5, scaff2 in enumerate(scaffold_list):
    current_list=''
    for count6, file3 in enumerate(split_file_string):
        current_list+=file3+'.'+scaff2+'.ThIsNeeDSRem0vinG '
    os.system('sort -k4,4n -m '+current_list+'> '+tmpdir+scaff2+output_name+'.sort.ThIsNeeDSRem0vinG')
    removal_list.append(tmpdir+scaff2+output_name+'.sort.ThIsNeeDSRem0vinG')
    os.system('cat '+tmpdir+scaff2+output_name+'.sort.ThIsNeeDSRem0vinG >> '+tmpdir+output_name+'.concat.txt')

#Remove temp files
for count7, file4 in enumerate(removal_list):
    os.system('rm '+file4)

#Testing command: python3 Sian_sort_for_ScanTools.py '/Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/dips.table.recode.txt /Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/tets.table.recode.txt' /Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/ British_dips_Vs_tets_WS10000_MS25
#Testing command: python3 Sian_sort_for_ScanTools.py '/Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/sort_tets.txt /Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/sort_dips.txt /Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/sort_trips.txt' /Users/brays/Desktop/Major_ScanTools_Issue/Test_data/Sort_test/ sort_dips_Vs_trips_Vs_tets_WS10000_MS25
