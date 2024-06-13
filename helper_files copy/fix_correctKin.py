import pandas as pd
import sys


path_to_file = sys.argv[1]

filename = sys.argv[2]

correctKin = pd.DataFrame(pd.read_csv(path_to_file + '/' + filename, delimiter="\t"))
correctKin_fix = correctKin.rename(columns={'ID1': 'a', 'ID2': 'b'})  

relate = pd.DataFrame(pd.read_csv('sim20/no_deam/cov10/results/relate/real_data_out', delimiter="\t" ))


final = pd.merge(correctKin_fix, relate[['a', 'b', 'kinship', 'degree' ]], on=['a', 'b'])



final.to_csv(path_to_file + '/out_angsdHaploCall.relatives.tsv' , sep="\t",  index=False)

