#### changes a column in original files with rmse value of caparison with ibds calculated from fasta files

import pandas as pd
import numpy as np
import sys

# names = sys.argv[1]
names='names.txt'
with open(names, "r") as file:
    contents = file.read()
individuals = contents.split("\n")[:-1]

true = pd.DataFrame(pd.read_csv('trueIBD', delimiter="\t"))

inpath=sys.argv[1] 
filename = sys.argv[2] 
outfilename= sys.argv[3]

observed = pd.DataFrame(pd.read_csv(inpath + filename, delimiter="\t"))


# path_to_file= '/projects/korneliuss&en/people/vlr112/final_simulations/results/lcmlkin/'

if filename == 'lcmlkin_results_out':
    observed.rename(columns={'k0_hat': 'J9'}, inplace=True)
    observed.rename(columns={'k1_hat': 'J8'}, inplace=True)
    observed.rename(columns={'k2_hat': 'J7'}, inplace=True)
    observed.rename(columns={'Ind1': 'a'}, inplace=True)
    observed.rename(columns={'Ind2': 'b'}, inplace=True)

if filename == 'plink_out.genome':
    observed.rename(columns={'Z0': 'J9'}, inplace=True)
    observed.rename(columns={'Z1': 'J8'}, inplace=True)
    observed.rename(columns={'Z2': 'J7'}, inplace=True)
    observed.rename(columns={'IID1': 'a'}, inplace=True)
    observed.rename(columns={'IID2': 'b'}, inplace=True)


## remove all possible trash
observed = observed[observed['a'].isin(individuals)]

# after trahs removal (if aplicable) be sure to set J columns to float type
observed['J9'] = observed['J9'].astype('float')
observed['J8'] = observed['J8'].astype('float')
observed['J7'] = observed['J7'].astype('float')


# # confirm that both individuals exist in both true and observed df
tmp = observed[['a', 'b']]
true_final = pd.merge(true, tmp, how='inner', on=['a', 'b'])


cols = ['J9','J8','J7']
observed['rmse'] = np.sqrt(np.sum(np.square(observed[cols] - true_final[cols]), axis=1)/3)

observed.to_csv(outfilename, sep="\t",  index=False)
