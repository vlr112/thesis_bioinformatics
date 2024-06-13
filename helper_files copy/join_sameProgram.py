# import sys
# import pandas as pd
# from fractions import Fraction
# import numpy as np
# import os

import pandas as pd
import glob
import sys
import os
from functools import reduce                # Import reduce function

# prefix = sys.argv[1]
file_list = sys.argv[1]
outifle = sys.argv[2]
with open(file_list, "r") as file:
    contents = file.read()
files = contents.split("\n")[:-1]

count =0

colslist = ['Ind1', 'Ind2', 'k0_hat', 'k1_hat', 'k2_hat', 'pi_HAT', 'nbSNP', 'cov', 'deamination', 'kinship', 'degree', 'program']

dfs = []
for i in range(len(files)):

    # with open(files[i]) as f:
    #     first_line = f.readline()
    #     print(first_line)
    # this serves to check if file is not empty. it'll be empty if fx, calculations for sim are completed 
    if os.stat(files[i]).st_size != 0:
        # print(files[i])
        df = pd.read_csv(files[i], sep="\t")
        # df = df[['Ind1', 'Ind2', 'k0_hat', 'k1_hat', 'k2_hat', 'pi_HAT', 'nbSNP', 'cov', 'deamination', 'kinship', 'degree', 'program']]
        # # print(len(df.columns))
        # # if ! colslist.isin.df.columns
        # if files[i] == './sim20/deam/cov8/results/lcmlkin/lcmlkin_results_update_out' :
        #     print(list(df.columns))
        #     print(df)
        #     print(df[df[df.columns].isin(colslist)])
        #         # print('buiaaaa')



        if len(df.columns) > 3:
        #     # print(files[i])
        #     # df1 = df.filter(['ID1', 'ID2', 'N_SNP', 'HetHet', 'IBS0', 'Kinship', 'cov', 'deamination', 'kinship', 'degree', 'program'], axis=1)
            dfs.append(df)
#             # print(df.columns)
#             # print(files[i])
# #             if len(df.columns)>15:
# #                 print(files[i])

# #                 count += 1
# # print(count)


# # merged = pd.concat(dfs, join="inner", axis=0)

merged = pd.concat(dfs, ignore_index=True)

# # Write the merged dataframe to a new file
merged.to_csv(outifle, sep="\t", index=False)






