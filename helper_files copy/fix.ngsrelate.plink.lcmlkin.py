import sys
import pandas as pd
from fractions import Fraction
import numpy as np
import os


# theoretical_dict = {'Relationship': ['PO', 'HS', 'FS', 'FC', 
#                          'DFC', 'SC', 'A', 'Offspring of sib-matings', 'Unrelated'],
#         'J9': [Fraction(0, 1), Fraction(1,2), Fraction(1, 4), Fraction(3, 4), Fraction(9, 16), Fraction(15, 16), Fraction(1, 2), Fraction(1, 16), Fraction(0, 1)],
#         'J8': [Fraction(1, 1), Fraction(1,2), Fraction(1, 2), Fraction(1, 4), Fraction(6, 16), Fraction(1, 16), Fraction(1, 2), Fraction(5, 16), Fraction(0, 1)],
#         'J7': [Fraction(0,1), Fraction(0,1), Fraction(1, 4), Fraction(0, 1), Fraction(1, 16), Fraction(0, 1), Fraction(0, 1), Fraction(7, 32), Fraction(0, 1)],
#         'J6': [ Fraction(0,1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 32), Fraction(0, 1)],
#         'J5': [ Fraction(0,1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 8), Fraction(0, 1)],
#         'J4': [ Fraction(0,1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 32), Fraction(0, 1)],
#         'J3': [ Fraction(0,1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1,8), Fraction(0, 1)],
#         'J2': [Fraction(0,1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 32), Fraction(0, 1)],
#         'J1': [Fraction(0,1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(0, 1), Fraction(1, 16), Fraction(0, 1)]
#         }

theoretical_dict = {'Relationship': ['PO', 'HS', 'FS', 'FC', 
                         'DFC', 'SC', 'A', 'Offspring of sib-matings', 'Unrelated'],
        'J9': [Fraction(0, 1), Fraction(1,2), Fraction(1, 4), Fraction(3, 4), Fraction(9, 16), Fraction(15, 16), Fraction(1, 2), Fraction(1, 16), Fraction(0, 1)],
        'J8': [Fraction(1, 1), Fraction(1,2), Fraction(1, 2), Fraction(1, 4), Fraction(6, 16), Fraction(1, 16), Fraction(1, 2), Fraction(5, 16), Fraction(0, 1)],
        'J7': [Fraction(0,1), Fraction(0,1), Fraction(1, 4), Fraction(0, 1), Fraction(1, 16), Fraction(0, 1), Fraction(0, 1), Fraction(7, 32), Fraction(0, 1)]
        }

theoretical= pd.DataFrame.from_dict(theoretical_dict, orient='columns')
theoretical.set_index("Relationship", inplace=True)
theoretical = theoretical.astype(float)
theoretical.reset_index(inplace=True)


data_file = sys.argv[1] 
path_to_file= sys.argv[2]
# open the input file for reading
with open(data_file, 'r') as f:
    # read the header row
    header = f.readline().strip()
    headers = header.strip().split('\t')

    # ## when data_file is lcmlkin:
    # if data_file == 'lcmlkin.txt':
    #     # observed.rename(columns={'J9': 'lol'}, inplace=True)

    # #     headers[2]='J9'
    # #     headers[3]='J8'
    # #     headers[2]='J7'

    # # ## when data_file is lcmlkin:
    # # if data_file == 'plink.txt':
    # #     headers[5]='J9'
    # #     headers[6]='J8'
    # #     headers[7]='J7'

    # initialize a dictionary to store the rows
    rows = {}
    # read each row and add it to the dictionary
    for line in f:
        # split the line into columns
        cols = line.strip().split('\t')
        if len(cols) < len(headers):
            continue
        # use the last two columns as the key
        key = tuple(cols[-3:])
        # add the row to the dictionary
        if key not in rows:
            # create a new list for this key
            rows[key] = [header, line]
        else:
            # add the row to the existing list for this key
            rows[key].append(line)

# write each group of rows to a separate file
for key, group in rows.items():

    # create the filename by joining the last two columns with '_'
    filename = '_'.join(key) + '.txt'

    if filename == 'kinship_degree_program.txt':
        continue
    # open the output file for writing
    with open(path_to_file+filename, 'w') as f:
        # write the header row
        f.write(group[0] + '\n')
        # write each row in the group
        for line in group[1:]:
            f.write(line)
        # f.close()


    observed = pd.DataFrame(pd.read_csv(path_to_file+filename, delimiter="\t"))
    ## when data_file is lcmlkin:
    if path_to_file == "/projects/korneliussen/people/vlr112/final_simulations/results/lcmlkin/":  
    # 'lcmlkin.txt':
    # list_lcmlkin=['A_2_lcmlkin.txt', 'FC_3_lcmlkin.txt', 'FS_1_lcmlkin.txt', 'GP_2_lcmlkin.txt', 'HC_4_lcmlkin.txt', 'HS_2_lcmlkin.txt', 'PO_1_lcmlkin.txt', 'Unrelated_Unrelated_lcmlkin.txt']
    # if filename in list_lcmlkin:
        observed.rename(columns={'k0_hat': 'J9'}, inplace=True)
        observed.rename(columns={'k1_hat': 'J8'}, inplace=True)
        observed.rename(columns={'k2_hat': 'J7'}, inplace=True)
        observed.rename(columns={'Ind1': 'a'}, inplace=True)
        observed.rename(columns={'Ind2': 'b'}, inplace=True)


    ## when data_file is lcmlkin:
    if path_to_file == "/projects/korneliussen/people/vlr112/final_simulations/results/plink/":
    # list_plink=['A_2_plink.txt', 'FC_3_plink.txt', 'FS_1_plink.txt', 'GP_2_plink.txt', 'HC_4_plink.txt', 'HS_2_plink.txt', 'PO_1_plink.txt', 'Unrelated_Unrelated_plink.txt']
    # if filename in list_plink:
    
    # if filename == 'plink.txt':
        observed.rename(columns={'Z0': 'J9'}, inplace=True)
        observed.rename(columns={'Z1': 'J8'}, inplace=True)
        observed.rename(columns={'Z2': 'J7'}, inplace=True)
        observed.rename(columns={'IID1': 'a'}, inplace=True)
        observed.rename(columns={'IID2': 'b'}, inplace=True)


    # get the kinship value from observed data
    kinship = observed['kinship'][0]

    # get the row with the kinship value from theoretical data
    row = theoretical.loc[theoretical['Relationship'] == kinship]

    if not row.empty:
        common_cols = list(set(observed.columns) & set(row.columns))
        common_cols.sort(reverse=True)
        df3 = observed.copy()

        coisa = row[common_cols].reset_index()
        coisa.drop('index', axis=1, inplace=True)

        rmse=[]
        for index, row in df3[common_cols].iterrows():
            tmp=df3[common_cols].iloc[index].reset_index()
            tmp.drop('index', axis=1, inplace=True)
            diff = np.sqrt(np.sum(np.square(np.array(tmp[index]) - np.array(coisa.T[0])))/3)
            rmse.append(diff)

        # observed.index+=1
        observed['RMSE'] = rmse
        

        observed.to_csv(path_to_file+filename, sep="\t",  index=False)

        f.close()

    else:
        continue

    



