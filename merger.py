import pandas as pd
import glob
import sys

# prefix = sys.argv[1]
dir_str = sys.argv[1:-1]
name = sys.argv[-1]
# print(dir_str)

# dir_list = dir_str.split(',')


# Get all files that start with prefix and end with .txt in subdirectories
# file_list = glob.glob("./results/{}/{}*.txt".format(dir,prefix))

# dfs = []
# for dir in dir_str:
#     # file= (glob.glob("./results/{}/{}*.txt".format(dir, dir)))
#     df = pd.read_csv("results/{}/{}.txt".format(dir, dir), sep="\t")
#     if dir == "lcmlkin":  

#         df.rename(columns={'k0_hat': 'J9'}, inplace=True)
#         df.rename(columns={'k1_hat': 'J8'}, inplace=True)
#         df.rename(columns={'k2_hat': 'J7'}, inplace=True)
#         df.rename(columns={'Ind1': 'a'}, inplace=True)
#         df.rename(columns={'Ind2': 'b'}, inplace=True)
#         df = df.dropna(axis=0, thresh=11)
#         # df['J9'] = df['J9'].astype(float)


#         # print(df.iloc[0])
#         # print('LCMLKIN')
#         # print(df.dtypes)

#     ## when data_file is lcmlkin:
#     if dir == "plink":

#         df.rename(columns={'Z0': 'J9'}, inplace=True)
#         df.rename(columns={'Z1': 'J8'}, inplace=True)
#         df.rename(columns={'Z2': 'J7'}, inplace=True)
#         df.rename(columns={'IID1': 'a'}, inplace=True)
#         df.rename(columns={'IID2': 'b'}, inplace=True)

#         # # print(df.iloc[0])
#         # print(df.dtypes)
#         # print('PLINK')

#     dfs.append(df)

# merged = pd.concat(dfs, join="inner", axis=0)

# # Write the merged dataframe to a new file
# merged.to_csv('results/merged_ALL_IBD.txt', sep="\t", index=False)



#####

# correctKin = pd.DataFrame(pd.read_csv('results/correctKin/correctKin.txt', delimiter="\t" ))
# correctKin['program'] = 'correctKin'
# relate = pd.DataFrame(pd.read_csv('results/relateREAL/relateREAL.txt', delimiter="\t" ))
# king = pd.DataFrame(pd.read_csv('results/king/king.txt', delimiter="\t" ))
# plink = pd.DataFrame(pd.read_csv('results/plink/plink.txt', delimiter="\t" ))
# lcmlkin = pd.DataFrame(pd.read_csv('results/lcmlkin/lcmlkin.txt', delimiter="\t" ))
# # read = pd.DataFrame(pd.read_csv('results/READ/READ_results', delimiter="\t" ))

# relate.rename(columns={'a': 'ID1', 'b': 'ID2'}, inplace=True)
# king.rename(columns={'Kinship': 'theta'}, inplace=True)
# lcmlkin.rename(columns={'Ind1': 'ID1', 'Ind2': 'ID2', 'pi_HAT': 'theta'}, inplace= True)
# plink.rename(columns={'IID1': 'ID1', 'IID2': 'ID2', 'PI_HAT': 'theta'}, inplace=True)
# correctKin.rename(columns={ 'a': 'ID1', 'b': 'ID2','corr. kin. coeff.': 'theta'}, inplace=True)

# dfs = [relate, king, lcmlkin, plink, correctKin]

# dfs_final = []
# for i in range(len(dfs)):
#     coisa = pd.DataFrame(dfs[i])[['ID1', 'ID2', 'theta', 'cov', 'deamination', 'kinship', 'degree',	'program']]
#     dfs_final.append(coisa)

# df = pd.concat(dfs_final)


# df['ID'] = df['ID1'] + '-' +  df['ID2']
# df['cov'].replace(1, 1.0, inplace=True)
# df['cov'].replace(2.0, 2, inplace=True)
# df['cov'].replace(4.0, 4, inplace=True)
# df['cov'].replace(8.0, 8, inplace=True)
# df['cov'].replace(10.0, 10, inplace=True)
# df['cov'].replace(20.0, 20, inplace=True)

# # print(df)

# df.to_csv('results/merged_ALL.txt', sep="\t", index=False)


##############################


relate = pd.DataFrame(pd.read_csv('results/relateREAL/relateREAL.txt', delimiter="\t" ))
plink = pd.DataFrame(pd.read_csv('results/plink/plink.txt', delimiter="\t" ))
lcmlkin = pd.DataFrame(pd.read_csv('results/lcmlkin/lcmlkin.txt', delimiter="\t" ))

# relate.rename(columns={'a': 'ID1', 'b': 'ID2', 'J9': 'k0', 'J8': 'k1', 'J7': 'k2' }, inplace=True)
# lcmlkin.rename(columns={'Ind1': 'ID1', 'Ind2': 'ID2', 'k0_hat': 'J9', 'k1_hat': 'J8', 'k2_hat': 'J7'}, inplace= True)
# plink.rename(columns={'IID1': 'ID1', 'IID2': 'ID2', 'Z0': 'k0', 'Z1': 'k1', 'Z2': 'k2'}, inplace=True)

# relate.rename(columns={'a': 'ID1', 'b': 'ID2' }, inplace=True)
lcmlkin.rename(columns={'Ind1': 'a', 'Ind2': 'b', 'k0_hat': 'J9', 'k1_hat': 'J8', 'k2_hat': 'J7'}, inplace= True)
plink.rename(columns={'IID1': 'a', 'IID2': 'b', 'Z0': 'J9', 'Z1': 'J8', 'Z2': 'J7'}, inplace=True)



dfs = [relate, lcmlkin, plink]

dfs_final = []
for i in range(len(dfs)):
    coisa = pd.DataFrame(dfs[i])[['a', 'b', 'J9', 'J8', 'J7', 'cov', 'deamination', 'kinship', 'degree',	'program']]
    dfs_final.append(coisa)

df = pd.concat(dfs_final)


df['ID'] = df['a'] + '-' +  df['b']
df['cov'].replace(1, 1.0, inplace=True)
df['cov'].replace(2.0, 2, inplace=True)
df['cov'].replace(4.0, 4, inplace=True)
df['cov'].replace(8.0, 8, inplace=True)
df['cov'].replace(10.0, 10, inplace=True)
df['cov'].replace(20.0, 20, inplace=True)

# print(df)

df.to_csv('results/merged_ALL_IBD.txt', sep="\t", index=False)
