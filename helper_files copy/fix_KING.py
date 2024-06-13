import pandas as pd
import sys

path_to_file = sys.argv[1]
relate = sys.argv[2]

filename = sys.argv[3]
df = pd.DataFrame(pd.read_csv(path_to_file + '/' + filename, delimiter="\t"))
correction = pd.DataFrame(pd.read_csv(relate, delimiter="\t"))

final = pd.merge(correction[['a','b']], df  , right_on=['ID1','ID2'], left_on=['a', 'b'], how='left').drop(columns=['a', 'b'])
final.iloc[0] = df.loc[0]

final.to_csv(path_to_file + '/' + filename + '2', sep="\t",  index=False)
