"""
Gather _result_table.tsv files from AmpliconSuite for further analysis

python gather_results.py [results directory] [output CSV]
"""

import sys
import os
import pandas as pd

results_dir = sys.argv[1]
outfn = sys.argv[2]
# results_dir = '/scratch/users/tbencomo/ranson-wgs/results'
samples = os.listdir(results_dir)

df_list = []
for sample in samples:
    fn = f'{sample}_result_table.tsv'
    fp = os.path.join(results_dir, f'{sample}', f'{sample}_classification', fn)
    if os.path.exists(fp):
        df = pd.read_csv(fp, sep='\t')
        df_list.append(df)
    else:
        print(f'No results file found for {sample}! Expecting {fp}')
df = pd.concat(df_list, ignore_index = True)
print(df.head())
df.to_csv(outfn, index=False)
# df.to_csv('AC_Results.csv', index = False)
