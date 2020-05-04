#!/usr/bin/env python
#Author: Alba Sanchis-Juan (as2635@cam.ac.uk)

import sys
import pandas as pd
import numpy as np

infile = sys.argv[ 1 ]

df = pd.read_csv(infile, sep = '\t', names = ['chr', 'pos', 'cov'])

mean = df['cov'].mean()
median = df['cov'].median()
max_cov = df['cov'].max()

df_count = df['cov'].value_counts(normalize=True).to_frame().reset_index()

df_count.columns = ['cov', 'percent']

cols = ['cov', 'min_percent']
df_percent = pd.DataFrame(columns = cols)

for v in range(1, int(max_cov)+1):
    col2 = df_count.loc[df_count['cov'] >= v, 'percent'].sum()
    data = pd.DataFrame([{cols[0] : v, cols[1]: col2}])
    df_percent = df_percent.append(data)

print( 'mean', mean )
print( 'median', median )
print( df_percent.to_csv(sep=' ', index=False, header=False) )
