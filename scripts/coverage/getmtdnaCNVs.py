import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Functions definition
def getwsavg(index,coverage, lastpos):
    diff = lastpos - index
    if diff >= ws:
        avg = np.mean(np.array(coverage[index:int(index+ws)]))
        return(avg)
    else:
        avg = np.mean(np.array(coverage[index:lastpos] + coverage[0:int(ws-diff)]))
        return(avg)

def getderivative(index, coverage, lastpos):
    index = int(index)
    if index == lastpos:
        deriv = coverage[lastpos] - coverage[0]
    else:
        deriv = coverage[index] - coverage[int(index + 1)]
    return(deriv)


# define thw Window size (ws)
ws = 50

# Read the coverage data
file = 'mt_CPCT02350007T_coverage.txt'
df = pd.read_csv(file, header = None, sep ='\t')
df.rename(columns=dict(zip([0,1,2],['CHROM','POS', 'coverage'])), inplace = True)


# Calculate the average coverage window size
df['avgws'] = df.apply(lambda x: getwsavg(int(x['POS']), list(df['coverage']), int(df['POS'].tail(1).values[0])), axis = 1)
# Calculate the coverage difference between two positions (the derivative)
df['deriv'] = df.apply(lambda x: getderivative(x['POS'], list(df['coverage']), df['POS'].tail(1).values[0]), axis = 1)
# Fid out those position passing a given threshold
thr = 15
df.loc[df['deriv'].abs() > thr]
