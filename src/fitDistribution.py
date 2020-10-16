import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import math
import scipy.stats as st
import statsmodels as sm
import warnings
import time

def main(file1, file2, filterBool, distributionNumber, binEnd, outputDir):
    if filterBool == "ERROR: INVALID BOOL SUBMITTED":
        print("ERROR: INVALID BOOL SUBMITTED")
        return

    tTotal = time.time()

    distributions = [st.betaprime, st.halfgennorm, st.pareto, st.lomax, st.genpareto, st.gamma, 
                    st.genexpon, st.expon, st.mielke, st.exponweib, st.loglaplace, st.chi, st.chi2,
                    st.nakagami, st.burr, st.ncx2, st.pearson3]

    distribution = distributions[distributionNumber]

    file1Path = Path(file1)
    file2Path = Path(file2)

    # Read in the data
    file1DF = pd.read_table(file1Path, header=None, sep="\s+")
    file2DF = pd.read_table(file2Path, header=None, sep="\s+")

    # Converting to a np array for faster functions later
    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=float)
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=float)

    plt.rcParams['agg.path.chunksize'] = 10000
    matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
    matplotlib.style.use('ggplot')

    file1SumArr = np.sum(file1Arr, axis=1)
    file2SumArr = np.sum(file2Arr, axis=1)

    distances = np.sum(np.square(file1Arr - file2Arr), axis=1)

    if filterBool:
        idx = [i for i in range(file1SumArr.shape[0]) if file1SumArr[i] > 1 or file2SumArr[i] > 1]
        data = pd.Series(distances[idx])
    else:
        data = pd.Series(distances)

    if binEnd == "Max" or binEnd == "max":
        y, x = np.histogram(data, bins=100, range=(0, np.amax(data)), density=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0
    else:
        y, x = np.histogram(data, bins=100, range=(0, float(binEnd)), density=True)
        x = (x + np.roll(x, -1))[:-1] / 2.0

    # Fit the data
    params = distribution.fit(data)

    # Separate parts of parameters
    distArgs = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Calculate SSE
    pdf = distribution.pdf(x, loc=loc, scale=scale, *distArgs)
    sse = np.sum(np.power(y - pdf, 2.0))

    distName = distribution.name

    param_names = (distribution.shapes + ', loc, scale').split(', ') if distribution.shapes else ['loc', 'scale']
    param_str = ', '.join(['{}={:0.2f}'.format(k,v) for k,v in zip(param_names, params)])
    dist_str = '{}({})'.format(distName, param_str)

    print()
    print("File 1:", file1)
    print("File 2:", file2)
    print(dist_str)
    print("SSE:", sse)

    allSSEPath = Path(outputDir) / "allSSE.txt"
    with open(allSSEPath, 'a') as allSSE:
        allSSE.write("{}\t{}\n".format(distName, sse))

    allParamsPath = Path(outputDir) / "allParams.txt"
    with open(allParamsPath, 'a') as allParams:
        allParams.write("{}\n".format(dist_str))

    print()
    print("    Time Elapsed:", time.time() - tTotal)

def strToBool(string):
    if string in ["True", "true", "T", "t", "y", "Y", "yes", "Yes"]:
        return True
    elif string in ["False", "false", "F", "f", "y", "Y", "yes", "Yes"]:
        return False
    else:
        return "ERROR: INVALID BOOL SUBMITTED"

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], strToBool(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6])