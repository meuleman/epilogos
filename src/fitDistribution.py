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

def main(file1, file2, observationFile, filterBool, fullBool, distributionNumber, binEnd, numStates, outputDir):
    if filterBool == "ERROR: INVALID BOOL SUBMITTED":
        print("ERROR: INVALID BOOL SUBMITTED")
        return

    tTotal = time.time()

    if fullBool:
        print("Fitting to Full Distribution")
        distributions = [st.cauchy, st.exponnorm, st.t, st.genlogistic, st.gennorm, st.gumbel_r, 
                    st.gumbel_l, st.gausshyper, st.hypsecant, st.johnsonsu, st.loglaplace, 
                    st.laplace, st.logistic, st.foldnorm, st.norm, st.norminvgauss, 
                    st.powerlognorm, st.powernorm, st.lognorm, st.skewnorm]
    else:
        print("Fitting to Full Distribution")
        distributions = [st.betaprime, st.halfgennorm, st.pareto, st.lomax, st.genpareto, st.gamma, 
                        st.genexpon, st.expon, st.mielke, st.exponweib, st.loglaplace, st.chi, st.chi2,
                        st.nakagami, st.burr, st.ncx2, st.pearson3]


    distribution = distributions[distributionNumber]

    file1Path = Path(file1)
    file2Path = Path(file2)
    observationPath = Path(observationFile)

    names = ["chr", "binStart", "binEnd"]
    for i in range(1, numStates + 1):
        names.append("s{}".format(i))
        
    chrOrder = []
    for i in range(1, 23):
        chrOrder.append("chr{}".format(i))
    chrOrder.append("chrX")

    # Read in the data
    file1DF = pd.read_table(file1Path, header=None, sep="\s+", names=names)
    file2DF = pd.read_table(file2Path, header=None, sep="\s+", names=names)
    observationDF = pd.read_table(observationPath, header=None, sep="\s+", names=["chr", "binStart", "binEnd", "maxDiffState", "maxDiff", "sign", "totalDiff"])

    file1DF["chr"] = pd.Categorical(file1DF["chr"], categories=chrOrder, ordered=True)
    file2DF["chr"] = pd.Categorical(file2DF["chr"], categories=chrOrder, ordered=True)
    observationDF["chr"] = pd.Categorical(observationDF["chr"], categories=chrOrder, ordered=True)

    file1DF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)
    file2DF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)
    observationDF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)


    # Converting to a np array for faster functions later
    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=float)
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=float)
    observationArr = observationDF.iloc[:,3:].to_numpy(dtype=float)

    plt.rcParams['agg.path.chunksize'] = 10000

    if fullBool:
        distances = np.sqrt(np.sum(np.square(file1Arr - file2Arr), axis=1)) * observationArr[:,2]
    else:
        distances = np.sqrt(np.sum(np.square(file1Arr - file2Arr), axis=1))

    if filterBool:
        quiescentVal1 = round(file1Arr[0][-1], 5)
        quiescentVal2 = round(file2Arr[0][-1], 5)
        idx = [i for i in range(file1Arr.shape[0]) if round(file1Arr[i][-1], 5) != quiescentVal1 or round(file2Arr[i][-1], 5) != quiescentVal2]
        # data = pd.Series(distances[idx])
        
        distancesPositive = distances[idx][np.where(distances[idx] >= 0)[0]]
        distancesNegative = -distances[idx][np.where(distances[idx] <= 0)[0]]
        dataPositive = pd.Series(distancesPositive)
        dataNegative = pd.Series(distancesNegative)
    else:
        data = pd.Series(distances)

    if fullBool:
        if binEnd == "Max" or binEnd == "max":
            # y, x = np.histogram(data.values, bins=100, range=(np.amin(data), np.amax(data)), density=True)
            # x = (x + np.roll(x, -1))[:-1] / 2.0

            yPos, xPos = np.histogram(dataPositive.values, bins=100, range=(np.amin(dataPositive), np.amax(dataPositive)), density=True)
            xPos = (xPos + np.roll(xPos, -1))[:-1] / 2.0

            yNeg, xNeg = np.histogram(dataNegative.values, bins=100, range=(np.amin(dataNegative), np.amax(dataNegative)), density=True)
            xNeg = (xNeg + np.roll(xNeg, -1))[:-1] / 2.0
        else:
            # y, x = np.histogram(data.values, bins=100, range=(-float(binEnd), float(binEnd)), density=True)
            # x = (x + np.roll(x, -1))[:-1] / 2.0

            yPos, xPos = np.histogram(dataPositive.values, bins=100, range=(-float(binEnd), float(binEnd)), density=True)
            xPos = (xPos + np.roll(xPos, -1))[:-1] / 2.0

            yNeg, xNeg = np.histogram(dataNegative.values, bins=100, range=(-float(binEnd), float(binEnd)), density=True)
            xNeg = (xNeg + np.roll(xNeg, -1))[:-1] / 2.0
    else:
        if binEnd == "Max" or binEnd == "max":
            # y, x = np.histogram(data.values, bins=100, range=(0, np.amax(data)), density=True)
            # x = (x + np.roll(x, -1))[:-1] / 2.0

            yPos, xPos = np.histogram(dataPositive.values, bins=100, range=(0, np.amax(dataPositive)), density=True)
            xPos = (xPos + np.roll(xPos, -1))[:-1] / 2.0

            yNeg, xNeg = np.histogram(dataNegative.values, bins=100, range=(0, np.amax(dataNegative)), density=True)
            xNeg = (xNeg + np.roll(xNeg, -1))[:-1] / 2.0
        else:
            # y, x = np.histogram(data.values, bins=100, range=(0, float(binEnd)), density=True)
            # x = (x + np.roll(x, -1))[:-1] / 2.0

            yPos, xPos = np.histogram(dataPositive.values, bins=100, range=(0, float(binEnd)), density=True)
            xPos = (xPos + np.roll(xPos, -1))[:-1] / 2.0

            yNeg, xNeg = np.histogram(dataNegative.values, bins=100, range=(0, float(binEnd)), density=True)
            xNeg = (xNeg + np.roll(xNeg, -1))[:-1] / 2.0


    #####################################
    #                                   #
    #      Positive Side Fitting        #
    #                                   #
    #####################################
    params, sse, mle = fit(distribution, dataPositive, xPos, yPos)

    distArgs = params[:-2]
    loc = params[-2]
    scale = params[-1]
    distName = distribution.name

    param_names = (distribution.shapes + ', loc, scale').split(', ') if distribution.shapes else ['loc', 'scale']
    param_str = ', '.join(['{}={:0.5f}'.format(k,v) for k,v in zip(param_names, params)])
    dist_str = '{}({})'.format(distName, param_str)
    dist_str_comma = ", ".join("{:0.5f}".format(v) for v in params)

    print()
    print("File 1:", file1)
    print("File 2:", file2)
    print("POSITIVE")
    print("Dist Args:", distArgs)
    print("loc:", loc)
    print("scale:", scale)
    print(dist_str)
    print(dist_str_comma)
    print("SSE:", sse)
    print("MLE:", mle)

    allSSEPath = Path(outputDir) / "allPositiveSSE.txt"
    with open(allSSEPath, 'a') as allSSE:
        allSSE.write("{}    {}\n".format(distName, sse))

    allParamsPath = Path(outputDir) / "allPositiveParams.txt"
    with open(allParamsPath, 'a') as allParams:
        allParams.write("{}\n".format(dist_str))

    allMLEPath = Path(outputDir) / "allPositiveMLE.txt"
    with open(allMLEPath, 'a') as allMLE:
        allMLE.write("{}    {}\n".format(distName, mle))

    allParamsCommaPath = Path(outputDir) / "allPositiveParamsComma.txt"
    with open(allParamsCommaPath, 'a') as allParamsComma:
        allParamsComma.write("{}: {}\n".format(distName, dist_str_comma))


    #####################################
    #                                   #
    #      Negative Side Fitting        #
    #                                   #
    #####################################
    params, sse, mle = fit(distribution, dataNegative, xNeg, yNeg)

    distArgs = params[:-2]
    loc = params[-2]
    scale = params[-1]
    distName = distribution.name

    param_names = (distribution.shapes + ', loc, scale').split(', ') if distribution.shapes else ['loc', 'scale']
    param_str = ', '.join(['{}={:0.5f}'.format(k,v) for k,v in zip(param_names, params)])
    dist_str = '{}({})'.format(distName, param_str)
    dist_str_comma = ", ".join("{:0.5f}".format(v) for v in params)

    print()
    print("File 1:", file1)
    print("File 2:", file2)
    print("NEGATIVE")
    print("Dist Args:", distArgs)
    print("loc:", loc)
    print("scale:", scale)
    print(dist_str)
    print(dist_str_comma)
    print("SSE:", sse)
    print("MLE:", mle)

    allSSEPath = Path(outputDir) / "allNegativeSSE.txt"
    with open(allSSEPath, 'a') as allSSE:
        allSSE.write("{}    {}\n".format(distName, sse))

    allParamsPath = Path(outputDir) / "allNegativeParams.txt"
    with open(allParamsPath, 'a') as allParams:
        allParams.write("{}\n".format(dist_str))

    allMLEPath = Path(outputDir) / "allNegativeMLE.txt"
    with open(allMLEPath, 'a') as allMLE:
        allMLE.write("{}    {}\n".format(distName, mle))

    allParamsCommaPath = Path(outputDir) / "allNegativeParamsComma.txt"
    with open(allParamsCommaPath, 'a') as allParamsComma:
        allParamsComma.write("{}: {}\n".format(distName, dist_str_comma))


    print()
    print("    Time Elapsed:", time.time() - tTotal)

def fit(distribution, data, x, y):
    # ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Fit the data
        params = distribution.fit(data)

        # Separate parts of parameters
        distArgs = params[:-2]
        loc = params[-2]
        scale = params[-1]

        # Calculate SSE and MLE
        pdf = distribution.pdf(x, loc=loc, scale=scale, *distArgs)
        sse = np.sum(np.power(y - pdf, 2.0))
        mle = distribution.nnlf(params, data)

        return params, sse, mle


def strToBool(string):
    if string in ["True", "true", "T", "t", "y", "Y", "yes", "Yes"]:
        return True
    elif string in ["False", "false", "F", "f", "y", "Y", "yes", "Yes"]:
        return False
    else:
        return "ERROR: INVALID BOOL SUBMITTED"

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], strToBool(sys.argv[4]), strToBool(sys.argv[5]), int(sys.argv[6]), sys.argv[7], int(sys.argv[8]), sys.argv[9])