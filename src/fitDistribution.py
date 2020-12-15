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

def main(file1, file2, observationFile, filterBool, fullBool, distanceMetric, distributionNumber, numStates, outputDir):
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
        print("Fitting to Positive and Negative Distributions")
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

    if distanceMetric == "sqrt":
        distances = np.sqrt(np.sum(np.square(file1Arr - file2Arr), axis=1)) * observationArr[:,2]
    elif distanceMetric == "nonsqrt":
        distances = np.sum(np.square(file1Arr - file2Arr), axis=1) * observationArr[:,2]
    elif distanceMetric == "cosine":
        distances = (np.einsum('ij, ij->i', file1Arr, file2Arr) / (np.linalg.norm(file1Arr, axis=1) * np.linalg.norm(file2Arr, axis=1))) * observationArr[:,2]

    if fullBool:
        fittingToFull(distances, filterBool, file1Arr, file2Arr, outputDir, distribution, file1, file2)
    else:
        fittingPosNeg(distances, filterBool, file1Arr, file2Arr, outputDir, distribution, file1, file2)

    print()
    print("    Time Elapsed:", time.time() - tTotal)

def fittingToFull(distances, filterBool, file1Arr, file2Arr, outputDir, distribution, file1, file2):
    if filterBool:
        quiescentVal1 = round(file1Arr[0][-1], 5)
        quiescentVal2 = round(file2Arr[0][-1], 5)
        idx = [i for i in range(file1Arr.shape[0]) if round(file1Arr[i][-1], 5) != quiescentVal1 or round(file2Arr[i][-1], 5) != quiescentVal2]
        data = pd.Series(distances[idx])
    else:
        data = pd.Series(distances)

    y, x = np.histogram(data.values, bins=100, range=(np.amin(data), np.amax(data)), density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    params, sse, mle = fit(distribution, data, x, y)

    writeOut(outputDir, "", file1, file2, distribution, params, sse, mle)

def fittingPosNeg(distances, filterBool, file1Arr, file2Arr, outputDir, distribution, file1, file2):
    if filterBool:
        quiescentVal1 = round(file1Arr[0][-1], 5)
        quiescentVal2 = round(file2Arr[0][-1], 5)
        idx = [i for i in range(file1Arr.shape[0]) if round(file1Arr[i][-1], 5) != quiescentVal1 or round(file2Arr[i][-1], 5) != quiescentVal2]
        
        distancesPositive = distances[idx][np.where(distances[idx] >= 0)[0]]
        distancesNegative = -distances[idx][np.where(distances[idx] <= 0)[0]]
        dataPositive = pd.Series(distancesPositive)
        dataNegative = pd.Series(distancesNegative)

    else:
        distancesPositive = distances[np.where(distances >= 0)[0]]
        distancesNegative = -distances[np.where(distances <= 0)[0]]
        dataPositive = pd.Series(distancesPositive)
        dataNegative = pd.Series(distancesNegative)

    #  Positive Side Fitting
    yPos, xPos = np.histogram(dataPositive.values, bins=100, range=(0, np.amax(dataPositive)), density=True)
    xPos = (xPos + np.roll(xPos, -1))[:-1] / 2.0
    params, sse, mle = fit(distribution, dataPositive, xPos, yPos)
    writeOut(outputDir, "Positive", file1, file2, distribution, params, sse, mle)

    # Negative Side Fitting
    yNeg, xNeg = np.histogram(dataNegative.values, bins=100, range=(0, np.amax(dataNegative)), density=True)
    xNeg = (xNeg + np.roll(xNeg, -1))[:-1] / 2.0
    params, sse, mle = fit(distribution, dataNegative, xNeg, yNeg)
    writeOut(outputDir, "Negative", file1, file2, distribution, params, sse, mle)


def writeOut(outputDir, sign, file1, file2, distribution, params, sse, mle):
    distArgs = params[:-2]
    loc = params[-2]
    scale = params[-1]
    distName = distribution.name

    param_names = (distribution.shapes + ', loc, scale').split(', ') if distribution.shapes else ['loc', 'scale']
    param_str = ', '.join(['{}={:0.5f}'.format(k,v) for k,v in zip(param_names, params)])
    dist_str = '{}({})'.format(distName, param_str)
    dist_str_comma = ", ".join(str(v) for v in params)

    print()
    print("File 1:", file1)
    print("File 2:", file2)
    print(sign)
    print("Dist Args:", distArgs)
    print("loc:", loc)
    print("scale:", scale)
    print(dist_str)
    print(dist_str_comma)
    print("SSE:", sse)
    print("MLE:", mle)

    allSSEPath = Path(outputDir) / "all{}SSE.txt".format(sign)
    with open(allSSEPath, 'a') as allSSE:
        allSSE.write("{}    {}\n".format(distName, sse))

    allParamsPath = Path(outputDir) / "all{}Params.txt".format(sign)
    with open(allParamsPath, 'a') as allParams:
        allParams.write("{}\n".format(dist_str))

    allMLEPath = Path(outputDir) / "all{}MLE.txt".format(sign)
    with open(allMLEPath, 'a') as allMLE:
        allMLE.write("{}    {}\n".format(distName, mle))

    allParamsCommaPath = Path(outputDir) / "all{}ParamsComma.txt".format(sign)
    with open(allParamsCommaPath, 'a') as allParamsComma:
        allParamsComma.write("{}: {}\n".format(distName, dist_str_comma))


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
    main(sys.argv[1], sys.argv[2], sys.argv[3], strToBool(sys.argv[4]), strToBool(sys.argv[5]), sys.argv[6], int(sys.argv[7]), int(sys.argv[8]), sys.argv[9])