import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time
import operator as op
from functools import reduce
import multiprocessing
import itertools

def main(filename1, filename2, numStates, saliency, outputDirPath):
    tTotal = time.time()

    file1Path = Path(filename1)
    file2Path = Path(filename2)
    outputDirPath = Path(outputDirPath)

    # For distinguishing the chunkwise output files
    fileTag = "{}_{}_{}".format(file1Path.parent.name, file2Path.parent.name, file1Path.name.split(".")[0])

    # Read in the data
    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, header=None, sep="\t")
    print("    Time: ", time.time() - tRead1)

    print("Reading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, header=None, sep="\t")
    print("    Time: ", time.time() - tRead2)

    # Converting to a np array for faster functions later
    print("Converting to numpy arrays...")
    tConvert = time.time()
    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int) - 1
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int) - 1
    locationArr = file1DF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    # Combining the arrays to calculate a combined background
    dataArr = np.concatenate((file1Arr, file2Arr), axis=1)
    dataDF = pd.concat((file1DF, file2DF.iloc[:,3:]), axis=1, ignore_index=True)

    if saliency == 1:
        s1Exp(dataDF, dataArr, numStates, outputDirPath, fileTag)
    elif saliency == 2:
        s2Exp(dataDF, dataArr, numStates, outputDirPath, fileTag)
    elif saliency == 3:
        print("A saliency metric of 3 is unsupported for pairwise comparisons")
        return
        # s3Exp(dataDF, dataArr, numStates, outputDirPath, fileTag)
    else:
        print("Inputed saliency value not supported")
        return

    print("Total Time:", time.time() - tTotal)

# Function that calculates the expected frequencies for the S1 metric
def s1Exp(dataDF, dataArr, numStates, outputDirPath, fileTag):
    numRows, numCols = dataArr.shape

    # Calculate the expected frequencies of each state
    stateIndices = list(range(1, numStates + 1))
    expFreqSeries = pd.Series(np.zeros(numStates), index=stateIndices)
    dfSize = numRows * numCols
    for i in range(3, numCols + 3):
        stateCounts = dataDF[i].value_counts()
        for state, count in stateCounts.items():
            expFreqSeries.loc[state] += count / dfSize
    expFreqArr = expFreqSeries.to_numpy()

    storeExpArray(expFreqArr, outputDirPath, fileTag)

# Function that calculates the expected frequencies for the S2 metric
def s2Exp(dataDF, dataArr, numStates, outputDirPath, fileTag):
    numRows, numCols = dataArr.shape

    expFreqArr = np.zeros((numStates, numStates))

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    if sys.version_info < (3, 8):
        print("\nFor maximum efficiency please update python to version 3.8 or later")
        print("NOTE: The code will still run in a lower version, but will be slightly slower\n")
        combinations = ncr(numCols, 2)
        for row in range(numRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i in range(len(uniqueStates)):
                for j in range(len(uniqueStates)):
                    if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                        expFreqArr[uniqueStates[i], uniqueStates[j]] += stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                    elif uniqueStates[i] == uniqueStates[j]:
                        expFreqArr[uniqueStates[i], uniqueStates[j]] += ncr(stateCounts[i], 2) / combinations
    else:
        combinations = math.comb(numCols, 2)
        for row in range(numRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True) 
            for i in range(len(uniqueStates)):
                for j in range(len(uniqueStates)):
                    if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                        expFreqArr[uniqueStates[i], uniqueStates[j]] += stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                    elif uniqueStates[i] == uniqueStates[j]:
                        expFreqArr[uniqueStates[i], uniqueStates[j]] += math.comb(stateCounts[i], 2) / combinations

    expFreqArr = expFreqArr / numRows

    storeExpArray(expFreqArr, outputDirPath, fileTag)

# Helper to store the expected frequency arrays
def storeExpArray(expFreqArr, outputDirPath, fileTag):
    # Creating a file path
    expFreqFilename = "temp_exp_freq_{}.npy".format(fileTag)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)

# Helper to calculate combinations
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer // denom

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])