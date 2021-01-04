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

def main(filename, filename2, numStates, saliency, outputDirPath, fileTag):
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirPath)
    if filename2 == "null":
        # Reading in the data
        print("\nReading data from file...")
        tRead = time.time()
        dataDF = pd.read_table(dataFilePath, header=None, sep="\t")
        print("    Time: ", time.time() - tRead)

        # Converting to a np array for faster functions later
        print("Converting to numpy array...")
        tConvert = time.time()
        dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
        print("    Time: ", time.time() - tConvert)
    else:
        pairwisePath = Path(filename2)
        # Read in the data
        print("\nReading data from file 1...")
        tRead1 = time.time()
        file1DF = pd.read_table(dataFilePath, header=None, sep="\t")
        print("    Time: ", time.time() - tRead1)

        print("\nReading data from file 2...")
        tRead2 = time.time()
        file2DF = pd.read_table(pairwisePath, header=None, sep="\t")
        print("    Time: ", time.time() - tRead2)

        # Converting to a np array for faster functions later
        print("Converting to numpy arrays...")
        tConvert = time.time()
        file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int)
        file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int)
        locationArr = file1DF.iloc[:,0:3].to_numpy(dtype=str)
        print("    Time: ", time.time() - tConvert)

        # Combining the arrays for per row shuffling
        dataArr = np.concatenate((file1Arr, file2Arr), axis=1)
        dataDF = pd.concat((file1DF, file2DF.iloc[:,3:]), axis=1, ignore_index=True)

        print("File 1 DF Shape:", file1DF.shape)
        print("File 2 DF Shape:", file2DF.shape)
        print("Combined DF Shape:", dataDF.shape)
        print("File 1 Arr Shape:", file1Arr.shape)
        print("File 2 Arr Shape:", file2Arr.shape)
        print("Combined Arr Shape:", dataArr.shape)

    if saliency == 1:
        s1Exp(dataDF, dataArr, numStates, outputDirPath, fileTag)
    elif saliency == 2:
        s2Exp(dataDF, dataArr, numStates, outputDirPath, fileTag)
    elif saliency == 3:
        s3Exp(dataDF, dataArr, numStates, outputDirPath, fileTag)
    else:
        print("Inputed saliency value not supported")
        return

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

    storeExpArray(dataDF, expFreqArr, outputDirPath, fileTag)

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

    storeExpArray(dataDF, expFreqArr, outputDirPath, fileTag)

# Function that calculates the expected frequencies for the S3 metric
def s3Exp(dataDF, dataArr, numStates, outputDirPath, fileTag):
    numRows, numCols = dataArr.shape
    numProcesses = multiprocessing.cpu_count()

    # Use multiprocessing to speed up expected frequency calculation time
    # Calculate expected frequencies

    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2))).T

    # Initializing needed variables
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))
    expQueue = multiprocessing.Queue()
    expProcesses = []

    print("Num Processes:", numProcesses)

    # Creating the expected frequency processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRows // numProcesses, (i+1) * numRows // numProcesses)
        p = multiprocessing.Process(target=s3ExpMulti, args=(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, expQueue))
        expProcesses.append(p)
        p.start()

    print("Queue Size:", expQueue.qsize())

    # Combine all the calculated expvalue arrays into one
    for process in expProcesses:
        expFreqArr += expQueue.get()

    # Shut down all the processes
    for process in expProcesses:
        process.join()

    # Normalize the array
    expFreqArr /= numRows * numCols * (numCols - 1)

    storeExpArray(dataDF, expFreqArr, outputDirPath, fileTag)

# Helper function for the multiprocessing
def s3ExpMulti(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, queue):
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))
    print("Rows in Multiprocess:", rowsToCalculate)
    for row in rowsToCalculate:
        expFreqArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]] += np.ones(basePermutationArr.shape[1])
    queue.put(expFreqArr)

# Helper to store the expected frequency arrays
def storeExpArray(dataDF, expFreqArr, outputDirPath, fileTag):
    # Creating a file path
    locationTag = "{}_{}_{}".format(dataDF.iloc[0, 0], dataDF.iloc[0,1], dataDF.iloc[0,2])
    expFreqFilename = "temp_exp_freq_{}_{}.npy".format(fileTag, locationTag)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)

# Helper to calculate combinations
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer // denom

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6])