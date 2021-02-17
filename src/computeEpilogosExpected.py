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
from contextlib import closing
import gzip

def main(filename, numStates, saliency, outputDirPath, fileTag, numProcesses, verbose):
    if verbose: tTotal = time.time()

    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirPath)

    filename = dataFilePath.name.split(".")[0]

    if not verbose: print("    {}\t".format(filename), end="", flush=True)

    if dataFilePath.name.endswith("gz"):
        with gzip.open(dataFilePath, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))
    else:
        with open(dataFilePath, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count()

    # Split the rows up according to the number of cores we have available
    rowList = []
    for i in range(numProcesses):
        rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
        rowList.append(rowsToCalculate)

    determineSaliency(saliency, dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose)

    print("Total Time:", time.time() - tTotal, flush=True) if verbose else print("\t[Done]", flush=True)


def determineSaliency(saliency, dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose):
    if saliency == 1:
        s1Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose)
    elif saliency == 2:
        s2Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose)
    elif saliency == 3:
        s3Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose)
    else:
        raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")


# Function that calculates the expected frequencies for the S1 metric
def s1Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses)
    
    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s1Calc, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(numStates), itertools.repeat(verbose)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis=0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, filename)


# Function that reads in data and calculates the expected frequencies for the S2 metric over a chunk of the data
def s1Calc(dataFilePath, rowsToCalculate, numStates, verbose):
    # Reading in the data
    if verbose and rowsToCalculate[0] == 0: print("Reading data from {}...".format(dataFilePath), flush=True); tRead = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    dataArr = pd.read_table(dataFilePath, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1 
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time.time() - tRead, flush=True)

    expFreqArr = np.zeros(numStates, dtype=np.int32)

    # Simply count all states across out our subset of data
    uniqueStates, stateCounts = np.unique(dataArr, return_counts=True)
    
    if verbose and rowsToCalculate[0] == 0:
        print("Calculating expected frequency...", flush=True)
        tExp = time.time()
        percentDone = 0
    printCheckmarks = [int(len(uniqueStates) * float(i / 10)) for i in range(1, 10)]

    if rowsToCalculate[0] == 0: print("Unique States:", uniqueStates)
    if rowsToCalculate[0] == 0: print("printCheckmarks:", printCheckmarks)

    for i, state in enumerate(uniqueStates):

        if verbose and rowsToCalculate[0] == 0 and i in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and i in printCheckmarks: print(".", end="", flush=True)

        expFreqArr[state] += stateCounts[i]

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time.time() - tExp, flush=True)

    return expFreqArr

# Function that deploys the processes used to calculate the expected frequencies for the s2 metric. Also calls function to store expected frequency
def s2Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s2Calc, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(numStates), itertools.repeat(verbose)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis = 0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, filename)


# Function that reads in data and calculates the expected frequencies for the S2 metric over a chunk of the data
def s2Calc(dataFilePath, rowsToCalculate, numStates, verbose):
    # Reading in the data
    if verbose and rowsToCalculate[0] == 0: print("Reading data from {}...".format(dataFilePath), flush=True); tRead = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    dataArr = pd.read_table(dataFilePath, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time.time() - tRead, flush=True)

    multiprocessRows = dataArr.shape[0]

    expFreqArr = np.zeros((numStates, numStates), dtype=np.int32)

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tExp = time.time()
        percentDone = 0
    printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]

    # SumOverRows: Within a row, how many ways can you choose x and y to be together (will normalize later)
    # Can choose x and y to be together n*m ways if n != m and n(n-1) ways if n == m (where n and m are the number of times that x and y show up respectively)
    for row in range(multiprocessRows):

        if verbose and rowsToCalculate[0] == 0 and row in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and row in printCheckmarks: print(".", end="", flush=True)

        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True) 
        for i, state1 in enumerate(uniqueStates):
            for j, state2 in enumerate(uniqueStates):
                if state1 == state2:
                    expFreqArr[state1, state2] += stateCounts[i] * (stateCounts[i] - 1)
                else: # state1 > state2 or state1 < state2
                    expFreqArr[state1, state2] += stateCounts[i] * stateCounts[j]

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time.time() - tExp, flush=True)

    return expFreqArr

# Function that deploys the processes used to calculate the expected frequencies for the s3 metric. Also call function to store expected frequency
def s3Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s3Calc, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(numStates), itertools.repeat(verbose)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing
    expFreqArr = np.sum(results, axis = 0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, filename)

# Function that calculates the expected frequencies for the S3 metric over a chunk of the data
def s3Calc(dataFilePath, rowsToCalculate, numStates, verbose):
    # Reading in the data
    if verbose and rowsToCalculate[0] == 0: print("Reading data from {}...".format(dataFilePath), flush=True); tRead = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    dataArr = pd.read_table(dataFilePath, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time.time() - tRead, flush=True)

    multiprocessRows, numCols = dataArr.shape

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2))).T

    expFreqArr = np.zeros((numCols, numCols, numStates, numStates), dtype=np.int32)
    
    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tExp = time.time()
        percentDone = 0
    printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]

    # We tally a one for all the state/column combinations we observe (e.g. for state 18 in column 2 and state 15 in column 6 we would add one to index [5, 1, 17, 14])
    for row in range(multiprocessRows):

        if verbose and rowsToCalculate[0] == 0 and row in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and row in printCheckmarks: print(".", end="", flush=True)

        expFreqArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]] += 1

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time.time() - tExp, flush=True)

    return expFreqArr

# Helper to store the expected frequency arrays
def storeExpArray(expFreqArr, outputDirPath, fileTag, filename):
    # Creating a file path
    expFreqFilename = "temp_exp_freq_{}_{}.npy".format(fileTag, filename)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)


# Helper for reading number of lines in input file
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

# Helper for slurm to send boolean values
def strToBool(string):
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        raise ValueError("Invalid boolean string")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], int(sys.argv[6]), strToBool(sys.argv[7]))