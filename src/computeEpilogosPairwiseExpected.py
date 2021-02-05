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

def main(filename1, filename2, numStates, saliency, outputDirPath, fileTag, numProcesses):
    tTotal = time.time()

    file1Path = Path(filename1)
    file2Path = Path(filename2)
    outputDirPath = Path(outputDirPath)

    ###################################################
    ##
    ##   chrName is hacky, find a better way later
    ##
    ####################################################
    chrName  = file1Path.name.split("_")[-1].split(".")[0]

    # Get the genome file
    genomeFileList = list(file1Path.parents[0].glob("*.genome"))
    if len(genomeFileList) > 1:
        raise IOError("Too many '.genome' files provided")
    elif len(genomeFileList) < 1:
        raise IOError("No '.genome' file provided")

    # Read line by line until we are on correct chromosome, then read out number of bp
    for genomeFile in genomeFileList:
        with open(genomeFile, "r") as gf:
            line = gf.readline()
            while chrName not in line:
                line = gf.readline()
            basePairs = int(line.split()[1])

    totalRows = math.ceil(basePairs / 200)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count()

    # Split the rows up according to the number of cores we have available
    rowList = []
    for i in range(numProcesses):
        rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
        rowList.append(rowsToCalculate)

    if saliency == 1:
        s1Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, chrName, numProcesses)
    elif saliency == 2:
        s2Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, chrName, numProcesses)
    elif saliency == 3:
        raise ValueError("A saliency metric of 3 is unsupported for pairwise comparisons")
    else:
        raise ValueError("Inputed saliency value not supported")

    print("Total Time:", time.time() - tTotal)


# Function that deploys the processes used to calculate the expected frequencies for the s1 metric. Also calls function to store expected frequency
def s1Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)
    
    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s1Calc, zip(itertools.repeat(file1Path), itertools.repeat(file2Path), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis=0) / totalRows

    storeExpArray(expFreqArr, outputDirPath, fileTag, chrName)

# Function that reads in data and calculates the expected frequencies for the S1 metric over a chunk of the data
def s1Calc(file1Path, file2Path, rowsToCalculate, numStates):
    # Read in the data
    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead1)

    print("Reading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead2)

    # Converting to a np array for faster functions
    # Also combine the arrays to calculate a combined background
    print("Converting to numpy arrays and combining...")
    tConvert = time.time()
    dataArr = np.concatenate((file1DF.iloc[:,3:].to_numpy(dtype=int) - 1, file2DF.iloc[:,3:].to_numpy(dtype=int) - 1), axis=1)
    print("    Time: ", time.time() - tConvert)

    numCols = dataArr.shape[1]

    expFreqArr = np.zeros(numStates, dtype=np.float32)

    # Simply count all states across out our subset of data
    uniqueStates, stateCounts = np.unique(dataArr, return_counts=True)
    for i, state in enumerate(uniqueStates):
        expFreqArr[state] += stateCounts[i] / numCols

    return expFreqArr

# Function that deploys the processes used to calculate the expected frequencies for the s2 metric. Also calls function to store expected frequency
def s2Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    if (sys.version_info < (3, 8)):
        print("\nFor maximum efficiency please update python to version 3.8 or later")
        print("NOTE: The code will still run in a lower version, but will be slightly slower\n")

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s2Calc, zip(itertools.repeat(file1Path), itertools.repeat(file2Path), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis = 0) / totalRows

    storeExpArray(expFreqArr, outputDirPath, fileTag, chrName)


# Function that reads in data and calculates the expected frequencies for the S2 metric over a chunk of the data
def s2Calc(file1Path, file2Path, rowsToCalculate, numStates):
    # Read in the data
    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead1)

    print("Reading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead2)

    # Converting to a np array for faster functions
    # Also combine the arrays to calculate a combined background
    print("Converting to numpy arrays and combining...")
    tConvert = time.time()
    dataArr = np.concatenate((file1DF.iloc[:,3:].to_numpy(dtype=int) - 1, file2DF.iloc[:,3:].to_numpy(dtype=int) - 1), axis=1)
    print("    Time: ", time.time() - tConvert)    

    multiprocessRows, numCols = dataArr.shape

    expFreqArr = np.zeros((numStates, numStates), dtype=np.float32)

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    if (sys.version_info < (3, 8)):
        combinations = ncr(numCols, 2)
        for row in range(multiprocessRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i, state1 in enumerate(uniqueStates):
                for j, state2 in enumerate(uniqueStates):
                    if state1 == state2:
                        expFreqArr[state1, state2] += ncr(stateCounts[i], 2) / combinations
                    else: # state1 > state2 or state1 < state2
                        expFreqArr[state1, state2] += stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
    else:
        combinations = math.comb(numCols, 2)
        for row in range(multiprocessRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True) 
            for i, state1 in enumerate(uniqueStates):
                for j, state2 in enumerate(uniqueStates):
                    if state1 == state2:
                        expFreqArr[state1, state2] += math.comb(stateCounts[i], 2) / combinations
                    else: # state1 > state2 or state1 < state2
                        expFreqArr[state1, state2] += stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix

    return expFreqArr

# Helper to store the expected frequency arrays
def storeExpArray(expFreqArr, outputDirPath, fileTag, chrName):
    # Creating a file path
    expFreqFilename = "temp_exp_freq_{}_{}.npy".format(fileTag, chrName)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)

# Helper to calculate combinations
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer // denom

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], int(sys.argv[7]))
