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

def main(filename, numStates, saliency, outputDirPath, fileTag, numProcesses):
    tTotal = time.time()

    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirPath)

    ###################################################
    ##
    ##   chrName is hacky, find a better way later
    ##
    ####################################################
    # chrName = dataFilePath.name.split("_")[-1].split(".")[0]

    # # Get the genome file
    # genomeFileList = list(dataFilePath.parents[0].glob("*.genome"))
    # if len(genomeFileList) > 1:
    #     raise IOError("Too many '.genome' files provided")
    # elif len(genomeFileList) < 1:
    #     raise IOError("No '.genome' file provided")

    # # Read line by line until we are on correct chromosome, then read out number of bp
    # for genomeFile in genomeFileList:
    #     with open(genomeFile, "r") as gf:
    #         line = gf.readline()
    #         while chrName not in line:
    #             line = gf.readline()
    #         basePairs = int(line.split()[1])

    # totalRows = math.ceil(basePairs / 200)

    # Chrname is actually filename
    chrName = dataFilePath.name.split(".")[0]

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

    determineSaliency(saliency, dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses)

    print("Total Time:", time.time() - tTotal)


def determineSaliency(saliency, dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses):
    if saliency == 1:
        s1Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses)
    elif saliency == 2:
        s2Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses)
    elif saliency == 3:
        s3Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses)
    else:
        print("Inputed saliency value not supported")
        return


# Function that calculates the expected frequencies for the S1 metric
def s1Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)
    
    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s1Calc, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis=0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, chrName)


# Function that reads in data and calculates the expected frequencies for the S2 metric over a chunk of the data
def s1Calc(dataFilePath, rowsToCalculate, numStates):
    # Reading in the data
    if rowsToCalculate[0] == 0:
        print("\nReading data from file...")
        tRead = time.time()
    dataDF = pd.read_table(dataFilePath, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions
    if rowsToCalculate[0] == 0:
        print("Converting to numpy array...")
        tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tConvert)

    expFreqArr = np.zeros(numStates, dtype=np.int32)

    # Simply count all states across out our subset of data
    uniqueStates, stateCounts = np.unique(dataArr, return_counts=True)
    
    if rowsToCalculate[0] == 0:
        print("Calculating expected frequency...")
        tExp = time.time()
        printCheckmarks = [int(len(uniqueStates) * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    for i, state in enumerate(uniqueStates):

        if rowsToCalculate[0] == 0 and i in printCheckmarks:
            percentDone += 10
            print("{}% Completed".format(percentDone))

        expFreqArr[state] += stateCounts[i]

    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tExp)

    return expFreqArr

# Function that deploys the processes used to calculate the expected frequencies for the s2 metric. Also calls function to store expected frequency
def s2Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s2Calc, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis = 0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, chrName)


# Function that reads in data and calculates the expected frequencies for the S2 metric over a chunk of the data
def s2Calc(dataFilePath, rowsToCalculate, numStates):
    # Reading in the data
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    print("    Time: ", time.time() - tConvert)

    multiprocessRows = dataArr.shape[0]

    expFreqArr = np.zeros((numStates, numStates), dtype=np.int32)

    if rowsToCalculate[0] == 0:
        print("\nCalculating Scores...")
        tExp = time.time()
        printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]
        percentDone = 0
    
    # SumOverRows: Within a row, how many ways can you choose x and y to be together (will normalize later)
    # Can choose x and y to be together n*m ways if n != m and n(n-1) ways if n == m (where n and m are the number of times that x and y show up respectively)
    for row in range(multiprocessRows):

        if rowsToCalculate[0] == 0 and row in printCheckmarks:
            percentDone += 10
            print("{}% Completed".format(percentDone))

        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True) 
        for i, state1 in enumerate(uniqueStates):
            for j, state2 in enumerate(uniqueStates):
                if state1 == state2:
                    expFreqArr[state1, state2] += stateCounts[i] * (stateCounts[i] - 1)
                else: # state1 > state2 or state1 < state2
                    expFreqArr[state1, state2] += stateCounts[i] * stateCounts[j]

    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tExp)

    return expFreqArr

# Function that deploys the processes used to calculate the expected frequencies for the s3 metric. Also call function to store expected frequency
def s3Exp(dataFilePath, rowList, numStates, outputDirPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s3Calc, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing
    expFreqArr = np.sum(results, axis = 0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, chrName)

# Function that calculates the expected frequencies for the S3 metric over a chunk of the data
def s3Calc(dataFilePath, rowsToCalculate, numStates):
    # Reading in the data
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    print("    Time: ", time.time() - tConvert)

    multiprocessRows, numCols = dataArr.shape

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2))).T

    expFreqArr = np.zeros((numCols, numCols, numStates, numStates), dtype=np.int32)
    
    if rowsToCalculate[0] == 0:
        print("\nCalculating Scores...")
        tExp = time.time()
        printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # We tally a one for all the state/column combinations we observe (e.g. for state 18 in column 2 and state 15 in column 6 we would add one to index [5, 1, 17, 14])
    for row in range(multiprocessRows):

        if rowsToCalculate[0] == 0 and row in printCheckmarks:
            percentDone += 10
            print("{}% Completed".format(percentDone))

        expFreqArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]] += 1

    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tExp)

    return expFreqArr

# Helper to store the expected frequency arrays
def storeExpArray(expFreqArr, outputDirPath, fileTag, chrName):
    # Creating a file path
    expFreqFilename = "temp_exp_freq_{}_{}.npy".format(fileTag, chrName)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)


# Helper for reading number of lines in input file
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], int(sys.argv[6]))