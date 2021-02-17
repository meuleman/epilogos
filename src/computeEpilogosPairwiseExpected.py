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
    # chrName  = file1Path.name.split("_")[-1].split(".")[0]

    # # Get the genome file
    # genomeFileList = list(file1Path.parents[0].glob("*.genome"))
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
    filename = file1Path.name.split(".")[0]

    if file1Path.name.endswith("gz"):
        with gzip.open(file1Path, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))
    else:
        with open(file1Path, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count()

    # Split the rows up according to the number of cores we have available
    rowList = []
    for i in range(numProcesses):
        rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
        rowList.append(rowsToCalculate)

    if saliency == 1:
        s1Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, filename, numProcesses)
    elif saliency == 2:
        s2Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, filename, numProcesses)
    elif saliency == 3:
        raise ValueError("A saliency metric of 3 is unsupported for pairwise comparisons")
    else:
        raise ValueError("Inputed saliency value not supported")

    print("Total Time:", time.time() - tTotal)


# Function that deploys the processes used to calculate the expected frequencies for the s1 metric. Also calls function to store expected frequency
def s1Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, filename, numProcesses):
    print("\nNumber of Processes:", numProcesses)
    
    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s1Calc, zip(itertools.repeat(file1Path), itertools.repeat(file2Path), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis=0) / totalRows

    storeExpArray(expFreqArr, outputDirPath, fileTag, filename)

# Function that reads in data and calculates the expected frequencies for the S1 metric over a chunk of the data
def s1Calc(file1Path, file2Path, rowsToCalculate, numStates):
    # Read in the data
    if rowsToCalculate[0] == 0:
        print("Reading data from file 1...")
        tRead1 = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file1Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file1Arr = pd.read_table(file1Path, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead1)

    if rowsToCalculate[0] == 0:
        print("Reading data from file 2...")
        tRead2 = time.time()
    cols = range(3, pd.read_table(file2Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file2Arr = pd.read_table(file2Path, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead2)

    # Combine the arrays to calculate a combined background
    if rowsToCalculate[0] == 0:
        print("Combining input matrices...")
        tConvert = time.time()
    dataArr = np.concatenate((file1Arr, file2Arr), axis=1)
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
            print("    {}% Completed".format(percentDone))

        expFreqArr[state] += stateCounts[i]

    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tExp)

    return expFreqArr

# Function that deploys the processes used to calculate the expected frequencies for the s2 metric. Also calls function to store expected frequency
def s2Exp(file1Path, file2Path, rowList, totalRows, numStates, outputDirPath, fileTag, filename, numProcesses):
    print("\nNumber of Processes:", numProcesses)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses)) as pool:
        results = pool.starmap(s2Calc, zip(itertools.repeat(file1Path), itertools.repeat(file2Path), rowList, itertools.repeat(numStates)))
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis = 0) / totalRows

    storeExpArray(expFreqArr, outputDirPath, fileTag, filename)


# Function that reads in data and calculates the expected frequencies for the S2 metric over a chunk of the data
def s2Calc(file1Path, file2Path, rowsToCalculate, numStates):
    # Read in the data
    if rowsToCalculate[0] == 0:
        print("Reading data from file 1...")
        tRead1 = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file1Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file1Arr = pd.read_table(file1Path, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead1)

    if rowsToCalculate[0] == 0:
        print("Reading data from file 2...")
        tRead2 = time.time()
    cols = range(3, pd.read_table(file2Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file2Arr = pd.read_table(file2Path, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead2)

    # Converting to a np array for faster functions
    # Also combine the arrays to calculate a combined background
    if rowsToCalculate[0] == 0:
        print("Combining input matrices...")
        tConvert = time.time()
    dataArr = np.concatenate((file1Arr, file2Arr), axis=1)
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tConvert)    

    multiprocessRows = dataArr.shape[0]

    expFreqArr = np.zeros((numStates, numStates), dtype=np.int32)

    if rowsToCalculate[0] == 0:
        print("Calculating Scores...")
        tExp = time.time()
        printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # SumOverRows: Within a row, how many ways can you choose x and y to be together (will normalize later)
    # Can choose x and y to be together n*m ways if n != m and n(n-1) ways if n == m (where n and m are the number of times that x and y show up respectively)
    for row in range(multiprocessRows):

        if rowsToCalculate[0] == 0 and row in printCheckmarks:
            percentDone += 10
            print("    {}% Completed".format(percentDone))

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

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], int(sys.argv[7]))
