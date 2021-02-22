import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time
import numpy.ma as ma
import operator as op
from functools import reduce
import multiprocessing
import itertools
import ctypes as c
from contextlib import closing
import gzip

def main(file, numStates, saliency, outputDirPath, expFreqPath, fileTag, numProcesses, verbose):
    # try:
    if verbose: tTotal = time.time()
    
    dataFilePath = Path(file)
    outputDirPath = Path(outputDirPath)

    filename = dataFilePath.name.split(".")[0]

    if not verbose: print("    {}\t".format(filename), end="", flush=True)

    try:
        if dataFilePath.name.endswith("gz"):
            with gzip.open(dataFilePath, "rb") as f:
                totalRows = sum(bl.count(b'\n') for bl in blocks(f))
        else:
            with open(dataFilePath, "rb") as f:
                totalRows = sum(bl.count(b'\n') for bl in blocks(f))
    except:
        print(sys.exc_info()[0], flush=True)
        print("IN line number determination", flush=True)

    try:
        # If user doesn't want to choose number of cores use as many as available
        if numProcesses == 0:
            numProcesses = multiprocessing.cpu_count()
    except:
        print(sys.exc_info()[0], flush=True)
        print("In numprocesses", flush=True)

    try:
        # Split the rows up according to the number of cores we have available
        rowList = []
        for i in range(numProcesses):
            rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
            rowList.append(rowsToCalculate)
    except:
        print(sys.exc_info()[0], flush=True)
        print("In rowList", flush=True)

    determineSaliency(saliency, dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose)

    print("Total Time:", time.time() - tTotal, flush=True) if verbose else print("\t[Done]", flush=True)
    # except OSError as err:
    #     if err.errno == 16:
    #         print("Warning: OSError 16 thrown. In testing we have found this does not effect the program output, but please check output file to be certain")
    #     else:
    #         print(err)
    # except:
    #     print(sys.exc_info()[0], flush=True)

def determineSaliency(saliency, dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose):
    try:
        if saliency == 1:
            s1Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose)
        elif saliency == 2:
            s2Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose)
        elif saliency == 3:
            s3Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose)
        else:
            raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")
    except:
        print(sys.exc_info()[0], flush=True)
        print("In determine saliency", flush=True)


# Helper for unflattening a shared array into a 2d numpy array
def sharedToNumpy(sharedArr, numRows, numStates):
    try:
        return np.frombuffer(sharedArr.get_obj(), dtype=np.float32).reshape((numRows, numStates))
    except:
        print(sys.exc_info()[0], flush=True)
        print("In shared to numpy", flush=True)

# initiliazer for multiprocessing
def _init(sharedArr_):
    try:
        global sharedArr
        sharedArr = sharedArr_
    except:
        print(sys.exc_info()[0], flush=True)
        print("In init", flush=True)

# Function that deploys the processes used to calculate the scores for the s1 metric. Also call function to store scores
def s1Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses)

    try:
        sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    except:
        print(sys.exc_info()[0], flush=True)
        print("In shared array", flush=True)

    # Start the processes
    try:
        with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
            pool.starmap(s1Score, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(expFreqPath), itertools.repeat(verbose)))
        pool.join()
    except:
        print(sys.exc_info()[0], flush=True)
        print("In multiprocessing", flush=True)

    try:
        chrName = pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").iloc[0, 0]
    except:
        print(sys.exc_info()[0], flush=True)
        print("in chrName determinations", flush=True)

    storeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputDirPath, fileTag, filename, chrName)


# Calculates the scores for the s1 metric over a given range of rows
def s1Score(dataFilePath, rowsToCalculate, expFreqPath, verbose):
    try:
        # Read in the data
        if verbose and rowsToCalculate[0] == 0: print("Reading data from file...", flush=True); tRead = time.time()
        # Dont want to read in locations
        cols = range(3, pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").shape[1])
        # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
        dataArr = pd.read_table(dataFilePath, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
        if verbose and rowsToCalculate[0] == 0: print("    Time: ", time.time() - tRead, flush=True)
    except:
        print(sys.exc_info()[0], flush=True)
        print("In reading data", flush=True)

    try:
        # Loading the expected frequency array
        expFreqArr = np.load(expFreqPath, allow_pickle=False)
    except:
        print(sys.exc_info()[0], flush=True)
        print("In reading in expFreqArr", flush=True)

    numCols = dataArr.shape[1]

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tScore = time.time()
        percentDone = 0
    printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]
    
    try:
        scoreArr = sharedToNumpy(*sharedArr)
        # Calculate the observed frequencies and final scores for the designated rows
        for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):
            
            if verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
            if not verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: print(".", end="", flush=True)

            uniqueStates, stateCounts = np.unique(dataArr[obsRow], return_counts=True)
            for i, state in enumerate(uniqueStates):
                # Function input is obsFreq and expFreq
                scoreArr[scoreRow, state] = klScore(stateCounts[i] / numCols, expFreqArr[state])
    except:
        print(sys.exc_info()[0], flush=True)
        print("In score calculation", flush=True)

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time.time() - tScore, flush=True)


# Function that deploys the processes used to calculate the scores for the s2 metric. Also call function to store scores
def s2Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
        pool.starmap(s2Score, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(expFreqPath), itertools.repeat(verbose)))
    pool.join()

    chrName = pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").iloc[0, 0]

    storeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputDirPath, fileTag, filename, chrName)


# Calculates the scores for the s2 metric over a given range of rows
def s2Score(dataFilePath, rowsToCalculate, expFreqPath, verbose):
    # Read in the data
    if verbose and rowsToCalculate[0] == 0: print("Reading data from file...", flush=True); tRead = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    dataArr = pd.read_table(dataFilePath, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time.time() - tRead, flush=True)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    numCols = dataArr.shape[1]
    numStates = sharedArr[2]

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tScore = time.time()
        percentDone = 0
    printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]

    # Calculte the scores and store them in the shared array
    scoreArr = sharedToNumpy(*sharedArr)
    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ordered ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1) ways if same (where n is the number of times that x/y shows up)
    rowObsArr = np.zeros((numStates, numStates))
    permutations = numCols * (numCols - 1)
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: print(".", end="", flush=True)

        uniqueStates, stateCounts = np.unique(dataArr[obsRow], return_counts=True)
        for i, state1 in enumerate(uniqueStates):
            for j, state2 in enumerate(uniqueStates):
                if state1 == state2:
                    rowObsArr[state1, state2] = stateCounts[i] * (stateCounts[i] - 1) / permutations # Equates to statecounts[i] permute 2 / permutations
                else: # state1 > state2 or state1 < state2
                    rowObsArr[state1, state2] = stateCounts[i] * stateCounts[j] / permutations
        
        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        scoreArr[scoreRow] = klScoreND(rowObsArr, expFreqArr).sum(axis=0)

        # Reset the array so it doesn't carry over values
        rowObsArr.fill(0)
    
    if verbose and rowsToCalculate[0] == 0: print("    Time:", time.time() - tScore, flush=True)
    

# Function that deploys the processes used to calculate the scores for the s3 metric. Also call function to store scores
def s3Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start all the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
        pool.starmap(s3Score, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(expFreqPath), itertools.repeat(verbose)))
    pool.join()

    chrName = pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").iloc[0, 0]

    storeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputDirPath, fileTag, filename, chrName)

# Helper for the multiprocessing implemented in s3
def s3Score(dataFilePath, rowsToCalculate, expFreqPath, verbose):
    # Read in the data
    if verbose and rowsToCalculate[0] == 0: print("Reading data from file...", flush=True); tRead = time.time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(dataFilePath, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    dataArr = pd.read_table(dataFilePath, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time.time() - tRead, flush=True) 

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    multiprocessRows, numCols = dataArr.shape
    numStates = sharedArr[2]

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2)), dtype=np.int16).T

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates), dtype=np.float32) / (numCols * (numCols - 1)), expFreqArr)

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tScore = time.time()
        percentDone = 0
    printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]

    # Calculte the scores and store them in the shared array
    scoreArr = sharedToNumpy(*sharedArr)
    rowScoreArr = np.zeros(numStates, dtype=np.float32)
    for dataRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if verbose and rowsToCalculate[0] == 0 and dataRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and dataRow in printCheckmarks: print(".", end="", flush=True)

        if dataRow < dataArr.shape[0]:
            # Pull the scores from the precalculated score array add them to the correct index in the rowScoreArr
            np.add.at(rowScoreArr, dataArr[dataRow, basePermutationArr[1]], scoreArrOnes[basePermutationArr[0], basePermutationArr[1], dataArr[dataRow, basePermutationArr[0]], dataArr[dataRow, basePermutationArr[1]]])

            # Store the scores in the shared score array
            scoreArr[scoreRow] = rowScoreArr

            # Reset the array so it doesn't carry over scores from other rows
            rowScoreArr.fill(0)

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time.time() - tScore, flush=True)

# Helper to store the score arrays combined with the location arrays
def storeScores(scoreArr, outputDirPath, fileTag, filename, chrName):
    # # Create a location array
    # numRows = scoreArr.shape[0]
    # locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(numRows)])

    # scoreFilename = "temp_scores_{}_{}.npz".format(fileTag, filename)
    # scoreFilePath = outputDirPath / scoreFilename

    # # Savez saves space allowing location to be stored as string and scoreArr as float
    # np.savez_compressed(scoreFilePath, locationArr=locationArr, scoreArr=scoreArr)



# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)


# Helper to calculate KL-score for Nd arrays (cleans up the code)
def klScoreND(obs, exp):
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)


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
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]), strToBool(sys.argv[8]))