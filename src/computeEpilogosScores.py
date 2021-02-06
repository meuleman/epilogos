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

def main(file, numStates, saliency, outputDirPath, expFreqPath, fileTag, numProcesses):
    tTotal = time.time()
    dataFilePath = Path(file)
    outputDirPath = Path(outputDirPath)

    ###################################################
    ##
    ##   chrName is hacky, find a better way later
    ##
    ####################################################
    chrName  = dataFilePath.name.split("_")[-1].split(".")[0]

    # Get the genome file
    genomeFileList = list(dataFilePath.parents[0].glob("*.genome"))
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

    determineSaliency(saliency, dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses)

    print("Total Time:", time.time() - tTotal)


def determineSaliency(saliency, dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses):
    if saliency == 1:
        s1Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses)
    elif saliency == 2:
        s2Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses)
    elif saliency == 3:
        s3Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses)
    else:
        print("Inputed saliency value not supported")
        return


# Helper for unflattening a shared array into a 2d numpy array
def sharedToNumpy(sharedArr, numRows, numStates):
    return np.frombuffer(sharedArr.get_obj(), dtype=np.float32).reshape((numRows, numStates))

# initiliazer for multiprocessing
def _init(sharedArr_):
    global sharedArr
    sharedArr = sharedArr_

# Function that deploys the processes used to calculate the scores for the s1 metric. Also call function to store scores
def s1Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
        pool.starmap(s1Score, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(expFreqPath)))
    pool.join()

    storeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputDirPath, fileTag, chrName)

# Calculates the scores for the s1 metric over a given range of rows
def s1Score(dataFilePath, rowsToCalculate, expFreqPath):
    # Read in the data
    if rowsToCalculate[0] == 0:
        print("\nReading data from file...")
        tRead = time.time()

    dataDF = pd.read_table(dataFilePath, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    if rowsToCalculate[0] == 0:
        print("Converting to numpy array...")
        tConvert = time.time()

    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 

    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tConvert)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    numCols = dataArr.shape[1]

    if rowsToCalculate[0] == 0:
        print("\nCalculating Scores...")
        tScore = time.time()
        printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]
        percentDone = 0
    
    scoreArr = sharedToNumpy(*sharedArr)
    # Calculate the observed frequencies and final scores for the designated rows
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):
        
        if rowsToCalculate[0] == 0 and obsRow in printCheckmarks:
                print("{}% Completed".format(percentDone))
                percentDone += 10

        if obsRow < dataArr.shape[0]:
            uniqueStates, stateCounts = np.unique(dataArr[obsRow], return_counts=True)
            for i, state in enumerate(uniqueStates):
                # Function input is obsFreq and expFreq
                scoreArr[scoreRow, state] = klScore(stateCounts[i] / numCols, expFreqArr[state])
    
    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tScore)


# Function that deploys the processes used to calculate the scores for the s2 metric. Also call function to store scores
def s2Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    if (sys.version_info < (3, 8)):
        print("\nFor maximum efficiency please update python to version 3.8 or later")
        print("NOTE: The code will still run in a lower version, but will be slightly slower\n")

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
        pool.starmap(s2Score, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(expFreqPath)))
    pool.join()

    storeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputDirPath, fileTag, chrName)


# Calculates the scores for the s2 metric over a given range of rows
def s2Score(dataFilePath, rowsToCalculate, expFreqPath):
    # Read in the data
    if rowsToCalculate[0] == 0:
        print("\nReading data from file...")
        tRead = time.time()

    dataDF = pd.read_table(dataFilePath, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    if rowsToCalculate[0] == 0:
        print("Converting to numpy array...")
        tConvert = time.time()

    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tConvert)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    multiprocessRows, numCols = dataArr.shape
    numStates = sharedArr[2]

    # Calculate the observed frequencies
    obsFreqArr = np.zeros((multiprocessRows, numStates, numStates))

    if rowsToCalculate[0] == 0:
        print("\nCalculating Observed Frequencies...")
        tObs = time.time()
        printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    if (sys.version_info < (3, 8)):
        combinations = ncr(numCols, 2)
        for row in range(multiprocessRows):

            if rowsToCalculate[0] == 0 and row in printCheckmarks:
                print("{}% Completed".format(percentDone))
                percentDone += 10

            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i, state1 in enumerate(uniqueStates):
                for j, state2 in enumerate(uniqueStates):
                    if state1 == state2:
                        obsFreqArr[row, state1, state2]  = ncr(stateCounts[i], 2) / combinations
                    else: # state1 > state2 or state1 < state2
                        obsFreqArr[row, state1, state2]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
    else:
        combinations = math.comb(numCols, 2)
        for row in range(multiprocessRows):

            if rowsToCalculate[0] == 0 and row in printCheckmarks:
                print("{}% Completed".format(percentDone))
                percentDone += 10

            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i, state1 in enumerate(uniqueStates):
                for j, state2 in enumerate(uniqueStates):
                    if state1 == state2:
                        obsFreqArr[row, state1, state2]  = math.comb(stateCounts[i], 2) / combinations
                    else: # state1 > state2 or state1 < state2
                        obsFreqArr[row, state1, state2]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
    
    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tObs)

    if rowsToCalculate[0] == 0:
        print("\nCalculating Scores...")
        tScore = time.time()
        printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # Calculte the scores and store them in the shared array
    scoreArr = sharedToNumpy(*sharedArr)
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if rowsToCalculate[0] == 0 and obsRow in printCheckmarks:
                print("{}% Completed".format(percentDone))
                percentDone += 10

        if obsRow < obsFreqArr.shape[0]:
            # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
            scoreArr[scoreRow] = klScoreND(obsFreqArr[obsRow], expFreqArr).sum(axis=0)

    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tScore)
    

# Function that deploys the processes used to calculate the scores for the s3 metric. Also call function to store scores
def s3Multi(dataFilePath, rowList, totalRows, numStates, outputDirPath, expFreqPath, fileTag, chrName, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start all the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
        pool.starmap(s3Score, zip(itertools.repeat(dataFilePath), rowList, itertools.repeat(expFreqPath)))
    pool.join()

    storeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputDirPath, fileTag, chrName)

# Helper for the multiprocessing implemented in s3
def s3Score(dataFilePath, rowsToCalculate, expFreqPath):
    # Read in the data
    if rowsToCalculate[0] == 0:
        print("\nReading data from file...")
        tRead = time.time()

    dataDF = pd.read_table(dataFilePath, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    if rowsToCalculate[0] == 0:
        print("Converting to numpy array...")
        tConvert = time.time()

    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    
    if rowsToCalculate[0] == 0:
        print("    Time: ", time.time() - tConvert)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    multiprocessRows, numCols = dataArr.shape
    numStates = sharedArr[2]

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2)), dtype=np.int16).T

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates), dtype=np.float32) / (numCols * (numCols - 1)), expFreqArr)

    if rowsToCalculate[0] == 0:
        print("\nCalculating Scores...")
        tScore = time.time()
        printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # Calculte the scores and store them in the shared array
    scoreArr = sharedToNumpy(*sharedArr)
    rowScoreArr = np.zeros(numStates, dtype=np.float32)
    for dataRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if rowsToCalculate[0] == 0 and dataRow in printCheckmarks:
                print("{}% Completed".format(percentDone))
                percentDone += 10

        if dataRow < dataArr.shape[0]:

            # Reset the array so it doesn't carry over scores from other rows
            rowScoreArr.fill(0)

            # Pull the scores from the precalculated score array add them to the correct index in the rowScoreArr
            np.add.at(rowScoreArr, dataArr[dataRow, basePermutationArr[1]], scoreArrOnes[basePermutationArr[0], basePermutationArr[1], dataArr[dataRow, basePermutationArr[0]], dataArr[dataRow, basePermutationArr[1]]])

            # Store the scores in the shared score array
            scoreArr[scoreRow] = rowScoreArr

    if rowsToCalculate[0] == 0:
        print("    Time:", time.time() - tScore)

# Helper to store the score arrays combined with the location arrays
def storeScores(scoreArr, outputDirPath, fileTag, chrName):
    # Create a location array
    numRows = scoreArr.shape[0]
    locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(numRows)])

    # Creating a file path
    scoreFilename = "temp_scores_{}_{}.npz".format(fileTag, chrName)
    scoreFilePath = outputDirPath / scoreFilename

    # Savez saves space allowing location to be stored as string and scoreArr as float
    np.savez_compressed(scoreFilePath, locationArr=locationArr, scoreArr=scoreArr)

# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to calculate KL-score for Nd arrays (cleans up the code)
def klScoreND(obs, exp):
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)

# Helper to calculate combinations
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer // denom

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]))