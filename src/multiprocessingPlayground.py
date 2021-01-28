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
    filePath = Path(file)
    outputDirPath = Path(outputDirPath)
    filename = filePath.name.split(".")[0]

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    print("FILE:", filename)

    # Read in the data
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(filePath, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    locationArr = dataDF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count()

    numRows, numCols = dataArr.shape

    # Split the rows up according to the number of cores we have available
    rowList = []
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRows // numProcesses, (i+1) * numRows // numProcesses)
        rowList.append(rowsToCalculate)

    # Split the data array up so we don't have to pass it all to everything
    dataArrList = []
    for rowsToCalculate in rowList:
        dataArrList.append(dataArr[rowsToCalculate])


    determineSaliency(saliency, dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses)

    print("Total Time:", time.time() - tTotal)

def determineSaliency(saliency, dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses):
    if saliency == 1:
        s1Multi(dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses)
    elif saliency == 2:
        s2Multi(dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses)
    elif saliency == 3:
        s3Multi(dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses)
    else:
        print("Inputed saliency value not supported")
        return


# Helper for unflattening a shared array into a 2d numpy array
def sharedToNumpy(sharedArr, numRows, numStates):
    return np.frombuffer(sharedArr.get_obj(), dtype=np.float32).reshape((numRows, numStates))

sharedArr=None
inputInfo=None

# initiliazer for multiprocessing
def _init(sharedArr_, inputInfo_):
    global sharedArr
    global inputInfo
    sharedArr = sharedArr_
    inputInfo = inputInfo_

# Function that deploys the processes used to calculate the scores for the s1 metric. Also call function to store scores
def s1Multi(dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), numRows * numStates)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, numRows, numStates), (expFreqArr, numCols)))) as pool:
        pool.starmap(s1Score, zip(rowList, dataArrList))
    pool.join()

    storeScores(sharedToNumpy(sharedArr, numRows, numStates), locationArr, outputDirPath, fileTag, filename)

# Calculates the scores for the s1 metric over a given range of rows
def s1Score(rowsToCalculate, dataArr):
    # inputInfo[0] is expFreqArr
    # inputInfo[1] is numCols

    scoreArr = sharedToNumpy(*sharedArr)
    # Calculate the observed frequencies and final scores for the designated rows
    for row in range(len(rowsToCalculate)):
        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
        for i in range(len(uniqueStates)):
            # Function input is obsFreq and expFreq
            scoreArr[row, uniqueStates[i]] = klScore(stateCounts[i] / (inputInfo[1]), inputInfo[0][uniqueStates[i]])


# Function that deploys the processes used to calculate the scores for the s2 metric. Also call function to store scores
def s2Multi(dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), numRows * numStates)

    if (sys.version_info < (3, 8)):
        print("\nFor maximum efficiency please update python to version 3.8 or later")
        print("NOTE: The code will still run in a lower version, but will be slightly slower\n")

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, numRows, numStates), (expFreqArr, numCols)))) as pool:
        pool.starmap(s2Score, zip(rowList, dataArrList))
    pool.join()

    storeScores(sharedToNumpy(sharedArr, numRows, numStates), locationArr, outputDirPath, fileTag, filename)


# Calculates the scores for the s2 metric over a given range of rows
def s2Score(rowsToCalculate, dataArr):
    # sharedArr[2] is numStates
    # inputInfo[0] is expFreqArr
    # inputInfo[1] is numCols

    # Calculate the observed frequencies
    obsFreqArr = np.zeros((len(rowsToCalculate), sharedArr[2], sharedArr[2]))

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    if (sys.version_info < (3, 8)):
        combinations = ncr(inputInfo[1], 2)
        for writeRow, readRow in enumerate(range(len(rowsToCalculate))):
            uniqueStates, stateCounts = np.unique(dataArr[readRow], return_counts=True)
            for i in range(len(uniqueStates)):
                for j in range(len(uniqueStates)):
                    if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                        obsFreqArr[writeRow, uniqueStates[i], uniqueStates[j]]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                    elif uniqueStates[i] == uniqueStates[j]:
                        obsFreqArr[writeRow, uniqueStates[i], uniqueStates[j]]  = ncr(stateCounts[i], 2) / combinations
    else:
        combinations = math.comb(inputInfo[1], 2)
        for writeRow, readRow in enumerate(range(len(rowsToCalculate))):
            uniqueStates, stateCounts = np.unique(dataArr[readRow], return_counts=True)
            for i in range(len(uniqueStates)):
                for j in range(len(uniqueStates)):
                    if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                        obsFreqArr[writeRow, uniqueStates[i], uniqueStates[j]]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                    elif uniqueStates[i] == uniqueStates[j]:
                        obsFreqArr[writeRow, uniqueStates[i], uniqueStates[j]]  = math.comb(stateCounts[i], 2) / combinations

    # Calculte the scores and store them in the shared array
    scoreArr = sharedToNumpy(*sharedArr)
    for obsRow, scoreRow in enumerate(rowsToCalculate):
        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        scoreArr[scoreRow] = klScoreND(obsFreqArr[obsRow], inputInfo[0]).sum(axis=0)
    

# Function that deploys the processes used to calculate the scores for the s3 metric. Also call function to store scores
def s3Multi(dataArrList, rowList, locationArr, numRows, numCols, numStates, outputDirPath, expFreqArr, fileTag, filename, numProcesses):
    print("NUM PROCESSES:", numProcesses)

    sharedArr = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), numRows * numStates)

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2))).T

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates)) / (numCols * (numCols - 1)), expFreqArr)

    # Start all the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=((sharedArr, numRows, numStates), (basePermutationArr, scoreArrOnes, numCols)))) as pool:
        pool.map(s3Score, rowList)
    pool.join()

    storeScores(sharedToNumpy(sharedArr, numRows, numStates), locationArr, outputDirPath, fileTag, filename)

# Helper for the multiprocessing implemented in s3
def s3Score(rowsToCalculate, dataArr):
    # sharedArr[2] is numStates
    # inputInfo[0] is basePermutationArr
    # inputInfo[1] is scoreArrOnes
    # inputInfo[2] is numCols

    scoreArr = sharedToNumpy(*sharedArr)

    rowScoreArr = np.zeros((inputInfo[2], inputInfo[2], sharedArr[1], sharedArr[1]), dtype=np.float32)
    for row in range(len(rowsToCalculate)):
        # Reset the array so it doesn't carry over scores from other rows
        rowScoreArr.fill(0)

        # Pull the scores from the precalculated score array and put them into the correct positions for the state combinations that we observe
        rowScoreArr[inputInfo[0][0], inputInfo[0][1], dataArr[row, inputInfo[0][0]], dataArr[row, inputInfo[0][1]]] = inputInfo[1][inputInfo[0][0], inputInfo[0][1], dataArr[row, inputInfo[0][0]], dataArr[row, inputInfo[0][1]]]

        # Flatten the scores and put them into the shared score array
        scoreArr[row] = rowScoreArr.sum(axis=(0,1,2))

# Helper to store the score arrays combined with the location arrays
def storeScores(scoreArr, locationArr, outputDirPath, fileTag, filename):
    # Creating a file path
    scoreFilename = "temp_scores_{}_{}.npz".format(fileTag, filename)
    scoreFilePath = outputDirPath / scoreFilename

    # Savez saves space allowing location to be stored as string and scoreArr as float
    np.savez_compressed(scoreFilePath, locationArr=locationArr, scoreArr=scoreArr)

# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to calculate KL-score for 2d arrays (cleans up the code)
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