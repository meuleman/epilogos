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
import gzip
from contextlib import closing

def main(file1, file2, numStates, saliency, outputDirPath, expFreqPath, fileTag, numProcesses, realOrNull):
    tTotal = time.time()

    file1Path = Path(file1)
    file2Path = Path(file2)
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
                line = gf.readline
            basePairs = int(line.split()[1])

    totalRows = math.ceil(basePairs / 200)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = multiprocessing.cpu_count()

    # Split the rows up according to the number of cores we have available
    rowList = []
    # for i in range(numProcesses):
    #     rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
    #     rowList.append(rowsToCalculate)
    rowsToCalculate = (0, totalRows // numProcesses)
    rowList.append(rowsToCalculate)

    print("Calculating Scores...")
    tScore = time.time()
    score1Arr, score2Arr = determineSaliency(saliency, file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, realOrNull)
    print("    Time: ", time.time() - tScore)

    print("Calculating Raw Differences...")
    tDiff = time.time()
    diffArr = score1Arr - score2Arr
    print("    Time:", time.time() - tDiff)

    # Only calculate the distances for the null data in this step
    if realOrNull.lower() == "null":
        print("Calculating Squared Euclidean Distance and Maximum Contributing Difference...")
        tDistance = time.time()
        diffSign = np.sign(np.sum(diffArr, axis=1))
        nullDistancesArr = np.sum(np.square(diffArr), axis=1) * diffSign
        print("    Time:", time.time() - tDistance)

    print("Writing output to disk...")
    tWrite = time.time()
    # If it's the real data, we will just write the delta and calculate metrics in computeEpilogosPairwiseVisual
    # If it's the null data, we will just write the signed squared euclidean distances
    if realOrNull.lower() == "real":
        writeReal(diffArr, outputDirPath, fileTag, chrName)
    elif realOrNull.lower() == "null":
        writeNull(nullDistancesArr, outputDirPath, fileTag, chrName)
    else:
        raise ValueError("Could not determine if writing real or null data")
    print("    Time:", time.time() - tWrite)

    print("Total Time:", time.time() - tTotal)


# Helper function to determine the pick the correct saliency function
def determineSaliency(saliency, file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, realOrNull):
    if saliency == 1:
        return s1Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, realOrNull)
    elif saliency == 2:
        return s2Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, realOrNull)
    elif saliency == 3:
        raise ValueError("A saliency value of 3 is not supported for pairwise comparison")
    else:
        raise ValueError("Inputed saliency value not supported")

def readInData(file1Path, file2Path, rowsToCalculate, realOrNull):
    # Read in the data
    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead1)

    print("Reading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t")
    print("    Time: ", time.time() - tRead2)

    if realOrNull.lower() == "real":
        # Converting to a np array for faster functions later
        print("Converting to numpy arrays...")
        tConvert = time.time()
        file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int) - 1
        file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int) - 1
        print("    Time: ", time.time() - tConvert)
    elif realOrNull.lower() == "null":
        # Converting to a np array for faster functions later
        print("Converting to numpy arrays...")
        tConvert = time.time()
        unshuffledFile1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int) - 1
        unshuffledFile2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int) - 1
        print("    Time: ", time.time() - tConvert)

        print("Shuffling input matrices...")
        tShuffle = time.time()
        # Combining the arrays for per row shuffling
        combinedArr = np.concatenate((unshuffledFile1Arr, unshuffledFile2Arr), axis=1)

        # Row independent vectorized shuffling of the 2 arrays
        randomIndices = np.argsort(np.random.rand(*combinedArr.shape), axis=1)
        shuffledCombinedArr = np.take_along_axis(combinedArr, randomIndices, axis=1)
        file1Arr = shuffledCombinedArr[:,:unshuffledFile1Arr.shape[1]]
        file2Arr = shuffledCombinedArr[:,unshuffledFile1Arr.shape[1]:]
        print("    Time:", time.time() - tShuffle)
    else:
        raise ValueError("Could not determine whether score calculation is for real or null data")

    return file1Arr, file2Arr

# Helper for unflattening a shared array into a 2d numpy array
def sharedToNumpy(sharedArr, numRows, numStates):
    return np.frombuffer(sharedArr.get_obj(), dtype=np.float32).reshape((numRows, numStates))

# initiliazer for multiprocessing
def _init(sharedArr1_, sharedArr2_, totalRows, numStates):
    global sharedArr1
    global sharedArr2
    sharedArr1 = (sharedArr1_, totalRows, numStates)
    sharedArr2 = (sharedArr2_, totalRows, numStates)

# Function that calculates the scores for the S1 metric
def s1Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, realOrNull):
    print("NUM PROCESSES:", numProcesses)

    sharedArr1 = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    sharedArr2 = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=(sharedArr1, sharedArr2, totalRows, numStates))) as pool:
        pool.starmap(s1Score, zip(itertools.repeat(file1Path), itertools.repeat(file2Path), rowList, itertools.repeat(expFreqPath), itertools.repeat(realOrNull)))
    pool.join()

    return sharedToNumpy(sharedArr1, totalRows, numStates), sharedToNumpy(sharedArr2, totalRows, numStates)

def s1Score(file1Path, file2Path, rowsToCalculate, expFreqPath, realOrNull):
    file1Arr, file2Arr = readInData(file1Path, file2Path, rowsToCalculate, realOrNull)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    numCols1 = file1Arr.shape[1]
    numCols2 = file2Arr.shape[1]

    scoreArr1 = sharedToNumpy(*sharedArr1)
    scoreArr2 = sharedToNumpy(*sharedArr2)

    print(rowsToCalculate)
    print(file1Arr.shape)
    print(file2Arr.shape)
    print(expFreqArr.shape)
    print(scoreArr1.shape)
    print(scoreArr2.shape)
    print(rowsToCalculate[0], rowsToCalculate[1])
    
    # Calculate the observed frequencies and final scores for the designated rows
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):
        uniqueStates, stateCounts = np.unique(file1Arr[obsRow], return_counts=True)
        for i, state in enumerate(uniqueStates):
            # Function input is obsFreq and expFreq
            scoreArr1[scoreRow, state] = klScore(stateCounts[i] / numCols1, expFreqArr[state])
        uniqueStates, stateCounts = np.unique(file2Arr[obsRow], return_counts=True)
        for i, state in enumerate(uniqueStates):
            # Function input is obsFreq and expFreq
            scoreArr2[scoreRow, state] = klScore(stateCounts[i] / numCols2, expFreqArr[state])
        

# Function that calculates the scores for the S2 metric
def s2Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, realOrNull):
    print("NUM PROCESSES:", numProcesses)

    sharedArr1 = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    sharedArr2 = multiprocessing.Array(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    if (sys.version_info < (3, 8)):
        print("\nFor maximum efficiency please update python to version 3.8 or later")
        print("NOTE: The code will still run in a lower version, but will be slightly slower\n")

    # Start the processes
    with closing(multiprocessing.Pool(numProcesses, initializer=_init, initargs=(sharedArr1, sharedArr2, totalRows, numStates))) as pool:
        pool.starmap(s2Score, zip(itertools.repeat(file1Path), itertools.repeat(file2Path), rowList, itertools.repeat(expFreqPath), itertools.repeat(realOrNull)))
    pool.join()

    return sharedToNumpy(sharedArr1, totalRows, numStates), sharedToNumpy(sharedArr2, totalRows, numStates)


def s2Score(file1Path, file2Path, rowsToCalculate, expFreqPath, realOrNull):
    file1Arr, file2Arr = readInData(file1Path, file2Path, rowsToCalculate, realOrNull)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    obsFreqArr1 = s2Obs(file1Arr, sharedArr1[2])
    obsFreqArr2 = s2Obs(file2Arr, sharedArr1[2])

    # Calculte the scores and store them in the shared array
    scoreArr1 = sharedToNumpy(*sharedArr1)
    scoreArr2 = sharedToNumpy(*sharedArr2)

    print(rowsToCalculate)
    print(file1Arr.shape)
    print(file2Arr.shape)
    print(expFreqArr.shape)
    print(obsFreqArr1.shape)
    print(obsFreqArr2.shape)
    print(scoreArr1.shape)
    print(scoreArr2.shape)
    print(rowsToCalculate[0], rowsToCalculate[1])

    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):
        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        scoreArr1[scoreRow] = klScoreND(obsFreqArr1[obsRow], expFreqArr).sum(axis=0)
        scoreArr2[scoreRow] = klScoreND(obsFreqArr2[obsRow], expFreqArr).sum(axis=0)

# Helper for calculating the observed frequencies in the s2 metric
def s2Obs(dataArr, numStates):
    multiprocessRows, numCols = dataArr.shape

    obsFreqArr = np.zeros((multiprocessRows, numStates, numStates))

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
                        obsFreqArr[row, state1, state2]  = ncr(stateCounts[i], 2) / combinations
                    else: # state1 > state2 or state1 < state2
                        obsFreqArr[row, state1, state2]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
    else:
        combinations = math.comb(numCols, 2)
        for row in range(multiprocessRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i, state1 in enumerate(uniqueStates):
                for j, state2 in enumerate(uniqueStates):
                    if state1 == state2:
                        obsFreqArr[row, state1, state2]  = math.comb(stateCounts[i], 2) / combinations
                    else: # state1 > state2 or state1 < state2
                        obsFreqArr[row, state1, state2]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix

    return obsFreqArr

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


# Helper for writing when we are working with real data
def writeReal(diffArr, outputDirPath, fileTag, chrName):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    deltaTxtPath = outputDirPath / "pairwiseDelta_{}_{}.txt.gz".format(fileTag, chrName)
    deltaTxt = gzip.open(deltaTxtPath, "wt")

    # Create a location array
    numRows = diffArr.shape[0]
    locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(numRows)])

    # Creating a string to write out the raw differences (faster than np.savetxt)
    deltaTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(diffArr.shape[1] - 1)) + "{1[%d]:.5f}\n" % (diffArr.shape[1] - 1)
    deltaStr = "".join(deltaTemplate.format(locationArr[i], diffArr[i]) for i in range(diffArr.shape[0]))

    deltaTxt.write(deltaStr)
    deltaTxt.close()


# Helper for writing when we are working with null data
def writeNull(nullDistancesArr, outputDirPath, fileTag, chrName):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    nullDistancesTxtPath = outputDirPath / "nullDistances_{}_{}.txt.gz".format(fileTag, chrName)
    nullDistancesTxt = gzip.open(nullDistancesTxtPath, "wt")

    # Create a location array
    numRows = nullDistancesArr.shape[0]
    locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(numRows)])

    # Creating a string to write out the nullDistancess array
    nullDistancesTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\n"
    nullDistancesStr = "".join(nullDistancesTemplate.format(locationArr[i], nullDistancesArr[i]) for i in range(len(nullDistancesArr)))

    nullDistancesTxt.write(nullDistancesStr)
    nullDistancesTxt.close()


# Helper to store the score arrays combined with the location arrays
def storeScores(dataArr, scoreArr, locationArr, outputDirPath, fileTag, filename):
    # Creating a file path
    scoreFilename = "temp_scores_{}_{}.npy".format(fileTag, filename)
    scoreFilePath = outputDirPath / scoreFilename

    # Concatenating the locationArr and dataArr into one helps writing later
    combinedArr = np.concatenate((locationArr, scoreArr), axis=1)

    np.save(scoreFilePath, combinedArr, allow_pickle=False)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7], int(sys.argv[8]), sys.argv[9])