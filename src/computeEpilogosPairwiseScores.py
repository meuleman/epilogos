import numpy as np
import sys
from pathlib import Path
from math import log2
import pandas as pd
from time import time
from multiprocessing import Pool, cpu_count, RawArray
from itertools import repeat
import gzip
from contextlib import closing

def main(file1, file2, numStates, saliency, outputDirPath, expFreqPath, fileTag, numProcesses):
    tTotal = time()

    file1Path = Path(file1)
    file2Path = Path(file2)
    outputDirPath = Path(outputDirPath)

    filename = file1Path.name.split(".")[0]

    if file1Path.name.endswith("gz"):
        with gzip.open(file1Path, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))
    else:
        with open(file1Path, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    # Split the rows up according to the number of cores we have available
    rowList = []
    for i in range(numProcesses):
        rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
        rowList.append(rowsToCalculate)

    # Calculate the scores
    realScoreArr1, realScoreArr2, nullScoreArr1, nullScoreArr2 = determineSaliency(saliency, file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses)
    # realScoreArr1, realScoreArr2 = determineSaliency(saliency, file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses, "real")

    print("Calculating Raw Differences...")
    tDiff = time()
    nullDiffArr = nullScoreArr1 - nullScoreArr2
    realDiffArr = realScoreArr1 - realScoreArr2
    print("    Time:", time() - tDiff)

    # Only calculate the distances for the null data in this step
    print("Calculating Squared Euclidean Distance and Maximum Contributing Difference...")
    tDistance = time()
    diffSign = np.sign(np.sum(nullDiffArr, axis=1))
    nullDistancesArr = np.sum(np.square(nullDiffArr), axis=1) * diffSign
    print("    Time:", time() - tDistance)

    print("Writing output to disk...")
    tWrite = time()
    chrName = pd.read_table(file1Path, nrows=1, header=None, sep="\t").iloc[0, 0]
    # If it's the real data, we will just write the delta and calculate metrics in computeEpilogosPairwiseVisual
    # If it's the null data, we will just write the signed squared euclidean distances
    writeReal(realDiffArr, outputDirPath, fileTag, filename, chrName)
    writeNull(nullDistancesArr, outputDirPath, fileTag, filename, chrName)
    print("    Time:", time() - tWrite)

    print("Total Time:", time() - tTotal)

# Helper function to determine the pick the correct saliency function
def determineSaliency(saliency, file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses):
    if saliency == 1:
        return s1Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses)
    elif saliency == 2:
        return s2Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses)
    elif saliency == 3:
        raise ValueError("A saliency value of 3 is not supported for pairwise comparison")
    else:
        raise ValueError("Inputed saliency value not supported")

def readInData(file1Path, file2Path, rowsToCalculate):
    # Read in the data
    if rowsToCalculate[0] == 0: print("Reading data from file 1..."); tRead1 = time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file1Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file1Arr = pd.read_table(file1Path, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if rowsToCalculate[0] == 0: print("    Time: ", time() - tRead1)

    if rowsToCalculate[0] == 0: print("Reading data from file 2..."); tRead2 = time()
    cols = range(3, pd.read_table(file2Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file2Arr = pd.read_table(file2Path, usecols=cols, skiprows=rowsToCalculate[0], nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if rowsToCalculate[0] == 0: print("    Time: ", time() - tRead2)

    if rowsToCalculate[0] == 0: print("Shuffling input matrices..."); tShuffle = time()
    # Combining the arrays for per row shuffling
    combinedArr = np.concatenate((file1Arr, file2Arr), axis=1)

    # Row independent vectorized shuffling of the 2 arrays
    randomIndices = np.argsort(np.random.rand(*combinedArr.shape), axis=1)
    shuffledCombinedArr = np.take_along_axis(combinedArr, randomIndices, axis=1)
    if rowsToCalculate[0] == 0: print("    Time:", time() - tShuffle)
    
    # shuffledCombinedArr is split by size of the original arrays
    return file1Arr, file2Arr, shuffledCombinedArr[:,:file1Arr.shape[1]], shuffledCombinedArr[:,file1Arr.shape[1]:]

# Helper for unflattening a shared array into a 2d numpy array
def sharedToNumpy(sharedArr, numRows, numStates):
    return np.frombuffer(sharedArr, dtype=np.float32).reshape((numRows, numStates))

# initiliazer for multiprocessing
def _init(sharedArr1_, sharedArr2_, shuffledSharedArr1_, shuffledSharedArr2_, totalRows, numStates):
    global sharedArr1
    global sharedArr2
    global shuffledSharedArr1
    global shuffledSharedArr2

    sharedArr1 = (sharedArr1_, totalRows, numStates)
    sharedArr2 = (sharedArr2_, totalRows, numStates)
    shuffledSharedArr1 = (shuffledSharedArr1_, totalRows, numStates)
    shuffledSharedArr2 = (shuffledSharedArr2_, totalRows, numStates)

# Function that calculates the scores for the S1 metric
def s1Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses):
    print("\nNumber of Processes:", numProcesses)

    # Shared arrays across the multiple processes for the scores
    # We avoid race conditions by writing to separate parts of the array in each process
    sharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    sharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(Pool(numProcesses, initializer=_init, initargs=(sharedArr1, sharedArr2, shuffledSharedArr1, shuffledSharedArr2, totalRows, numStates))) as pool:
        pool.starmap(s1Score, zip(repeat(file1Path), repeat(file2Path), rowList, repeat(expFreqPath)))
    pool.join()

    return sharedToNumpy(sharedArr1, totalRows, numStates), sharedToNumpy(sharedArr2, totalRows, numStates), sharedToNumpy(shuffledSharedArr1, totalRows, numStates), sharedToNumpy(shuffledSharedArr2, totalRows, numStates)

def s1Score(file1Path, file2Path, rowsToCalculate, expFreqPath):
    file1Arr, file2Arr, shuffledFile1Arr, shuffledFile2Arr = readInData(file1Path, file2Path, rowsToCalculate)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    numCols1 = file1Arr.shape[1]
    numCols2 = file2Arr.shape[1]
    numStates = sharedArr1[2]

    realScoreArr1 = sharedToNumpy(*sharedArr1)
    realScoreArr2 = sharedToNumpy(*sharedArr2)
    nullScoreArr1 = sharedToNumpy(*shuffledSharedArr1)
    nullScoreArr2 = sharedToNumpy(*shuffledSharedArr2)
    
    if rowsToCalculate[0] == 0:
        print("Calculating Scores...")
        tScore = time()
        printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # Calculate the observed frequencies and final scores for the designated rows
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if rowsToCalculate[0] == 0 and obsRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone))

        # klScoreND input is obsFreq and expFreq
        realScoreArr1[scoreRow] = klScoreND(rowObsS1(file1Arr, obsRow, numCols1, numStates), expFreqArr)
        realScoreArr2[scoreRow] = klScoreND(rowObsS1(file2Arr, obsRow, numCols2, numStates), expFreqArr)
        nullScoreArr1[scoreRow] = klScoreND(rowObsS1(shuffledFile1Arr, obsRow, numCols1, numStates), expFreqArr)
        nullScoreArr2[scoreRow] = klScoreND(rowObsS1(shuffledFile2Arr, obsRow, numCols2, numStates), expFreqArr)

        # uniqueStates, stateCounts = np.unique(file1Arr[obsRow], return_counts=True)
        # for i, state in enumerate(uniqueStates):
        #     # klScore input is obsFreq and expFreq
        #     realScoreArr1[scoreRow, state] = klScore(stateCounts[i] / numCols1, expFreqArr[state])

        # uniqueStates, stateCounts = np.unique(file2Arr[obsRow], return_counts=True)
        # for i, state in enumerate(uniqueStates):
        #     realScoreArr2[scoreRow, state] = klScore(stateCounts[i] / numCols2, expFreqArr[state])

        # uniqueStates, stateCounts = np.unique(shuffledFile1Arr[obsRow], return_counts=True)
        # for i, state in enumerate(uniqueStates):
        #     nullScoreArr1[scoreRow, state] = klScore(stateCounts[i] / numCols1, expFreqArr[state])

        # uniqueStates, stateCounts = np.unique(shuffledFile2Arr[obsRow], return_counts=True)
        # for i, state in enumerate(uniqueStates):
        #     nullScoreArr2[scoreRow, state] = klScore(stateCounts[i] / numCols2, expFreqArr[state])

    if rowsToCalculate[0] == 0: print("    Time:", time() - tScore)
        
# Helper for calculating observed in s1 case
def rowObsS1(dataArr, row, numCols, numStates):
    rowObsArr = np.zeros(numStates)
    uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
    for i, state in enumerate(uniqueStates):
        rowObsArr[state] = stateCounts[i] / numCols
    return rowObsArr

# Function that calculates the scores for the S2 metric
def s2Multi(file1Path, file2Path, rowList, totalRows, numStates, expFreqPath, numProcesses):
    print("\nNumber of Processes:", numProcesses)

    # Shared arrays across the multiple processes for the scores
    # We avoid race conditions by writing to separate parts of the array in each process
    sharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    sharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(Pool(numProcesses, initializer=_init, initargs=(sharedArr1, sharedArr2, shuffledSharedArr1, shuffledSharedArr2, totalRows, numStates))) as pool:
        pool.starmap(s2Score, zip(repeat(file1Path), repeat(file2Path), rowList, repeat(expFreqPath)))
    pool.join()

    return sharedToNumpy(sharedArr1, totalRows, numStates), sharedToNumpy(sharedArr2, totalRows, numStates), sharedToNumpy(shuffledSharedArr1, totalRows, numStates), sharedToNumpy(shuffledSharedArr2, totalRows, numStates)


def s2Score(file1Path, file2Path, rowsToCalculate, expFreqPath):
    file1Arr, file2Arr, shuffledFile1Arr, shuffledFile2Arr = readInData(file1Path, file2Path, rowsToCalculate)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    numCols1 = file1Arr.shape[1]
    numCols2 = file2Arr.shape[1]
    numStates = sharedArr1[2]

    if rowsToCalculate[0] == 0:
        print("Calculating Scores...")
        tScore = time()
        printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]
        percentDone = 0

    # Calculte the scores and store them in the shared array
    realScoreArr1 = sharedToNumpy(*sharedArr1)
    realScoreArr2 = sharedToNumpy(*sharedArr2)
    nullScoreArr1 = sharedToNumpy(*shuffledSharedArr1)
    nullScoreArr2 = sharedToNumpy(*shuffledSharedArr2)

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    permutations1 = numCols1 * (numCols1 - 1)
    permutations2 = numCols2 * (numCols2 - 1)
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if rowsToCalculate[0] == 0 and obsRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone))
        
        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        realScoreArr1[scoreRow] = klScoreND(rowObsS2(file1Arr, obsRow, permutations1, numStates), expFreqArr).sum(axis=0)
        realScoreArr2[scoreRow] = klScoreND(rowObsS2(file2Arr, obsRow, permutations2, numStates), expFreqArr).sum(axis=0)
        nullScoreArr1[scoreRow] = klScoreND(rowObsS2(shuffledFile1Arr, obsRow, permutations1, numStates), expFreqArr).sum(axis=0)
        nullScoreArr2[scoreRow] = klScoreND(rowObsS2(shuffledFile2Arr, obsRow, permutations2, numStates), expFreqArr).sum(axis=0)

    if rowsToCalculate[0] == 0: print("    Time:", time() - tScore)

# Helper for calculating observed in s2 case
def rowObsS2(dataArr, row, permutations, numStates):
    rowObsArr = np.zeros((numStates, numStates))
    uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
    for i, state1 in enumerate(uniqueStates):
        for j, state2 in enumerate(uniqueStates):
            if state1 == state2:
                rowObsArr[state1, state2] = stateCounts[i] * (stateCounts[i] - 1) / permutations
            else: # state1 > state2 or state1 < state2
                rowObsArr[state1, state2] = stateCounts[i] * stateCounts[j] / permutations
    return rowObsArr

# # Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
# def klScore(obs, exp):
#     if obs == 0.0:
#         return 0.0
#     else:
#         return obs * log2(obs / exp)


# Helper to calculate KL-score for Nd arrays (cleans up the code)
def klScoreND(obs, exp):
    return obs * np.ma.log2(np.ma.divide(obs, exp).filled(0)).filled(0)


# Helper for writing when we are working with real data
def writeReal(diffArr, outputDirPath, fileTag, filename, chrName):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    deltaTxtPath = outputDirPath / "pairwiseDelta_{}_{}.txt.gz".format(fileTag, filename)
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
def writeNull(nullDistancesArr, outputDirPath, fileTag, filename, chrName):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    nullDistancesPath = outputDirPath / "temp_nullDistances_{}_{}.npz".format(fileTag, filename)

    np.savez_compressed(nullDistancesPath, chrName=np.array([chrName]), nullDistances=nullDistancesArr)


# Helper for reading number of lines in input file
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7], int(sys.argv[8]))