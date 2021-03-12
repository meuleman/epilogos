import numpy as np
from sys import argv
from pathlib import Path
from math import log2
import pandas as pd
from time import time
import numpy.ma as ma
from multiprocessing import cpu_count, Pool, RawArray
from itertools import repeat, permutations
from contextlib import closing
from epilogosHelpers import strToBool, blocks, splitRows, readStates
import gzip

def main(file1, file2, numStates, saliency, outputDirPath, expFreqPath, fileTag, numProcesses, verbose):
    # try:
    if verbose: tTotal = time()
    
    file1Path = Path(file1)
    file2Path = Path(file2)
    outputDirPath = Path(outputDirPath)

    filename = file1Path.name.split(".")[0]

    if not verbose: print("    {}\t".format(filename), end="", flush=True)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    # Determine which rows to assign to each core
    rowList = splitRows(file1Path, numProcesses)

    if file2 == "null":
        calculateScores(saliency, file1Path, rowList, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose)
    else:
        calculateScoresPairwise(saliency, file1Path, file2Path, rowList, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose)
    
    print("Total Time:", time() - tTotal, flush=True) if verbose else print("\t[Done]", flush=True)


# Helper for unflattening a shared array into a 2d numpy array
def sharedToNumpy(sharedArr, numRows, numStates):
    # return np.frombuffer(sharedArr.get_obj(), dtype=np.float32).reshape((numRows, numStates))
    return np.frombuffer(sharedArr, dtype=np.float32).reshape((numRows, numStates))
    
# initiliazer for multiprocessing
def _init(sharedArr_):
    global sharedArr
    sharedArr = sharedArr_

# initiliazer for multiprocessing
def _initPairwise(sharedArr1_, sharedArr2_, shuffledSharedArr1_, shuffledSharedArr2_, totalRows, numStates):
    global sharedArr1
    global sharedArr2
    global shuffledSharedArr1
    global shuffledSharedArr2

    sharedArr1 = (sharedArr1_, totalRows, numStates)
    sharedArr2 = (sharedArr2_, totalRows, numStates)
    shuffledSharedArr1 = (shuffledSharedArr1_, totalRows, numStates)
    shuffledSharedArr2 = (shuffledSharedArr2_, totalRows, numStates)

# Function that deploys the processes used to calculate the scores for the s1 metric. Also call function to store scores
def calculateScores(saliency, file1Path, rowList, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses, flush=True)

    totalRows = rowList[-1][-1]

    # Shared arrays across the multiple processes for the scores
    # We avoid race conditions by writing to separate parts of the array in each process
    sharedArr = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(Pool(numProcesses, initializer=_init, initargs=((sharedArr, totalRows, numStates), ))) as pool:
        if saliency == 1:
            pool.starmap(s1Score, zip(repeat(file1Path), repeat(Path("null")), rowList, repeat(expFreqPath), repeat(verbose)))
        elif saliency == 2:
            pool.starmap(s2Score, zip(repeat(file1Path), repeat(Path("null")), rowList, repeat(expFreqPath), repeat(verbose)))
        elif saliency == 3:
            pool.starmap(s3Score, zip(repeat(file1Path), rowList, repeat(expFreqPath), repeat(verbose)))
        else:
            raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")
    pool.join()

    chrName = pd.read_table(file1Path, nrows=1, header=None, sep="\t").iloc[0, 0]

    outputTxtPath = outputDirPath / "scores_{}_{}.txt.gz".format(fileTag, filename)
    writeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputTxtPath, chrName)

def calculateScoresPairwise(saliency, file1Path, file2Path, rowList, numStates, outputDirPath, expFreqPath, fileTag, filename, numProcesses, verbose):
    if verbose: print("\nNumber of Processes:", numProcesses, flush=True)

    totalRows = rowList[-1][-1]

    # Shared arrays across the multiple processes for the scores
    # We avoid race conditions by writing to separate parts of the array in each process
    sharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    sharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(Pool(numProcesses, initializer=_initPairwise, initargs=(sharedArr1, sharedArr2, shuffledSharedArr1, shuffledSharedArr2, totalRows, numStates))) as pool:
        if saliency == 1:
            pool.starmap(s1Score, zip(repeat(file1Path), repeat(file2Path), rowList, repeat(expFreqPath), repeat(verbose)))
        elif saliency == 2:
            pool.starmap(s2Score, zip(repeat(file1Path), repeat(file2Path), rowList, repeat(expFreqPath), repeat(verbose)))
        else:
            raise ValueError("Please ensure that saliency metric is either 1 or 2 for Pairwise Epilogos")
    pool.join()

    # Calculate the differences between array 1 and 2 in both the real and null case
    if verbose: print("Calculating Raw Differences...", flush=True); tDiff = time()
    realDiffArr = sharedToNumpy(sharedArr1, totalRows, numStates) - sharedToNumpy(sharedArr2, totalRows, numStates)
    nullDiffArr = sharedToNumpy(shuffledSharedArr1, totalRows, numStates) - sharedToNumpy(shuffledSharedArr2, totalRows, numStates)
    if verbose: print("    Time:", time() - tDiff, flush=True)

    # Only calculate the distances for the null data in this step
    if verbose: print("Calculating Squared Euclidean Distance and Maximum Contributing Difference...", flush=True); tDistance = time()
    diffSign = np.sign(np.sum(nullDiffArr, axis=1))
    nullDistancesArr = np.sum(np.square(nullDiffArr), axis=1) * diffSign
    if verbose: print("    Time:", time() - tDistance, flush=True)

    # If it's the real data, we will just write the delta and calculate metrics in computeEpilogosPairwiseVisual
    # If it's the null data, we will just save the signed squared euclidean distances as a temporary npz file as we don't care about saving the data
    if verbose: print("Writing output to disk...", flush=True); tWrite = time()
    chrName = pd.read_table(file1Path, nrows=1, header=None, sep="\t").iloc[0, 0]
    realOutputPath = outputDirPath / "pairwiseDelta_{}_{}.txt.gz".format(fileTag, filename)
    writeScores(realDiffArr, realOutputPath, chrName)
    nullOutputPath = outputDirPath / "temp_nullDistances_{}_{}.npz".format(fileTag, filename)
    np.savez_compressed(nullOutputPath, chrName=np.array([chrName]), nullDistances=nullDistancesArr)
    if verbose: print("    Time:", time() - tWrite, flush=True)

# Calculates the scores for the s1 metric over a given range of rows
def s1Score(file1Path, file2Path, rowsToCalculate, expFreqPath, verbose):
    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    # Loading the data and creating the shared arrays for the scores
    if str(file2Path) == "null":
        dataArr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalculate=rowsToCalculate, expBool=False, verbose=verbose)

        numCols = dataArr.shape[1]
        numStates = sharedArr[2]

        scoreArr = sharedToNumpy(*sharedArr)
    else:
        file1Arr, file2Arr, shuffledFile1Arr, shuffledFile2Arr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalculate=rowsToCalculate, expBool=False, verbose=verbose)

        numCols1 = file1Arr.shape[1]
        numCols2 = file2Arr.shape[1]
        numStates = sharedArr1[2]

        realScoreArr1 = sharedToNumpy(*sharedArr1)
        realScoreArr2 = sharedToNumpy(*sharedArr2)
        nullScoreArr1 = sharedToNumpy(*shuffledSharedArr1)
        nullScoreArr2 = sharedToNumpy(*shuffledSharedArr2)

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tScore = time()
        percentDone = 0
    printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]

    # Calculate the observed frequencies and final scores for the designated rows
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):
        
        if verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: print(".", end="", flush=True)

        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        if str(file2Path) == "null":
            scoreArr[scoreRow] = klScoreND(rowObsS1(dataArr, obsRow, numCols, numStates), expFreqArr)
        else:
            realScoreArr1[scoreRow] = klScoreND(rowObsS1(file1Arr, obsRow, numCols1, numStates), expFreqArr)
            realScoreArr2[scoreRow] = klScoreND(rowObsS1(file2Arr, obsRow, numCols2, numStates), expFreqArr)
            nullScoreArr1[scoreRow] = klScoreND(rowObsS1(shuffledFile1Arr, obsRow, numCols1, numStates), expFreqArr)
            nullScoreArr2[scoreRow] = klScoreND(rowObsS1(shuffledFile2Arr, obsRow, numCols2, numStates), expFreqArr)

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time() - tScore, flush=True)

# Helper for calculating observed in s1 case
def rowObsS1(dataArr, row, numCols, numStates):
    rowObsArr = np.zeros(numStates)
    # Count number of each state in row and return array
    uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
    for i, state in enumerate(uniqueStates):
        rowObsArr[state] = stateCounts[i] / numCols
    return rowObsArr

# Calculates the scores for the s2 metric over a given range of rows
def s2Score(file1Path, file2Path, rowsToCalculate, expFreqPath, verbose):
    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    # Loading the data and creating the shared arrays for the scores
    if str(file2Path) == "null":
        dataArr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalculate=rowsToCalculate, expBool=False, verbose=verbose)

        numCols = dataArr.shape[1]
        numStates = sharedArr[2]

        scoreArr = sharedToNumpy(*sharedArr)

        # Need the permuations to effective count state pairs (see rowObsS2() for theory)
        permutations = numCols * (numCols - 1)
    else:
        file1Arr, file2Arr, shuffledFile1Arr, shuffledFile2Arr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalculate=rowsToCalculate, expBool=False, verbose=verbose)

        numCols1 = file1Arr.shape[1]
        numCols2 = file2Arr.shape[1]
        numStates = sharedArr1[2]

        realScoreArr1 = sharedToNumpy(*sharedArr1)
        realScoreArr2 = sharedToNumpy(*sharedArr2)
        nullScoreArr1 = sharedToNumpy(*shuffledSharedArr1)
        nullScoreArr2 = sharedToNumpy(*shuffledSharedArr2)

        # Need the permuations to effective count state pairs (see rowObsS2() for theory)
        permutations1 = numCols1 * (numCols1 - 1)
        permutations2 = numCols2 * (numCols2 - 1)

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tScore = time()
        percentDone = 0
    printCheckmarks = [int(rowsToCalculate[1] * float(i / 10)) for i in range(1, 10)]

    # Find scores for each row that the core is responsible for
    for obsRow, scoreRow in enumerate(range(rowsToCalculate[0], rowsToCalculate[1])):

        if verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: percentDone += 10; print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalculate[0] == 0 and obsRow in printCheckmarks: print(".", end="", flush=True)
        
        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        if str(file2Path) == "null":
            scoreArr[scoreRow] = klScoreND(rowObsS2(dataArr, obsRow, permutations, numStates), expFreqArr).sum(axis=0)
        else:
            realScoreArr1[scoreRow] = klScoreND(rowObsS2(file1Arr, obsRow, permutations1, numStates), expFreqArr).sum(axis=0)
            realScoreArr2[scoreRow] = klScoreND(rowObsS2(file2Arr, obsRow, permutations2, numStates), expFreqArr).sum(axis=0)
            nullScoreArr1[scoreRow] = klScoreND(rowObsS2(shuffledFile1Arr, obsRow, permutations1, numStates), expFreqArr).sum(axis=0)
            nullScoreArr2[scoreRow] = klScoreND(rowObsS2(shuffledFile2Arr, obsRow, permutations2, numStates), expFreqArr).sum(axis=0)

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time() - tScore, flush=True)
    
# Helper for calculating observed in s2 case
def rowObsS2(dataArr, row, permutations, numStates):
    # (Within a row, how many ways can you choose x and y to be together) / (how many ordered ways can you choose 2 states) = Prob of choosing x and y
    # Can choose x and y to be together x*y ways if different and n(n-1) ways if same (where n is the number of times that x/y shows up)
    rowObsArr = np.zeros((numStates, numStates))
    uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
    for i, state1 in enumerate(uniqueStates):
        for j, state2 in enumerate(uniqueStates):
            if state1 == state2:
                rowObsArr[state1, state2] = stateCounts[i] * (stateCounts[i] - 1) / permutations # Equates to statecounts[i] permute 2 / permutations
            else: # state1 > state2 or state1 < state2
                rowObsArr[state1, state2] = stateCounts[i] * stateCounts[j] / permutations
    return rowObsArr

# Helper for the multiprocessing implemented in s3
def s3Score(file1Path, rowsToCalculate, expFreqPath, verbose):
    dataArr = readStates(file1Path=file1Path, rowsToCalculate=rowsToCalculate, verbose=verbose)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    multiprocessRows, numCols = dataArr.shape
    numStates = sharedArr[2]

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(permutations(range(numCols), 2)), dtype=np.int16).T

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates), dtype=np.float32) / (numCols * (numCols - 1)), expFreqArr)

    if verbose and rowsToCalculate[0] == 0:
        print("Calculating Scores...", flush=True)
        tScore = time()
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

    if verbose and rowsToCalculate[0] == 0: print("    Time:", time() - tScore, flush=True)

# Helper for writing when we are working with real data
def writeScores(dataArr, outputTxtPath, chrName):
    outputTxt = gzip.open(outputTxtPath, "wt")

    # Create a location array
    numRows = dataArr.shape[0]
    numStates = dataArr.shape[1]
    locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(numRows)])

    # Creating a string to write out the data (faster than np.savetxt)
    outputTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates - 1)) + "{1[%d]:.5f}\n" % (numStates - 1)
    outputStr = "".join(outputTemplate.format(locationArr[i], dataArr[i]) for i in range(numRows))

    # Write out the string
    outputTxt.write(outputStr)
    outputTxt.close()

# Helper to calculate KL-score for Nd arrays
def klScoreND(obs, exp):
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)

if __name__ == "__main__":
    main(argv[1], argv[2], int(argv[3]), int(argv[4]), argv[5], argv[6], argv[7], int(argv[8]), strToBool(argv[9]))