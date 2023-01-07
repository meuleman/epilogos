import numpy as np
from sys import argv
from pathlib import Path
import pandas as pd
from time import time
import numpy.ma as ma
from multiprocessing import cpu_count, Pool, RawArray
from itertools import repeat, permutations
from contextlib import closing
from epilogos.helpers import strToBool, countRows, splitRows, readStates, sharedToNumpy
import gzip


def main(file1, file2, numStates, saliency, outputDir, expFreqPath, fileTag, numProcesses, quiescentState, groupSize,
         verbose):
    """
    Wrapper function which prepares inputs for the score calculation

    Input:
    file1          -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2          -- The path of the second file to read states from (paired epilogos)
    numStates      -- The number of states in the state model
    saliency       -- The saliency metric being used in the epilogos run
    outputDir      -- The path of the output directory
    expFreqPath    -- The path of the expected frequency array
    fileTag        -- A string which helps ensure outputed files are named similarly within an epilogos run
    numProcesses   -- The number of cores to run on
    quiescentState -- The state used to filter out quiescent bins
    groupSize      -- Size of the outputted null array
    verbose        -- Boolean which if True, causes much more detailed prints
    """
    if verbose: tTotal = time()

    file1Path = Path(file1)
    file2Path = Path(file2)
    outputDirPath = Path(outputDir)

    filename = file1Path.name.split(".")[0]

    if not verbose: print("    {}\t".format(filename), end="", flush=True)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    # Determine which rows to assign to each core
    rowList = splitRows(countRows(file1Path), numProcesses)

    if file2 == "null":
        calculateScores(saliency, file1Path, rowList, numStates, outputDirPath, expFreqPath, fileTag, filename,
                        numProcesses, verbose)
    else:
        calculateScoresPairwise(saliency, file1Path, file2Path, rowList, numStates, outputDirPath, expFreqPath, fileTag,
                                filename, numProcesses, quiescentState, groupSize, verbose)

    print("Total Time:", time() - tTotal, flush=True) if verbose else print("\t[Done]", flush=True)


def _init(sharedArr_, expFreqPath_, verbose_):
    """
    Initializes global variables for multiprocessing in the single epilogos case

    Input:
    sharedArr_   -- A tuple containing relevant information about the shared score array
    expFreqPath_ -- A pathlib path to the expected frequency array
    verbose_     -- A boolean which tells us the amount we need to print
    """
    global sharedArr
    global expFreqPath
    global verbose

    sharedArr = sharedArr_
    expFreqPath = expFreqPath_
    verbose = verbose_


def _initPairwise(sharedArr1_, sharedArr2_, shuffledSharedArr1_, shuffledSharedArr2_, quiescenceSharedArr_, totalRows,
                  numStates, quiescentState_, expFreqPath_, groupSize_, verbose_):
    """
    Initializes global variables for multiprocessing in the paired epilogos case

    Input:
    sharedArr1_          -- The first shared score array
    sharedArr2_          -- The second shared score array
    shuffledSharedArr1_  -- The first shared null score array
    shuffledSharedArr2_  -- The second shared null score array
    quiescenceSharedArr_ -- Shared array containing T/F values for whether a bin is quiescent
    totalRows            -- The number of rows of the input files
    numStates            -- The number of states in the state model
    quiescentState_      -- The state used to filter out quiescent bins
    expFreqPath_         -- A pathlib path to the expected frequency array
    groupSize_           -- Size of outputed null array
    verbose_             -- A boolean which tells us the amount we need to print
    """
    global sharedArr1
    global sharedArr2
    global shuffledSharedArr1
    global shuffledSharedArr2
    global quiescenceSharedArr
    global quiescentState
    global expFreqPath
    global verbose
    global groupSize

    sharedArr1 = (sharedArr1_, totalRows, numStates)
    sharedArr2 = (sharedArr2_, totalRows, numStates)
    shuffledSharedArr1 = (shuffledSharedArr1_, totalRows, numStates)
    shuffledSharedArr2 = (shuffledSharedArr2_, totalRows, numStates)
    quiescenceSharedArr = quiescenceSharedArr_
    quiescentState = quiescentState_
    expFreqPath = expFreqPath_
    verbose = verbose_
    groupSize = groupSize_


def calculateScores(saliency, file1Path, rowList, numStates, outputDirPath, expFreqPath, fileTag, filename,
                    numProcesses, verbose):
    """
    Function responsible for deploying the processes used to calculate the scores in the single epilogos case

    Input:
    saliency      -- The saliency metric being used in the epilogos run
    file1Path     -- The path of the only file to read states from
    rowList       -- A list of tuples which contain the first and last rows for each core to use
    numStates     -- The number of states in the state model
    outputDirPath -- The path of the output directory
    expFreqPath   -- The path to the expected frequency array
    fileTag       -- A string which helps ensure outputed files are named similarly within an epilogos run
    filename      -- The name of the file for which we are calculating the expected frequencies
    numProcesses  -- The number of cores to run on
    verbose       -- Boolean which if True, causes much more detailed prints

    Output:
    scores*.txt.gz file containing the scores for the input file
    temp_scores*.npz file containing the scores as a npz file for fast reading later
    """

    if verbose: print("\nNumber of Processes:", numProcesses, flush=True)

    totalRows = rowList[-1][-1]

    # Shared arrays across the multiple processes for the scores
    # We avoid race conditions by writing to separate parts of the array in each process
    sharedArr = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)

    # Start the processes
    with closing(Pool(numProcesses, initializer=_init,
                      initargs=((sharedArr, totalRows, numStates), expFreqPath, verbose))) as pool:
        if saliency == 1:
            pool.starmap(s1Score, zip(repeat(file1Path), repeat(Path("null")), rowList))
        elif saliency == 2:
            pool.starmap(s2Score, zip(repeat(file1Path), repeat(Path("null")), rowList))
        elif saliency == 3:
            pool.starmap(s3Score, zip(repeat(file1Path), rowList))
        else:
            raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")
    pool.join()

    # chrName = pd.read_table(file1Path, nrows=1, header=None, sep="\t").iloc[0, 0]
    # locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(totalRows)])
    locationArr = pd.read_table(file1Path, header=None, sep="\t", usecols=[0, 1, 2]).to_numpy()

    outputTxtPath = outputDirPath / "scores_{}_{}.txt.gz".format(fileTag, filename)
    writeScores(sharedToNumpy(sharedArr, totalRows, numStates), outputTxtPath, locationArr)
    # Temporary scores for regions of interest (np files are faster to read in)
    chrName = pd.read_table(file1Path, nrows=1, header=None, usecols=[0], sep="\t").iloc[0, 0]
    tempScoresPath = outputDirPath / "temp_scores_{}_{}.npz".format(fileTag, filename)
    np.savez_compressed(tempScoresPath, chrName=np.array([chrName]),
                        scoreArr=sharedToNumpy(sharedArr, totalRows, numStates), locationArr=locationArr)


def calculateScoresPairwise(saliency, file1Path, file2Path, rowList, numStates, outputDirPath, expFreqPath, fileTag,
                            filename, numProcesses, quiescentState, groupSize, verbose):
    """
    Function responsible for deploying the processes used to calculate the scores in the paired epilogos case

    Input:
    saliency       -- The saliency metric being used in the epilogos run
    file1Path      -- The path of the first file to read states from
    file2Path      -- The path of the second file to read states from
    rowList        -- A list of tuples which contain the first and last rows for each core to use
    numStates      -- The number of states in the state model
    outputDirPath  -- The path of the output directory
    expFreqPath    -- The path to the expected frequency array
    fileTag        -- A string which helps ensure outputed files are named similarly within an epilogos run
    filename       -- The name of the file for which we are calculating the expected frequencies
    numProcesses   -- The number of cores to run on
    quiescentState -- The state used to filter out quiescent bins
    groupSize      -- Size of the outputted null array
    verbose        -- Boolean which if True, causes much more detailed prints

    Output:
    pairwiseDelta*.txt.gz file containing the differences between the scores of input files
    temp_nullDistances*.npz file containing the null distances as a npz file for fast reading later
    temp_quiescence*.npz file containing booleans identifying bins which are completely quiescent
    """
    if verbose: print("\nNumber of Processes:", numProcesses, flush=True)

    totalRows = rowList[-1][-1]

    # Shared arrays across the multiple processes for the scores
    # We avoid race conditions by writing to separate parts of the array in each process
    sharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    sharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr1 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    shuffledSharedArr2 = RawArray(np.ctypeslib.as_ctypes_type(np.float32), totalRows * numStates)
    quiescenceSharedArr = RawArray(np.ctypeslib.as_ctypes_type(np.bool_), totalRows)

    # Start the processes
    with closing(Pool(numProcesses, initializer=_initPairwise,
                      initargs=(sharedArr1, sharedArr2, shuffledSharedArr1, shuffledSharedArr2, quiescenceSharedArr,
                                totalRows, numStates, quiescentState, expFreqPath, groupSize, verbose))) as pool:
        if saliency == 1:
            pool.starmap(s1Score, zip(repeat(file1Path), repeat(file2Path), rowList))
        elif saliency == 2:
            pool.starmap(s2Score, zip(repeat(file1Path), repeat(file2Path), rowList))
        else:
            raise ValueError("Please ensure that saliency metric is either 1 or 2 for Pairwise Epilogos")
    pool.join()

    # Calculate the differences between array 1 and 2 in both the real and null case
    if verbose: print("Calculating Raw Differences...", flush=True); tDiff = time()
    realDiffArr = sharedToNumpy(sharedArr1, totalRows, numStates) - sharedToNumpy(sharedArr2, totalRows, numStates)
    nullDiffArr = sharedToNumpy(shuffledSharedArr1, totalRows, numStates)\
        - sharedToNumpy(shuffledSharedArr2, totalRows, numStates)
    if verbose: print("    Time:", time() - tDiff, flush=True)

    # Only calculate the distances for the null data in this step
    if verbose: print("Calculating Squared Euclidean Distance and Maximum Contributing Difference...", flush=True); \
        tDistance = time()
    diffSign = np.sign(np.sum(nullDiffArr, axis=1))
    nullDistancesArr = np.sum(np.square(nullDiffArr), axis=1) * diffSign
    if verbose: print("    Time:", time() - tDistance, flush=True)

    # If it's the real data, we will just write the delta and calculate metrics in computeEpilogosPairwiseVisual
    # If it's the null data, we will just save the signed squared euclidean distances as a temporary npz file as we
    # don't care about saving the data
    # We also want to write the locations of any quiescent bins (where all states are in the quiescent state).
    # This ensures that our eventual fit on the null data is more accurate.
    if verbose: print("Writing output to disk...", flush=True); tWrite = time()

    # chrName = pd.read_table(file1Path, nrows=1, header=None, sep="\t").iloc[0, 0]
    # locationArr = np.array([[chrName, 200*i, 200*i+200] for i in range(totalRows)])
    locationArr = pd.read_table(file1Path, header=None, sep="\t", usecols=[0, 1, 2]).to_numpy()

    realOutputPath = outputDirPath / "pairwiseDelta_{}_{}.txt.gz".format(fileTag, filename)
    writeScores(realDiffArr, realOutputPath, locationArr)

    chrName = pd.read_table(file1Path, nrows=1, header=None, usecols=[0], sep="\t").iloc[0, 0]

    nullOutputPath = outputDirPath / "temp_nullDistances_{}_{}.npz".format(fileTag, filename)
    np.savez_compressed(nullOutputPath, chrName=np.array([chrName]), nullDistances=nullDistancesArr)
    quiescentOutputPath = outputDirPath / "temp_quiescence_{}_{}.npz".format(fileTag, filename)
    np.savez_compressed(quiescentOutputPath, chrName=np.array([chrName]),
                        quiescenceArr=np.frombuffer(quiescenceSharedArr, dtype=np.bool_))
    if verbose: print("    Time:", time() - tWrite, flush=True)


def s1Score(file1Path, file2Path, rowsToCalc):
    """
    Function responsible for score calculation over a set of rows for a saliency metric of 1. Note that there are global
    variables used for the shared score array(s). There is no output as scores are put directly into the shared score
    array(s)

    Input:
    file1Path  -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2Path  -- The path of the second file to read states from (paired epilogos)
    rowsToCalc -- The rows to count expected frequencies from the files
    """
    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    # Loading the data and creating the shared arrays for the scores
    if str(file2Path) == "null":
        dataArr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalc=rowsToCalc, expBool=False,
                             verbose=verbose)

        numStates = sharedArr[2]

        scoreArr = sharedToNumpy(*sharedArr)
    else:
        file1Arr, file2Arr, shuffledFile1Arr, shuffledFile2Arr = readStates(file1Path=file1Path, file2Path=file2Path,
                                                                            rowsToCalc=rowsToCalc, expBool=False,
                                                                            verbose=verbose, groupSize=groupSize)

        numStates = sharedArr1[2]

        realScoreArr1 = sharedToNumpy(*sharedArr1)
        realScoreArr2 = sharedToNumpy(*sharedArr2)
        nullScoreArr1 = sharedToNumpy(*shuffledSharedArr1)
        nullScoreArr2 = sharedToNumpy(*shuffledSharedArr2)
        quiescenceArr = np.frombuffer(quiescenceSharedArr, dtype=np.bool_)

        if quiescentState != -1:
            # Figure out the quiescent bins
            sortedArr1 = np.sort(file1Arr, axis=1)
            sortedArr2 = np.sort(file2Arr, axis=1)
            # If all values in a bin are equal to the quiescentState in both file1Arr and file2Arr,
            # then the bin is quiescent
            quiescenceArr[rowsToCalc[0]:rowsToCalc[1]][np.where((sortedArr1[:, 0] == quiescentState)
                                                                & (sortedArr1[:, -1] == quiescentState)
                                                                & (sortedArr2[:, 0] == quiescentState)
                                                                & (sortedArr2[:, -1] == quiescentState))] = True

    if verbose and rowsToCalc[0] == 0: print("Calculating Scores...", flush=True); tScore = time(); percentDone = 0
    printCheckmarks = [int(rowsToCalc[1] * float(i / 10)) for i in range(1, 10)]

    # Calculate the observed frequencies and final scores for the designated rows
    for obsRow, scoreRow in enumerate(range(rowsToCalc[0], rowsToCalc[1])):

        if verbose and rowsToCalc[0] == 0 and obsRow in printCheckmarks: percentDone += 10; \
            print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalc[0] == 0 and obsRow in printCheckmarks: print(".", end="", flush=True)

        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        if str(file2Path) == "null":
            scoreArr[scoreRow] = klScoreND(rowObsS1(dataArr, obsRow, numStates), expFreqArr)
        else:
            realScoreArr1[scoreRow] = klScoreND(rowObsS1(file1Arr, obsRow, numStates), expFreqArr)
            realScoreArr2[scoreRow] = klScoreND(rowObsS1(file2Arr, obsRow, numStates), expFreqArr)
            nullScoreArr1[scoreRow] = klScoreND(rowObsS1(shuffledFile1Arr, obsRow, numStates), expFreqArr)
            nullScoreArr2[scoreRow] = klScoreND(rowObsS1(shuffledFile2Arr, obsRow, numStates), expFreqArr)

    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tScore, flush=True)


def rowObsS1(dataArr, row, numStates):
    """
    Calculates the observed counts for each state for a saliency metric of 1

    Input:
    dataArr   -- The numpy array which contains the states to count
    row       -- The row from which to calculate the observed counts
    numStates -- The number of states in the state model

    Output:
    rowObsArr -- Observed counts for each state in a row for a saliency metric of 1
    """
    rowObsArr = np.zeros(numStates)
    # Count number of each state in row and return array
    uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
    for i, state in enumerate(uniqueStates):
        rowObsArr[state] = stateCounts[i] / dataArr.shape[1]
    return rowObsArr


def s2Score(file1Path, file2Path, rowsToCalc):
    """
    Function responsible for score calculation over a set of rows for a saliency metric of 2. Note that there are global
    variables used for the shared score array(s). There is no output as scores are put directly into the shared score
    array(s)

    Input:
    file1Path  -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2Path  -- The path of the second file to read states from (paired epilogos)
    rowsToCalc -- The rows to count expected frequencies from the files
    """
    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    # Loading the data and creating the shared arrays for the scores
    if str(file2Path) == "null":
        dataArr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalc=rowsToCalc, expBool=False,
                             verbose=verbose)

        numStates = sharedArr[2]

        scoreArr = sharedToNumpy(*sharedArr)

        # Need the permuations to effective count state pairs (see rowObsS2() for theory)
        permutations = dataArr.shape[1] * (dataArr.shape[1] - 1)
    else:
        file1Arr, file2Arr, shuffledFile1Arr, shuffledFile2Arr = readStates(file1Path=file1Path, file2Path=file2Path,
                                                                            rowsToCalc=rowsToCalc, expBool=False,
                                                                            verbose=verbose, groupSize=groupSize)

        numStates = sharedArr1[2]

        realScoreArr1 = sharedToNumpy(*sharedArr1)
        realScoreArr2 = sharedToNumpy(*sharedArr2)
        nullScoreArr1 = sharedToNumpy(*shuffledSharedArr1)
        nullScoreArr2 = sharedToNumpy(*shuffledSharedArr2)
        quiescenceArr = np.frombuffer(quiescenceSharedArr, dtype=np.bool_)

        if quiescentState != -1:
            # Figure out the quiescent bins
            sortedArr1 = np.sort(file1Arr, axis=1)
            sortedArr2 = np.sort(file2Arr, axis=1)
            # If all values in a bin are equal to the quiescentState in both file1Arr and file2Arr, then the bin is
            # quiescent
            quiescenceArr[rowsToCalc[0]:rowsToCalc[1]][np.where((sortedArr1[:, 0] == quiescentState)
                                                                & (sortedArr1[:, -1] == quiescentState)
                                                                & (sortedArr2[:, 0] == quiescentState)
                                                                & (sortedArr2[:, -1] == quiescentState))] = True

        # Need the permuations to effective count state pairs (see rowObsS2() for theory)
        permutations1 = file1Arr.shape[1] * (file1Arr.shape[1] - 1)
        permutations2 = file2Arr.shape[1] * (file2Arr.shape[1] - 1)

    if verbose and rowsToCalc[0] == 0: print("Calculating Scores...", flush=True); tScore = time(); percentDone = 0
    printCheckmarks = [int(rowsToCalc[1] * float(i / 10)) for i in range(1, 10)]

    # Find scores for each row that the core is responsible for
    for obsRow, scoreRow in enumerate(range(rowsToCalc[0], rowsToCalc[1])):

        if verbose and rowsToCalc[0] == 0 and obsRow in printCheckmarks: percentDone += 10; \
            print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalc[0] == 0 and obsRow in printCheckmarks: print(".", end="", flush=True)

        # Inputs to klScoreND are obsFreqArr and expFreqArr respectively
        if str(file2Path) == "null":
            scoreArr[scoreRow] = klScoreND(rowObsS2(dataArr, obsRow, permutations, numStates), expFreqArr).sum(axis=0)
        else:
            realScoreArr1[scoreRow] = klScoreND(rowObsS2(file1Arr, obsRow, permutations1, numStates),
                                                expFreqArr).sum(axis=0)
            realScoreArr2[scoreRow] = klScoreND(rowObsS2(file2Arr, obsRow, permutations2, numStates),
                                                expFreqArr).sum(axis=0)
            nullScoreArr1[scoreRow] = klScoreND(rowObsS2(shuffledFile1Arr, obsRow, permutations1, numStates),
                                                expFreqArr).sum(axis=0)
            nullScoreArr2[scoreRow] = klScoreND(rowObsS2(shuffledFile2Arr, obsRow, permutations2, numStates),
                                                expFreqArr).sum(axis=0)

    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tScore, flush=True)


def rowObsS2(dataArr, row, permutations, numStates):
    """
    Calculates the observed counts for each pair of states for a saliency metric of 2

    Input:
    dataArr      -- The numpy array which contains the states to count
    row          -- The row from which to calculate the observed counts
    permutations -- The number of permutations of epigenomes possible within the row of the dataArr
    numStates    -- The number of states in the state model

    Output:
    rowObsArr -- Observed counts for each pair of states within a bin for a saliency metric of 2
    """
    # (Within a row, # of ways can you choose x and y to be together) / (# of ordered ways can you choose 2 states)
    #  = Prob of choosing x and y
    # Can choose x and y to be together x*y ways if different and n(n-1) ways if same
    # (where n is the number of times that x/y shows up)
    rowObsArr = np.zeros((numStates, numStates))
    uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
    for i, state1 in enumerate(uniqueStates):
        for j, state2 in enumerate(uniqueStates):
            if state1 == state2:
                # Equates to statecounts[i] permute 2 / permutations
                rowObsArr[state1, state2] = stateCounts[i] * (stateCounts[j] - 1) / permutations
            else:  # state1 > state2 or state1 < state2
                rowObsArr[state1, state2] = stateCounts[i] * stateCounts[j] / permutations
    return rowObsArr


def s3Score(file1Path, rowsToCalc):
    """
    Function responsible for score calculation over a set of rows for a saliency metric of 3. Note that there are global
    variables used for the shared score array(s). There is no output as scores are put directly into the shared score
    array(s)

    Input:
    file1Path  -- The path of the only file to read states from
    rowsToCalc -- The rows to count expected frequencies from the files
    """
    dataArr = readStates(file1Path=file1Path, rowsToCalc=rowsToCalc, verbose=verbose)

    # Loading the expected frequency array
    expFreqArr = np.load(expFreqPath, allow_pickle=False)

    numCols = dataArr.shape[1]
    numStates = sharedArr[2]

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(permutations(range(numCols), 2)), dtype=np.int16).T

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the
    # scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates), dtype=np.float32)
                             / (numCols * (numCols - 1)), expFreqArr)

    if verbose and rowsToCalc[0] == 0: print("Calculating Scores...", flush=True); tScore = time(); percentDone = 0
    printCheckmarks = [int(rowsToCalc[1] * float(i / 10)) for i in range(1, 10)]

    # Calculte the scores and store them in the shared array
    scoreArr = sharedToNumpy(*sharedArr)
    rowScoreArr = np.zeros(numStates, dtype=np.float32)
    for dataRow, scoreRow in enumerate(range(rowsToCalc[0], rowsToCalc[1])):

        if verbose and rowsToCalc[0] == 0 and dataRow in printCheckmarks: percentDone += 10; \
            print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalc[0] == 0 and dataRow in printCheckmarks: print(".", end="", flush=True)

        if dataRow < dataArr.shape[0]:
            # Pull the scores from the precalculated score array add them to the correct index in the rowScoreArr
            np.add.at(rowScoreArr, dataArr[dataRow, basePermutationArr[1]],
                      scoreArrOnes[basePermutationArr[0], basePermutationArr[1],
                                   dataArr[dataRow, basePermutationArr[0]], dataArr[dataRow, basePermutationArr[1]]])

            # Store the scores in the shared score array
            scoreArr[scoreRow] = rowScoreArr

            # Reset the array so it doesn't carry over scores from other rows
            rowScoreArr.fill(0)

    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tScore, flush=True)


def writeScores(dataArr, outputTxtPath, locationArr):
    """
    Writes the scores out to a gziped text file

    Input:
    dataArr -- Numpy array containing the scores to write
    outputTextPath -- The path of the file to write to
    locationArr -- Numpy array containing the genomic locations of the scores

    Output:
    scores*.txt.gz file containing the scores for the input file
    OR
    pairwiseDelta*.txt.gz file containing the differences between the scores of input files
    """
    outputTxt = gzip.open(outputTxtPath, "wt")

    # Create a location array
    numRows = dataArr.shape[0]
    numStates = dataArr.shape[1]

    # Creating a string to write out the data (faster than np.savetxt)
    outputTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates - 1))\
                     + "{1[%d]:.5f}\n" % (numStates - 1)
    outputStr = "".join(outputTemplate.format(locationArr[i], dataArr[i]) for i in range(numRows))

    # Write out the string
    outputTxt.write(outputStr)
    outputTxt.close()


def klScoreND(obs, exp):
    """
    Calculates the kullback leibler scores using an observed and expected numpy array

    Input:
    obs -- Nd numpy array containing observed frequencies
    exp -- Nd numpy array containing expected frequencies

    Output:
    Nd numpy array with the kullback leibler scores
    """
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)


if __name__ == "__main__":
    main(argv[1], argv[2], int(argv[3]), int(argv[4]), argv[5], argv[6], argv[7], int(argv[8]), int(argv[9]),
         int(argv[10]), strToBool(argv[11]))
