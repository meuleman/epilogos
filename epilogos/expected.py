import numpy as np
from sys import argv
from pathlib import Path
from time import time
from multiprocessing import cpu_count, Pool
from itertools import repeat, permutations
from contextlib import closing
from epilogos.helpers import strToBool, countRows, splitRows, readStates


def main(file1, file2, numStates, saliency, outputDir, fileTag, numProcesses, verbose):
    """
    Wrapper function which prepares inputs for the expected frequency calculation

    Input:
    file1        -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2        -- The path of the second file to read states from (paired epilogos)
    numStates    -- The number of states in the state model
    saliency     -- The saliency metric being used in the epilogos run
    outputDir    -- The path of the output directory
    fileTag      -- A string which helps ensure outputed files are named similarly within an epilogos run
    numProcesses -- The number of cores to run on
    verbose      -- Boolean which if True, causes much more detailed prints
    """
    if verbose: tTotal = time()

    file1Path = Path(file1)
    file2Path = Path(file2)
    outputDirPath = Path(outputDir)

    filename = file1Path.name.split(".")[0]

    if not verbose: print("    {}\t".format(filename), end="", flush=True)

    # If user doesn't want to choose number of cores, use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    # Determine which rows to assign to each core
    rowList = splitRows(countRows(file1Path), numProcesses)

    calculateExpected(saliency, file1Path, file2Path, rowList, numStates, outputDirPath, fileTag, filename,
                      numProcesses, verbose)

    print("Total Time:", time() - tTotal, flush=True) if verbose else print("\t[Done]", flush=True)


def calculateExpected(saliency, file1Path, file2Path, rowList, numStates, outputDirPath, fileTag, filename,
                      numProcesses, verbose):
    """
    Function responsible for deploying the processes used to calculate the expected frequencies

    Input:
    saliency      -- The saliency metric being used in the epilogos run
    file1Path     -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2Path     -- The path of the second file to read states from (paired epilogos)
    rowList       -- A list of tuples which contain the first and last rows for each core to use
    numStates     -- The number of states in the state model
    outputDirPath -- The path of the output directory
    fileTag       -- A string which helps ensure outputed files are named similarly within an epilogos run
    filename      -- The name of the file for which we are calculating the expected frequencies
    numProcesses  -- The number of cores to run on
    verbose       -- Boolean which if True, causes much more detailed prints

    Output:
    A .npy array containing the expected frequencies of the input file
    """
    if verbose: print("\nNumber of Processes:", numProcesses, flush=True)

    # Start the processes
    with closing(Pool(numProcesses)) as pool:
        if saliency == 1:
            results = pool.starmap(s1Calc, zip(repeat(file1Path), repeat(file2Path), rowList, repeat(numStates),
                                               repeat(verbose)))
        elif saliency == 2:
            results = pool.starmap(s2Calc, zip(repeat(file1Path), repeat(file2Path), rowList, repeat(numStates),
                                               repeat(verbose)))
        elif saliency == 3:
            results = pool.starmap(s3Calc, zip(repeat(file1Path), rowList, repeat(numStates), repeat(verbose)))
        else:
            raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")
    pool.join()

    # Sum all the expected frequency arrays from the seperate processes and normalize by dividing by numRows
    expFreqArr = np.sum(results, axis=0)

    storeExpArray(expFreqArr, outputDirPath, fileTag, filename)


def s1Calc(file1Path, file2Path, rowsToCalc, numStates, verbose):
    """
    Function responsible for expected frequency calculation over a set of rows for a saliency metric of 1

    Input:
    file1Path  -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2Path  -- The path of the second file to read states from (paired epilogos)
    rowsToCalc -- The rows to count expected frequencies from the files
    numStates  -- The number of states in the state model
    verbose    -- Boolean which if True, causes much more detailed prints

    Output:
    A numpy array containing the counts of each state within the specified rows of the file
    """
    dataArr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalc=rowsToCalc, verbose=verbose)

    expFreqArr = np.zeros(numStates, dtype=np.int64)

    if verbose and rowsToCalc[0] == 0: print("Calculating expected frequencies...", flush=True); tExp = time()

    # Simply count all states across out our subset of data
    uniqueStates, stateCounts = np.unique(dataArr, return_counts=True)
    for i, state in enumerate(uniqueStates):
        expFreqArr[state] += stateCounts[i]
    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tExp, flush=True)

    return expFreqArr


def s2Calc(file1Path, file2Path, rowsToCalc, numStates, verbose):
    """
    Function responsible for expected frequency calculation over a set of rows for a saliency metric of 2

    Input:
    file1Path  -- The path of the only (single epilogos) or first (paired epilogos) file to read states from
    file2Path  -- The path of the second file to read states from (paired epilogos)
    rowsToCalc -- The rows to count expected frequencies from the files
    numStates  -- The number of states in the state model
    verbose    -- Boolean which if True, causes much more detailed prints

    Output:
    A numpy array containing the counts of each pair of states within the specified rows of the file
    """
    dataArr = readStates(file1Path=file1Path, file2Path=file2Path, rowsToCalc=rowsToCalc, verbose=verbose)

    multiprocessRows = dataArr.shape[0]

    expFreqArr = np.zeros((numStates, numStates), dtype=np.int64)

    if verbose and rowsToCalc[0] == 0:
        print("Calculating expected frequencies...", flush=True); tExp = time(); percentDone = 0
    printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]

    # SumOverRows: Within a row, how many ways can you choose x and y to be together (will normalize later)
    # Can choose x and y to be together n*m ways if n != m and n(n-1) ways if n == m
    # (where n and m are the number of times that x and y show up respectively)
    for row in range(multiprocessRows):

        if verbose and rowsToCalc[0] == 0 and row in printCheckmarks: percentDone += 10; \
            print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalc[0] == 0 and row in printCheckmarks: print(".", end="", flush=True)

        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
        for i, state1 in enumerate(uniqueStates):
            for j, state2 in enumerate(uniqueStates):
                if state1 == state2:
                    expFreqArr[state1, state2] += stateCounts[i] * (stateCounts[j] - 1)
                else:  # state1 > state2 or state1 < state2
                    expFreqArr[state1, state2] += stateCounts[i] * stateCounts[j]

    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tExp, flush=True)

    return expFreqArr


def s3Calc(file1Path, rowsToCalc, numStates, verbose):
    """
    Function responsible for expected frequency calculation over a set of rows for a saliency metric of 3

    Input:
    file1Path  -- The path of the only file to read states from
    rowsToCalc -- The first and last rows to count expected frequencies from the file
    numStates  -- The number of states in the state model
    verbose    -- Boolean which if True, causes much more detailed prints

    Output:
    A numpy array containing the counts of each grouping of states and epigenomes within the specified rows of the file
    """
    dataArr = readStates(file1Path=file1Path, rowsToCalc=rowsToCalc, verbose=verbose)

    multiprocessRows, numCols = dataArr.shape

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(permutations(range(numCols), 2))).T

    expFreqArr = np.zeros((numCols, numCols, numStates, numStates), dtype=np.int32)

    if verbose and rowsToCalc[0] == 0:
        print("Calculating expected frequencies...", flush=True); tExp = time(); percentDone = 0
    printCheckmarks = [int(multiprocessRows * float(i / 10)) for i in range(1, 10)]

    # We tally a one for all the state/column combinations we observe
    # (e.g. for state 18 in column 2 and state 15 in column 6 we would add one to index [5, 1, 17, 14])
    for row in range(multiprocessRows):

        if verbose and rowsToCalc[0] == 0 and row in printCheckmarks: percentDone += 10; \
            print("    {}% Completed".format(percentDone), flush=True)
        if not verbose and rowsToCalc[0] == 0 and row in printCheckmarks: print(".", end="", flush=True)

        expFreqArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]],
                   dataArr[row, basePermutationArr[1]]] += 1

    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tExp, flush=True)

    return expFreqArr


def storeExpArray(expFreqArr, outputDirPath, fileTag, filename):
    """
    Stores the expected frequency array of the file as a temporary .npy file

    Input:
    expFreqArr    -- Numpy array of the expected frequencies
    outputDirPath -- The output directory for epilogos
    fileTag       -- A string which helps ensure outputed files are named similarly within an epilogos run
    filename      -- The name of the file for which we calculated the expected frequencies

    Output:
    A .npy array containing the expected frequencies of the input file
    """
    expFreqFilename = "temp_exp_freq_{}_{}.npy".format(fileTag, filename)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)


if __name__ == "__main__":
    main(argv[1], argv[2], int(argv[3]), int(argv[4]), argv[5], argv[6], int(argv[7]), strToBool(argv[8]))
