import gzip
from time import time
import pandas as pd
import numpy as np
from pathlib import Path


def getNumStates(stateFile):
    """
    Input:
    stateFile -- tsv containing information on all the states in the state model

    Output:
    Number of states
    """
    return pd.read_table(Path(stateFile), header=0, sep="\t").shape[0]


def getStateNames(stateFile):
    """
    Input:
    stateFile -- tsv containing information on all the states in the state model

    Output:
    Shorthand names for all the states
    """
    return pd.read_table(Path(stateFile), header=0, sep="\t")['short_name'].values


def getStateColorsRGB(stateFile):
    """
    Input:
    stateFile -- tsv containing information on all the states in the state model

    Output:
    n x 3 numpy array of rgb values of all the states
    """
    series = pd.read_table(Path(stateFile), header=0, sep="\t")['rgba']
    rgbList = []
    for rgba in series:
        rgba = rgba.split("rgba(")[-1].split(",")
        rgbList.append((int(rgba[0]) / 255, int(rgba[1]) / 255, int(rgba[2]) / 255))
    return np.array(rgbList)


def strToBool(string):
    """
    Input:
    string -- a string with the value either 'True' or 'False'

    Output:
    A boolean equivalent to the string
    """
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        raise ValueError("Invalid boolean string")


def blocks(file, size=65536):
    """
    Helper which speeds determining the number of lines in a file by splitting it into blocks

    Input:
    file -- the file to block
    size -- the size of the blocks [default=65536]

    Output:
    Each of the blocks of the file
    """
    while True:
        b = file.read(size)
        if not b: break
        yield b


def splitRows(dataFilePath, numProcesses):
    """
    Determines what rows to assign to each core

    Input:
    dataFilePath -- the data file which we want to chunk up and send to seperate cores
    numProcesses -- the number of cores which we want to work over

    Output:
    A list of tuples which contain the first and last rows for each core to use
    """
    # Calculate the number of rows in the input file
    if dataFilePath.name.endswith("gz"):
        with gzip.open(dataFilePath, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))
    else:
        with open(dataFilePath, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))

    # Split the rows up according to the number of cores we have available
    # Row list contain tuples of the first and last rows to be assigned to each core
    rowList = []
    for i in range(numProcesses):
        rowsToCalc = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
        rowList.append(rowsToCalc)

    return rowList


def readStates(file1Path=Path("null"), file2Path=Path("null"), rowsToCalc=(0, 0), expBool=True, verbose=True, groupSize=-1):
    """
    Reads the states from the relevant rows of the inputed data file(s)

    Input:
    file1Path -- The path to the first file to read from (used in both single and paired epilogos) [default=Path("null")]
    file2Path -- The path to the second file to read from (used only in paired epilogos) [default=Path("null")]
    rowsToCalc -- The first and last rows to read from the files [default=(0, 0)]
    expBool -- Tells us if we are calculating the expected frequencies or the scores [default=True]
    verbose -- If True, we print out updates [default=True]
    groupSize -- When returning null output, return 2 evenly sized groups of this size

    Output:
        Single Epilogos:
            file1Arr -- 2d numpy array of the state info for each epigenome in the data file over the relevant rows
        Paired Epilogos Expected Frequency:
            combinedArr -- 2d numpy array of the state info for each epigenome of both data files over the relevant rows
        Paired Epilogos Scores:
            file1Arr -- 2d numpy array of the state info for each epigenome in the first data file over the relevant rows
            file2Arr -- 2d numpy array of the state info for each epigenome in the second data file over the relevant rows
            shuffledCombinedArr[:, :file1Arr.shape[1]] -- 2d numpy array of the state info shuffled between first and second
                                                          data files, with the same width as file1Arr
            shuffledCombinedArr[:, file1Arr.shape[1]:] -- 2d numpy array of the state info shuffled between first and second
                                                          data files, with the same width as file2Arr
    """
    # Read in the data from file 1
    if verbose and rowsToCalc[0] == 0: print("Reading data from file 1...", flush=True); tRead1 = time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file1Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file1Arr = pd.read_table(file1Path, usecols=cols, skiprows=rowsToCalc[0],
                             nrows=rowsToCalc[1]-rowsToCalc[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalc[0] == 0: print("    Time: ", time() - tRead1, flush=True)

    # If we are doing single group epilogos, just return the data array, otherwise, we continute to the second file
    if str(file2Path) == "null":
        return file1Arr

    # Read in the data from file 2
    if verbose and rowsToCalc[0] == 0: print("Reading data from file 2...", flush=True); tRead2 = time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file2Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file2Arr = pd.read_table(file2Path, usecols=cols, skiprows=rowsToCalc[0],
                             nrows=rowsToCalc[1]-rowsToCalc[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalc[0] == 0: print("    Time: ", time() - tRead2, flush=True)

    # Combining the arrays for per row shuffling
    if verbose and rowsToCalc[0] == 0: print("Combining input matrices..."); tCombine = time()
    combinedArr = np.concatenate((file1Arr, file2Arr), axis=1)
    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tCombine, flush=True)

    # If we are calculating the expected frequencies for pairwise epilogos,
    # we can just return the combined array of file 1 and file 2
    if expBool:
        return combinedArr

    # Row independent vectorized shuffling of the 2 arrays
    if verbose and rowsToCalc[0] == 0: print("Shuffling input matrices...", flush=True); tShuffle = time()
    randomIndices = np.argsort(np.random.rand(*combinedArr.shape), axis=1)
    shuffledCombinedArr = np.take_along_axis(combinedArr, randomIndices, axis=1)
    if verbose and rowsToCalc[0] == 0: print("    Time:", time() - tShuffle, flush=True)

    # In the case of calculating the scores for pairwise epilogos,
    # we need the original file 1 and file 2 arrays as well as their shuffled counterparts
    # shuffledCombinedArr is split by size of the original arrays
    if groupSize == -1:
        return file1Arr, file2Arr, shuffledCombinedArr[:, :file1Arr.shape[1]], shuffledCombinedArr[:, file1Arr.shape[1]:]
    else:
        return file1Arr, file2Arr, shuffledCombinedArr[:, :groupSize], shuffledCombinedArr[:, groupSize:2*groupSize]