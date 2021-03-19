import gzip
from time import time
import pandas as pd
import numpy as np
from pathlib import Path


def getNumStates(stateFile):
    return pd.read_table(Path(stateFile), header=0, sep="\t").shape[0]


def getStateNames(stateFile):
    return pd.read_table(Path(stateFile), header=0, sep="\t")['short_name'].values


def getStateColorsRGB(stateFile):
    series = pd.read_table(Path(stateFile), header=0, sep="\t")['rgba']
    rgbList = []
    for rgba in series:
        rgba = rgba.split("rgba(")[-1].split(",")
        rgbList.append((int(rgba[0]) / 255, int(rgba[1]) / 255, int(rgba[2]) / 255))
    return np.array(rgbList)


# Helper for slurm to send boolean values
def strToBool(string):
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        raise ValueError("Invalid boolean string")


# Helper for reading number of lines in input file
def blocks(files, size=65536):
    while True:
        b = files.read(size)
        if not b: break
        yield b


# Helper to determine which rows to assign to each core
def splitRows(dataFilePath, numProcesses):
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
        rowsToCalculate = (i * totalRows // numProcesses, (i+1) * totalRows // numProcesses)
        rowList.append(rowsToCalculate)

    return rowList


# Helper to read in the states from the datafiles
def readStates(file1Path=Path("null"), file2Path=Path("null"), rowsToCalculate=(0, 0), expBool=True, verbose=True):
    # Read in the data from file 1
    if verbose and rowsToCalculate[0] == 0: print("Reading data from file 1...", flush=True); tRead1 = time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file1Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file1Arr = pd.read_table(file1Path, usecols=cols, skiprows=rowsToCalculate[0], \
        nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time() - tRead1, flush=True)

    # If we are doing single group epilogos, just return the data array, otherwise, we continute to the second file
    if str(file2Path) == "null":
        return file1Arr

    # Read in the data from file 2
    if verbose and rowsToCalculate[0] == 0: print("Reading data from file 2...", flush=True); tRead2 = time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file2Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file2Arr = pd.read_table(file2Path, usecols=cols, skiprows=rowsToCalculate[0], \
        nrows=rowsToCalculate[1]-rowsToCalculate[0], header=None, sep="\t").to_numpy(dtype=int) - 1
    if verbose and rowsToCalculate[0] == 0: print("    Time: ", time() - tRead2, flush=True)

    # Combining the arrays for per row shuffling
    if verbose and rowsToCalculate[0] == 0: print("Combining input matrices..."); tCombine = time()
    combinedArr = np.concatenate((file1Arr, file2Arr), axis=1)
    if verbose and rowsToCalculate[0] == 0: print("    Time:", time() - tCombine, flush=True)

    # If we are calculating the expected frequencies for pairwise epilogos, 
    # we can just return the combined array of file 1 and file 2
    if expBool:
        return combinedArr

    # Row independent vectorized shuffling of the 2 arrays
    if verbose and rowsToCalculate[0] == 0: print("Shuffling input matrices...", flush=True); tShuffle = time()
    randomIndices = np.argsort(np.random.rand(*combinedArr.shape), axis=1)
    shuffledCombinedArr = np.take_along_axis(combinedArr, randomIndices, axis=1)
    if verbose and rowsToCalculate[0] == 0: print("    Time:", time() - tShuffle, flush=True)
    
    # In the case of calculating the scores for pairwise epilogos, 
    # we need the original file 1 and file 2 arrays as well as their shuffled counterparts
    # shuffledCombinedArr is split by size of the original arrays
    return file1Arr, file2Arr, shuffledCombinedArr[:, :file1Arr.shape[1]], shuffledCombinedArr[:, file1Arr.shape[1]:]
