import gzip
from time import time
import pandas as pd
import numpy as np
from pathlib import Path
import re
from epilogos.filter_regions import Filter

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


def countRows(dataFilePath):
    """
    Determines the number of rows in a file

    Input:
    dataFilePath -- the data file for which we want to count rows

    Output:
    totalRows -- the number of rows in in the datafile
    """

    # Calculate the number of rows in the input file
    if dataFilePath.name.endswith("gz"):
        with gzip.open(dataFilePath, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))
    else:
        with open(dataFilePath, "rb") as f:
            totalRows = sum(bl.count(b'\n') for bl in blocks(f))

    return totalRows


def splitRows(totalRows, numProcesses):
    """
    Determines what rows of a datafile to assign to each core

    Input:
    totalRows    -- the number of rows we want to split between cores
    numProcesses -- the number of cores which we want to work over

    Output:
    A list of tuples which contain the first and last rows for each core to use
    """
    # Split the rows up according to the number of cores we have available
    # Row list contain tuples of the first and last rows to be assigned to each core
    rowList = []
    for i in range(numProcesses):
        rowsToCalc = (i * totalRows // numProcesses, (i + 1) * totalRows // numProcesses)
        rowList.append(rowsToCalc)

    return rowList


def readStates(file1Path=Path("null"), file2Path=Path("null"), rowsToCalc=(0, 0), expBool=True, verbose=True,
               groupSize=-1):
    """
    Reads the states from the relevant rows of the inputed data file(s)

    Input:
    file1Path  -- The path to the first file to read from (used in both single and paired) [default=Path("null")]
    file2Path  -- The path to the second file to read from (used only in paired epilogos) [default=Path("null")]
    rowsToCalc -- The first and last rows to read from the files [default=(0, 0)]
    expBool    -- Tells us if we are calculating the expected frequencies or the scores [default=True]
    verbose    -- If True, we print out updates [default=True]
    groupSize  -- When returning null output, return 2 evenly sized groups of this size

    Output:
        Single Epilogos:
            file1Arr    -- 2d numpy array of the state info for each epigenome in the data file over the relevant rows
        Paired Epilogos Expected Frequency:
            combinedArr -- 2d numpy array of the state info for each epigenome of both data files over the relevant rows
        Paired Epilogos Scores:
            file1Arr    -- 2d numpy array of the state info for each epigenome in the 1st data file over relevant rows
            file2Arr    -- 2d numpy array of the state info for each epigenome in the 2nd data file over relevant rows
            shuffledCombinedArr[:, :file1Arr.shape[1]] -- 2d numpy array of the state info shuffled between 1st and 2nd
                                                          data files, with the same width as file1Arr
            shuffledCombinedArr[:, file1Arr.shape[1]:] -- 2d numpy array of the state info shuffled between 1st and 2nd
                                                          data files, with the same width as file2Arr
    """
    # Read in the data from file 1
    if verbose and rowsToCalc[0] == 0: print("Reading data from file 1...", flush=True); tRead1 = time()
    # Dont want to read in locations
    cols = range(3, pd.read_table(file1Path, nrows=1, header=None, sep="\t").shape[1])
    # Read using pd.read_table and convert to numpy array for faster calculation (faster than np.genfromtext())
    file1Arr = pd.read_table(file1Path, usecols=cols, skiprows=rowsToCalc[0],
                             nrows=rowsToCalc[1] - rowsToCalc[0], header=None, sep="\t").to_numpy(dtype=int) - 1
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
                             nrows=rowsToCalc[1] - rowsToCalc[0], header=None, sep="\t").to_numpy(dtype=int) - 1
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
        return file1Arr, file2Arr, shuffledCombinedArr[:, :file1Arr.shape[1]],\
            shuffledCombinedArr[:, file1Arr.shape[1]:]
    else:
        return file1Arr, file2Arr, shuffledCombinedArr[:, :groupSize], shuffledCombinedArr[:, groupSize:2 * groupSize]


def generateRegionArr(query):
    """
    Takes in one or more query regions and builds an 2 numpy array containing their coordinates.
    Used to handle both single queries and files containing multiple queries

    Input:
    query -- Query coordinates formatted chr:start-end OR
             Path to bed file containing query coordinates on each line, formatted chr[tab]start[tab]end

    Output:
    2D Numpy array with the first dimension being query regions and the second dimension being coordinates
    """

    # if its a single query build a one coord array
    if re.fullmatch("chr[a-zA-z\\d]+:[\\d]+-[\\d]+", query):
        chr = query.split(":")[0]
        start = query.split(":")[1].split("-")[0]
        end = query.split(":")[1].split("-")[1]
        return np.array([[chr, int(start), int(end)]], dtype=object)
    # if it is a path to a file, build a numpy array with all coords
    elif Path(query).is_file():
        return pd.read_table(Path(query), sep="\t", header=None, usecols=[0, 1, 2]).to_numpy(dtype=object)
    else:
        raise ValueError("Please input valid query (region formatted as chr:start-end"
                         + "or path to bed file containing query regions)")


def orderChromosomes(chromosomes):
    """
    Takes in an unordered iterable of chromosomes and orders them first by number then alphabetically
    e.g. For humans: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y

    Input:
    chromosomes -- An iterable of strings containing chromosome names formatted chr[] ([] is replaced with the chr name)

    Output:
    chrOrder -- The ordered list of chromosomes
    """

    # Orders chromosomes so numbers count up, followed by alphabetical
    rawChrNamesInts = []
    rawChrNamesStrs = []
    for chromosome in chromosomes:
        try:
            rawChrNamesInts.append(int(chromosome.split("chr")[-1]))
        except ValueError:
            rawChrNamesStrs.append(chromosome.split("chr")[-1])
    rawChrNamesInts.sort()
    rawChrNamesStrs.sort()
    chrOrder = rawChrNamesInts + rawChrNamesStrs
    for i in range(len(chrOrder)):
        chrOrder[i] = "chr" + str(chrOrder[i])

    return chrOrder


def maxMean(inputArr, roiWidth, maxRegions):
    """
    Wrapper function for the filter-regions maxmean function. Generates 100 results and sorts them in descending order

    Input:
    inputArr   -- A numpy array containing coordinates in the first 3 columns and score metric in the 4th column
    roiWidth   -- The desired width of the regions of interest in bins
    maxReginos -- The number of regions for the maxmean algorithm to find

    Output:
    rois      -- Pandas dataframe containing the regions chosen by maxmean
    indices   -- The original indices within the inputArr for the center of each of the chosen regions
    """
    f = Filter(method='maxmean', input=inputArr, input_type='bedgraph', aggregation_method='max',
                  window_bins=roiWidth, max_elements=maxRegions, preserve_cols=True, quiet=False)
    f.read()
    f.filter()
    rois = f.output_df
    rois = rois.sort_values(by=["RollingMax", "RollingMean", "Score"], ascending=False).reset_index(drop=True)
    indices = rois["OriginalIdx"].to_numpy()

    return rois, indices


def generateROIIndicesArr(indices, roiWidth):
    """
    Creates a 2d numpy array containing the full indices [start...end] for each region chosen by maxmean

    Input:
    indices  -- The index of the centerpoint for each region of interest
    roiWidth -- The width of each region of interest in bins

    Output:
    np.array(upperLower, dtype=np.int32) -- 2D numpy array with the 1st dimension being regions
                                            and the 2nd dimension being all the indices within each region
    """
    upperLower = []
    for idx in indices:
        # If the region of interest is an odd number of bins we have to add 1 to the upper index
        lowerIdx = idx - roiWidth // 2
        upperIdx = idx + roiWidth // 2
        if roiWidth % 2: upperIdx += 1
        upperLower.append(np.arange(lowerIdx, upperIdx, dtype=np.int32))
    return np.array(upperLower, dtype=np.int32)


def findSign(x):
    """
    Returns a string containing the sign of the inputted number

    Input:
    x -- Any number

    Output:
    "+" or "-"
    """
    if (x >= 0):
        return "+"
    else:
        return "-"


def sharedToNumpy(sharedArr, numRows, numStates):
    """
    Helper for unflattening a shared array into a 2d numpy array

    Input:
    sharedArr -- The shared array to shape
    numRows   -- The number of rows for the numpy array
    numStates -- The number of columns for the numpy array

    Output:
    numpy array generated from a buffered shared array
    """
    return np.frombuffer(sharedArr, dtype=np.float32).reshape((numRows, numStates))
