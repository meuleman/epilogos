from sys import argv
from pathlib import Path
import numpy as np
from time import time
from os import remove
import pandas as pd
from epilogos.helpers import maxMean, generateROIIndicesArr, orderChromosomes, strToBool, getStateNames, findSign


def main(outputDir, stateInfo, fileTag, expFreqPath, roiWidth, verbose):
    """
    Finds the top scoring regions across all epilogos score files and puts them into a txt file

    Input:
    outputDir   -- The path of the output directory
    stateInfo   -- State model tab seperated information file
    fileTag     -- A string which helps ensure outputed files are named similarly within an epilogos run
    expFreqPath -- The location of the stored expected frequency array
    roiWidth    -- Size of regions of interest in bins
    verbose     -- Boolean which if True, causes much more detailed prints
    """
    outputDirPath = Path(outputDir)

    stateNameList = getStateNames(stateInfo)

    if verbose: print("\nReading in score files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, scoreArr = readInData(outputDirPath)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    if verbose: print("\nFinding regions of interest...", flush=True); tRoi = time()
    else: print("    Regions of interest txt\t", end="", flush=True)
    roiPath = outputDirPath / "regionsOfInterest_{}.txt".format(fileTag)
    createTopScoresTxt(roiPath, locationArr, scoreArr, stateNameList, roiWidth)
    if verbose: print("    Time:", time() - tRoi, flush=True)
    else: print("\t[Done]", flush=True)

    # Removing the expected frequency array
    remove(Path(expFreqPath))


def readInData(outputDirPath):
    """
    Reads all the epilogos score files in and combines them into a numpy array ordered by location

    Input:
    outputDirPath -- Path to the epilogos output directory (this contains the score files)

    Output:
    locationArr -- Numpy array containing the genomic locations for all the scores
    scoreArr    -- Numpy array containing the KL-Scores for each genomic bin
    """
    results = map(unpackNPZ, outputDirPath.glob("temp_scores_*.npz"))

    # Split up read results into tuples of chromosomes, scores, and locations
    dataChunks = list(zip(*results))

    # Figuring out chromosome order
    chrOrder = orderChromosomes(dataChunks[0])

    # Sorting the dataframes by chromosomal location
    # Creating array of scores and locations ordered by chromosome based on the read in chunks
    index = dataChunks[0].index(chrOrder[0])
    scoreArr = dataChunks[1][index]
    locationArr = dataChunks[2][index]
    for chrName in chrOrder[1:]:
        index = dataChunks[0].index(chrName)
        scoreArr = np.concatenate((scoreArr, dataChunks[1][index]))
        locationArr = np.concatenate((locationArr, dataChunks[2][index]))

    # Cleaning up the temp files after we've read them
    for file in outputDirPath.glob("temp_scores_*.npz"):
        remove(file)

    return locationArr, scoreArr


def unpackNPZ(file):
    """
    Takes in an .npz file containing scores and returns the individual numpy arrays

    Input:
    file -- The .npz file to unpack

    Output:
    npzFile['chrName'][0]  -- The name of the chromosome which the scores are of
    npzFile['scoreArr']    -- Numpy array containing the kullback leibler scores
    npzFile['locationArr'] -- Numpy array containing the genomic locations for all the scores
    """
    npzFile = np.load(Path(file), allow_pickle=True)
    return npzFile['chrName'][0], npzFile['scoreArr'], npzFile['locationArr']


def createTopScoresTxt(filePath, locationArr, scoreArr, nameArr, roiWidth):
    """
    Finds the top 100 regions of interest using filter-region's maxmean algorithm then outputs them in a txt with some
    information about each (chromosome, bin start, bin end, state name, absolute value of score, sign of score)

    Input:
    filePath    -- The path of the file to write to
    locationArr -- Numpy array containing the genomic locations of all the bins
    scoreArr    -- Numpy array containing the scores within each bin
    nameArr     -- Numpy array containing the names of all the states
    roiWidth    -- size of regions of interest in bins

    Output:
    regionsOfInterest*.txt of top 100 regions of interest formatted as follows:
        Column 1: Chromosome
        Column 2: Start coordinate
        Column 3: End coordinate
        Column 4: Name of the largest scoring state
        Column 5: Sum of the Kullback-Leibler scores
        Column 6: Sign of the sum of Kullback-Leibler scores
    """
    with open(filePath, 'w') as outFile:
        # Concatenate the location and score arrays so that filter regions outputs the region coords as well as scores
        rois, indices = maxMean(np.concatenate((locationArr, scoreArr.sum(axis=1).reshape(len(scoreArr), 1)), axis=1),
                                roiWidth, 100)

        # Generating array of all indices for vectorized calculations
        roiIndicesArr = generateROIIndicesArr(indices, roiWidth)

        # Calculate the maximum contributing state in each region
        # In the case of a tie, the higher number state wins (e.g. last state wins if all states are 0)
        # Flip makes it so tie leads to the higher number state
        # Calculate the maximum value for each state in and then from there the argmax for the max state
        # Subtract from the shape of the array to reverse the effect of flip and get the state number
        maxStates = scoreArr.shape[1] - np.argmax(np.max(np.flip(scoreArr[roiIndicesArr], axis=2), axis=1), axis=1)

        # Build pandas dataframe for writing
        locations = pd.DataFrame(rois.loc[:, ["Chromosome", "Start", "End", "Score"]])\
            .astype({"Chromosome": str, "Start": np.int64, "End": np.int64, "Score": np.float32})
        locations["MaxScoreState"] = maxStates.astype(np.int32)

        # Write all the locations to the file
        outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3}\n"
        outString = "".join(outTemplate.format(locations.iloc[i], nameArr[int(float(locations.iloc[i, 4])) - 1],
                            abs(float(locations.iloc[i, 3])), findSign(float(locations.iloc[i, 3])))
                            for i in range(locations.shape[0]))

        outFile.write(outString)


if __name__ == "__main__":
    main(argv[1], argv[2], argv[3], argv[4], int(argv[5]), strToBool(argv[6]))
