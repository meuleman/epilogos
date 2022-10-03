from sys import argv
from pathlib import Path
import numpy as np
from time import time
from os import remove
import pandas as pd
from epilogos.pairwiseVisual import findSign
from epilogos.helpers import strToBool, getStateNames
import filter_regions as fr

def main(outputDir, stateInfo, fileTag, expFreqPath, exemplarWidth, verbose):
    """
    Finds the top scoring regions across all epilogos score files and puts them into a txt file

    Input:
    outputDir -- The path of the output directory
    stateInfo -- State model tab seperated information file
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    expFreqPath -- The location of the stored expected frequency array
    exemplarWidth -- 2*exemplarWidth+1 = size of exemplar regions
    verbose -- Boolean which if True, causes much more detailed prints
    """
    outputDirPath = Path(outputDir)

    stateNameList = getStateNames(stateInfo)

    if verbose: print("\nReading in score files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, scoreArr = readInData(outputDirPath)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    if verbose: print("\nFinding greatest hits...", flush=True); tHits = time()
    else: print("    Greatest hits txt\t", end="", flush=True)
    greatestHitsPath = outputDirPath / "greatestHits_{}.txt".format(fileTag)
    createTopScoresTxt(greatestHitsPath, locationArr, scoreArr, stateNameList, exemplarWidth)
    if verbose: print("    Time:", time() - tHits, flush=True)
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
    scoreArr.sum(axis=1) -- Numpy array containing sums of the KL-Scores for each genomic bin
    maxScoreArr -- Numpy array containing the state which had the highest KL-score in each bin
    """
    results = map(unpackNPZ, outputDirPath.glob("temp_scores_*.npz"))

    # Split up read results into tuples of chromosomes, scores, and locations
    dataChunks = list(zip(*results))

    # Figuring out chromosome order
    rawChrNamesInts = []
    rawChrNamesStrs = []
    for chromosome in dataChunks[0]:
        try:
            rawChrNamesInts.append(int(chromosome.split("chr")[-1]))
        except ValueError:
            rawChrNamesStrs.append(chromosome.split("chr")[-1])
    rawChrNamesInts.sort()
    rawChrNamesStrs.sort()
    chrOrder = rawChrNamesInts + rawChrNamesStrs
    for i in range(len(chrOrder)):
        chrOrder[i] = "chr" + str(chrOrder[i])

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

    return locationArr, scoreArr.sum(axis=1)


def unpackNPZ(file):
    """
    Takes in an .npz file containing scores and returns the individual numpy arrays

    Input:
    file -- The .npz file to unpack

    Output:
    npzFile['chrName'][0] -- The name of the chromosome which the scores are of
    npzFile['scoreArr'] -- Numpy array containing the kullback leibler scores
    npzFile['locationArr'] -- Numpy array containing the genomic locations for all the scores
    """
    npzFile = np.load(Path(file))
    return npzFile['chrName'][0], npzFile['scoreArr'], npzFile['locationArr']


def createTopScoresTxt(filePath, locationArr, scoreArr, nameArr, exemplarWidth):
    """
    Finds the top 1000 scoring bins and merges adjacent bins, then outputs a txt containing these top scoring regions and
    some information about each (chromosome, bin start, bin end, state name, absolute value of score, sign of score)

    Input:
    filePath -- The path of the file to write to
    locationArr -- Numpy array containing the genomic locations of all the bins
    scoreArr -- Numpy array containing the sum of the scores within each bin
    nameArr -- Numpy array containing the names of all the states
    exemplarWidth -- 2*exemplarWidth+1 = size of exemplar regions
    """
    with open(filePath, 'w') as outFile:
        # Concatenate the location and score arrays so that filter regions outputs the region coords as well as scores
        f = fr.Filter(method='maxmean', input=np.concatenate((locationArr, scoreArr.reshape(len(scoreArr), 1)), axis=1), input_type='bedgraph', aggregation_method='max', window_bins=2 * exemplarWidth + 1, max_elements=100, preserve_cols=True, quiet=False)
        f.read()
        f.filter()
        exemplars = f.output_df

        exemplars = exemplars.sort_values(by=["RollingMax", "RollingMean", "Score"], ascending=False).reset_index(drop=True)

        indices = exemplars["OriginalIdx"].to_numpy()

        # Calculate the maximum contributing state in each region
        # In the case of a tie, the higher number state wins (e.g. last state wins if all states are 0)
        maxStates = np.zeros((len(indices), 1))
        for i, idx in enumerate(indices):
            # Flip makes it so tie leads to the higher number state
            # Calculate the maximum value for each state in and then from there the argmax for the max state
            # Subtract from the shape of the array to reverse the effect of flip and get the state number
            maxStates[i, 0] = scoreArr.shape[1] - np.argmax(np.max(np.flip(scoreArr[idx - exemplarWidth:idx + exemplarWidth + 1], axis=1), axis=0))

        locations = pd.DataFrame(exemplars.loc[:, ["Chromosome", "Start", "End", "Score"]]).astype({"Chromosome": str, "Start": np.int32, "End": np.int32, "Score": np.float32})
        locations["MaxScoreState"] = maxStates.astype(np.int32)

        # Write all the locations to the file
        outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3}\n"
        outString = "".join(outTemplate.format(locations.iloc[i], nameArr[int(float(locations.iloc[i, 4])) - 1],
                            abs(float(locations.iloc[i, 3])), findSign(float(locations.iloc[i, 3])))
                            for i in range(locations.shape[0]))

        outFile.write(outString)


if __name__ == "__main__":
    main(argv[1], argv[2], argv[3], argv[4], int(argv[5]), strToBool(argv[6]))