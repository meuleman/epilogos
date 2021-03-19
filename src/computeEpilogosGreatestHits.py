from sys import argv
from computeEpilogosExpectedCombination import strToBool
from computeEpilogosPairwiseVisual import hasAdjacent, mergeAdjacent, findSign
from epilogosHelpers import strToBool, getStateNames
from pathlib import Path
import numpy as np
from time import time

def main(outputDir, stateInfo, fileTag, verbose):
    outputDirPath = Path(outputDir)

    stateNameList = getStateNames(stateInfo)

    if verbose: print("\nReading in score files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, scoreArr, maxScoreArr = readInData(outputDirPath)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    if verbose: print("\nFinding greatest hits...", flush=True); tHits = time()
    else: print("    Greatest hits txt\t", end="", flush=True)
    greatestHitsPath = outputDirPath / "greatestHits_{}.txt".format(fileTag)
    createTopScoresTxt(greatestHitsPath, locationArr, scoreArr, maxScoreArr, stateNameList)
    if verbose: print("    Time:", time() - tHits, flush=True)
    else: print("\t[Done]", flush=True)


def readInData(outputDirPath):
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

    maxScoreArr = np.abs(np.argmax(np.abs(np.flip(scoreArr, axis=1)), axis=1) - scoreArr.shape[1]).astype(int)

    return locationArr, scoreArr.sum(axis=1), maxScoreArr


def unpackNPZ(file):
    npzFile = np.load(Path(file))
    return npzFile['chrName'][0], npzFile['scoreArr'], npzFile['locationArr']


def createTopScoresTxt(filePath, locationArr, scoreArr, maxScoreArr, nameArr):
    with open(filePath, 'w') as f:
        # Sort the values
        indices = (-np.abs(scoreArr)).argsort()[:1000]

        locations = np.concatenate((locationArr[indices], scoreArr[indices].reshape(len(indices), 1), maxScoreArr[indices].reshape(len(indices), 1)), axis=1)

        # Iterate until all is merged
        while(hasAdjacent(locations)):
            locations = mergeAdjacent(locations)
            
        # Write all the locations to the file
        outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3}\n"
        outString = "".join(outTemplate.format(locations[i], nameArr[int(float(locations[i, 4])) - 1], abs(float(locations[i, 3])), findSign(float(locations[i, 3]))) for i in range(locations.shape[0]))
        f.write(outString)


if __name__ == "__main__":
    main(argv[1], argv[2], argv[3], strToBool(argv[4]))