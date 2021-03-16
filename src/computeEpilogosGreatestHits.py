from sys import argv
from computeEpilogosExpectedCombination import strToBool
from computeEpilogosPairwiseVisual import hasAdjacent, mergeAdjacent, findSign
from epilogosHelpers import strToBool
from pathlib import Path
from contextlib import closing
from multiprocessing import Pool, cpu_count
import pandas as pd
from itertools import repeat
import numpy as np
from time import time


def main(outputDir, numStates, fileTag, numProcesses, verbose):
    outputDirPath = Path(outputDir)

    if numStates == 18:
        stateNameList = np.array(["TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk", "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "ZNF/Rpts", "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"])
    elif numStates == 15:
        stateNameList = np.array(["TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"])
    else:
        raise ValueError("State model not supported for plotting")

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    if verbose: print("\nReading in score files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, scoreArr, maxScoreArr = readInData(outputDirPath, numProcesses, numStates)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    if verbose: print("\nFinding greatest hits...", flush=True); tHits = time()
    else: print("    Greatest hits txt\t", end="", flush=True)
    greatestHitsPath = outputDirPath / "greatestHits_{}.txt".format(fileTag)
    createTopScoresTxt(greatestHitsPath, locationArr, scoreArr, maxScoreArr, stateNameList)
    if verbose: print("    Time:", time() - tHits, flush=True)
    else: print("\t[Done]", flush=True)


def readInData(outputDirPath, numProcesses, numStates):
    # For keeping the data arrays organized correctly
    names = ["chr", "binStart", "binEnd"] + ["s{}".format(i) for i in range(1, numStates + 1)]

    # Data frame to dump inputed data into
    diffDF = pd.DataFrame(columns=names)

    with closing(Pool(numProcesses)) as pool:
        # results = pool.starmap(readTableMulti, zip(outputDirPath.glob("scores_*.txt.gz"), repeat(names)))
        results = pool.map(readTableMulti, outputDirPath.glob("temp_scores_*.npz"))
    pool.join()

    # # Concatenating all chunks to the real differences dataframe
    # for diffDFChunk in results:
    #     diffDF = pd.concat((diffDF, diffDFChunk), axis=0, ignore_index=True)

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


    # # Sorting the dataframes by chromosomal location
    # diffDF["chr"] = pd.Categorical(diffDF["chr"], categories=chrOrder, ordered=True)
    # diffDF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)
    # Creating array of null distances ordered by chromosome based on the read in chunks
    print(dataChunks[0])
    print(chrOrder)
    print(chrOrder[0])
    print(type(dataChunks[0][0]))
    print(type(chrOrder[0][0]))
    index = dataChunks[0].index(chrOrder[0])
    scoreArr = dataChunks[1][index]
    locationArr = dataChunks[2][index]
    for chrName in chrOrder[1:]:
        index = dataChunks[0].index(chrName)
        scoreArr = np.concatenate((scoreArr, dataChunks[1][index]))
        locationArr = np.concatenate((locationArr, dataChunks[2][index]))

    # Convert dataframes to np arrays for easier manipulation
    # locationArr     = diffDF.iloc[:,0:3].to_numpy(dtype=str)
    # scoreArr        = diffDF.iloc[:,3:].to_numpy(dtype=float)

    maxScoreArr = np.abs(np.argmax(np.abs(np.flip(scoreArr, axis=1)), axis=1) - scoreArr.shape[1]).astype(int)

    return locationArr, scoreArr.sum(axis=1), maxScoreArr


# def readTableMulti(file, names):
def readTableMulti(file):
    # diffDFChunk = pd.read_table(Path(file), header=None, sep="\t", names=names)
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
        outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2}\t{3}\n"
        outString = "".join(outTemplate.format(locations[i], nameArr[int(float(locations[i, 4])) - 1], abs(float(locations[i, 3])), findSign(float(locations[i, 3]))) for i in range(locations.shape[0]))
        f.write(outString)


if __name__ == "__main__":
    main(argv[1], int(argv[2]), argv[3], int(argv[4]), strToBool(argv[5]))