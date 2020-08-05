import gzip
import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time
import click

# @click.command()
# @click.option("-d", "--dataset",        type=str,   required=True, help="Source publication or dataset (ROADMAP, ADSERA, or GORKIN)")
# @click.option("-a", "--assembly",       type=str,   required=True, help="Genomic assembly (hg19, hg38, or mm10)")
# @click.option("-m", "--state-model",    type=int,   required=True, help="State model (15, 18, or 25 for ROADMAP; 15 or 18 for ADSERA; 15 for GORKIN)")
# @click.option("-g", "--group",          type=str,   required=True, help="Individual dataset group name (using \"new\" naming scheme, ref. /net/seq/data/projects/Epilogos/epilogos-by-sample-group)")
# @click.option("-l", "--saliency-level", type=str,   required=True, help="Saliency level (S1, S2, or S3)")
# @click.option("-c", "--chromosome",     type=str,   required=True, help="Query chromosome")
# @click.option("-s", "--start",          type=int,   required=True, help="Query start position")
# @click.option("-e", "--end",            type=int,   required=True, help="Query end position")

def main(filename, numStates, saliency, outputDirectory, columnSpecification = "All"):
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    print("Reading data from file...")

    #Determine the number of columns in the data
    nColsData = len(np.loadtxt(fname = dataFilePath, dtype="str", max_rows= 1)) - 3

    # Read in the data
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, header = None, sep="\t", usecols=columnSpecificationAsTuple(columnSpecification, nColsData))
    numRows, numCols = dataDF.shape
    dataArr = dataDF.to_numpy()
    print("    Time: ", time.time() - tRead)

    # Calculate the expected frequencies of each state
    print("Calculating expected frequencies...")
    tExp = time.time()
    stateIndices = list(range(1, 16))
    expFreqSeries = pd.Series(np.zeros(15), index=stateIndices)
    # Raw counts
    for i in range(3, numCols):
        stateCounts = dataDF[i].value_counts()
        for state, count in stateCounts.items():
            expFreqSeries.loc[state] += count
    # Computing the frequencies
    dfSize = expFreqSeries.sum()
    for state, count in expFreqSeries.items():
        expFreqSeries.loc[state] = count / dfSize

    print("    Time: ", time.time() - tExp)

    # Calculate the observed frequencies and final scores in one loop
    print("Calculating scores...")
    tScore = time.time()
    expFreqArr = expFreqSeries.to_numpy()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(numRows):
        uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
        for i in range(len(uniqueStates)):
            # Function input is obsFreq and expFreq
            scoreArr[row, int(uniqueStates[i]) - 1] = klScore(stateCounts[i] / nColsData, expFreqArr[int(uniqueStates[i]) - 1])
    print("    Time: ", time.time() - tScore)

    # Writing the data to the files
    print("Writing to files...")
    tWrite = time.time()

    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    observationsTxtPath = outputDirPath / "observations.txt"
    scoresTxtPath = outputDirPath / "scores.txt"

    observationsTxt = open(observationsTxtPath, "w")
    scoresTxt = open(scoresTxtPath, "w")

    # Write each row in both observations and scores
    for i in range(numRows):
        # Write in the coordinates
        locationData = dataArr[i, :3]
        for location in locationData:
            observationsTxt.write("%s\t" % location)
            scoresTxt.write("%s\t" % location)
        
        # Write to observations
        maxContribution = np.amax(scoreArr[i])
        maxContributionLoc = np.argmax(scoreArr[i]) + 1
        totalScore = np.sum(scoreArr[i])

        observationsTxt.write("%d\t" % maxContributionLoc)
        observationsTxt.write("%.5f\t" % maxContribution)
        observationsTxt.write("1\t")
        observationsTxt.write("%.5f\t\n" % totalScore)

        # Write to scores
        for j in range(numStates):
            scoresTxt.write("%.5f\t" % scoreArr[i, j])
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

    print("    Time: ", time.time() - tWrite)

# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to create a sequence from the inputed column specification
def columnSpecificationAsTuple(columnSpecification, ncols):
    columnListExpanded = [0, 1, 2]
    if columnSpecification == "All":
        # Include all columns except the first 3
        columnListExpanded = list(range(0, ncols + 3))
    else:
        columnList = columnSpecification.split(',')
        for str in columnList:
            # Expand the hyphenated parts into a full list (e.g. "33-36" -> [33, 34, 35, 36])
            # Add 3 because of 3 not data columns in file
            if "-" in str:
                strList = str.split("-")
                columnListExpanded.append(list(range(int(strList[0]) + 3, int(strList[1] + 1 + 3))))
            else:
                columnListExpanded.append(int(str) + 3)
    return tuple(columnListExpanded)

if __name__ == "__main__":
    # Checking that the arguments are all correct
    if len(sys.argv) - 1 < 4:
        # Argument info if wrong number
        print("You must provide at least 4 arguments:\n")
        print("   1. Filename to read from\n")
        print("   2. Number of states in chromatin state model\n")
        print("   3. Saliency metric (1-3)\n")
        print("   4. Output directory\n")
        print("   5. (Optional) Specification of columns to read (all by default)\n")
    elif int(sys.argv[2]) != 15:
        print("We currently only offer support for a 15-state chromatin state model")
    elif int(sys.argv[3]) != 1:
        print("We currently only offer support for a saliency of one")
    elif len(sys.argv) - 1 == 4:
        # if the last argument is not used use the default
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    else:
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5])