import gzip
import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time
import numpy.ma as ma
# import click

# @click.command()
# @click.option("-d", "--dataset",        type=str,   required=True, help="Source publication or dataset (ROADMAP, ADSERA, or GORKIN)")
# @click.option("-a", "--assembly",       type=str,   required=True, help="Genomic assembly (hg19, hg38, or mm10)")
# @click.option("-m", "--state-model",    type=int,   required=True, help="State model (15, 18, or 25 for ROADMAP; 15 or 18 for ADSERA; 15 for GORKIN)")
# @click.option("-g", "--group",          type=str,   required=True, help="Individual dataset group name (using \"new\" naming scheme, ref. /net/seq/data/projects/Epilogos/epilogos-by-sample-group)")
# @click.option("-l", "--saliency-level", type=str,   required=True, help="Saliency level (S1, S2, or S3)")
# @click.option("-c", "--chromosome",     type=str,   required=True, help="Query chromosome")
# @click.option("-s", "--start",          type=int,   required=True, help="Query start position")
# @click.option("-e", "--end",            type=int,   required=True, help="Query end position")

def main(filename, numStates, saliency, outputDirectory):
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    # Read in the data
    print("Reading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, header = None, sep="\t")
    print("    Time: ", time.time() - tRead)

    if saliency == 1:
        s1Score(dataDF, numStates, outputDirPath)
    elif saliency == 2:
        s2Score(dataDF, numStates, outputDirPath)
    elif saliency == 3:
        s3Score(dataDF, numStates, outputDirPath)
    else:
        print("Inputed saliency value not supported")
        return

# Function that calculates the scores for the S1 metric
def s1Score(dataDF, numStates, outputDirPath):
    numRows, numCols = dataDF.shape
    dataArr = dataDF.to_numpy(dtype = str)

    # Calculate the expected frequencies of each state
    print("Calculating expected frequencies...")
    tExp = time.time()
    stateIndices = list(range(1, numStates + 1))
    expFreqSeries = pd.Series(np.zeros(15), index=stateIndices)
    dfSize = numRows * (numCols - 3)
    for i in range(3, numCols):
        stateCounts = dataDF[i].value_counts()
        for state, count in stateCounts.items():
            expFreqSeries.loc[state] += count / dfSize
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
            column = int(uniqueStates[i]) - 1
            scoreArr[row, column] = klScore(stateCounts[i] / (numCols - 3), expFreqArr[column])
    print("    Time: ", time.time() - tScore)

    # Writing the scores to the files
    print("Writing to files...")
    tWrite = time.time()
    writeScores(dataArr, scoreArr, outputDirPath, numRows, numStates, 1)
    print("    Time: ", time.time() - tWrite)

# Function that calculates the scores for the S2 metric
def s2Score(dataDF, numStates, outputDirPath):
    numRows, numCols = dataDF.shape
    dataArr = dataDF.to_numpy(dtype = str)

    # Calculate the expected frequencies and observed frequencies in the same loop
    print("Calculating expected and observed frequencies...")
    tExp = time.time()
    expFreqArr = np.zeros((numStates, numStates))
    obsFreqArr = np.zeros((numRows, numStates, numStates))

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states) / numRows / numStates^2
    # SumOverRows: (Prob of choosing x and y) / numRows / numStates^2
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    for row in range(numRows):
        uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
        for i in range(len(uniqueStates)):
            for j in range(len(uniqueStates)):
                if int(uniqueStates[i]) > int(uniqueStates[j]) or int(uniqueStates[i]) < int(uniqueStates[j]):
                    expFreqArr[int(uniqueStates[i]) - 1, int(uniqueStates[j]) - 1] += stateCounts[i] * stateCounts[j] / math.comb(numCols - 3, 2) / numRows / numStates**2
                    obsFreqArr[row, int(uniqueStates[i]) - 1, int(uniqueStates[j]) - 1] += stateCounts[i] * stateCounts[j] / math.comb(numCols - 3, 2)
                elif int(uniqueStates[i]) == int(uniqueStates[j]):
                    expFreqArr[int(uniqueStates[i]) - 1, int(uniqueStates[j]) - 1] += stateCounts[i] * (stateCounts[i] - 1) / 2 / math.comb(numCols - 3, 2) / numRows / numStates**2
                    obsFreqArr[row, int(uniqueStates[i]) - 1, int(uniqueStates[j]) - 1] += stateCounts[i] * (stateCounts[i] - 1) / 2 / math.comb(numCols - 3, 2)
    print("    Time: ", time.time() - tExp)

    # print("Calculating scores...")
    # # Calculate the KL Scores
    # tScore = time.time()
    # scoreArr = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     scoreArr[row] = klScore2D(obsFreqArr[row], expFreqArr).sum(axis=0)
    # print("    Time: ", time.time() - tScore)

    print("Calculating scores...")
    # Calculate the KL Scores
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates, numStates))
    for row in range(numRows):
        scoreArr[row] = klScore2D(obsFreqArr[row], expFreqArr)
    print("    Time: ", time.time() - tScore)

    # Writing the scores to the files
    print("Writing to files...")
    tWrite = time.time()
    writeScores(dataArr, scoreArr, outputDirPath, numRows, numStates, 2)
    print("    Time: ", time.time() - tWrite)
    

# Function that calculates the scores for the S3 metric
def s3Score(dataDF, numStates, outputDirPath):
    return

# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to calculate KL-score for 2d arrays (cleans up the code)
def klScore2D(obs, exp):
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)

# Helper to write the final scores to files
def writeScores(dataArr, scoreArr, outputDirPath, numRows, numStates, saliency):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    observationsTxtPath = outputDirPath / "observations2.txt"
    scoresTxtPath = outputDirPath / "scores2.txt"

    observationsTxt = open(observationsTxtPath, "w")
    scoresTxt = open(scoresTxtPath, "w")

    # Write each row in both observations and scores
    for i in range(numRows):
        # Write in the coordinates
        locationData = dataArr[i, :3]
        for location in locationData:
            observationsTxt.write("%s\t" % location)
            scoresTxt.write("%s\t" % location)
        
        if saliency == 1:
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
        elif saliency == 2 or saliency == 3:

            #################################################################
            #    maxContributionS1 and maxContributionLocS1 need to fixed   #
            #################################################################

            maxContributionS1 = np.amax(scoreArr[i].sum(axis=0))
            maxContributionLocS1 = np.argmax(scoreArr[i].sum(axis=0)) + 1
            maxContributionS2 = np.amax(scoreArr[i])

            coords = np.where(scoreArr[i] == np.amax(scoreArr[i]))
            maxContributionLocS2 = (coords[0][0] + 1, coords[1][0] + 1)

            totalScore = np.sum(scoreArr[i])

            observationsTxt.write("%d\t" % maxContributionLocS1)
            observationsTxt.write("%.5f\t" % maxContributionS1)
            observationsTxt.write("1\t")
            observationsTxt.write("{}\t".format(maxContributionLocS2))
            observationsTxt.write("%.5f\t" % maxContributionS2)
            observationsTxt.write("1\t")
            observationsTxt.write("%.5f\t\n" % totalScore)

            # Write to scores
            #################################################################
            #    Scores need to fixed (its just sums not s1 scores)         #
            #################################################################
            for j in range(numStates):
                scoresTxt.write("%.5f\t" % scoreArr[i].sum(axis=0)[j])
            scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    # Checking that the arguments are all correct
    if len(sys.argv) - 1 < 4:
        # Argument info if wrong number
        print("You must provide at least 4 arguments:\n")
        print("   1. Filename to read from\n")
        print("   2. Number of states in chromatin state model\n")
        print("   3. Saliency metric (1-3)\n")
        print("   4. Output directory\n")
    elif int(sys.argv[2]) != 15:
        print("We currently only offer support for a 15-state chromatin state model")
    elif (not int(sys.argv[3]) == 1) and (not int(sys.argv[3]) == 2):
        print("We currently only offer support for a saliency of 1 or 2")
    else:
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])