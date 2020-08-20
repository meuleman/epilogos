import gzip
import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time
import numpy.ma as ma
import operator as op
from functools import reduce
import itertools
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
    tTotal = time.time()
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    colNames = ["Chromosome", "Start", "End"] + list(range(127))

    # Read in the data
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, header=None, names=colNames, sep="\t")
    dataDF.iloc[:,3:] -= 1
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int)
    locationArr = dataDF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    if saliency == 1:
        scoreArr = s1Score(dataDF, dataArr, numStates, outputDirPath)
    elif saliency == 2:
        scoreArr =s2Score(dataDF, dataArr, numStates, outputDirPath)
    elif saliency == 3:
        scoreArr = s3Score(dataDF, dataArr, numStates, outputDirPath)
    else:
        print("Inputed saliency value not supported")
        return

    # Writing the scores to the files
    print("Writing to files...")
    tWrite = time.time()
    writeScores(locationArr, scoreArr, outputDirPath, numStates)
    print("    Time: ", time.time() - tWrite)

    print("Total Time: ", time.time() - tTotal)

# Function that calculates the scores for the S1 metric
def s1Score(dataDF, dataArr, numStates, outputDirPath):
    numRows, numCols = dataArr.shape

    # Calculate the expected frequencies of each state
    print("Calculating expected frequencies...")
    tExp = time.time()
    stateIndices = list(range(1, numStates + 1))
    expFreqSeries = pd.Series(np.zeros(numStates), index=stateIndices)
    dfSize = numRows * numCols
    for i in range(3, numCols + 3):
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
        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
        for i in range(len(uniqueStates)):
            # Function input is obsFreq and expFreq
            scoreArr[row, uniqueStates[i]] = klScore(stateCounts[i] / (numCols), expFreqArr[uniqueStates[i]])
    print("    Time: ", time.time() - tScore)

    return scoreArr

# Function that calculates the scores for the S2 metric
def s2Score(dataDF, dataArr, numStates, outputDirPath):
    numRows, numCols = dataArr.shape

    # Calculate the observed frequencies
    print("Calculating expected and observed frequencies...")
    tExp = time.time()
    # expFreqArr = np.zeros((numStates, numStates))
    obsFreqArr = np.zeros((numRows, numStates, numStates))

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    if sys.version_info < (3, 8):
        print("\nFor maximum efficiency please update python to version 3.8 or later")
        print("NOTE: The code will still run in a lower version, but will be slightly slower\n")
        combinations = ncr(numCols, 2)
        for row in range(numRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i in range(len(uniqueStates)):
                for j in range(len(uniqueStates)):
                    if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                        obsFreqArr[row, uniqueStates[i], uniqueStates[j]]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                    elif uniqueStates[i] == uniqueStates[j]:
                        obsFreqArr[row, uniqueStates[i], uniqueStates[j]]  = ncr(stateCounts[i], 2) / combinations
    else:
        combinations = math.comb(numCols, 2)
        for row in range(numRows):
            uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
            for i in range(len(uniqueStates)):
                for j in range(len(uniqueStates)):
                    if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                        obsFreqArr[row, uniqueStates[i], uniqueStates[j]]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                    elif uniqueStates[i] == uniqueStates[j]:
                        obsFreqArr[row, uniqueStates[i], uniqueStates[j]]  = math.comb(stateCounts[i], 2) / combinations

    # Calculate the expected frequencies by summing the observed frequencies for each row
    expFreqArr = obsFreqArr.sum(axis=0) / numRows
    print("    Time: ", time.time() - tExp)

    print("Calculating scores...")
    # Calculate the KL Scores
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(numRows):
        scoreArr[row] = klScoreND(obsFreqArr[row], expFreqArr).sum(axis=0)
    print("    Time: ", time.time() - tScore)

    return scoreArr
    
# Function that calculates the scores for the S3 metric
def s3Score(dataDF, dataArr, numStates, outputDirPath):
    numRows, numCols = dataArr.shape


    # FOR TESTING
    rowsToCalculate = range(500000, 500100)

    # # Calculate expected frequencies
    # print("Calculating Expected Frequencies...")
    # tExp = time.time()
    # expFreqArr = np.zeros((numCols, numCols, numStates, numStates))

    # # s1 = state 1, s2 = state 2
    # for row in rowsToCalculate:
    #     it1 = np.nditer(dataArr[row], flags=["c_index"])
    #     for s1 in it1:
    #         it2 = np.nditer(dataArr[row], flags=["c_index"])
    #         for s2 in it2:
    #             if it1.index != it2.index: # POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
    #                 expFreqArr[it1.index, it2.index, int(s1), int(s2)] += 1

    # # Normalize the array
    # expFreqArr /= len(rowsToCalculate) * numCols * (numCols - 1)
    # print("    Time: ", numRows * (time.time() - tExp) / len(rowsToCalculate))

    # # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # # This saves a lot of time in the loop as we are just looking up references and not calculating
    # scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates)) / (numCols * (numCols - 1)), expFreqArr)

    # print("Calculating observed frequencies and scores...")
    # tScore = time.time()
    # scoreArr = np.zeros((numRows, numStates))
    # for row in rowsToCalculate:
    #     tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
    #     it1 = np.nditer(dataArr[row], flags=["c_index"])
    #     for s1 in it1:
    #         it2 = np.nditer(dataArr[row], flags=["c_index"])
    #         for s2 in it2:
    #             if it1.index != it2.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
    #                 tempScoreArr[it1.index, it2.index, int(s1), int(s2)] = scoreArrOnes[it1.index, it2.index, int(s1), int(s2)]
    #     scoreArr[row] = tempScoreArr.sum(axis=(0,1,2))
    # print("    Time: ", numRows * (time.time() - tScore) / len(rowsToCalculate))


    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr = np.zeros((numCols, numStates, numCols, numStates))

    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: # POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    expFreqArr[it1.index, int(s1), it2.index, int(s2)] += 1

    # Normalize the array
    expFreqArr /= len(rowsToCalculate) * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / len(rowsToCalculate))

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numStates, numCols, numStates)) / (numCols * (numCols - 1)), expFreqArr)

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates))
    for row in rowsToCalculate:
        tempScoreArr = np.zeros((numCols, numStates, numCols, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    tempScoreArr[it1.index, int(s1), it2.index, int(s2)] = scoreArrOnes[it1.index, int(s1), it2.index, int(s2)]
        scoreArr[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore) / len(rowsToCalculate))


    print()
    print("Vectorized Format Permutations:")

    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr2 = np.zeros((numCols, numStates, numCols, numStates))

    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        permutationArr = np.array(list(itertools.permutations(dataDF.iloc[row, 3:].items(), 2))).reshape(-1, 4).T
        expFreqArr2[permutationArr[0], permutationArr[1], permutationArr[2], permutationArr[3]] += np.ones(permutationArr.shape[1])

    # Normalize the array
    expFreqArr2 /= len(rowsToCalculate) * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / len(rowsToCalculate))

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numStates, numCols, numStates)) / (numCols * (numCols - 1)), expFreqArr2)

    scoreArr2 = np.zeros((numRows, numStates))
    for row in rowsToCalculate:
        # Creates an array of all the permutations of columns and states
        permutationArr = np.array(list(itertools.permutations(dataDF.iloc[row, 3:].items(), 2))).reshape(-1, 4).T
        
        # Pull the scores from the precalculated score array
        rowScoreArr = np.zeros((numCols, numStates, numCols, numStates))
        rowScoreArr[permutationArr[0], permutationArr[1], permutationArr[2], permutationArr[3]] = scoreArrOnes[permutationArr[0], permutationArr[1], permutationArr[2], permutationArr[3]]

        scoreArr2[row] = rowScoreArr.sum(axis=(0,1,2)) 
    print("    Time: ", numRows * (time.time() - tScore) / len(rowsToCalculate))


    print()
    print("Vectorized Addition Combinations:")
    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr3 = np.zeros((numCols, numStates, numCols, numStates))

    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        combinationArr = np.array(list(itertools.combinations(dataDF.iloc[row, 3:].items(), 2))).reshape(-1, 4).T
        expFreqArr3[combinationArr[0], combinationArr[1], combinationArr[2], combinationArr[3]] += np.ones(combinationArr.shape[1])
        expFreqArr3[combinationArr[2], combinationArr[3], combinationArr[0], combinationArr[1]] += np.ones(combinationArr.shape[1])

    # Normalize the array
    expFreqArr3 /= len(rowsToCalculate) * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / len(rowsToCalculate))

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numStates, numCols, numStates)) / (numCols * (numCols - 1)), expFreqArr3)

    scoreArr3 = np.zeros((numRows, numStates))
    for row in rowsToCalculate:
        # Creates an array of all the permutations of columns and states
        combinationArr = np.array(list(itertools.combinations(dataDF.iloc[row, 3:].items(), 2))).reshape(-1, 4).T
        
        # Pull the scores from the precalculated score array
        rowScoreArr = np.zeros((numCols, numStates, numCols, numStates))
        indices1 = (combinationArr[0], combinationArr[1], combinationArr[2], combinationArr[3])
        indices2 = (combinationArr[2], combinationArr[3], combinationArr[0], combinationArr[1])
        rowScoreArr[indices1] = scoreArrOnes[indices1]
        rowScoreArr[indices2] = scoreArrOnes[indices2]


        scoreArr3[row] = rowScoreArr.sum(axis=(0,1,2)) 
    print("    Time: ", numRows * (time.time() - tScore) / len(rowsToCalculate))

    print()
    print("Fully Vectorized Implementation Combinations:")
    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr4 = np.zeros((numCols, numStates, numCols, numStates))

    baseCombinationArr = np.array(list(itertools.combinations(range(127), 2))).T

    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        fullCombinationArr = np.array([baseCombinationArr[0], dataArr[row, baseCombinationArr[0]], baseCombinationArr[1], dataArr[row, baseCombinationArr[1]]])
        expFreqArr4[fullCombinationArr[0], fullCombinationArr[1], fullCombinationArr[2], fullCombinationArr[3]] += np.ones(fullCombinationArr.shape[1])
        expFreqArr4[fullCombinationArr[2], fullCombinationArr[3], fullCombinationArr[0], fullCombinationArr[1]] += np.ones(fullCombinationArr.shape[1])

    # Normalize the array
    expFreqArr4 /= len(rowsToCalculate) * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / len(rowsToCalculate))

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numStates, numCols, numStates)) / (numCols * (numCols - 1)), expFreqArr4)

    scoreArr4 = np.zeros((numRows, numStates))
    for row in rowsToCalculate:
        # Creates an array of all the permutations of columns and states
        fullCombinationArr = np.array([baseCombinationArr[0], dataArr[row, baseCombinationArr[0]], baseCombinationArr[1], dataArr[row, baseCombinationArr[1]]])
        
        # Pull the scores from the precalculated score array
        rowScoreArr = np.zeros((numCols, numStates, numCols, numStates))
        indices1 = (fullCombinationArr[0], fullCombinationArr[1], fullCombinationArr[2], fullCombinationArr[3])
        indices2 = (fullCombinationArr[2], fullCombinationArr[3], fullCombinationArr[0], fullCombinationArr[1])
        rowScoreArr[indices1] = scoreArrOnes[indices1]
        rowScoreArr[indices2] = scoreArrOnes[indices2]

        scoreArr4[row] = rowScoreArr.sum(axis=(0,1,2))

    print("    Time: ", numRows * (time.time() - tScore) / len(rowsToCalculate))

    print()
    print("Fully Vectorized Implementation Permutations:")
    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr5 = np.zeros((numCols, numStates, numCols, numStates))

    basePermutationArr = np.array(list(itertools.permutations(range(127), 2))).T

    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        fullPermutationArr = np.array([basePermutationArr[0], dataArr[row, basePermutationArr[0]], basePermutationArr[1], dataArr[row, basePermutationArr[1]]])
        expFreqArr5[fullPermutationArr[0], fullPermutationArr[1], fullPermutationArr[2], fullPermutationArr[3]] += np.ones(fullPermutationArr.shape[1])

    # Normalize the array
    expFreqArr5 /= len(rowsToCalculate) * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / len(rowsToCalculate))

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numStates, numCols, numStates)) / (numCols * (numCols - 1)), expFreqArr5)

    scoreArr5 = np.zeros((numRows, numStates))
    for row in rowsToCalculate:
        # Creates an array of all the permutations of columns and states
        fullPermutationArr = np.array([basePermutationArr[0], dataArr[row, basePermutationArr[0]], basePermutationArr[1], dataArr[row, basePermutationArr[1]]])
        
        # Pull the scores from the precalculated score array
        rowScoreArr = np.zeros((numCols, numStates, numCols, numStates))
        rowScoreArr[fullPermutationArr[0], fullPermutationArr[1], fullPermutationArr[2], fullPermutationArr[3]] = scoreArrOnes[fullPermutationArr[0], fullPermutationArr[1], fullPermutationArr[2], fullPermutationArr[3]]

        scoreArr5[row] = rowScoreArr.sum(axis=(0,1,2))

    print("    Time: ", numRows * (time.time() - tScore) / len(rowsToCalculate))





    print()
    print("Expected Frequency Arrays Equal: ", (expFreqArr == expFreqArr4).all())
    print("Score Arrays Equal: ", (scoreArr == scoreArr4).all())
    print()
    print("Expected Frequency Arrays Equal (permutations and combinations): ", (expFreqArr5 == expFreqArr4).all())
    print("Score Arrays Equal (permutations and combinations): ", (scoreArr5 == scoreArr4).all())

    return scoreArr

# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to calculate KL-score for 2d arrays (cleans up the code)
def klScoreND(obs, exp):
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)

# Helper to write the final scores to files
def writeScores(locationArr, scoreArr, outputDirPath, numStates):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    observationsTxtPath = outputDirPath / "observationsV.txt.gz"
    scoresTxtPath = outputDirPath / "scoresV.txt.gz"

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    # Write each row in both observations and scores
    for i in range(locationArr.shape[0]):
        # Write in the coordinates
        for location in locationArr[i]:
            observationsTxt.write("{}\t".format(location))
            scoresTxt.write("{}\t".format(location))
        
        # Write to observations
        maxContribution = np.amax(scoreArr[i])
        maxContributionLoc = np.argmax(scoreArr[i]) + 1
        totalScore = np.sum(scoreArr[i])

        observationsTxt.write("{}\t".format(maxContributionLoc))
        observationsTxt.write("{0:.5f}\t".format(maxContribution))
        observationsTxt.write("1\t")
        observationsTxt.write("{0:.5f}\t\n".format(totalScore))

        # Write to scores
        for j in range(numStates):
            scoresTxt.write("{0:.5f}\t".format(scoreArr[i, j]))
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

# Helper to calculate combinations
def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n - r, -1), 1)
    denom = reduce(op.mul, range(1, r + 1), 1)
    return numer // denom

if __name__ == "__main__":
    # Checking that the arguments are all correct
    if len(sys.argv) - 1 < 4:
        # Argument info if wrong number
        print("\nYou must provide at least 4 arguments:\n")
        print("   1. Filename to read from\n")
        print("   2. Number of states in chromatin state model (only supports up to 127)\n")
        print("   3. Saliency metric (1-3)\n")
        print("   4. Output directory\n")
        print("NOTE: Please make sure you are have Python 3.8 or later installed for maximum efficiency (Python 3.6 is the oldest possible version)")
    elif int(sys.argv[2]) != 15:
        print("We currently only offer support for a 15-state chromatin state model")
    elif (not int(sys.argv[3]) == 1) and (not int(sys.argv[3]) == 2) and (not int(sys.argv[3]) == 3):
        print("We currently only offer support for a saliency of 1, 2, or 3")
    else:
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])