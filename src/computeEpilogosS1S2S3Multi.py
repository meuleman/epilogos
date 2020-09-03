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
import multiprocessing
import itertools
import click

@click.command()
@click.option("-f", "--filename", type=str, required=True, help="Path to file to read from")
@click.option("-m", "--state-model", "numStates", type=int, required=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency-level", "saliency", type=int, required=True, help="Saliency level (1, 2, or 3)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, help="Output Directory")
@click.option("-e", "--store-expected", "storeExpected", is_flag=True, help="[Flag] Store the expected frequency array for later calculations (Must be used in conjunction with '-d' and cannot be used in conjunction with '-u')")
@click.option("-u", "--use-expected", "useStoredExpected", is_flag=True, help="[Flag] Use previously stored expected frequency array (Must be used in conjunction with '-d' and cannot be used in conjunction with '-e')")
@click.option("-d", "--expected-directory", "expFreqDir", type=str, default="null", help="Path to the stored expected frequency array (Must be used in conjunction with either '-e' or '-u')")
def main(filename, numStates, saliency, outputDirectory, storeExpected, useStoredExpected, expFreqDir):
    tTotal = time.time()
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    # Check that paths are valid before doing anything
    if not dataFilePath.exists():
        print("ERROR: Given file path does not exist")
        return

    if not outputDirPath.is_dir():
        print("ERROR: Output directory is not a directory")
        return

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)


    # Read in the data
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    locationArr = dataDF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    # Variables for opting in to storing/retrieving the expected frequency array
    global storedExpDir
    if expFreqDir != "null":
        expFreqFilename = "exp_freq_" + str(dataArr.shape[1]) + "_" + str(numStates) + "_s" + str(saliency) + ".npy"
        storedExpDir = Path(expFreqDir) / expFreqFilename

    global storeExp
    storeExp = storeExpected

    global useStoredExp
    useStoredExp = useStoredExpected

    if saliency == 1:
        scoreArr = s1Score(dataDF, dataArr, numStates, outputDirPath)
    elif saliency == 2:
        scoreArr =s2Score(dataArr, numStates, outputDirPath)
    elif saliency == 3:
        scoreArr = s3Score(dataArr, numStates, outputDirPath)
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

    # If user desires, use the stored expected frequency array
    if useStoredExp:
        try:
            expFreqArr = np.load(storedExpDir, allow_pickle=False)
        except IOError:
            print("ERROR: Could not load stored expected value array.\n\tPlease check that the directory is correct or that the file exits")
    else:
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
        expFreqArr = expFreqSeries.to_numpy()

        # If user desires, store away the expected frequency array
        if storeExp:
            np.save(storedExpDir, expFreqArr, allow_pickle=False)

        print("    Time: ", time.time() - tExp)

    # Calculate the observed frequencies and final scores in one loop
    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(numRows):
        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
        for i in range(len(uniqueStates)):
            # Function input is obsFreq and expFreq
            scoreArr[row, uniqueStates[i]] = klScore(stateCounts[i] / (numCols), expFreqArr[uniqueStates[i]])
    print("    Time: ", time.time() - tScore)

    return scoreArr

# Function that calculates the scores for the S2 metric
def s2Score(dataArr, numStates, outputDirPath):
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

    # If the user desires, store the expected frequency array
    if storeExp:
        np.save(storedExpDir, expFreqArr, allow_pickle=False)

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
def s3Score(dataArr, numStates, outputDirPath):
    numRows, numCols = dataArr.shape
    numProcesses = multiprocessing.cpu_count()
    print("CPU COUNT: ", numProcesses)

    # FOR TESTING
    numRowsToCalculate = 100
    # FOR TESTING

    # Use multiprocessing to speed up expected frequency calculation time
    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()

    # If the user desires, use the stored expected frequency array
    if useStoredExp:
        expFreqArr = np.load(storedExpDir, allow_pickle=False)
    else:
        basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2))).T

        # Initializing needed variables
        expFreqArr = np.zeros((numCols, numCols, numStates, numStates))
        expQueue = multiprocessing.Queue()
        expProcesses = []

        # Creating the expected frequency processes and starting them
        for i in range(numProcesses):
            rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
            p = multiprocessing.Process(target=s3Exp, args=(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, expQueue))
            expProcesses.append(p)
            p.start()

        # Combine all the calculated expvalue arrays into one
        for process in expProcesses:
            expFreqArr += expQueue.get()

        # Shut down all the processes
        for process in expProcesses:
            process.join()

        # Normalize the array
        expFreqArr /= numRowsToCalculate * numCols * (numCols - 1)

        # If the user desires, store the expected frequency array
        if storeExp:
            np.save(storedExpDir, expFreqArr, allow_pickle=False)

    print("    Time: ", numRows * (time.time() - tExp) / numRowsToCalculate)

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates)) / (numCols * (numCols - 1)), expFreqArr)

    # Initializing necessary variables
    scoreArr = np.zeros((numRows, numStates))
    obsQueue = multiprocessing.Queue()
    obsProcesses = []

    tCalc = time.time()
    # Creating the observed frequency/score processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
        p = multiprocessing.Process(target=s3Obs, args=(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, scoreArrOnes, obsQueue))
        obsProcesses.append(p)
        p.start()

    print("    Calculation Time: ", numRows * (time.time() - tCalc) / numRowsToCalculate)

    tStore = time.time()
    # Move all the scores from the queue to the score array
    for i in range(numRowsToCalculate):
        scoreRow = obsQueue.get()
        scoreArr[scoreRow[0]] = scoreRow[1]

    print("    Storing Time: ", numRows * (time.time() - tStore) / numRowsToCalculate)

    # Shut down all the processes
    for process in obsProcesses:
        process.join()

    print("    Time: ", numRows * (time.time() - tScore) / numRowsToCalculate)
    return scoreArr

def s3Exp(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, queue):
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))
    for row in rowsToCalculate:
        expFreqArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]] += np.ones(basePermutationArr.shape[1])
    queue.put(expFreqArr)

def s3Obs(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, scoreArrOnes, queue):
    for row in rowsToCalculate:
        # Pull the scores from the precalculated score array
        rowScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        rowScoreArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]] = scoreArrOnes[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]]

        queue.put((row, rowScoreArr.sum(axis=(0,1,2))))

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
    observationsTxtPath = outputDirPath / "observationsM.txt.gz"
    scoresTxtPath = outputDirPath / "scoresM.txt.gz"

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
    main()