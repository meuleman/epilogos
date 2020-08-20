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
import ctypes
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

    # # Writing the scores to the files
    # print("Writing to files...")
    # tWrite = time.time()
    # writeScores(locationArr, scoreArr, outputDirPath, numStates)
    # print("    Time: ", time.time() - tWrite)

    # print("Total Time: ", time.time() - tTotal)

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
    print("Calculating observed frequencies and scores...")
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
    numProcesses = 32

    # FOR TESTING
    numRowsToCalculate = 10000
    # FOR TESTING

    # Use multiprocessing to speed up expected frequency calculation time
    print()
    print("Multiprocessing Queue")
    print("Calculating Expected Frequencies...")
    tExp = time.time()

    # Initializing needed variables
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))
    expQueue = multiprocessing.Queue()
    expProcesses = []

    # Creating the expected frequency processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
        p = multiprocessing.Process(target=s3ExpMulti, args=(dataArr, numCols, numStates, rowsToCalculate, expQueue))
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
    print("    Time: ", numRows * (time.time() - tExp) / numRowsToCalculate)

    # Use multiprocessing to speed up the observed frequency and score calculation time
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
        p = multiprocessing.Process(target=s3ObsMulti, args=(dataArr, numCols, numStates, rowsToCalculate, scoreArrOnes, obsQueue))
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




    print()
    print("Fully Vectorized Implementation Permutations:")
    # Calculate expected frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr = np.zeros((numCols, numStates, numCols, numStates))

    basePermutationArr = np.array(list(itertools.permutations(range(127), 2))).T

    # Initializing needed variables
    expFreqArr2 = np.zeros((numCols, numCols, numStates, numStates))
    expQueue2 = multiprocessing.Queue()
    expProcesses2 = []

    # Creating the expected frequency processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
        p = multiprocessing.Process(target=s3ExpMultiVectorized, args=(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, expQueue2))
        expProcesses2.append(p)
        p.start()

    # Combine all the calculated expvalue arrays into one
    for process in expProcesses2:
        expFreqArr2 += expQueue2.get()

    # Shut down all the processes
    for process in expProcesses2:
        process.join()

    # Normalize the array
    expFreqArr2 /= numRowsToCalculate * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / numRowsToCalculate)

    print("Calculating observed frequencies and scores...")
    tScore = time.time()
    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes2 = klScoreND(np.ones((numCols, numCols, numStates, numStates)) / (numCols * (numCols - 1)), expFreqArr2)

    # Initializing necessary variables
    scoreArr2 = np.zeros((numRows, numStates))
    obsQueue2 = multiprocessing.Queue()
    obsProcesses2 = []

    tCalc = time.time()
    # Creating the observed frequency/score processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
        p = multiprocessing.Process(target=s3ObsMultiVectorized, args=(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, scoreArrOnes2, obsQueue2))
        obsProcesses2.append(p)
        p.start()

    print("    Calculation Time: ", numRows * (time.time() - tCalc) / numRowsToCalculate)

    tStore = time.time()
    # Move all the scores from the queue to the score array
    for i in range(numRowsToCalculate):
        scoreRow = obsQueue2.get()
        scoreArr2[scoreRow[0]] = scoreRow[1]

    print("    Storing Time: ", numRows * (time.time() - tStore) / numRowsToCalculate)

    # Shut down all the processes
    for process in obsProcesses2:
        process.join()

    print("    Time: ", numRows * (time.time() - tScore) / numRowsToCalculate)




    # Use multiprocessing to speed up expected frequency calculation time
    print()
    print("Calculating Expected Frequencies Shared mp.array...")
    tExp = time.time()

    # Initializing needed variables
    expProcesses = []

    # Creating a shared expFreqArr
    expFreqArrShared2 = multiprocessing.Array(ctypes.c_int, numCols*numCols*numStates*numStates)

    # Creating the expected frequency processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
        p = multiprocessing.Process(target=s3ExpMultiShared2, args=(dataArr, numCols, numStates, rowsToCalculate, expFreqArrShared2))
        expProcesses.append(p)
        p.start()

    # Shut down all the processes
    for process in expProcesses:
        process.join()

    # Normalize the array
    print("    Time: ", numRows * (time.time() - tExp) / numRowsToCalculate)

    # Use multiprocessing to speed up the observed frequency and score calculation time
    print("Calculating observed frequencies and scores...")
    tScore = time.time()

    # Creating a shared expFreqArr
    scoreArrShared2 = multiprocessing.Array(ctypes.c_int, numRows*numStates)

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND((np.ones((numCols, numCols, numStates, numStates)) / (numCols * (numCols - 1))).flatten(), expFreqArrShared2).reshape(numCols, numCols, numStates, numStates)

    # Initializing necessary variables
    obsProcesses = []

    # Creating the observed frequency/score processes and starting them
    for i in range(numProcesses):
        rowsToCalculate = range(i * numRowsToCalculate // numProcesses, (i+1) * numRowsToCalculate // numProcesses)
        p = multiprocessing.Process(target=s3ObsMultiShared2, args=(dataArr, numRows, numCols, numStates, rowsToCalculate, scoreArrOnes, scoreArrShared2))
        obsProcesses.append(p)
        p.start()

    # Shut down all the processes
    for process in obsProcesses:
        process.join()
    
    print("    Time: ", numRows * (time.time() - tScore) / numRowsToCalculate)


    # print("Original multiprocess = shared np.ndarray: ", (scoreArr == scoreArrShared).all())
    # print("Original multiprocess = shared mp.array:   ", (scoreArr == scoreArrShared2).all())
    # print("shared mp.array       = shared np.ndarray: ", (scoreArrShared2 == scoreArrShared).all())


    return scoreArr

# Helper function to allow for multiprocessing of the expected frequency calculation
def s3ExpMulti(dataArr, numCols, numStates, rowsToCalculate, queue):
    # initialize the return array
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))

    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: # POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    expFreqArr[it1.index, it2.index, int(s1), int(s2)] += 1

    queue.put(expFreqArr)

# Helper function to allow for multiprocessing of the observed frequency and score calculation
def s3ObsMulti(dataArr, numCols, numStates, rowsToCalculate, scoreArrOnes, queue):
    for row in rowsToCalculate:
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    tempScoreArr[it1.index, it2.index, int(s1), int(s2)] = scoreArrOnes[it1.index, it2.index, int(s1), int(s2)]
        queue.put((row, tempScoreArr.sum(axis=(0,1,2))))


# Helper function to allow for multiprocessing of the expected frequency calculation
def s3ExpMultiShared(dataArr, numCols, numStates, rowsToCalculate, expFreqArr):
    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: # POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    expFreqArr[it1.index, it2.index, int(s1), int(s2)] += 1

# Helper function to allow for multiprocessing of the observed frequency and score calculation
def s3ObsMultiShared(dataArr, numCols, numStates, rowsToCalculate, scoreArrOnes, scoreArr):
    for row in rowsToCalculate:
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    tempScoreArr[it1.index, it2.index, int(s1), int(s2)] = scoreArrOnes[it1.index, it2.index, int(s1), int(s2)]
        scoreArr[row] = tempScoreArr.sum(axis=(0,1,2))


# Helper function to allow for multiprocessing of the expected frequency calculation
def s3ExpMultiShared2(dataArr, numCols, numStates, rowsToCalculate, expFreqArr):
    sharedArrFlat = np.frombuffer(expFreqArr.get_obj(), dtype=ctypes.c_int)
    sharedArr = sharedArrFlat.reshape(numCols, numCols, numStates, numStates)
    # s1 = state 1, s2 = state 2
    for row in rowsToCalculate:
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: # POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    sharedArr[it1.index, it2.index, int(s1), int(s2)] += 1 / (100 * numCols * (numCols - 1))

# Helper function to allow for multiprocessing of the observed frequency and score calculation
def s3ObsMultiShared2(dataArr, numRows, numCols, numStates, rowsToCalculate, scoreArrOnes, scoreArr):
    sharedArrFlat = np.frombuffer(scoreArr.get_obj(), dtype=ctypes.c_int)
    sharedArr = sharedArrFlat.reshape(numRows, numStates)
    for row in rowsToCalculate:
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                    tempScoreArr[it1.index, it2.index, int(s1), int(s2)] = scoreArrOnes[it1.index, it2.index, int(s1), int(s2)]
        sharedArr[row] = tempScoreArr.sum(axis=(0,1,2))

def s3ExpMultiVectorized(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, queue):
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))
    for row in rowsToCalculate:
        fullPermutationArr = np.array([basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]])
        expFreqArr[fullPermutationArr[0], fullPermutationArr[1], fullPermutationArr[2], fullPermutationArr[3]] += np.ones(fullPermutationArr.shape[1])
    queue.put(expFreqArr)

def s3ObsMultiVectorized(dataArr, numCols, numStates, rowsToCalculate, basePermutationArr, scoreArrOnes, queue):
    for row in rowsToCalculate:
        # Creates an array of all the permutations of columns and states
        fullPermutationArr = np.array([basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]])
        
        # Pull the scores from the precalculated score array
        rowScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        rowScoreArr[fullPermutationArr[0], fullPermutationArr[1], fullPermutationArr[2], fullPermutationArr[3]] = scoreArrOnes[fullPermutationArr[0], fullPermutationArr[1], fullPermutationArr[2], fullPermutationArr[3]]

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
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

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