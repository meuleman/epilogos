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
    tTotal = time.time()
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    # Read in the data
    print("Reading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, header = None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=np.int8)
    locationArr = dataDF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    print(dataArr.shape)
    print(locationArr.shape)
    print(dataArr[0][1], ", ", type(dataArr[0][1]))

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
    numRows, numCols = dataDF.shape
    numCols -= 3

    # Calculate the expected frequencies of each state
    print("Calculating expected frequencies...")
    tExp = time.time()
    stateIndices = list(range(1, numStates + 1))
    expFreqSeries = pd.Series(np.zeros(15), index=stateIndices)
    dfSize = numRows * numCols
    for i in range(3, numCols + 3):
        stateCounts = dataDF[i].value_counts()
        for state, count in stateCounts.items():
            expFreqSeries.loc[state] += count / dfSize
    print("    Time: ", time.time() - tExp)

    print(expFreqSeries)

    # Calculate the observed frequencies and final scores in one loop
    print("Calculating scores...")
    tScore = time.time()
    expFreqArr = expFreqSeries.to_numpy()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(numRows):
        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
        for i in range(len(uniqueStates)):
            # Function input is obsFreq and expFreq
            column = uniqueStates[i] - 1
            scoreArr[row, column] = klScore(stateCounts[i] / (numCols), expFreqArr[column])
    print("    Time: ", time.time() - tScore)

    return scoreArr

# Function that calculates the scores for the S2 metric
def s2Score(dataDF, dataArr, numStates, outputDirPath):
    numRows, numCols = dataDF.shape
    numCols -= 3

    # Calculate the observed frequencies
    print("Calculating expected and observed frequencies...")
    tExp = time.time()
    # expFreqArr = np.zeros((numStates, numStates))
    obsFreqArr = np.zeros((numRows, numStates, numStates))

    # SumOverRows: (Within a row, how many ways can you choose x and y to be together) / (how many ways can you choose 2 states)
    # SumOverRows: (Prob of choosing x and y)
    # Can choose x and y to be together x*y ways if different and n(n-1)/2 ways if same (where n is the number of times that x/y shows up)
    combinations = math.comb(numCols, 2)
    for row in range(numRows):
        uniqueStates, stateCounts = np.unique(dataArr[row], return_counts=True)
        for i in range(len(uniqueStates)):
            for j in range(len(uniqueStates)):
                if uniqueStates[i] > uniqueStates[j] or uniqueStates[i] < uniqueStates[j]:
                    obsFreqArr[row, uniqueStates[i] - 1, uniqueStates[j] - 1]  = stateCounts[i] * stateCounts[j] / combinations / 2 # Extra 2 is to account for the symmetric matrix
                elif uniqueStates[i] == uniqueStates[j]:
                    obsFreqArr[row, uniqueStates[i] - 1, uniqueStates[j] - 1]  = math.comb(stateCounts[i], 2) / combinations

    # Calculate the expected frequencies by summing the observed frequencies for each row
    expFreqArr = obsFreqArr.sum(axis=0) / numRows
    print("    Time: ", time.time() - tExp)

    print("Calculating scores...")
    # Calculate the KL Scores
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(numRows):
        scoreArr[row] = klScore2D(obsFreqArr[row], expFreqArr).sum(axis=0)
    print("    Time: ", time.time() - tScore)

    return scoreArr
    

# Function that calculates the scores for the S3 metric
def s3Score(dataDF, dataArr, numStates, outputDirPath):
    numRows, numCols = dataDF.shape
    numCols -= 3

    # Calculate observed frequencies
    print("Calculating Expected Frequencies...")
    tExp = time.time()
    expFreqArr = np.zeros((numCols, numCols, numStates, numStates))

    # s1 = state 1, s2 = state 2
    # should we count the epigenome paired with itself?
    # only iterate over the data (not location)
    dataIter = np.nditer(dataArr[0:50], flags=["multi_index", "refs_ok"])
    for s1 in dataIter:
        rowIter = np.nditer(dataArr[dataIter.multi_index[0]], flags=["c_index", "refs_ok"])
        for s2 in rowIter:
            if dataIter.multi_index[1] != rowIter.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
                expFreqArr[dataIter.multi_index[1], rowIter.index, int(s1) - 1, int(s2) - 1] += 1

    # Normalize the array
    expFreqArr /= 50 * numCols * (numCols - 1)
    print("    Time: ", numRows * (time.time() - tExp) / 50)
    print(expFreqArr.sum())

    # print("Calculating scores...")
    # # Calculate the KL Scores
    # tScore = time.time()
    # scoreArr = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     scoreArr[row] = klScore2D(obsFreqArr[row], expFreqArr).sum(axis=(0,1,2))
    # print("    Time: ", time.time() - tScore)

    # Calculate the obsFreq and scores
    # print("Calculating observed frequencies and scores...")
    # tScore = time.time()
    # scoreArr = np.zeros((numRows, numStates))
    # for row in range(50):
    #     obsFreqArr = np.zeros((numCols, numCols, numStates, numStates))
    #     it1 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
    #     for s1 in it1:
    #         it2 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
    #         for s2 in it2:
    #             if it1.index != it2.index: #POTENTIALLY DELETABLE IF WE DONT MIND COUNTING EPIGENOME WITH ITSELF
    #                 obsFreqArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] += 1
    #     obsFreqArr /= numCols * (numCols - 1)
    #     scoreArr[row] = klScore2D(obsFreqArr, expFreqArr).sum(axis=(0,1,2))
    # print("    Time: ", numRows * (time.time() - tScore) / 50)


    scoreArrOnes = klScore2D(np.ones((numCols, numCols, numStates, numStates)), expFreqArr)
    print("Calculating observed frequencies and scores (refs | no cast | sum axis)")
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore) / 50)
    print()

    print("Calculating observed frequencies and scores (refs | no cast | incremental sum)")
    tScore2 = time.time()
    scoreArr2 = np.zeros((numRows, numStates))
    for row in range(50):
        it1 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
            for s2 in it2:
                if it1.index != it2.index:
                    scoreArr2[row, s1 - 1] += scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1] / 2
                    scoreArr2[row, s2 - 1] += scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1] / 2
    print("    Time: ", numRows * (time.time() - tScore2) / 50)
    print()


    print("Calculating observed frequencies and scores (refs | cast | sum axis)")
    tScore3 = time.time()
    scoreArr3 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr3[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore3) / 50)
    print()

    print("Calculating observed frequencies and scores (refs | no cast | incremental sum)")
    tScore4 = time.time()
    scoreArr4 = np.zeros((numRows, numStates))
    for row in range(50):
        it1 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "refs_ok"])
            for s2 in it2:
                if it1.index != it2.index:
                    scoreArr4[row, int(s1) - 1] += scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1] / 2
                    scoreArr4[row, int(s2) - 1] += scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1] / 2
    print("    Time: ", numRows * (time.time() - tScore4) / 50)
    print()

    print("Calculating observed frequencies and scores (raw loop)")
    tScore5 = time.time()
    scoreArr5 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        for s1 in range(numCols):
            for s2 in range(numCols):
                if s1 != s2:
                    tempScoreArr[s1, s2, dataArr[row, s1] - 1, dataArr[row, s2] - 1] = scoreArrOnes[s1, s2, dataArr[row, s1] - 1, dataArr[row, s2] - 1]
        scoreArr5[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore5) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | cast | sum axis)")
    tScore6 = time.time()
    scoreArr6 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr6[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore6) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | no cast | sum axis)")
    tScore7 = time.time()
    scoreArr7 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr7[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore7) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | no cast | sum axis | buffered int32)")
    tScore8 = time.time()
    scoreArr8 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int32'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int32'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr8[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore8) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | cast | sum axis | buffered int32)")
    tScore9 = time.time()
    scoreArr9 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int32'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int32'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr9[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore9) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | no cast | sum axis | buffered int16)")
    tScore10 = time.time()
    scoreArr10 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int16'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int16'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr10[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore10) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | cast | sum axis | buffered int16)")
    tScore11 = time.time()
    scoreArr11 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int16'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int16'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr11[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore11) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | no cast | sum axis | buffered int8)")
    tScore12 = time.time()
    scoreArr12 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int8'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int8'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr12[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore12) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | cast | sum axis | buffered int8)")
    tScore13 = time.time()
    scoreArr13 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int8'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"], op_dtypes=['int8'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr13[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore13) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | no cast | sum axis | buffered)")
    tScore14 = time.time()
    scoreArr14 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr14[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore14) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | cast | sum axis | buffered)")
    tScore15 = time.time()
    scoreArr15 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index", "buffered"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index", "buffered"])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr15[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore15) / 50)
    print()
    
    print("Calculating observed frequencies and scores (no refs | no cast | sum axis | copy int8)")
    tScore16 = time.time()
    scoreArr16 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"], op_flags=["readonly", "copy"], op_dtypes=['int8'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"], op_flags=["readonly", "copy"], op_dtypes=['int8'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, s1 - 1, s2 - 1] = scoreArrOnes[it1.index, it2.index, s1 - 1, s2 - 1]
        scoreArr16[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore16) / 50)
    print()

    print("Calculating observed frequencies and scores (no refs | cast | sum axis | copy int8)")
    tScore17 = time.time()
    scoreArr17 = np.zeros((numRows, numStates))
    for row in range(50):
        tempScoreArr = np.zeros((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"], op_flags=["readonly", "copy"], op_dtypes=['int8'])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"], op_flags=["readonly", "copy"], op_dtypes=['int8'])
            for s2 in it2:
                if it1.index != it2.index:
                    tempScoreArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = scoreArrOnes[it1.index, it2.index, int(s1) - 1, int(s2) - 1]
        scoreArr17[row] = tempScoreArr.sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore17) / 50)
    print()


    print("Calculating observed frequencies and scores (no refs | cast | sum axis | masked)")
    tScore = time.time()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(50):
        booleanArr = np.ones((numCols, numCols, numStates, numStates))
        it1 = np.nditer(dataArr[row], flags=["c_index"])
        for s1 in it1:
            it2 = np.nditer(dataArr[row], flags=["c_index"])
            for s2 in it2:
                if it1.index != it2.index:
                    booleanArr[it1.index, it2.index, int(s1) - 1, int(s2) - 1] = 0
        scoreArr[row] = ma.masked_array(scoreArrOnes, booleanArr).filled(0).sum(axis=(0,1,2))
    print("    Time: ", numRows * (time.time() - tScore) / 50)
    print()

    print((scoreArr == scoreArr6).all())








    return scoreArr

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
def writeScores(locationArr, scoreArr, outputDirPath, numStates):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    observationsTxtPath = outputDirPath / "observations3.txt"
    scoresTxtPath = outputDirPath / "scores3.txt"

    observationsTxt = open(observationsTxtPath, "w")
    scoresTxt = open(scoresTxtPath, "w")

    # Write each row in both observations and scores
    for i in range(locationArr.shape[0]):
        # Write in the coordinates
        for location in locationArr[i]:
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

if __name__ == "__main__":
    # Checking that the arguments are all correct
    if len(sys.argv) - 1 < 4:
        # Argument info if wrong number
        print("You must provide at least 4 arguments:\n")
        print("   1. Filename to read from\n")
        print("   2. Number of states in chromatin state model (only supports up to 127)\n")
        print("   3. Saliency metric (1-3)\n")
        print("   4. Output directory\n")
    elif int(sys.argv[2]) != 15:
        print("We currently only offer support for a 15-state chromatin state model")
    elif (not int(sys.argv[3]) == 1) and (not int(sys.argv[3]) == 2) and (not int(sys.argv[3]) == 3):
        print("We currently only offer support for a saliency of 1, 2, or 3")
    else:
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])