import math
import time
import operator as op
import numpy as np
from functools import reduce
import multiprocessing
import ctypes as c
import itertools
import random
from pathlib import Path
import glob
import pandas as pd
import os
import sys
import subprocess
from pathlib import PurePath
import glob
import gzip
import numpy.ma as ma



def main():
    dataFilePath = Path("/home/jquon/epilogos/data/pyData/male/epilogos_matrix_chr1.txt.gz")
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, nrows=500, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    print("    Time: ", time.time() - tConvert)

    expFreqArr = np.load("/home/jquon/epilogosTesting_01122021/output/minimal/exp_freq_male_saliency3.npy", allow_pickle=False)

    numRows, numCols = dataArr.shape
    numStates = 18

    # Gives us everyway to combine the column numbers in numpy indexing form
    basePermutationArr = np.array(list(itertools.permutations(range(numCols), 2)), dtype=np.int16).T

    # Because each epigenome, epigenome, state, state combination only occurs once per row, we can precalculate all the scores assuming a frequency of 1/(numCols*(numCols-1))
    # This saves a lot of time in the loop as we are just looking up references and not calculating
    scoreArrOnes = klScoreND(np.ones((numCols, numCols, numStates, numStates), dtype=np.float32) / (numCols * (numCols - 1)), expFreqArr)


    print("Current way...")
    tCurrent = time.time()
    scoreArr = np.zeros((numRows, numStates), dtype=np.float32)
    rowScoreArr = np.zeros((numCols, numCols, numStates, numStates), dtype=np.float32)
    for row in range(numRows):
        # Reset the array so it doesn't carry over scores from other rows
        rowScoreArr.fill(0)

        # Pull the scores from the precalculated score array and put them into the correct positions for the state combinations that we observe
        rowScoreArr[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]] = scoreArrOnes[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]]

        # Flatten the scores and put them into the shared score array
        scoreArr[row] = rowScoreArr.sum(axis=(0,1,2))
    print("    Time:", time.time() - tCurrent)


    print("Changed Axes way...")
    scoreArrOnes.reshape((numStates, numStates, numCols, numCols))
    tAxes = time.time()
    scoreArr = np.zeros((numRows, numStates), dtype=np.float32)
    rowScoreArr = np.zeros((numStates, numStates, numCols, numCols), dtype=np.float32)
    for row in range(numRows):
        # Reset the array so it doesn't carry over scores from other rows
        rowScoreArr.fill(0)

        # Pull the scores from the precalculated score array and put them into the correct positions for the state combinations that we observe
        rowScoreArr[dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]], basePermutationArr[0], basePermutationArr[1]] = scoreArrOnes[dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]], basePermutationArr[0], basePermutationArr[1]]

        # Flatten the scores and put them into the shared score array
        scoreArr[row] = rowScoreArr.sum(axis=(1,2,3))
    print("    Time:", time.time() - tAxes)

    print("ravel way...")
    tRavel = time.time()
    scoreArr = np.zeros((numRows, numStates), dtype=np.float32)
    rowScoreArr = np.zeros((numCols, numCols, numStates, numStates), dtype=np.float32)
    for row in range(numRows):
        # Reset the array so it doesn't carry over scores from other rows
        rowScoreArr.fill(0)

        # Pull the scores from the precalculated score array and put them into the correct positions for the state combinations that we observe
        rowScoreArr[basePermutationArr[0] * rowScoreArr.shape[1] * rowScoreArr.shape[2] * rowScoreArr.shape[3] + basePermutationArr[1] * rowScoreArr.shape[1] * rowScoreArr.shape[2] + dataArr[row, basePermutationArr[1]] * rowScoreArr.shape[1] + dataArr[row, basePermutationArr[0]]] = scoreArrOnes[basePermutationArr[0] * scoreArrOnes.shape[1] * scoreArrOnes.shape[2] * scoreArrOnes.shape[3] + basePermutationArr[1] * scoreArrOnes.shape[1] * scoreArrOnes.shape[2] + dataArr[row, basePermutationArr[1]] * scoreArrOnes.shape[1] + dataArr[row, basePermutationArr[0]]]

        # Flatten the scores and put them into the shared score array
        scoreArr[row] = rowScoreArr.sum(axis=(0,1,2))
    print("    Time:", time.time() - tRavel)


    print("Less indexing way way...")
    tCurrent = time.time()
    scoreArr = np.zeros((numRows, numStates), dtype=np.float32)
    rowScoreArr = np.zeros(numStates, dtype=np.float32)
    for row in range(numRows):
        # Reset the array so it doesn't carry over scores from other rows
        rowScoreArr.fill(0)

        # Pull the scores from the precalculated score array and put them into the correct positions for the state combinations that we observe
        np.add.at(rowScoreArr, dataArr[row, basePermutationArr[1]], scoreArrOnes[basePermutationArr[0], basePermutationArr[1], dataArr[row, basePermutationArr[0]], dataArr[row, basePermutationArr[1]]])

        # Flatten the scores and put them into the shared score array
        scoreArr[row] = rowScoreArr
    print("    Time:", time.time() - tCurrent)


# Helper to calculate KL-score for Nd arrays (cleans up the code)
def klScoreND(obs, exp):
    return obs * ma.log2(ma.divide(obs, exp).filled(0)).filled(0)

if __name__ == "__main__":
    main()