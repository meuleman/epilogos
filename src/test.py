import math
import time
import operator as op
import numpy as np
from functools import reduce
import multiprocessing
import ctypes
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

def main():
    # tsvPath = Path("C:/Users/User/Downloads/adsera_dataset.tsv")
    # comparisonPath = Path("C:/Users/user/Desktop/epilogos/adsera_matrix_columns.txt")

    # # Reading in the data
    # colList = [i for i in range(15)]
    # names = ["id", "type"]

    # dataDF = pd.read_table(tsvPath, header=0, sep="\t", usecols=colList)
    # comparisonDF = pd.read_table(comparisonPath, header=None, sep="\s[|]\s", names=names, dtype=str)

    # print()
    # print(dataDF.shape)
    # print()

    # dataDF.drop_duplicates(subset=['id'], keep='first', inplace=True, ignore_index=True)

    # print()
    # print(dataDF.shape)
    # print()

    # cut = []

    # for i in range(len(dataDF)):
    #     if dataDF['id'][i] not in comparisonDF['id'].values:
    #         cut.append(dataDF['id'][i])
    #         dataDF.drop(i, inplace=True)
            

    # print()
    # print(dataDF.shape)
    # print()
    # print(cut)

    # idOrder = comparisonDF['id'].values

    # dataDF['id'] = pd.Categorical(dataDF['id'], categories=idOrder, ordered=True)

    # dataDF.sort_values(by=['id'], inplace=True)

    # dataDF.reset_index(drop=True, inplace=True)

    # outPath = Path("C:/Users/user/Desktop/epilogos/new_adsera_columns.txt")

    # print(dataDF.head(10))

    # dataDF.to_csv(outPath, sep="\t", index=False)

    # dataFilePath = Path("/home/jquon/EpilogosInput_AdseraMF/male/split/epilogos_matrix_chr10.txt.gz")
    # pairwisePath = Path("/home/jquon/EpilogosInput_AdseraMF/female/split/epilogos_matrix_chr10.txt.gz")


    # print("\nReading data from file 1...")
    # tRead1 = time.time()
    # file1DF = pd.read_table(dataFilePath, header=None, sep="\t")
    # print("    Time: ", time.time() - tRead1)

    # print("\nReading data from file 2...")
    # tRead2 = time.time()
    # file2DF = pd.read_table(pairwisePath, header=None, sep="\t")
    # print("    Time: ", time.time() - tRead2)

    # # Converting to a np array for faster functions later
    # print("Converting to numpy arrays...")
    # tConvert = time.time()
    # file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int)
    # file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int)
    # locationArr = file1DF.iloc[:,0:3].to_numpy(dtype=str)
    # print("    Time: ", time.time() - tConvert)

    # # Combining the arrays for per row shuffling
    # dataArr = np.concatenate((file1Arr, file2Arr), axis=1)
    # dataDF = pd.concat((file1DF, file2DF.iloc[:,3:]), axis=1, ignore_index=True)

    # print("File 1 DF Shape:", file1DF.shape)
    # print("File 2 DF Shape:", file2DF.shape)
    # print("Combined DF Shape:", dataDF.shape)
    # print("File 1 Arr Shape:", file1Arr.shape)
    # print("File 2 Arr Shape:", file2Arr.shape)
    # print("Combined Arr Shape:", dataArr.shape)

    # numStates=18

    # numRows, numCols = dataArr.shape

    # print(list(dataDF.columns))

    # print("numRows", numRows)
    # print("numCols", numCols)

    # # Calculate the expected frequencies of each state
    # stateIndices = list(range(1, numStates + 1))
    # expFreqSeries = pd.Series(np.zeros(numStates), index=stateIndices)
    # dfSize = numRows * numCols
    # for i in range(3, numCols + 3):
    #     stateCounts = dataDF[i].value_counts()
    #     for state, count in stateCounts.items():
    #         expFreqSeries.loc[state] += count / dfSize

    arr = np.array(np.arange(4, 837)).tolist()
    print(arr)
    random.shuffle(arr)
    size = int(len(arr)/2)
    print("{}\t1,2,3,{}".format("ARR 1", str(arr[:size]).strip('[]').replace(" ", "")))
    print()
    print("{}\t1,2,3,{}".format("ARR 2", str(arr[size:2*size]).strip('[]').replace(" ", "")))



if __name__ == "__main__":
    main()