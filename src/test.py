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

    dataFilePath = Path("/home/jquon/EpilogosInput_AdseraMF/male/split/epilogos_matrix_chr10.txt.gz")

    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr = dataDF.iloc[:,3:].to_numpy(dtype=int) - 1 
    print("    Time: ", time.time() - tConvert)

    numRows, numCols = dataArr.shape

    numStates=18

    print("numRows", numRows)
    print("numCols", numCols)

    # Calculate the expected frequencies of each state
    stateIndices = list(range(1, numStates + 1))
    expFreqSeries = pd.Series(np.zeros(numStates), index=stateIndices)
    dfSize = numRows * numCols
    for i in range(3, numCols + 3):
        stateCounts = dataDF[i].value_counts()
        print("Len(dataDF[i]): ", len(dataDF[i]))
        for state, count in stateCounts.items():
            expFreqSeries.loc[state] += count / dfSize


if __name__ == "__main__":
    main()