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

    

    argsortTimes = []
    partTimes = []

    for i in range(10):
        rand = np.random.rand(1000000, 721)
        tSort = time.time()
        rand.argsort(axis=1)
        argsortTimes.append(time.time() - tSort)

    for i in range(10):
        rand = np.random.rand(1000000, 721)
        tPart = time.time()
        np.argpartition(rand,0,axis=1)
        partTimes.append(time.time() - tPart)

    print("ArgSort Mean: {}\t|\tAll: {}".format(sum(argsortTimes)/len(argsortTimes), argsortTimes))
    print("Partition Mean: {}\t|\tAll: {}".format(sum(partTimes)/len(partTimes), partTimes))


if __name__ == "__main__":
    main()