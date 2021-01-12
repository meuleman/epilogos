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
import scipy.stats as st


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

    # arr = np.array(np.arange(4, 837)).tolist()
    # print(arr)
    # random.shuffle(arr)
    # size = int(len(arr)/2)
    # print("{}\t1,2,3,{}".format("ARR 1", str(arr[:size]).strip('[]').replace(" ", "")))
    # print()
    # print("{}\t1,2,3,{}".format("ARR 2", str(arr[size:2*size]).strip('[]').replace(" ", "")))

    # print("Removing Directory")
    # b = os.system("rm -rf /home/jquon/osSystemTest/")
    # print("done", b)
    # print()

    # print("Making new dir")
    # x = os.system("mkdir -p /home/jquon/osSystemTest/")
    # print("done", x)
    # print()

    # print("Cp to new dir")
    # y = os.system("cp /home/jquon/EpilogosOutput_Adsera/bloodTCell/scores_home_jquon_AdseraStateByGroup_bloodTCell_split_epilogos_matrix_chr* /home/jquon/osSystemTest/")
    # print("done", y)
    # print()

    # print("combine in new dir")
    # z = os.system("gzip -cd /home/jquon/osSystemTest/scores_home_jquon_AdseraStateByGroup_bloodTCell_split_epilogos_matrix_chr* > /home/jquon/osSystemTest/combined.txt")
    # print("done", z)
    # print()

    # print("zip in new dir")
    # a = os.system("gzip /home/jquon/osSystemTest/combined.txt")
    # print("done", a)
    # print()


    # path = Path("C:/Users/User/Desktop/epilogos/new_adsera_columns.txt")

    # bad = ['BSS00205', 'BSS00694', 'BSS00359', 'BSS00148', 'BSS01629', 'BSS00433', 'BSS01546', 'BSS00474', 'BSS01632', 'BSS01248', 'BSS01395', 'BSS01172', 'BSS01063', 'BSS01242', 'BSS01223', 'BSS00358', 'BSS00382', 'BSS01425', 'BSS00434', 'BSS01679', 'BSS00456', 'BSS00357', 'BSS00380', 'BSS00364', 'BSS00506', 'BSS01177', 'BSS00460', 'BSS01452', 'BSS00436', 'BSS01678', 'BSS01670', 'BSS00366', 'BSS00435', 'BSS00430', 'BSS00360', 'BSS01872', 'BSS00432', 'BSS01246', 'BSS01681', 'BSS01457', 'BSS00467', 'BSS00356', 'BSS00431', 'BSS01671', 'BSS00362', 'BSS01161', 'BSS00313', 'BSS00327', 'BSS01163', 'BSS00695', 'BSS00469', 'BSS01396', 'BSS00320', 'BSS01627', 'BSS01880', 'BSS01060', 'BSS00744', 'BSS00324', 'BSS00429', 'BSS00365', 'BSS01180', 'BSS01702', 'BSS01247', 'BSS00457', 'BSS01839', 'BSS01473', 'BSS01061', 'BSS01062', 'BSS00463', 'BSS00437', 'BSS01883', 'BSS01652', 'BSS00383', 'BSS00323', 'BSS01199', 'BSS00458', 'BSS00363', 'BSS00406', 'BSS01183', 'BSS00465', 'BSS00408', 'BSS01408', 'BSS01093']

    # good = []

    # containing = []


    # with open(path, 'r') as f:
    #     line = f.readline()
    #     while line:
    #         good.append(line.split()[0])
    #         line = f.readline()

    # for item in bad:
    #     if item in good:
    #         containing.append(item)

    # print(good)
    # print(bad)
    # print()
    # print()
    # print(containing)

    yticks = []
    ytickLabels = ["$10^{-16}$", "$10^{-15}$", "$10^{-14}$", "$10^{-13}$", "$10^{-12}$", "$10^{-11}$", "$10^{-10}$", "$10^{-9}$", "$10^{-8}$", "$10^{-7}$", "$10^{-6}$", "$10^{-5}$", "$10^{-4}$", "$1$", "$10^{-4}$", "$10^{-5}$", "$10^{-6}$", "$10^{-7}$", "$10^{-8}$", "$10^{-9}$", "$10^{-10}$", "$10^{-11}$", "$10^{-12}$", "$10^{-13}$", "$10^{-14}$", "$10^{-15}$", "$10^{-16}$"]

    yticks.append(-st.gennorm.isf(10**-4/2, .2, loc=0, scale=.2))
    yticksFinal = []
    ytickLabelsFinal = []
    yticksFinal.append(yticks[0])

    print(yticks)


    for i in range(4):
        ytickLabelsFinal.append(ytickLabels[i])

    print(ytickLabelsFinal)

if __name__ == "__main__":
    main()