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
import click

def main(file1, file2):
    
    file1Path = Path(file1)
    file2Path = Path(file2)
    
    colNames = ["chr", "binStart", "binEnd"] + ["s{}".format(i) for i in range(1, 16)]
    chrOrder = []
    for i in range(1, 23):
        chrOrder.append("chr{}".format(i))
    chrOrder.append("chrX")
    chrOrder.append("chrY")


    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, header=None, sep="\t")
    print("    Time: ", time.time() - tRead1)

    print("Reading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, header=None, sep="\t")
    print("    Time: ", time.time() - tRead2)

    # Sorting the dataframes by chromosomal location
    file1DF["chr"] = pd.Categorical(file1DF["chr"], categories=chrOrder, ordered=True)
    file2DF["chr"] = pd.Categorical(file2DF["chr"], categories=chrOrder, ordered=True)
    file1DF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)
    file2DF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)

    print(file1DF.shape)
    print(file2DF.shape)


    print("Converting to numpy arrays...")
    tConvert = time.time()
    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int) - 1
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int) - 1
    locationArr = file1DF.iloc[:,:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    percentDiffArr = (file1Arr - file2Arr) / file1Arr * 100

    fivePercentDiff = np.where(np.any(percentDiffArr >= 5, axis=1))[0]
    tenPercentDiff = np.where(np.any(percentDiffArr >= 10, axis=1))[0]

    print(len(fivePercentDiff))
    print(len(tenPercentDiff))

    print(locationArr[fivePercentDiff])

    print(locationArr[tenPercentDiff])
    

if __name__ == "__main__":
    main()