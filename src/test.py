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

def main(file1, file2, numStates):
    
    file1Path = Path(file1)
    file2Path = Path(file2)
    
    colNames = ["chr", "binStart", "binEnd"] + ["s{}".format(i) for i in range(1, numStates+1)]
    chrOrder = []
    for i in range(1, 23):
        chrOrder.append("chr{}".format(i))
    chrOrder.append("chrX")

    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, header=None, sep="\t", names=colNames)
    print("    Time: ", time.time() - tRead1)

    if "all127" in file1.split("/"):
        for chr in chrOrder:
            index = file1DF["chr"].where(file1DF["chr"] == chr).last_valid_index()
            file1DF.drop(index, inplace=True)

    print("Reading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, header=None, sep="\t", names=colNames)
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
    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int)
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int)
    print("    Time: ", time.time() - tConvert)

    print("Calculating percent difference...")
    tDiff = time.time()
    maxStateArr1 = np.argmax(file1Arr, axis=1) + 1
    maxStateArr2 = np.argmax(file2Arr, axis=1) + 1
    error = np.mean(maxStateArr1 != maxStateArr2)
    print("    Time: ", time.time() - tDiff)
    
    print("Percent Difference is:", error)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])