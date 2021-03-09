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
import functools
import errno
import logging
from multiprocessing import cpu_count, Pool
from contextlib import closing
import scipy.stats as st

def main(file1, file2):
    # filePath = Path(file)

    # names = ["beta", "loc", "scale", "mle"]

    # bootstrapDF = pd.read_table(filePath, header=None, sep="\t", names=names)


    # bootstrapDF.sort_values(by=["mle"], inplace=True)

    # print(bootstrapDF.shape)

    # print("Top bootstrap")
    # print("MLE:", bootstrapDF.iloc[0, 3])
    # print(bootstrapDF.head(10))

    # print()
    # print()
    # print("Bottom bootstrap")
    # print("MLE:", bootstrapDF.iloc[-1, 3])
    # print(bootstrapDF.tail(10))

    # print()
    # print()
    # print("Median of bootstrap")
    # print("MLE:", bootstrapDF.iloc[int((bootstrapDF.shape[0] - 1) / 2), 3])
    # print("{},{},{}".format(bootstrapDF.iloc[int((bootstrapDF.shape[0] - 1) / 2), 0], bootstrapDF.iloc[int((bootstrapDF.shape[0] - 1) / 2), 1], bootstrapDF.iloc[int((bootstrapDF.shape[0] - 1) / 2), 2]))        
    filePath1 = Path(file1)
    filePath2 = Path(file2)

    df1 = pd.read_table(filePath1, header=None, sep="\t", usecols=[0,1,2])
    df2 = pd.read_table(filePath2, header=None, sep="\t", usecols=[0,1,2])

    boolean = True

    for i in range(df1.shape[0]):
        if df1.iloc[i, 0] != df2.iloc[i, 0] or df1.iloc[i, 1] != df2.iloc[i, 1] or df1.iloc[i, 2] != df2.iloc[i, 2]:
            boolean = False

    print(boolean)



if __name__ == "__main__":
    main(sys.argv[1])