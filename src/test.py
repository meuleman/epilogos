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

def main(file):
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

    numProcesses = cpu_count()
    # np.random.seed(7062016)
    pairwiseMetrics = Path(file)
    # pairwiseMetrics   = Path("/home/jquon/bootstrappingTest/100000/pairwiseMetrics_mmo_fmo_s1.txt.gz")
    names = ['chr', 'binStart', 'binEnd', 'distance']
    distanceDF = pd.read_table(pairwiseMetrics, usecols=[0, 1, 2, 3], header=None, sep="\t", names=names)
    # distanceArr = np.array(distanceDF['distance'])

    # idx = [i for i in range(len(distanceArr)) if round(distanceArr[i], 5) != 0]
    # dataNull = pd.Series(distanceArr[idx])

    # # Multiprocess the reading
    # with closing(Pool(numProcesses)) as pool:
    #     results = pool.starmap(test, zip(itertools.repeat(dataNull, 1001), itertools.repeat(1000, 1001)))
    # pool.join()

    # index = [i for i in range(len(results))]
    # columns = [i for i in range(1000)]

    # randomDF = pd.DataFrame(index=index, columns=columns)

    # for i in range(len(results)):
    #     randomDF.iloc[i] = list(results[i])

    print("unique distances", len(distanceDF['distance'].unique()))
    # print("Unique values", len(randomDF.iloc[:,0].unique()))

def test(x, size):
    np.random.seed()
    return np.random.choice(x, size=size, replace=False)

if __name__ == "__main__":
    main(sys.argv[1])