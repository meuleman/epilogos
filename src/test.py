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
    filePath = Path(file)

    names = ["beta", "loc", "scale", "mle"]

    bootstrapDF = pd.read_table(filePath, header=None, sep="\t", names=names)


    bootstrapDF.sort_values(by=["mle"], inplace=True)

    print(bootstrapDF.shape)

    print("Top bootstrap")
    print(bootstrapDF.head(10))

    print()
    print()
    print("Bottom bootstrap")
    print(bootstrapDF.tail(10))

    print()
    print()
    print("Median of bootstrap")
    print(bootstrapDF.iloc[int((bootstrapDF.shape[0] - 1) / 2)])


if __name__ == "__main__":
    main(sys.argv[1])