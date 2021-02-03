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
    dataFilePath = Path("/home/jquon/AdseraStateByGroup/male/split/epilogos_matrix_chr22.txt.gz")
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(dataFilePath, skiprows=100000, nrows=500000, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    print(dataDF.shape)

if __name__ == "__main__":
    main()