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


def main():

    numProcesses = multiprocessing.cpu_count()

    pool = multiprocessing.Pool(processes=numProcesses)

    rangeList = [range(2, 4), range(1, 8), range(5, 9)]

    y = 10
    z = 20

    results = [pool.apply(dummy, args=(y, x, z)) for x in rangeList]

    print(results)

def dummy(y, x, z):
    lst = []
    for i in x:
        lst.append(i+y+z)
    return lst

if __name__ == "__main__":
    main()