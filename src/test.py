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

def main():

    path = Path("C:/Users/Jake/Desktop/epilogos/epilogos/data/")

    test = path.glob("chr*.gz")

    for x in test:
        numEpigenomes = pd.read_table(x, header=None, sep="\t", nrows=1).shape[1] - 3
        break
    

    print(numEpigenomes)
    

if __name__ == "__main__":
    main()