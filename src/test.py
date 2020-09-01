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

def main():

    strArr = np.array([["chr1", "0", "200"],
                       ["chr1", "200", "400"],
                       ["chr1", "400", "600"]], dtype=str)

    intArr = np.array([[1.7, 1.3, .9],
                       [2.0, .21, .08],
                       [7.3, 2.1, 33]], dtype=str)

    print(intArr)

    print(np.argmax(intArr.astype(float)))

    print(os.getcwd())
    print(Path.cwd())
    print(Path(__file__).is_absolute())
    for file in (Path.cwd() / Path(__file__).parents[0]).glob("*"):
        print(file.name.split(".")[0])


    print()

    print("{0}{0}{0}{1}".format("test", "success"))


    

if __name__ == "__main__":
    main()