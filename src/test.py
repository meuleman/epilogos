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
                       [7.3, 2.1, 33]], dtype=float)

    print(intArr)

    for i in range(10):
        arr = np.array(np.arange(i, i+9).reshape(3,3), dtype=float)
        np.save(Path.cwd() / "intArrTest{}.npy".format(i), arr, allow_pickle=False)

    count = 0
    for file in Path.cwd().glob("intArrTest*.npy"):
        if count == 0:
            expFreqArr = np.load(file, allow_pickle=False)
        else:
            expFreqArr += np.load(file, allow_pickle=False)
        count += 1
        # Delete file after we're done with it
        os.remove(file)

    print()
    print(expFreqArr)
    print(count)
    expFreqArr /= count

    print()
    print(expFreqArr)
    print()

    global x
    x= "fantabulous"

    test()


    print(os.getcwd())
    print(Path.cwd())
    print(Path(__file__).is_absolute())
    for file in (Path(__file__).parents[0]).glob("*"):
        print(file)


    print()

    print("{0}{0}{0}{1}".format("test", "success"))


def test():
    print(x)

if __name__ == "__main__":
    main()