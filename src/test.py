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
import subprocess
from pathlib import PurePath

def main():
    
    test = np.array([[1.578658568658, 2.46546, 3.35345345],
                    [3.45345435, 4.4353454, 5.4535435345],
                    [5.45435, 6, 7]])

    print(np.amax(test, axis=1).reshape((test.shape[0], 1)))
    print(np.sum(test, axis=1).reshape((test.shape[0], 1)))

    print(np.around(test, decimals=5))

    scoreFMT = "%s\t%s\t%s"
    for i in range(2):
        scoreFMT += "\t%.5f"

    print(scoreFMT)


if __name__ == "__main__":
    main()