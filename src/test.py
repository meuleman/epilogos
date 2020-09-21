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
    
    test = np.array([[1, 2, 3],
                    [3, 4, 5],
                    [5, 6, 7]])

    print(np.amax(test, axis=1).reshape((test.shape[0], 1)))
    print(np.sum(test, axis=1).reshape((test.shape[0], 1)))
    print(np.sum(test, axis=0))

if __name__ == "__main__":
    main()