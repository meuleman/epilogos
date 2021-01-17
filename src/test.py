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

    scoreArr = multiprocessing.Array(c.c_float, 10)

    print(np.all(scoreArr == 0))

    for i in range(10):
        print(scoreArr[i])

if __name__ == "__main__":
    main()