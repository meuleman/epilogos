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


def main():
    t = time.time()

    with gzip.open(Path("/home/jquon/AdseraStateByGroup/male/matrix.txt.gz"), 'rb') as f:
        for i, l in enumerate(f):
            pass
    print("File {1} contain {0} lines".format(i, "matrix"))

    print(time.time() - t)


    t = time.time()
    number_lines = sum(1 for line in open(Path("/home/jquon/AdseraStateByGroup/male/matrix.txt.gz")))
    print(number_lines)
    print(time.time() - t)

if __name__ == "__main__":
    main()