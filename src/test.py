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
    tRead = time.time()
    # Taking in the the score array
    filePath = Path("/home/jquon/chromHMMConversionTest/18state_b'chr3'.npz")
    npzFile = np.load(filePath)
    scoreArr = npzFile['scoreArr']
    locationArr = npzFile['locationArr']

    print(scoreArr[:10])
    print(locationArr[:10])

    print("Read Time:", time.time() - tRead)

if __name__ == "__main__":
    main()