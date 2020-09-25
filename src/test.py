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
    
    test = np.array(np.arange(30)).reshape(10,3)

    testlist = [1, 3, 5, 9]

    print(test)
    print(test[testlist])

    print(np.sum(test, axis=1))

    print( test / 2)

if __name__ == "__main__":
    main()