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
    file = Path("C:/Users/User/Desktop/epilogos/yout.txt")

    with open(file, "r") as f:
        lines = f.readlines()
        for line in lines:
            if int(line.split()[-1]) > 10:
                print(line)

if __name__ == "__main__":
    main()