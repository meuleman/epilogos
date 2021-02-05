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
import numpy.ma as ma
import click

def main():
    bedOne = Path("/home/jquon/epilogosTesting_01122021/output/pairwiseMF/largestDistanceLoci_split_split.bed")
    bedTwo = Path("/home/jquon/newPairwiseTest/maleFemale/largestDistanceLoci_split_split.bed")

    with open(bedOne, "r") as f1:
        with open(bedTwo, "r") as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
            for i in range(len(lines1)):
                line1Split = lines1[i].split()
                line2Split = lines2[i].split()

                if line1Split[-2] != line2Split[-2]:
                    print(i, line1Split[-2], line2Split[-2])

if __name__ == "__main__":
    main()