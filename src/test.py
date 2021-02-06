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
    dataFilePath = Path("/home/jquon/RoadmapStateByGroup/male/split/")
    
    for file in dataFilePath.glob("*"):
        print(file)
        # Skip over ".genome" files
        if file.name.split(".")[1] == "genome":
            print("GENOME")

if __name__ == "__main__":
    main()