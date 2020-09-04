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

def main(fileDir, outputDir):

    dataFilePath = Path(fileDir)
    outputPath = Path(outputDir)

    # Finding the location of the .py files that must be run
    if Path(__file__).is_absolute:
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    for file in dataFilePath.glob("*"):
        if not file.is_dir():
            filename = file.name.split(".")[0]
            jobName = "exp_freq_calc_{}".format(filename)
            file = dataFilePath / file
            computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

            pythonCommand = "python {} {} {} {} {}".format(computeExpectedPy, file, 15, 1, outputPath)
            slurmCommand = "sbatch --uid=jquon --job-name={0}.job --output=.out/{0}.out --error=.out/{0}.err --nodes=1 --ntasks=1 --wrap='{1}'".format(jobName, pythonCommand)

            process = subprocess.run(slurmCommand, shell=True, universal_newlines=True)

            print("STDOUT", process.stdout)

            # slurmCheck = subprocess.run("")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])