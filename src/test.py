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

def main(fileDir, outputDir):
    dataFilePath = Path(fileDir)
    outputDirPath = Path(outputDir)

    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath
        print("OUTPUTPATH: ", outputDirPath)

    if not PurePath(dataFilePath).is_absolute():
        dataFilePath = Path.cwd() / dataFilePath
        print("FILE PATH: ", dataFilePath)

    # For slurm output and error later
    (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)

    # Finding the location of the .py files that must be run
    if Path(__file__).is_absolute:
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    for file in dataFilePath.glob("*"):
        filename = file.name.split(".")[0]
        jobName = "exp_freq_calc_{}".format(filename)
        jobOutPath = outputDirPath / (".out/" + jobName + ".out")
        jobErrPath = outputDirPath / (".err/" + jobName + ".err")

        # Creating the out and err files for the batch job
        try:
            jout = open(jobOutPath, 'x')
            jout.close()
            jerr = open(jobErrPath, 'x')
            jerr.close()
        except FileExistsError:
            print("ERROR: sbatch '.out' or '.err' file already exists")
        
        computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

        pythonCommand = "python {} {} {} {} {}".format(computeExpectedPy, file, 15, 1, outputDirPath)
        slurmCommand = "sbatch --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)

        print(pythonCommand)
        print()
        print(slurmCommand)
        print()

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True)

        if not sp.stdout.startswith("Submitted batch"):
            print("ERROR: sbatch not submitted correctly")
        
        print("JOBID: ", int(sp.stdout.split()[-1]))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])