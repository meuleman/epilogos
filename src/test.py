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

    test = "hiya"

    print(Path(__file__) / (".out/" + test + ".out"))

    print(Path.cwd())

    print(test[2])

    dataFilePath = Path(fileDir)
    outputDirPath = Path(outputDir)

    # Finding the location of the .py files that must be run
    if Path(__file__).is_absolute:
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    for file in dataFilePath.glob("*"):
        print(file)
        # if not file.is_dir():
        #     filename = file.name.split(".")[0]
        #     jobName = "exp_freq_calc_{}".format(filename)
        #     computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

        #     # pythonCommand = "python {} {} {} {} {}".format(computeExpectedPy, file, 15, 1, outputDirPath)
        #     # slurmCommand = "sbatch --job-name={0}.job --output=.out/{0}.out --error=.out/{0}.err --nodes=1 --ntasks=1 --wrap='{1}'".format(jobName, pythonCommand)

        #     # process = subprocess.run(slurmCommand, shell=True, universal_newlines=True)

        #     # print("STDOUT", process.stdout)

        #     jobPath = outputDirPath / "{}.sh".format(jobName)
        #     with open(jobPath, "w") as script:
        #         script.writelines("#!/bin/bash\n")
        #         script.writelines("#SBATCH --job-name={}.job\n".format(jobName))
        #         script.writelines("#SBATCH --output=.out/{}.out\n".format(jobName))
        #         script.writelines("#SBATCH --error=.out/{}.err\n".format(jobName))
        #         script.writelines("#SBATCH --nodes=1\n".format())
        #         script.writelines("#SBATCH --ntasks=1\n".format())
        #         script.writelines("python {} {} {} {} {}".format(computeExpectedPy, file, 15, 1, outputDirPath))

        #     process = subprocess.run("sbatch {}".format(jobPath), shell=True, universal_newlines=True)

        #     print("STDOUT: ", process.stdout)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])