import sys
from pathlib import Path
from pathlib import PurePath
import os
import subprocess
import scipy.stats as st
import pandas as pd
import numpy as np

def main(file1, file2, outputDir, a, b, loc, scale):
    print("Submitting Slurm Jobs....")

    outputDirPath = Path(outputDir)
    file1Path = Path(file1)
    file2Path = Path(file2)

    (outputDirPath / "out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / "err/").mkdir(parents=True, exist_ok=True)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    dataLength = pd.read_table(file1Path, header=None, sep="\s+", usecols=[0]).shape[0]

    jobIDStr = ""
    # for i in np.arange(0, dataLength, 500000):
    for i in range(1):
        jobIDArr = []
        print(jobIDStr)
        # for j in np.arange(i, i+500000, 5000):
        for j in range(1):
            jobName = "pval_{}".format(j)
            jobOutPath = outputDirPath / ("out/" + jobName + ".out")
            jobErrPath = outputDirPath / ("err/" + jobName + ".err")

            # Creating the out and err files for the batch job
            if jobOutPath.exists():
                os.remove(jobOutPath)
            if jobErrPath.exists():
                os.remove(jobErrPath)
            try:
                jout = open(jobOutPath, 'x')
                jout.close()
                jerr = open(jobErrPath, 'x')
                jerr.close()
            except FileExistsError:
                # This error should never occur because we are deleting the files first
                print("ERROR: sbatch '.out' or '.err' file already exists")

            computePvalsPy = pythonFilesDir / "computePvals.py"

            pythonCommand = "python {} {} {} {} {} {} {} {} {}".format(computePvalsPy, file1, file2, outputDir, a, b, loc, scale, j)

            if i == 0:
                slurmCommand = "sbatch --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=40000 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)
            else:
                slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=40000 --wrap='{}'".format(jobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

            sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

            if not sp.stdout.startswith("Submitted batch"):
                print("ERROR: sbatch not submitted correctly")

            jobIDArr.append(int(sp.stdout.split()[-1]))
        jobIDStr = str(jobIDArr).strip('[]').replace(" ", "")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])