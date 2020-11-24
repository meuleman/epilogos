import sys
import time
from pathlib import Path
from pathlib import PurePath
import os
import subprocess

def main(inputDir1, inputDir2, outputDir1, outputDir2):
    inputDir1Path = Path(inputDir1)
    inputDir2Path = Path(inputDir2)
    outputDir1Path = Path(outputDir1)
    outputDir2Path = Path(outputDir2)
    
    # out and error can only go to one file so make it the first given output directory
    (outputDir1Path / "out/").mkdir(parents=True, exist_ok=True)
    (outputDir1Path / "err/").mkdir(parents=True, exist_ok=True)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    print()
    print("Submitting Slurm Jobs....")
    jobIDArr = []
    for file1 in inputDir1Path.glob("*"):
        # Only want the matching part of the to genome
        filename = file1.name.split(".")[0]
        for file2 in inputDir2Path.glob(filename + "*"):
            jobName = "shuffle_{}".format(filename)
            jobOutPath = outputDir1Path / ("out/" + jobName + ".out")
            jobErrPath = outputDir1Path / ("err/" + jobName + ".err")

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

            fitDistributionPy = pythonFilesDir / "fitDistribution.py"

            pythonCommand = "python {} {} {} {} {} {} {} {} {} {}".format(fitDistributionPy, file1, file2, observationFile, filterBool, fullBool, distNum, binEnd, numStates, outputDir)

            slurmCommand = "sbatch --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=16000 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)

            sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

            if not sp.stdout.startswith("Submitted batch"):
                print("ERROR: sbatch not submitted correctly")
                
                jobIDArr.append(int(sp.stdout.split()[-1]))

    jobIDStr = str(jobIDArr).strip('[]').replace(" ", "")

    print("    JobIDs:", jobIDStr)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
