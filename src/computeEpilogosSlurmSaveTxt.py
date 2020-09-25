import gzip
import numpy as np
from pathlib import Path
import pandas as pd
import time
import click
import os
import subprocess
from pathlib import PurePath

@click.command()
@click.option("-f", "--file-directory", "fileDirectory", type=str, required=True, help="Path to directory that contains files to read from (All files in this directory will be read in)")
@click.option("-s", "--state-model", "numStates", type=int, required=True, help="Number of states in chromatin state model")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-l", "--saliency-level", "saliency", type=int, default=1, show_default=True, help="Desired saliency level (1, 2, or 3)")
@click.option("-m", "--mode-of-operation", "modeOfOperation", type=click.Choice(["bg", "s", "both"]), default="both", show_default=True, help="bg for background, s for scores, both for both")
@click.option("-b", "--background-directory", "expFreqDir", type=str, default="null", help="Path to where the background frequency array is read from (-m s) or written to (-m bg, -m both) [default: output-directory]") # default output directory
def main(fileDirectory, numStates, saliency, outputDirectory, modeOfOperation, expFreqDir):
    """
    This script computes scores for chromatin states across the genome.
    """
    tTotal = time.time()
    dataFilePath = Path(fileDirectory)
    outputDirPath = Path(outputDirectory)

    fileTag = "_".join(str(dataFilePath).split("/")[-5:])

    print()
    print("Input Directory =", dataFilePath)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)
    print("Mode of Operation =", modeOfOperation)
    print("Background Directory =", expFreqDir)

    if expFreqDir == "null":
        expFreqDir = outputDirectory

    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath

    if not PurePath(dataFilePath).is_absolute():
        dataFilePath = Path.cwd() / dataFilePath

    # Check that paths are valid before doing anything
    if not dataFilePath.exists() or not dataFilePath.is_dir():
        print()
        print("ERROR: Given file path does not exist or is not a directory")
        print()
        return

    if not list(dataFilePath.glob("*")):
        print()
        print("ERROR: Ensure that file directory is not empty")
        print()
        return

    if saliency != 1 and saliency != 2 and saliency != 3:
        print()
        print("ERROR: Ensure that saliency metric is either 1, 2, or 3")
        print()
        return

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)
    
    # For slurm output and error later
    (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)

    if not outputDirPath.is_dir():
        print("ERROR: Output directory is not a directory")
        return

    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to path of the input file directory
    expFreqFilename = "exp_freq_{}.npy".format(fileTag)
    storedExpPath = Path(expFreqDir) / expFreqFilename
    print()
    print("Background Frequency Array Location:", storedExpPath)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    # Check if user wants to calculate it
    if modeOfOperation == "s":
        try:
            expFreqArr = np.load(storedExpPath, allow_pickle=False)
        except IOError:
            print("ERROR: Could not load stored expected value array.\n\tPlease check that the directory is correct or that the file exists")
    else:     
        expJobIDArr = []   
        print()
        print("Submitting Slurm Jobs for Per Datafile Background Frequency Calculation....")
        for file in dataFilePath.glob("*"):
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "exp_freq_calc_{}_{}".format(fileTag, filename)
                jobOutPath = outputDirPath / (".out/" + jobName + ".out")
                jobErrPath = outputDirPath / (".err/" + jobName + ".err")

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

                computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

                pythonCommand = "python {} {} {} {} {} {}".format(computeExpectedPy, file, numStates, saliency, outputDirPath, fileTag)

                if saliency == 1:
                    slurmCommand = "sbatch --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)
                elif saliency == 2:
                    slurmCommand = "sbatch --job-name=S2_{}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=2000 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)
                elif saliency == 3:
                    slurmCommand = "sbatch --job-name=S3_{}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=64000 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    print("ERROR: sbatch not submitted correctly")
                
                expJobIDArr.append(int(sp.stdout.split()[-1]))

        # Combining all the different chromosome expected frequency arrays into one
        # create a string for slurm dependency to work
        jobIDStrComb = str(expJobIDArr).strip('[]').replace(" ", "")

        print("    JobIDs:", jobIDStrComb)

        print()
        print("Submitting Slurm Job for Combining Background Frequency Arrays....")

        jobName = "exp_freq_comb_{}".format(fileTag)
        jobOutPath = outputDirPath / (".out/" + jobName + ".out")
        jobErrPath = outputDirPath / (".err/" + jobName + ".err")

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

        computeExpectedCombinationPy = pythonFilesDir / "computeEpilogosExpectedCombination.py"

        pythonCommand = "python {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, fileTag, storedExpPath)

        if saliency == 1:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --wrap='{}'".format(jobIDStrComb, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 2:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(jobIDStrComb, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 3:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name=S3_{}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=64000 --wrap='{}'".format(jobIDStrComb, jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            print("ERROR: sbatch not submitted correctly")
        
        combinationJobID = int(sp.stdout.split()[-1])

        print("    JobID:", combinationJobID)

    if modeOfOperation == "s" or modeOfOperation == "both":
        # Calculate the observed frequencies and scores
        print()
        print("Submitting Slurm Jobs for Score Calculation....")
        scoreJobIDArr = []
        for file in dataFilePath.glob("*"):
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "score_calc_{}_{}".format(fileTag, filename)
                jobOutPath = outputDirPath / (".out/" + jobName + ".out")
                jobErrPath = outputDirPath / (".err/" + jobName + ".err")

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
                
                computeScorePy = pythonFilesDir / "computeEpilogosScores.py"

                pythonCommand = "python {} {} {} {} {} {} {}".format(computeScorePy, file, numStates, saliency, outputDirPath, storedExpPath, fileTag)

                if saliency == 1:
                    slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, pythonCommand)
                elif saliency == 2:
                    slurmCommand = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, pythonCommand)
                elif saliency == 3:
                    slurmCommand = "sbatch --dependency=afterok:{} --job-name=S3_{}.job --output={} --error={} --nodes=1 --ntasks=1 --cpus-per-task=4 --mem-per-cpu=32000 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    print("ERROR: sbatch not submitted correctly")
                
                scoreJobIDArr.append(int(sp.stdout.split()[-1]))

        # create a string for slurm dependency to work
        jobIDStrWrite = str(scoreJobIDArr).strip('[]').replace(" ", "")
        
        print("    JobIDs:", jobIDStrWrite)

        # WRITING TO SCORE FILE
        print()
        print("Submitting Slurm Jobs for Writing to Score Files....")

        jobName = "writeSaveTxt_{}".format(fileTag)
        jobOutPath = outputDirPath / (".out/" + jobName + ".out")
        jobErrPath = outputDirPath / (".err/" + jobName + ".err")

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

        computeExpectedWritePy = pythonFilesDir / "computeEpilogosWriteSaveTxt.py"

        pythonCommand = "python {} {} {} {}".format(computeExpectedWritePy, fileTag, outputDirPath, numStates)

        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --mail-type=END --mail-user=jquon@altius.org --mem-per-cpu=32000 --wrap='{}'".format(jobIDStrWrite, jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            print("ERROR: sbatch not submitted correctly")

        print("    JobID:", int(sp.stdout.split()[-1]))

        print()
        print("All JobIDs: {}".format(jobIDStrComb + "," + combinationJobID + "," + jobIDStrWrite + "," + sp.stdout.split()[-1]))

if __name__ == "__main__":
    main()