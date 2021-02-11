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
@click.option("-f", "--file-directory", "fileDirectory", type=str, required=True, multiple=True, help="Path to directory that contains files to read from (All files in this directory will be read in)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, multiple=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-s", "--state-model", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-l", "--saliency-level", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1, 2, or 3)")
@click.option("-m", "--mode-of-operation", "modeOfOperation", type=click.Choice(["bg", "s", "both"]), default=["both"], show_default=True, multiple=True, help="bg for background, s for scores, both for both")
@click.option("-b", "--background-directory", "expFreqDir", type=str, default=["null"], multiple=True, help="Path to where the background frequency array is read from (-m s) or written to (-m bg, -m both) [default: output-directory]")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-x", "--exit-when-complete", "exitBool", is_flag=True, multiple=True, help="If flag is enabled program will exit upon completion of all processes rather than SLURM submission of all processes")
def main(fileDirectory, outputDirectory, numStates, saliency, modeOfOperation, expFreqDir, numProcesses, exitBool):
    """
    This script computes scores for chromatin states across the genome.
    """

    # Handling case if user inputs flag multiples times
    if len(fileDirectory) > 1:
        raise ValueError("Too many [-f, --file-directory] arguments provided")
    elif len(outputDirectory) > 1:
        raise ValueError("Too many [-o, --output-directory] arguments provided")
    elif len(numStates) > 1:
        raise ValueError("Too many [-s, --state-model] arguments provided")
    elif len(saliency) > 1:
        raise ValueError("Too many [-l, --saliency-level] arguments provided")
    elif len(modeOfOperation) > 1:
        raise ValueError("Too many [-m, --mode-of-operation] arguments provided")
    elif len(expFreqDir) > 1:
        raise ValueError("Too many [-b, --background-directory] arguments provided")
    elif len(numProcesses) > 1:
        raise ValueError("Too many [-c, --num-cores] arguments provided")
    elif len(exitBool) > 1:
        raise ValueError("Too many [-x, --exit-when-complete] arguments provided")
    fileDirectory = fileDirectory[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    modeOfOperation = modeOfOperation[0]
    expFreqDir = expFreqDir[0]
    numProcesses = numProcesses[0]
    if exitBool:
        exitBool = exitBool[0]


    dataFilePath = Path(fileDirectory)
    outputDirPath = Path(outputDirectory)

    print()
    print("Input Directory =", dataFilePath)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)
    print("Mode of Operation =", modeOfOperation)
    print("Background Directory =", expFreqDir)

    # If user does not specificy a directory to look for expected frequencies default the output directory
    if expFreqDir == "null":
        expFreqDir = outputDirectory

    # Making paths absolute
    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath
    if not PurePath(dataFilePath).is_absolute():
        dataFilePath = Path.cwd() / dataFilePath

    # For making sure all files are consistently named
    fileTag = "{}_saliency{}".format(dataFilePath.name, saliency)

    if saliency != 1 and saliency != 2 and saliency != 3:
        print("\nERROR: Ensure that saliency metric is either 1, 2, or 3\n")
        return

    # Check that paths are valid before doing anything
    if not dataFilePath.exists() or not dataFilePath.is_dir():
        print("\nERROR: Given file path does not exist or is not a directory\n")
        return

    if not list(dataFilePath.glob("*")):
        print("\nERROR: Ensure that file directory is not empty\n")
        return

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)
    
    if not outputDirPath.is_dir():
        print("\nERROR: Output directory is not a directory\n")
        return

    if numProcesses < 0:
        print("\nERROR: Number of cores must be positive or zero (0 means use all cores)\n")
        return
    elif numProcesses == 0:
        numTasks = "--exclusive"
    else:
        numTasks = "--ntasks={}".format(numProcesses)

    # For slurm output and error later
    (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)

    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to path of the input file directory
    storedExpPath = Path(expFreqDir) / "exp_freq_{}.npy".format(fileTag)
    print("\nBackground Frequency Array Location:", storedExpPath)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    # Only calculate the expected frequencies if user asks for it, otherwise just load from where the user said
    if modeOfOperation == "s":
        try:
            expFreqArr = np.load(storedExpPath, allow_pickle=False)
        except IOError:
            print("\nERROR: Could not load stored expected value array.\n\tPlease check that the directory is correct or that the file exists\n")
            return
    else:     
        expJobIDArr = []   
        print("\nSubmitting Slurm Jobs for Per Datafile Background Frequency Calculation....")
        for file in dataFilePath.glob("*"):
            # Skip over ".genome" files
            if file.name.split(".")[-1] == "genome":
                continue
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "exp_calc_{}_{}".format(fileTag, filename)
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
                    print("\nERROR: sbatch '.out' or '.err' file already exists\n")
                    return

                computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

                pythonCommand = "python {} {} {} {} {} {} {}".format(computeExpectedPy, file, numStates, saliency, outputDirPath, fileTag, numProcesses)

                if saliency == 1:
                    slurmCommand = "sbatch --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                elif saliency == 2:
                    slurmCommand = "sbatch --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                elif saliency == 3:
                    slurmCommand = "sbatch --job-name=S3_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    print("\nERROR: sbatch not submitted correctly\n")
                    return
                
                expJobIDArr.append(int(sp.stdout.split()[-1]))

        # Combining all the different chromosome expected frequency arrays into one
        # create a string for slurm dependency to work
        expJobIDStr = str(expJobIDArr).strip('[]').replace(" ", "")

        print("    JobIDs:", expJobIDStr)

        print("\nSubmitting Slurm Job for Combining Background Frequency Arrays....")

        jobName = "exp_comb_{}".format(fileTag)
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
            print("\nERROR: sbatch '.out' or '.err' file already exists\n")
            return

        computeExpectedCombinationPy = pythonFilesDir / "computeEpilogosExpectedCombination.py"

        pythonCommand = "python {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, storedExpPath, fileTag)

        if saliency == 1:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name=S1_{}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 2:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 3:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name=S3_{}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            print("\nERROR: sbatch not submitted correctly\n")
            return
        
        combinationJobID = int(sp.stdout.split()[-1])

        print("    JobID:", combinationJobID)

    if modeOfOperation == "s" or modeOfOperation == "both":
        # Calculate the observed frequencies and scores
        print("\nSubmitting Slurm Jobs for Score Calculation....")
        scoreJobIDArr = []
        for file in dataFilePath.glob("*"):
            # Skip over ".genome" files
            if file.name.split(".")[-1] == "genome":
                continue
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "score_{}_{}".format(fileTag, filename)
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
                    print("\nERROR: sbatch '.out' or '.err' file already exists\n")
                    return
                
                computeScorePy = pythonFilesDir / "computeEpilogosScores.py"

                pythonCommand = "python {} {} {} {} {} {} {} {}".format(computeScorePy, file, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses)

                if modeOfOperation == "s":
                    if saliency == 1:
                        slurmCommand = "sbatch --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 2:
                        slurmCommand = "sbatch --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 3:
                        slurmCommand = "sbatch --job-name=S3_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                else:
                    if saliency == 1:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 2:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 3:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name=S3_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    print("\nERROR: sbatch not submitted correctly\n")
                    return
                
                scoreJobIDArr.append(int(sp.stdout.split()[-1]))

        # Create a string for slurm dependency to work
        scoreJobIDStr = str(scoreJobIDArr).strip('[]').replace(" ", "")
        
        print("    JobIDs:", scoreJobIDStr)

        # Write scores out to gzipped text files
        print("\nSubmitting Slurm Jobs for Writing to Score Files....")
        writeJobIDArr = []
        for file in dataFilePath.glob("*"):
            # Skip over ".genome" files
            if file.name.split(".")[-1] == "genome":
                continue
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "write_{}_{}".format(fileTag, filename)
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
                    print("\nERROR: sbatch '.out' or '.err' file already exists\n")
                    return

                computeExpectedWritePy = pythonFilesDir / "computeEpilogosWrite.py"

                pythonCommand = "python {} {} {} {} {}".format(computeExpectedWritePy, file, numStates, outputDirPath, fileTag)

                slurmCommand = "sbatch --dependency=afterok:{} --job-name=S{}_{}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=64000 --wrap='{}'".format(scoreJobIDStr, saliency, jobName, jobOutPath, jobErrPath, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    print("\nERROR: sbatch not submitted correctly\n")
                    return

                writeJobIDArr.append(int(sp.stdout.split()[-1]))

        writeJobIDStr = str(writeJobIDArr).strip('[]').replace(" ", "")
        print("    JobIDs:", writeJobIDStr)
        if modeOfOperation == "both":
            print("\nAll JobIDs: {},{},{},{}".format(expJobIDStr, combinationJobID, scoreJobIDStr, writeJobIDStr))
        elif modeOfOperation == "s":
            print("\nAll JobIDs: {},{}".format(scoreJobIDStr, writeJobIDStr))
        elif modeOfOperation == "bg":
            print("\nAll JobIDs: {},{}".format(expJobIDStr, combinationJobID))

    # If the user wants to exit upon job completion rather than submission
    if exitBool:
        # The last job is different depending on mode of operation
        if modeOfOperation == "s" or modeOfOperation == "both":
            lastJob = writeJobIDStr
        else:
            lastJob = combinationJobID

        lastJobCheckStr = "sacct --format=State --jobs {}".format(lastJob)

        # Every ten seconds check if the final job is done, if it is exit the program
        while True:
            time.sleep(10)
            sp = subprocess.run(lastJobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
            if not (sp.stdout.strip().contains("RUNNING") or sp.stdout.strip().contains("PENDING")):
                break
                
if __name__ == "__main__":
    main()