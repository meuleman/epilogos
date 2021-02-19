import numpy as np
from pathlib import Path
import pandas as pd
import time
import click
import os
import subprocess
from pathlib import PurePath
import errno

from pandas.core.algorithms import value_counts

@click.command()
@click.option("-i", "--input-directory", "inputDirectory", type=str, required=True, multiple=True, help="Path to directory that contains files to read from (All files in this directory will be read in)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, multiple=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-n", "--number-of-states", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency-level", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1, 2, or 3)")
@click.option("-m", "--mode-of-operation", "modeOfOperation", type=click.Choice(["bg", "s", "both"]), default=["both"], show_default=True, multiple=True, help="bg for background, s for scores, both for both")
@click.option("-b", "--background-directory", "expFreqDir", type=str, default=["null"], multiple=True, help="Path to where the background frequency array is read from (-m s) or written to (-m bg, -m both) [default: output-directory]")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-x", "--exit-when-complete", "exitBool", is_flag=True, multiple=True, help="If flag is enabled program will exit upon submission of all SLURM processes rather than completion of all processes")
def main(inputDirectory, outputDirectory, numStates, saliency, modeOfOperation, expFreqDir, numProcesses, exitBool):
    """
    This script computes scores for chromatin states across the genome.
    """

    # Handling case if user inputs flag multiples times
    if len(inputDirectory) > 1:
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
    inputDirectory = inputDirectory[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    modeOfOperation = modeOfOperation[0]
    expFreqDir = expFreqDir[0]
    numProcesses = numProcesses[0]
    # if exitBool:
    #     exitBool = not exitBool[0]
    verbose=True

    inputDirPath = Path(inputDirectory)
    outputDirPath = Path(outputDirectory)

    print()
    print("Input Directory =", inputDirPath)
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
    if not PurePath(inputDirPath).is_absolute():
        inputDirPath = Path.cwd() / inputDirPath

    # For making sure all files are consistently named
    fileTag = "{}_s{}".format(inputDirPath.name, saliency)

    if saliency != 1 and saliency != 2 and saliency != 3:
        raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")

    # Check that paths are valid before doing anything
    if not inputDirPath.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath)))
    if not inputDirPath.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath)))
    if not list(inputDirPath.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty:", str(inputDirPath))

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)
    
    if not outputDirPath.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDirPath)))

    if numProcesses < 0:
        raise ValueError("Number of cores must be positive or zero (0 means use all cores)")
    elif numProcesses == 0:
        numTasks = "--exclusive"
    else:
        numTasks = "--ntasks={}".format(numProcesses)

    # For slurm output and error later
    (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)
    print("Slurm .out log files are located at: {}".format(outputDirPath / ".out/"))
    print("Slurm .err log files are located at: {}".format(outputDirPath / ".err/"))


    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to path of the input file directory
    storedExpPath = Path(expFreqDir) / "exp_freq_{}.npy".format(fileTag)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    # Only calculate the expected frequencies if user asks for it, otherwise just load from where the user said
    if modeOfOperation == "s":
        try:
            expFreqArr = np.load(storedExpPath, allow_pickle=False)
        except IOError as err:
            print(err)
            return
    else:     
        expJobIDArr = []   
        print("\nSubmitting Slurm Jobs for Per Datafile Background Frequency Calculation....")
        for file in inputDirPath.glob("*"):
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
                except FileExistsError as err:
                    # This error should never occur because we are deleting the files first
                    print(err)
                    return

                computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

                pythonCommand = "python {} {} {} {} {} {} {} {}".format(computeExpectedPy, file, numStates, saliency, outputDirPath, fileTag, numProcesses, verbose)

                if saliency == 1:
                    slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                elif saliency == 2:
                    slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                elif saliency == 3:
                    slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    raise ChildProcessError("SlurmError: sbatch not submitted correctly")
                
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
        except FileExistsError as err:
            # This error should never occur because we are deleting the files first
            print(err)
            return

        computeExpectedCombinationPy = pythonFilesDir / "computeEpilogosExpectedCombination.py"

        pythonCommand = "python {} {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, storedExpPath, fileTag, verbose)

        if saliency == 1:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 2:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 3:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            raise ChildProcessError("SlurmError: sbatch not submitted correctly")
        
        combinationJobID = int(sp.stdout.split()[-1])

        print("    JobID:", combinationJobID)

    if modeOfOperation == "s" or modeOfOperation == "both":
        # Calculate the observed frequencies and scores
        print("\nSubmitting Slurm Jobs for Score Calculation....")
        scoreJobIDArr = []
        for file in inputDirPath.glob("*"):
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
                except FileExistsError as err:
                    # This error should never occur because we are deleting the files first
                    print(err)
                    return
                
                computeScorePy = pythonFilesDir / "computeEpilogosScores.py"

                pythonCommand = "python {} {} {} {} {} {} {} {} {}".format(computeScorePy, file, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses, verbose)

                if modeOfOperation == "s":
                    if saliency == 1:
                        slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 2:
                        slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 3:
                        slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                else:
                    if saliency == 1:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 2:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 3:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    raise ChildProcessError("SlurmError: sbatch not submitted correctly")
                
                scoreJobIDArr.append(int(sp.stdout.split()[-1]))

        # Create a string for slurm dependency to work
        scoreJobIDStr = str(scoreJobIDArr).strip('[]').replace(" ", "")
        
        print("    JobIDs:", scoreJobIDStr)

        # Write scores out to gzipped text files
        print("\nSubmitting Slurm Jobs for Writing to Score Files....")
        writeJobIDArr = []
        for file in inputDirPath.glob("*"):
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
                except FileExistsError as err:
                    # This error should never occur because we are deleting the files first
                    print(err)
                    return

                computeExpectedWritePy = pythonFilesDir / "computeEpilogosWrite.py"

                pythonCommand = "python {} {} {} {} {} {}".format(computeExpectedWritePy, file, numStates, outputDirPath, fileTag, verbose)

                slurmCommand = "sbatch --dependency=afterok:{} --job-name=S{}_{}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=64000 --wrap='{}'".format(scoreJobIDStr, saliency, jobName, jobOutPath, jobErrPath, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    raise ChildProcessError("SlurmError: sbatch not submitted correctly")

                writeJobIDArr.append(int(sp.stdout.split()[-1]))

        writeJobIDStr = str(writeJobIDArr).strip('[]').replace(" ", "")
        print("    JobIDs:", writeJobIDStr)
        if modeOfOperation == "both":
            allJobIDs = "{},{},{},{}".format(expJobIDStr, combinationJobID, scoreJobIDStr, writeJobIDStr)
            print("\nAll JobIDs:", allJobIDs)
        elif modeOfOperation == "s":
            allJobIDs = "{},{}".format(scoreJobIDStr, writeJobIDStr)
            print("\nAll JobIDs:", allJobIDs)
        elif modeOfOperation == "bg":
            allJobIDs = "{},{}".format(expJobIDStr, combinationJobID)
            print("\nAll JobIDs:", allJobIDs)

    # If the user wants to exit upon job completion rather than submission
    # If a job fails, it cancels all other jobs
    if not exitBool:
        jobCheckStr = "sacct --format=JobID,JobName,State --jobs {}".format(allJobIDs)

        # Run the job check once before the while loop to get the info lines
        sp = subprocess.run(jobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
        spLines = sp.stdout.split("\n")
        # Print out the info lines
        print(spLines[0])
        print(spLines[1])
        completedJobs = []

        # Every ten seconds check if the final job is done, if it is exit the program
        while True:
            time.sleep(10)
            # Check if there was an error, if so cancel everything and exit the program
            if "FAILED" in sp.stdout or "CANCELLED" in sp.stdout:
                print("\nERROR RUNNING JOBS: CANCELLING ALL REMAINING JOBS\n")
                print("Please check error logs in: {}/.err/\n".format(outputDirPath))
                subprocess.run("scancel {}".format(allJobIDs), shell=True)
                break
            # If final job is done, exit the program
            if not ("RUNNING" in sp.stdout or "PENDING" in sp.stdout):
                break
            # Print out jobs when they are completed
            for line in spLines[2:]:
                if "COMPLETED" in line:
                    jobID = line.split()[0]
                    # Don't want to print if we have already printed
                    if jobID not in completedJobs:
                        completedJobs.append(jobID)
                        print(line)

            sp = subprocess.run(jobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)


            

                
if __name__ == "__main__":
    main()