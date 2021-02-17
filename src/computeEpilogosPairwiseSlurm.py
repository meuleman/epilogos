import gzip
from multiprocessing import Value
import numpy as np
from pathlib import Path
import pandas as pd
import time
import click
import os
import subprocess
from pathlib import PurePath
import errno

@click.command()
@click.option("-a", "--directory-one", "inputDirectory1", type=str, required=True, multiple=True, help="Path to directory that contains files to read from (All files in this directory will be read in)")
@click.option("-b", "--directory-two", "inputDirectory2", type=str, required=True, multiple=True, help="")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, multiple=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-n", "--number-of-states", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency-level", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1 or 2)")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-x", "--exit-when-complete", "exitBool", is_flag=True, multiple=True, help="If flag is enabled program will exit upon completion of all processes rather than SLURM submission of all processes")
def main(inputDirectory1, inputDirectory2, outputDirectory, numStates, saliency, numProcesses, exitBool):
    """
    This script computes and visualizes differences between chromatin states across epigenomes
    """
    # Handling case if user inputs flag multiples times
    if len(inputDirectory1) > 1:
        raise ValueError("Too many [-i, --directory-one] arguments provided")
    elif len(inputDirectory2) > 1:
        raise ValueError("Too many [-n, --directory-two] arguments provided")
    elif len(outputDirectory) > 1:
        raise ValueError("Too many [-o, --output-directory] arguments provided")
    elif len(numStates) > 1:
        raise ValueError("Too many [-s, --state-model] arguments provided")
    elif len(saliency) > 1:
        raise ValueError("Too many [-l, --saliency-level] arguments provided")
    elif len(numProcesses) > 1:
        raise ValueError("Too many [-c, --num-cores] arguments provided")
    elif len(exitBool) > 1:
        raise ValueError("Too many [-x, --exit-when-complete] arguments provided")
    inputDirectory1 = inputDirectory1[0]
    inputDirectory2 = inputDirectory2[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    numProcesses = numProcesses[0]
    if exitBool:
        exitBool = exitBool[0]

    inputDirPath1 = Path(inputDirectory1)
    inputDirPath2 = Path(inputDirectory2)
    outputDirPath = Path(outputDirectory)

    print()
    print("File Directory 1 =", inputDirPath1)
    print("File Directory 2 =", inputDirPath2)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)

    if saliency != 1 and saliency != 2 and saliency != 3:
        raise ValueError("Saliency Metric Invalid: {} Please ensure that saliency metric is either 1 or 2 (Saliency of 3 is unsupported for pairwise comparison".format(saliency))

    # Making paths absolute
    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath
    if not PurePath(inputDirPath1).is_absolute():
        inputDirPath1 = Path.cwd() / inputDirPath1
    if not PurePath(inputDirPath2).is_absolute():
        inputDirPath2 = Path.cwd() / inputDirPath2

    # For making sure all files are consistently named
    fileTag = "{}_{}_s{}".format(inputDirPath1.name, inputDirPath2.name, saliency)

    # Check that paths are valid before doing anything
    if not inputDirPath1.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath1)))
    if not inputDirPath1.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath1)))
    if not list(inputDirPath1.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty:", str(inputDirPath1))
    if not inputDirPath2.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath2)))
    if not inputDirPath2.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath2)))
    if not list(inputDirPath2.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty:", str(inputDirPath2))

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

    # Path for storing/retrieving the expected frequency array
    storedExpPath = outputDirPath / "exp_freq_{}.npy".format(fileTag)
    print("\nBackground Frequency Array Location:", storedExpPath)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    expJobIDArr = []   
    print("\nSubmitting Slurm Jobs for Per Datafile Background Frequency Calculation....")
    for file1 in inputDirPath1.glob("*"):
        # Skip over ".genome" files
        if file1.name.split(".")[-1] == "genome":
                continue
        # Find matching file in other directory
        if not list(inputDirPath2.glob(file1.name)):
            raise FileNotFoundError("File not found: {} Please ensure corresponding files within input directories directories 1 and 2 have the same name".format(str(inputDirPath2 / file1.name)))
        else:
            file2 = next(inputDirPath2.glob(file1.name))

        if not file1.is_dir() and not file2.is_dir():
            filename = file1.name.split(".")[0]
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

            # Create a string for the python command
            computeExpectedPy = pythonFilesDir / "computeEpilogosPairwiseExpected.py"
            pythonCommand = "python {} {} {} {} {} {} {} {}".format(computeExpectedPy, file1, file2, numStates, saliency, outputDirPath, fileTag, numProcesses)

            # Create a string for the slurm command
            if saliency == 1:
                slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
            elif saliency == 2:
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

    # Create a string for the python commmand
    computeExpectedCombinationPy = pythonFilesDir / "computeEpilogosExpectedCombination.py"
    pythonCommand = "python {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, storedExpPath, fileTag)

    # Create a string for the slurm command
    if saliency == 1:
        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
    elif saliency == 2:
        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

    sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    if not sp.stdout.startswith("Submitted batch"):
        raise ChildProcessError("SlurmError: sbatch not submitted correctly")

    
    combinationJobID = int(sp.stdout.split()[-1])
    print("    JobID:", combinationJobID)


    # Calculate the observed frequencies and scores
    print("\nSubmitting Slurm Jobs for Score Calculation....")
    scoreRealJobIDArr = []
    scoreNullJobIDArr = []
    for file1 in inputDirPath1.glob("*"):
        # Skip over ".genome" files
        if file1.name.split(".")[-1] == "genome":
                continue
        # Find matching file in other directory
        if not list(inputDirPath2.glob(file1.name)):
            raise FileNotFoundError("File not found: {} Please ensure corresponding files within input directories directories 1 and 2 have the same name".format(str(inputDirPath2 / file1.name)))
        else:
            file2 = next(inputDirPath2.glob(file1.name))

        if not file1.is_dir() and not file2.is_dir():
            filename = file1.name.split(".")[0]
            jobNameReal = "score_real_{}_{}".format(fileTag, filename)
            jobNameNull = "score_null_{}_{}".format(fileTag, filename)
            jobOutPathReal = outputDirPath / (".out/" + jobNameReal + ".out")
            jobErrPathReal = outputDirPath / (".err/" + jobNameReal + ".err")
            jobOutPathNull = outputDirPath / (".out/" + jobNameNull + ".out")
            jobErrPathNull = outputDirPath / (".err/" + jobNameNull + ".err")

            # Creating the out and err files for the batch job
            if jobOutPathReal.exists():
                os.remove(jobOutPathReal)
            if jobErrPathReal.exists():
                os.remove(jobErrPathReal)
            if jobOutPathNull.exists():
                os.remove(jobOutPathNull)
            if jobErrPathNull.exists():
                os.remove(jobErrPathNull)
            try:
                jout = open(jobOutPathReal, 'x')
                jout.close()
                jerr = open(jobErrPathReal, 'x')
                jerr.close()
                jout = open(jobOutPathNull, 'x')
                jout.close()
                jerr = open(jobErrPathNull, 'x')
                jerr.close()
            except FileExistsError as err:
                # This error should never occur because we are deleting the files first
                print(err)
                return
            
            # Create a string for the python commands
            computeScorePy = pythonFilesDir / "computeEpilogosPairwiseScores.py"
            pythonCommandReal = "python {} {} {} {} {} {} {} {} {} real".format(computeScorePy, file1, file2, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses)
            pythonCommandNull = "python {} {} {} {} {} {} {} {} {} null".format(computeScorePy, file1, file2, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses)

            # Create a string for the slurm commands
            if saliency == 1:
                slurmCommandReal = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameReal, jobOutPathReal, jobErrPathReal, numTasks, pythonCommandReal)
                slurmCommandNull = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameNull, jobOutPathNull, jobErrPathNull, numTasks, pythonCommandNull)
            elif saliency == 2:
                slurmCommandReal = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameReal, jobOutPathReal, jobErrPathReal, numTasks, pythonCommandReal)
                slurmCommandNull = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameNull, jobOutPathNull, jobErrPathNull, numTasks, pythonCommandNull)

            spReal = subprocess.run(slurmCommandReal, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
            spNull = subprocess.run(slurmCommandNull, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

            if not spReal.stdout.startswith("Submitted batch") or not spNull.stdout.startswith("Submitted batch"):
                raise ChildProcessError("SlurmError: sbatch not submitted correctly")


            scoreRealJobIDArr.append(int(spReal.stdout.split()[-1]))
            scoreNullJobIDArr.append(int(spNull.stdout.split()[-1]))

    # create a string for slurm dependency to work
    scoreRealJobIDStr = str(scoreRealJobIDArr).strip('[]').replace(" ", "")
    scoreNullJobIDStr = str(scoreNullJobIDArr).strip('[]').replace(" ", "")
    
    print("    Data JobIDs:", scoreRealJobIDStr)
    print("    Null JobIDs:", scoreNullJobIDStr)

    # Fitting, calculating p-values, and visualizing pairiwse differences
    print("\nSubmitting Slurm Jobs for data visualization....")

    jobName = "visual_{}".format(fileTag)
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

    # Create a string for the python commmand
    computeVisualPy = pythonFilesDir / "computeEpilogosPairwiseVisual.py"
    pythonCommand = "python {} {} {} {} {} {} {}".format(computeVisualPy, inputDirPath1.name, inputDirPath2.name, numStates, outputDirPath, fileTag, numProcesses)

    # Create a string for the slurm command
    if saliency == 1:
        slurmCommand = "sbatch --dependency=afterok:{},{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(scoreRealJobIDStr, scoreNullJobIDStr, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
    elif saliency == 2:
        slurmCommand = "sbatch --dependency=afterok:{},{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(scoreRealJobIDStr, scoreNullJobIDStr, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

    sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    if not sp.stdout.startswith("Submitted batch"):
        raise ChildProcessError("SlurmError: sbatch not submitted correctly")

    
    visualJobID = int(sp.stdout.split()[-1])
    print("    JobID:", visualJobID)

    allJobIDs = "{},{},{},{},{}".format(expJobIDStr, combinationJobID, scoreRealJobIDStr, scoreNullJobIDStr, visualJobID)
    print("\nAll JobIDs: {}".format(allJobIDs))

    # If the user wants to exit upon job completion rather than submission
    # If a job fails, it cancels all other jobs
    if exitBool:
        lastJobCheckStr = "sacct --format=State --jobs {}".format(allJobIDs)

        # Every ten seconds check if the final job is done, if it is exit the program
        while True:
            time.sleep(10)
            sp = subprocess.run(lastJobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
            if not ("RUNNING" in sp.stdout or "PENDING" in sp.stdout):
                break
            if "FAILED" in sp.stdout or "CANCELLED" in sp.stdout:
                print("\nERROR RUNNING JOBS: CANCELLING ALL REMAINING JOBS\n")
                print("Please check error logs in: {}/.err/\n".format(outputDirPath))
                subprocess.run("scancel {}".format(allJobIDs), shell=True)
                break

if __name__ == "__main__":
    main()