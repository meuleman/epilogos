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
@click.option("-i", "--directory-one", "fileDirectory1", type=str, required=True, multiple=True, help="Path to directory that contains files to read from (All files in this directory will be read in)")
@click.option("-n", "--directory-two", "fileDirectory2", type=str, required=True, multiple=True, help="")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, multiple=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-s", "--state-model", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-l", "--saliency-level", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1 or 2)")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-x", "--exit-when-complete", "exitBool", is_flag=True, multiple=True, help="If flag is enabled program will exit upon completion of all processes rather than SLURM submission of all processes")
def main(fileDirectory1, fileDirectory2, outputDirectory, numStates, saliency, numProcesses, exitBool):
    """
    This script computes and visualizes differences between chromatin states across epigenomes
    """
    # Handling case if user inputs flag multiples times
    if len(fileDirectory1) > 1:
        raise ValueError("Too many [-i, --directory-one] arguments provided")
    elif len(fileDirectory2) > 1:
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
    fileDirectory1 = fileDirectory1[0]
    fileDirectory2 = fileDirectory2[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    numProcesses = numProcesses[0]
    if exitBool:
        exitBool = exitBool[0]

    file1Path = Path(fileDirectory1)
    file2Path = Path(fileDirectory2)
    outputDirPath = Path(outputDirectory)

    print()
    print("File Directory 1 =", file1Path)
    print("File Directory 2 =", file2Path)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)

    if saliency != 1 and saliency != 2 and saliency != 3:
        print("\nERROR: Ensure that saliency metric is either 1 or 2 (Saliency of 3 is unsupported for pairwise comparison\n")
        return

    # Making paths absolute
    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath
    if not PurePath(file1Path).is_absolute():
        file1Path = Path.cwd() / file1Path
    if not PurePath(file2Path).is_absolute():
        file2Path = Path.cwd() / file2Path

    # For making sure all files are consistently named
    fileTag = "{}_{}_saliency{}".format(file1Path.name, file2Path.name, saliency)

    # Check that paths are valid before doing anything
    if not file1Path.exists() or not file1Path.is_dir():
        print("\nERROR: Given file path does not exist or is not a directory\n")
        return
    if not list(file1Path.glob("*")):
        print("\nERROR: Ensure that file directory 1 is not empty\n")
        return
    if not file2Path.exists() or not file2Path.is_dir():
        print("\nERROR: Given file path does not exist or is not a directory\n")
        return
    if not list(file2Path.glob("*")):
        print("\nERROR: Ensure that file directory 2 is not empty\n")
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
    storedExpPath = outputDirPath / "exp_freq_{}_{}.npy".format(file1Path.name, file2Path.name)
    print("\nBackground Frequency Array Location:", storedExpPath)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    expJobIDArr = []   
    print("\nSubmitting Slurm Jobs for Per Datafile Background Frequency Calculation....")
    for file1 in file1Path.glob("*"):
        # Skip over ".genome" files
        if file1.name.split(".")[1] == "genome":
                continue
        # Find matching file in other directory
        if not list(file2Path.glob(file1.name)):
            print("\nERROR: File names do not match in input directory and pairwise directory")
            print("\tNo match for {} in pairwise directory. Ensure corresponding files within directories 1 and 2 have the same name\n".format(file1.name))
            return
        else:
            file2 = next(file2Path.glob(file1.name))

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
            except FileExistsError:
                # This error should never occur because we are deleting the files first
                print("\nERROR: sbatch '.out' or '.err' file already exists\n")
                return

            # Create a string for the python command
            computeExpectedPy = pythonFilesDir / "computeEpilogosPairwiseExpected.py"
            pythonCommand = "python {} {} {} {} {} {} {} {}".format(computeExpectedPy, file1, file2, numStates, saliency, outputDirPath, fileTag, numProcesses)

            # Create a string for the slurm command
            if saliency == 1:
                slurmCommand = "sbatch --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
            elif saliency == 2:
                slurmCommand = "sbatch --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

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

    # Create a string for the python commmand
    computeExpectedCombinationPy = pythonFilesDir / "computeEpilogosExpectedCombination.py"
    pythonCommand = "python {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, storedExpPath, fileTag)

    # Create a string for the slurm command
    if saliency == 1:
        slurmCommand = "sbatch --dependency=afterok:{} --job-name=S1_{}.job --output={} --partition=queue1 --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
    elif saliency == 2:
        slurmCommand = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --partition=queue1 --error={} --nodes=1 --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

    sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    if not sp.stdout.startswith("Submitted batch"):
        print("\nERROR: sbatch not submitted correctly\n")
        return
    
    combinationJobID = int(sp.stdout.split()[-1])
    print("    JobID:", combinationJobID)


    # Calculate the observed frequencies and scores
    print("\nSubmitting Slurm Jobs for Score Calculation....")
    scoreRealJobIDArr = []
    scoreNullJobIDArr = []
    for file1 in file1Path.glob("*"):
        # Skip over ".genome" files
        if file1.name.split(".")[1] == "genome":
                continue
        # Find matching file in other directory
        if not list(file2Path.glob(file1.name)):
            print("\nERROR: File names do not match in input directory and pairwise directory")
            print("\tNo match for {} in pairwise directory. Ensure corresponding files within directories 1 and 2 have the same name\n".format(file1.name))
            return
        else:
            file2 = next(file2Path.glob(file1.name))

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
            except FileExistsError:
                # This error should never occur because we are deleting the files first
                print("\nERROR: sbatch '.out' or '.err' file already exists\n")
                return
            
            # Create a string for the python commands
            computeScorePy = pythonFilesDir / "computeEpilogosPairwiseScores.py"
            pythonCommandReal = "python {} {} {} {} {} {} {} {} {} real".format(computeScorePy, file1, file2, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses)
            pythonCommandNull = "python {} {} {} {} {} {} {} {} {} null".format(computeScorePy, file1, file2, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses)

            # Create a string for the slurm commands
            if saliency == 1:
                slurmCommandReal = "sbatch --dependency=afterok:{} --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameReal, jobOutPathReal, jobErrPathReal, numTasks, pythonCommandReal)
                slurmCommandNull = "sbatch --dependency=afterok:{} --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameNull, jobOutPathNull, jobErrPathNull, numTasks, pythonCommandNull)
            elif saliency == 2:
                slurmCommandReal = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameReal, jobOutPathReal, jobErrPathReal, numTasks, pythonCommandReal)
                slurmCommandNull = "sbatch --dependency=afterok:{} --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobNameNull, jobOutPathNull, jobErrPathNull, numTasks, pythonCommandNull)

            spReal = subprocess.run(slurmCommandReal, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
            spNull = subprocess.run(slurmCommandNull, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

            if not spReal.stdout.startswith("Submitted batch") or not spNull.stdout.startswith("Submitted batch"):
                print("\nERROR: sbatch not submitted correctly\n")
                return

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
    except FileExistsError:
        # This error should never occur because we are deleting the files first
        print("\nERROR: sbatch '.out' or '.err' file already exists\n")
        return

    # Create a string for the python commmand
    computeVisualPy = pythonFilesDir / "computeEpilogosPairwiseVisual.py"
    pythonCommand = "python {} {} {} {} {} {}".format(computeVisualPy, file1Path.name, file2Path.name, numStates, outputDirPath, numProcesses)

    # Create a string for the slurm command
    if saliency == 1:
        slurmCommand = "sbatch --dependency=afterok:{},{} --job-name=S1_{}.job --output={} --partition=queue1 --error={} {} --wrap='{}'".format(scoreRealJobIDStr, scoreNullJobIDStr, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
    elif saliency == 2:
        slurmCommand = "sbatch --dependency=afterok:{},{} --job-name=S2_{}.job --output={} --partition=queue1 --error={} {} --wrap='{}'".format(scoreRealJobIDStr, scoreNullJobIDStr, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

    sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    if not sp.stdout.startswith("Submitted batch"):
        print("\nERROR: sbatch not submitted correctly\n")
        return
    
    visualJobID = int(sp.stdout.split()[-1])
    print("    JobID:", visualJobID)

    print("\nAll JobIDs: {},{},{},{},{}".format(expJobIDStr, combinationJobID, scoreRealJobIDStr, scoreNullJobIDStr, visualJobID))

    # If the user wants to exit upon job completion rather than submission
    if exitBool:
        lastJobCheckStr = "sacct --format=State --jobs {}".format(visualJobID)
        # Every ten seconds check if the final job is done, if it is exit the program
        while True:
            time.sleep(10)
            sp = subprocess.run(lastJobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
            if not (sp.stdout.strip().endswith("RUNNING") or sp.stdout.strip().endswith("PENDING")):
                break

if __name__ == "__main__":
    main()