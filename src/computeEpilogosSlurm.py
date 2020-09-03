import gzip
import numpy as np
from pathlib import Path
import pandas as pd
import time
import click
import os
import subprocess

@click.command()
@click.option("-f", "--file-directory", "fileDirectory", type=str, required=True, help="Path to directory that contains files to read from (Please ensure that this directory contains only files you want to read from)")
@click.option("-m", "--state-model", "numStates", type=int, required=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency-level", "saliency", type=int, required=True, help="Saliency level (1, 2, or 3)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, help="Output Directory")
@click.option("-e", "--store-expected", "storeExp", is_flag=True, help="[Flag] Store the expected frequency array for later calculations (Must be used in conjunction with '-d' and cannot be used in conjunction with '-u')")
@click.option("-u", "--use-expected", "useStoredExp", is_flag=True, help="[Flag] Use previously stored expected frequency array (Must be used in conjunction with '-d' and cannot be used in conjunction with '-e')")
@click.option("-d", "--expected-directory", "expFreqDir", type=str, default="null", help="Path to the stored expected frequency array (Must be used in conjunction with either '-e' or '-u')")
def main(fileDirectory, numStates, saliency, outputDirectory, storeExp, useStoredExp, expFreqDir):
    tTotal = time.time()
    dataFilePath = Path(fileDirectory)
    outputDirPath = Path(outputDirectory)

    # Check that paths are valid before doing anything
    if not dataFilePath.exists() or not dataFilePath.is_dir():
        print("ERROR: Given file path does not exist or is not a directory")
        return

    if list(dataFilePath.glob("*")):
        print("ERROR: Ensure that file directory is not empty")

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    if not outputDirPath.is_dir():
        print("ERROR: Output directory is not a directory")
        return

    # Calculate the number of epigenomes (just read one line of one of the data files)
    for x in dataFilePath.glob("*"):
        numEpigenomes = pd.read_table(x, header=None, sep="\t", nrows=1).shape[1] - 3
        break

    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to the number of epigenomes, number of states, and saliency metric involved in the calculation
    if expFreqDir != "null":
        expFreqFilename = "exp_freq_{}_{}_s{}.npy".format(numEpigenomes, numStates, saliency)
        storedExpPath = Path(expFreqDir) / expFreqFilename

    # Finding the location of the .py files that must be run
    if Path(__file__).is_absolute:
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    # Calculate the expected frequencies
    tExp = time.time()
    print("Calculating expected frequencies...")

    # Check if its already been calculated before
    if useStoredExp:
        try:
            expFreqArr = np.load(storedExpPath, allow_pickle=False)
        except IOError:
            print("ERROR: Could not load stored expected value array.\n\tPlease check that the directory is correct or that the file exits")
    else:        
        print(numEpigenomes)
        for file in dataFilePath.glob("*"):
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "exp_freq_calc_{}".format(filename)
                file = dataFilePath / file
                computeExpectedPy = pythonFilesDir / "computeEpilogosExpected.py"

                pythonCommand = "python {} {} {} {} {}".format(computeExpectedPy, file, numStates, saliency, outputDirPath)
                slurmCommand = "sbatch --job-name={0}.job --output=.out/{0}.out --error=.out/{0}.err --nodes=1 --ntasks=1 --wrap='{1}'".format(jobName, pythonCommand)

                subprocess.call(slurmCommand, shell=True)

        # Loop over all the expected value arrays and add them up (normalize for number of chromosomes)
        count = 0
        expFreqArr = np.zeros((1,1))
        for file in outputDirPath.glob("temp_exp_freq_*.npy"):
            if count == 0:
                expFreqArr = np.load(file, allow_pickle=False)
            else:
                expFreqArr += np.load(file, allow_pickle=False)
            count += 1
            # Delete file after we're done with it
            os.remove(file)
        expFreqArr /= count

        # If user desires, store away the expected frequency array
        if storeExp:
            np.save(storedExpPath, expFreqArr, allow_pickle=False)

        # Save somewhere regardless for use in score calculation
        print(numEpigenomes)
        expFreqFilename = "temp_exp_freq_{}_{}_s{}.npy".format(numEpigenomes, numStates, saliency)
        expFreqPath = outputDirPath / expFreqFilename
        np.save(expFreqPath, expFreqArr, allow_pickle=False)

    print("    Time: ", time.time() - tExp)

    # Calculate the observed frequencies and scores
    for file in dataFilePath.glob("*"):
        if not file.is_dir():
            filename = file.name.split(".")[0]
            jobName = "score_calc_{}".format(filename)
            computeScorePy = pythonFilesDir / "computeEpilogosScores.py"

            pythonCommand = "python {} {} {} {} {} {}".format(computeScorePy, file, numStates, saliency, outputDirPath, expFreqPath)
            slurmCommand = "sbatch --job-name={0}.job --output=.out/{0}.out --error=.out/{0}.err --nodes=1 --ntasks=1 --wrap='{1}'".format(jobName, pythonCommand)

            subprocess.call(slurmCommand, shell=True)

    # Clean up the saved temporary expected frequency array
    for file in outputDirPath.glob("temp_exp_freq_*.npy"):
        os.remove(file)

    # Writing the scores to the files
    print("Writing to files...")
    tWrite = time.time()

    # Order matters to us when writing, so use sorted
    for file in sorted(outputDirPath.glob("temp_scores_*.npy")):
        scoreArr = np.load(file, allow_pickle=False)
        writeScores(scoreArr, outputDirPath, numStates)
        # Clean up
        os.remove(file)

    print("    Time: ", time.time() - tWrite)

    print("Total Time: ", time.time() - tTotal)

# Helper to write the final scores to files
def writeScores(scoreArr, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observationsSLURM.txt.gz"
    scoresTxtPath = outputDirPath / "scoresSLURM.txt.gz"

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    # Write each row in both observations and scores
    for i in range(scoreArr.shape[0]):
        # Write in the coordinates
        for location in scoreArr[i, 0:3]:
            observationsTxt.write("{}\t".format(location))
            scoresTxt.write("{}\t".format(location))
        
        # Write to observations
        maxContribution = np.amax(int(scoreArr[i, 3:].astype(float)))
        maxContributionLoc = np.argmax(scoreArr[i, 3:].astype(float)) + 1
        totalScore = np.sum(scoreArr[i, 3:].astype(float))

        observationsTxt.write("{}\t".format(maxContributionLoc))
        observationsTxt.write("{0:.5f}\t".format(maxContribution))
        observationsTxt.write("1\t")
        observationsTxt.write("{0:.5f}\t\n".format(totalScore))

        # Write to scores
        for j in range(numStates):
            scoresTxt.write("{0:.5f}\t".format(scoreArr[i, j + 3]))
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main()