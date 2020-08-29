import gzip
import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time
import numpy.ma as ma
import operator as op
from functools import reduce
import multiprocessing
import ctypes
import itertools
import click
import glob
import os

@click.command()
@click.option("-f", "--file-directory", "fileDirectory", type=str, required=True, help="Path to directory that contains files to read from")
@click.option("-m", "--state-model", "numStates", type=int, required=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency-level", "saliency", type=int, required=True, help="Saliency level (1, 2, or 3)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, help="Output Directory")
@click.option("-e", "--store-expected", "storeExp", is_flag=True, help="Flag: Store the expected frequency array for later calculations (Must be used in conjunction with '-d' and cannot be used in conjunction with '-u')")
@click.option("-u", "--use-expected", "useStoredExp", is_flag=True, help="Flag: Use previously stored expected frequency array (Must be used in conjunction with '-d' and cannot be used in conjunction with '-e')")
@click.option("-d", "--expected-directory", "expFreqDir", type=str, default="null", help="Path to the stored expected frequency array (Must be used in conjunction with either '-e' or '-u')")
def main(fileDirectory, numStates, saliency, outputDirectory, storeExp, useStoredExp, expFreqDir):
    tTotal = time.time()
    dataFilePath = Path(fileDirectory)
    outputDirPath = Path(outputDirectory)

    # Check that paths are valid before doing anything
    if not dataFilePath.exists() or not dataFilePath.is_dir():
        print("ERROR: Given file path does not exist or is not a directory")
        return

    if not outputDirPath.is_dir():
        print("ERROR: Output directory is not a directory")
        return

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    # Calculate the number of epigenomes
    for x in dataFilePath.glob("*"):
        numEpigenomes = pd.read_table(x, header=None, sep="\t", nrows=1).shape[1] - 3
        break

    # Path for storing/retrieving the expected frequency array
    if expFreqDir != "null":
        expFreqFilename = "exp_freq_" + str(numEpigenomes) + "_" + str(numStates) + "_s" + str(saliency) + ".npy"
        storedExpPath = Path(expFreqDir) / expFreqFilename

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
        for file in dataFilePath.glob("*"):

            ##############################################################
            #                                                            #
            #   Insert bash stuff here to run computeEpilogosExpected    #
            #                                                            #
            ##############################################################

            # py computeEpilogosExpected.py file numStates saliency outputDirPath

            print("lakes")

        # Loop over all the expected value arrays and add them up
        firstFile = True
        for file in outputDirPath.glob("temp_exp_freq_*.npy"):
            if firstFile:
                expFreqArr = np.load(file, allow_pickle=False)
                firstFile = False
            else:
                expFreqArr += np.load(file, allow_pickle=False)
            # Delete file after were done with it
            os.remove(file)

        # If user desires, store away the expected frequency array
        if storeExp:
            np.save(storedExpPath, expFreqArr, allow_pickle=False)

        # Save somewhere regardless for use in score calculation
        expFreqPath = "temp_exp_freq_" + str(numEpigenomes) + "_" + str(numStates) + "_s" + str(saliency) + ".npy"
        np.save(expFreqPath, expFreqArr, allow_pickle=False)

    print("    Time: ", time.time() - tExp)

    # Calculate the observed frequencies and scores
    for file in dataFilePath.glob("*"):

        ##############################################################
        #                                                            #
        #    Insert bash stuff here to run computeEpilogosScores     #
        #                                                            #
        ##############################################################

        # py computeEpilogosScores.py file numStates saliency outputDirPath expFreqPath

        print("lakes")


    # Writing the scores to the files
    print("Writing to files...")
    tWrite = time.time()
    writeScores(locationArr, outputDirPath, numStates)
    print("    Time: ", time.time() - tWrite)

    print("Total Time: ", time.time() - tTotal)

# Helper to write the final scores to files
def writeScores(locationArr, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observationsM.txt.gz"
    scoresTxtPath = outputDirPath / "scoresM.txt.gz"

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    # Write each row in both observations and scores
    for i in range(locationArr.shape[0]):
        # Write in the coordinates
        for location in locationArr[i]:
            observationsTxt.write("{}\t".format(location))
            scoresTxt.write("{}\t".format(location))
        
        # Write to observations
        maxContribution = np.amax(scoreArr[i])
        maxContributionLoc = np.argmax(scoreArr[i]) + 1
        totalScore = np.sum(scoreArr[i])

        observationsTxt.write("{}\t".format(maxContributionLoc))
        observationsTxt.write("{0:.5f}\t".format(maxContribution))
        observationsTxt.write("1\t")
        observationsTxt.write("{0:.5f}\t\n".format(totalScore))

        # Write to scores
        for j in range(numStates):
            scoresTxt.write("{0:.5f}\t".format(scoreArr[i, j]))
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main()