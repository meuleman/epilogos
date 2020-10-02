import os
import numpy as np
from pathlib import Path
import gzip
import sys
import time

def main(fileTag, outputDirectory, numStates, file):
    tTotal = time.time()
    outputDirPath = Path(outputDirectory)
    filePath = Path(file)

    writeScores(fileTag, outputDirPath, int(numStates), filePath)
    
    # Clean Up
    os.remove(filePath)

    print("Total Time:", time.time() - tTotal)

# Helper to write the final scores to files
def writeScores(fileTag, outputDirPath, numStates, filePath):
    locationTag = '_'.join(filePath.name.split(".")[-2].split("_")[-3:])

    # observationsTxtPath = outputDirPath / "observations_{}_{}.txt.gz".format(fileTag, locationTag)
    scoresTxtPath = outputDirPath / "scores_{}_{}.txt.gz".format(fileTag, locationTag)

    # observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    tCalc = time.time()
    # Taking in the the score array
    combinedArr = np.load(filePath, allow_pickle=False)
    scoreArr = combinedArr[:, 3:].astype(float)
    locationArr = combinedArr[:, 0:3]

    # # Calculating observation values
    # maxContributions = np.amax(scoreArr, axis=1).reshape((scoreArr.shape[0], 1))
    # maxContributionsLocs = np.argmax(scoreArr, axis=1).reshape((scoreArr.shape[0], 1)) + 1
    # totalScores = np.sum(scoreArr, axis=1).reshape((scoreArr.shape[0], 1))

    # # Splicing all the observation arrays together
    # observationArr = np.concatenate((maxContributionsLocs, maxContributions, totalScores), axis=1)

    print("Calc Time:", time.time() - tCalc)

    t1 = time.time()
    # observationStr = "".join("{0[0]}\t{0[1]}\t{0[2]}\t{1:d}\t{2[1]:.5f}\t1\t{2[2]:.5f}\t\n".format(locationArr[i], int(observationArr[i, 0]), observationArr[i]) for i in range(scoreArr.shape[0]))
    scoresTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates)) + "\n"
    scoreStr = "".join(scoresTemplate.format(locationArr[i], scoreArr[i]) for i in range(scoreArr.shape[0]))
    print("String creation time:", time.time() - t1)

    # tObs = time.time()
    # observationsTxt.write(observationStr)
    # print("Observation Write Time:", time.time() - tObs)

    tScore = time.time()
    scoresTxt.write(scoreStr)
    print("Score Write Time:", time.time() - tScore)
    
    # observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])