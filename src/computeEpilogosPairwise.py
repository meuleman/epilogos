import numpy as np
import time
import pandas as pd
from pathlib import Path
import sys
import gzip

def main(file1, file2, outputDirectory):
    tTotal = time.time()
    filepath1 = Path(file1)
    filepath2 = Path(file2)
    outputDirPath = Path(outputDirectory)

    # Read in the data
    print("\nReading data from files...")
    tRead = time.time()
    dataDF1 = pd.read_table(filepath1, header=None, sep="\t")
    dataDF2 = pd.read_table(filepath2, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)
    
    if dataDF1.shape != dataDF2.shape:
        print("ERROR: Epilogos score files do not have the same shape")
        return

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    dataArr1 = dataDF1.iloc[:,3:].to_numpy(dtype=int)
    dataArr2 = dataDF2.iloc[:,3:].to_numpy(dtype=int)
    locationArr = dataDF1.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    scoreArr = calculatePairwise(dataArr1, dataArr2)

    writeScores(scoreArr, locationArr, outputDirPath, dataArr1.shape[1])

    print("Total Time: ", time.time() - tTotal)


def calculatePairwise(dataArr1, dataArr2):
    numRows, numCols = dataArr1.shape

    # scoreArr = np.zeros(numRows, numCols)

    # for row in range(numRows):
    #     scoreArr[row] = dataArr1[row] - dataArr2[row]

    scoreArr = dataArr1 - dataArr2

    return scoreArr






# Helper to write the final scores to files
def writeScores(locationArr, scoreArr, outputDirPath, numStates):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    observationsTxtPath = outputDirPath / "observationsP.txt.gz"
    scoresTxtPath = outputDirPath / "scoresP.txt.gz"

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
    if len(sys.argv) - 1 < 3:
        print("\n You provide 2 arguments:\n")
        print("    1. Filepath for the first epilogos score file")
        print("    2. Filepath for the second epilogos score file")
        print("    3. Output directory")
    elif len(sys.argv) - 1 == 3:
        main(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        print("ERROR: incorrect number of arguments submitted")