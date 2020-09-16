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
    dataArr1 = dataDF1.iloc[:,3:].to_numpy(dtype=float)
    dataArr2 = dataDF2.iloc[:,3:].to_numpy(dtype=float)
    locationArr = dataDF1.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    # Calculate the difference between the two score files
    tScores = time.time()
    scoreArr = dataArr1 - dataArr2
    print("Score Time: ", time.time() - tScores)

    writeScores(locationArr, scoreArr, outputDirPath)

    print("Total Time: ", time.time() - tTotal)


# Helper to write the final scores to files
def writeScores(locationArr, scoreArr, outputDirPath):
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    observationsTxtPath = outputDirPath / "observationsPairwise.txt.gz"
    scoresTxtPath = outputDirPath / "scoresPairwise.txt.gz"

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    # Write each row in both observations and scores
    for i in range(locationArr.shape[0]):
        # Write in the coordinates
        for location in locationArr[i]:
            observationsTxt.write("{}\t".format(location))
            scoresTxt.write("{}\t".format(location))
        
        # Write to observations
        maxDifference = np.amax(np.absolute(scoreArr[i]))
        maxDifferenceLoc = np.argmax(np.absolute(scoreArr[i])) + 1
        differenceSign = np.sign(scoreArr[i, maxDifferenceLoc - 1])
        totalAbsoluteDifference = np.sum(np.absolute(scoreArr[i]))

        observationsTxt.write("{}\t{0:.5f}\t{}\t{0:.5f}\t\n".format(maxDifferenceLoc, maxDifference, differenceSign, totalAbsoluteDifference))

        # Write to scores
        for j in range(scoreArr.shape[1]):
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