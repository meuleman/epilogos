import os
import numpy as np
from pathlib import Path
import gzip
import sys

def main(fileTag, outputDirectory, numStates):
    outputDirPath = Path(outputDirectory)

    writeScores(fileTag, outputDirPath, numStates)
        
    # # Clean up
    # for file in outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag)):
    #     os.remove(file)

# Helper to write the final scores to files
def writeScores(fileTag, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observations_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scores_{}.txt.gz".format(fileTag)

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    first = True
    # for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
    for file in sorted(outputDirPath.glob("temp_scores*.npy")):
        combinedArr = np.load(file, allow_pickle=False)

        fileScoreArr = combinedArr[:, 3:].astype(float)
        fileLocationArr = combinedArr[:, 0:3]

        # Calculating observation values
        maxContributions = np.amax(fileScoreArr, axis=1).reshape((fileScoreArr.shape[0], 1))
        maxContributionsLocs = np.argmax(fileScoreArr, axis=1).reshape((fileScoreArr.shape[0], 1)) + 1
        totalScores = np.sum(fileScoreArr, axis=1).reshape((fileScoreArr.shape[0], 1))

        # Splicing all the observation arrays together
        fileObservationArr = np.concatenate((maxContributionsLocs, maxContributions, totalScores), axis=1)

        # Storing the per file arrays into the entire array
        if first:
            scoreArr = fileScoreArr
            locationArr = fileLocationArr
            observationArr = fileObservationArr
            first = False
        else:
            np.concatenate((scoreArr, fileScoreArr), axis=0)
            np.concatenate((locationArr, fileLocationArr), axis=0)
            np.concatenate((observationArr, fileObservationArr), axis=0)

    # Write each row in both observations and scores
    for i in range(scoreArr.shape[0]):
        # Write in the coordinates
        observationsTxt.write("{}\t{}\t{}\t".format(locationArr[i, 0], locationArr[i, 1], locationArr[i, 2]))
        scoresTxt.write("{}\t{}\t{}\t".format(locationArr[i, 0], locationArr[i, 1], locationArr[i, 2]))
        
        observationsTxt.write("{:d}\t{:.5f}\t1\t{:.5f}\t\n".format(int(observationArr[i, 0]), observationArr[i, 1], observationArr[i, 2]))

        # Write to scores
        for j in range(len(scoreArr[i])):
            scoresTxt.write("{:.5f}\t".format(scoreArr[i, j]))
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])