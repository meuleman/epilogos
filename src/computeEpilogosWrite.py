import os
import numpy as np
from pathlib import Path
import gzip
import sys

def main(fileTag, outputDirectory, numStates):
    outputDirPath = Path(outputDirectory)

    # Order matters to us when writing, so use sorted
    for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
        combinedArr = np.load(file, allow_pickle=False)
        writeScores(fileTag, combinedArr, outputDirPath, numStates)
        
    # Clean up
    for file in outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag)):
        os.remove(file)

# Helper to write the final scores to files
def writeScores(fileTag, combinedArr, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observations_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scores_{}.txt.gz".format(fileTag)

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    scoreArr = combinedArr[:, 3:]
    locationArr = combinedArr[:, 0:3]

    # Write each row in both observations and scores
    for i in range(scoreArr.shape[0]):
        # Write in the coordinates
        for location in locationArr:
            observationsTxt.write("{}\t".format(location))
            scoresTxt.write("{}\t".format(location))
        
        # Write to observations
        maxContribution = np.amax(scoreArr[i].astype(float))
        maxContributionLoc = np.argmax(scoreArr[i].astype(float)) + 1
        totalScore = np.sum(scoreArr[i].astype(float))

        observationsTxt.write("{}\t{:.5f}\t1\t{:.5f}\t\n".format(maxContributionLoc, maxContribution, totalScore))

        # Write to scores
        for j in range(len(scoreArr[i])):
            scoresTxt.write("{:.5f}\t".format(float(scoreArr[i, j])))
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])