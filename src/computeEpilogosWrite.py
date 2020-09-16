import os
import numpy as np
from pathlib import Path
import gzip
import sys

def main(fileTag, outputDirectory, numStates):
    outputDirPath = Path(outputDirectory)

    # Clean up the saved temporary expected frequency array
    for file in outputDirPath.glob("temp_exp_freq_{}_*.npy".format(fileTag)):
        os.remove(file)

    # Order matters to us when writing, so use sorted
    for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
        scoreArr = np.load(file, allow_pickle=False)
        writeScores(fileTag, scoreArr, outputDirPath, numStates)
        # Clean up
        os.remove(file)

# Helper to write the final scores to files
def writeScores(fileTag, scoreArr, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observations_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scores_{}.txt.gz".format(fileTag)

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

        observationsTxt.write("{}\t{0:.5f}\t1\t{0:.5f}\t\n".format(maxContributionLoc, maxContribution, totalScore))

        # Write to scores
        for j in range(scoreArr.shape[1]):
            scoresTxt.write("{0:.5f}\t".format(scoreArr[i, j]))
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])