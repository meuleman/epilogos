import os
import numpy as np
from pathlib import Path
import gzip
import sys

def main(fileDirectory, outputDirectory, numStates):
    outputDirPath = Path(outputDirectory)
    filePath = Path(fileDirectory)

    # formatting the name for the output files
    filename = "_".join(str(filePath).split("/")[-5:])

    # Clean up the saved temporary expected frequency array
    for file in outputDirPath.glob("temp_exp_freq_*.npy"):
        os.remove(file)

    # Order matters to us when writing, so use sorted
    for file in sorted(outputDirPath.glob("temp_scores_*.npy")):
        scoreArr = np.load(file, allow_pickle=False)
        writeScores(filename, scoreArr, outputDirPath, numStates)
        # Clean up
        os.remove(file)

# Helper to write the final scores to files
def writeScores(filename, scoreArr, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observations_{}.txt.gz".format(filename)
    scoresTxtPath = outputDirPath / "scores_{}.txt.gz".format(filename)

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
    main(sys.argv[1], sys.argv[2], sys.argv[3])