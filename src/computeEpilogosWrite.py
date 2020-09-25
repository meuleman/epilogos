import os
import numpy as np
from pathlib import Path
import gzip
import sys
import time

def main(fileTag, outputDirectory, numStates):
    outputDirPath = Path(outputDirectory)

    writeScores(fileTag, outputDirPath, numStates)
        
    # Clean up
    for file in outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag)):
        os.remove(file)

# Helper to write the final scores to files
def writeScores(fileTag, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observations_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scores_{}.txt.gz".format(fileTag)

    observationsTxt = gzip.open(observationsTxtPath, "wt")
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    tWrite = time.time()
    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
        combinedArr = np.load(file, allow_pickle=False)

        scoreArr = combinedArr[:, 3:].astype(float)
        locationArr = combinedArr[:, 0:3]

        # Write each row in both observations and scores
        for i in range(scoreArr.shape[0]):
            # Write in the coordinates
            for location in locationArr[i]:
                observationsTxt.write("{}\t".format(location))
                scoresTxt.write("{}\t".format(location))
            
            # Write to observations
            maxContribution = np.amax(scoreArr[i])
            maxContributionLoc = np.argmax(scoreArr[i]) + 1
            totalScore = np.sum(scoreArr[i])

            observationsTxt.write("{}\t{:.5f}\t1\t{:.5f}\t\n".format(maxContributionLoc, maxContribution, totalScore))

            # Write to scores
            for j in range(len(scoreArr[i])):
                scoresTxt.write("{:.5f}\t".format(scoreArr[i, j]))
            scoresTxt.write("\n")

    print("write time:", time.time() - tWrite)

    observationsTxt.close()
    scoresTxt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])