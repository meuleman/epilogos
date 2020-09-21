import os
import numpy as np
from pathlib import Path
import gzip
import sys
import shutil

def main(fileTag, outputDirectory, numStates):
    outputDirPath = Path(outputDirectory)

    writeScores(fileTag, outputDirPath, numStates)
        
    # # Clean up
    # for file in outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag)):
    #     os.remove(file)

# Helper to write the final scores to files
def writeScores(fileTag, outputDirPath, numStates):
    observationsTxtPath = outputDirPath / "observationsGZIPLAST_{}.txt".format(fileTag)
    scoresTxtPath = outputDirPath / "scoresGZIPLAST_{}.txt".format(fileTag)

    gzipObservationsTxtPath = outputDirPath / "observationsGZIPLAST_{}.txt.gz".format(fileTag)
    gzipScoresTxtPath = outputDirPath / "scoresGZIPLAST_{}.txt.gz".format(fileTag)

    observationsTxt = open(observationsTxtPath, "w")
    scoresTxt = open(scoresTxtPath, "w")

    gzipObservationsTxt = gzip.open(gzipObservationsTxtPath, "w")
    gzipScoresTxt = gzip.open(gzipScoresTxtPath, "w")

    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    # for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
    for file in sorted(outputDirPath.glob("temp_scores*.npy")):
        combinedArr = np.load(file, allow_pickle=False)

        scoreArr = combinedArr[:, 3:]
        locationArr = combinedArr[:, 0:3]

        # Write each row in both observations and scores
        for i in range(scoreArr.shape[0]):
            # Write in the coordinates
            for location in locationArr[i]:
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

    shutil.copyfileobj(observationsTxt, gzipObservationsTxt)
    shutil.copyfileobj(scoresTxt, gzipScoresTxt)

    observationsTxt.close()
    scoresTxt.close()
    gzipObservationsTxt.close()
    gzipScoresTxt.close()

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])