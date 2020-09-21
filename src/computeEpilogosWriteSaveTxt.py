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
    observationsTxtPath = outputDirPath / "observationsSAVETXT_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scoresSAVETXT_{}.txt.gz".format(fileTag)

    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    first = True
    # for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
    for file in sorted(outputDirPath.glob("temp_scores*.npy")):
        combinedArr = np.load(file, allow_pickle=False)

        fileScoreArr = combinedArr[:, 3:]
        fileLocationArr = combinedArr[:, 0:3]

        # Calculating observation values
        maxContributions = np.amax(fileScoreArr.astype(float), axis=1).reshape((fileScoreArr.shape[0], 1))
        maxContributionsLocs = np.argmax(fileScoreArr.astype(float), axis=1).reshape((fileScoreArr.shape[0], 1)) + 1
        totalScores = np.sum(fileScoreArr.astype(float), axis=1).reshape((fileScoreArr.shape[0], 1))

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

    scoreFMT = "%s\t%s\t%s"
    for i in numStates:
        scoreFMT += "\t%.5f"
    np.savetxt(scoresTxtPath, np.concatenate((locationArr, scoreArr), axis=1), fmt=scoreFMT)

    observationFMT = "%s\t%s\t%s\t%d\t%.5f\t1\t%.5f"
    np.savetxt(observationsTxtPath, np.concatenate((locationArr, observationArr), axis=1), fmt=observationFMT)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])