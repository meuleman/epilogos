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
    observationsTxtPath = outputDirPath / "observationsSAVETXT_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scoresSAVETXT_{}.txt.gz".format(fileTag)

    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    tLoop = time.time()
    first = True
    for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
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
            scoreArr = np.concatenate((scoreArr, fileScoreArr), axis=0)
            locationArr = np.concatenate((locationArr, fileLocationArr), axis=0)
            observationArr = np.concatenate((observationArr, fileObservationArr), axis=0)

    scoreArr = np.around(scoreArr, decimals=5)
    observationArr = np.around(observationArr, decimals=5)

    print("Observation Calculation Time:", time.time() - tLoop)

    tConcat = time.time()
    scoreConcatArr = np.concatenate((locationArr, scoreArr), axis=1)
    obsConcatArr = np.concatenate((locationArr, observationArr), axis=1)

    print("Concatenation Time:", time.time() - tConcat)


    tScore = time.time()
    np.savetxt(scoresTxtPath, scoreConcatArr, fmt="%s", delimiter="\t")
    print("Score SaveTxt Time:", time.time() - tScore)

    tObs = time.time()
    observationFMT = "%s\t%s\t%s\t%s\t%s\t1\t%s"
    np.savetxt(observationsTxtPath, obsConcatArr, fmt=observationFMT)
    print("Observation SaveTxt:", time.time() - tScore)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])