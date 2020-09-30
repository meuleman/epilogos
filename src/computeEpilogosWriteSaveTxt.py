import os
import numpy as np
from pathlib import Path
import gzip
import sys
import time

def main(fileTag, outputDirectory, numStates):
    tTotal = time.time()
    outputDirPath = Path(outputDirectory)

    writeScores(fileTag, outputDirPath, int(numStates))
        
    tRemove = time.time()
    # Clean up
    for file in outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag)):
        os.remove(file)
    print("Remove file time:", time.time() - tRemove)

    print("Total Time:", time.time() - tTotal)

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
            scoreArr = np.around(fileScoreArr, decimals=5).astype(str)
            locationArr = fileLocationArr
            observationArr = np.concatenate((fileObservationArr[:,0].reshape(fileObservationArr.shape[0], 1).astype(int).astype(str), np.around(fileObservationArr[:,1], decimals=5).reshape(fileObservationArr.shape[0], 1), np.ones((fileObservationArr.shape[0], 1), dtype=int), np.around(fileObservationArr[:,2], decimals=5).reshape(fileObservationArr.shape[0], 1)), axis=1)
            first = False
        else:
            formattedScore = np.around(fileScoreArr, decimals=5).astype(str)
            scoreArr = np.concatenate((scoreArr, fileScoreArr), axis=0)
            locationArr = np.concatenate((locationArr, fileLocationArr), axis=0)
            formattedObs = np.concatenate((fileObservationArr[:,0].reshape(fileObservationArr.shape[0], 1).astype(int).astype(str), np.around(fileObservationArr[:,1], decimals=5).reshape(fileObservationArr.shape[0], 1), np.ones((fileObservationArr.shape[0], 1), dtype=int), np.around(fileObservationArr[:,2], decimals=5).reshape(fileObservationArr.shape[0], 1)), axis=1)
            observationArr = np.concatenate((observationArr, formattedObs), axis=0)

    print("Observation Calculation Time:", time.time() - tLoop)

    tConcat = time.time()
    scoreConcatArr = np.concatenate((locationArr, scoreArr), axis=1)
    obsConcatArr = np.concatenate((locationArr, observationArr), axis=1)

    print("Concatenation Time:", time.time() - tConcat)


    tScore = time.time()
    np.savetxt(scoresTxtPath, scoreConcatArr, delimiter="\t")
    print("Score SaveTxt Time:", time.time() - tScore)

    tObs = time.time()
    np.savetxt(observationsTxtPath, obsConcatArr, delimiter="\t")
    print("Observation SaveTxt:", time.time() - tObs)


    # print()
    # print("___STRUCTURED ARRAY TESTS___")

    # observationsTxtPath2 = outputDirPath / "observationsSAVETXT2_{}.txt.gz".format(fileTag)
    # scoresTxtPath2 = outputDirPath / "scoresSAVETXT2_{}.txt.gz".format(fileTag)

    # tStruct = time.time()
    # scoreTypes = (("U10", "int", "int") + ("int" for i in range(numStates)))
    # scoreNames = (("chr", "binstart", "binend") + tuple("s{}".format(i) for i in range(numStates)))
    # scoreArrStructured = np.zeros((scoreArr.shape[0], scoreArr.shape[1] + 3), dtype={"names":scoreNames, "formats":scoreTypes})
    # scoreArrStructured["chr"] = locationArr[:0]
    # scoreArrStructured["binstart"] = locationArr[:1]
    # scoreArrStructured["binend"] = locationArr[:2]
    # for i in range(numStates):
    #     scoreArrStructured["s{}".format(i)] = scoreArr[:,i]

    # obsTypes = ("U10", "int", "int", "int", "float", "int", "float")
    # obsNames = ("chr", "binstart", "binend", "maxloc", "maxval", "one", "totalscore")
    # obsArrStructured = np.zeros((observationArr.shape[0], 7), dtype={"names":obsNames, "formats":obsTypes})
    # obsArrStructured["chr"] = locationArr[:0]
    # obsArrStructured["binstart"] = locationArr[:1]
    # obsArrStructured["binend"] = locationArr[:2]
    # obsArrStructured["maxloc"] = observationArr[:,0]
    # obsArrStructured["maxval"] = observationArr[:,1]
    # obsArrStructured["one"] = np.ones((observationArr.shape[0], 1), dtype=int)
    # obsArrStructured["totalscore"] = observationArr[:,2]
    # print("Structured Creation Time:", time.time() - tStruct)

    # tScoreStruct = time.time()
    # np.savetxt(scoresTxtPath2, scoreArrStructured, delimiter="\t")
    # print("Score Structured SaveTxt Time:", time.time() - tScoreStruct)

    # tObsStruct = time.time()
    # np.savetxt(observationsTxtPath2, obsArrStructured, delimiter="\t")
    # print("Observation Structured SaveTxt Time:", time.time() - tObsStruct)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])