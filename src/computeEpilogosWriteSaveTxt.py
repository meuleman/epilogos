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

    # # Order matters to us when writing, so use sorted
    # # Loop over all score files and write them all to scores and observations txt
    # tLoop = time.time()
    # first = True
    # for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
    #     combinedArr = np.load(file, allow_pickle=False)

    #     fileScoreArr = combinedArr[:, 3:].astype(float)
    #     fileLocationArr = combinedArr[:, 0:3]

    #     # Calculating observation values
    #     maxContributions = np.amax(fileScoreArr, axis=1).reshape((fileScoreArr.shape[0], 1))
    #     maxContributionsLocs = np.argmax(fileScoreArr, axis=1).reshape((fileScoreArr.shape[0], 1)) + 1
    #     totalScores = np.sum(fileScoreArr, axis=1).reshape((fileScoreArr.shape[0], 1))

    #     # Splicing all the observation arrays together
    #     fileObservationArr = np.concatenate((maxContributionsLocs.astype(str), np.around(maxContributions, decimals=5), np.ones((maxContributions.shape[0], 1), dtype=int), np.around(totalScores, decimals=5)), axis=1)

    #     # Storing the per file arrays into the entire array
    #     if first:
    #         scoreArr = np.around(fileScoreArr, decimals=5).astype(str)
    #         locationArr = fileLocationArr
    #         observationArr = fileObservationArr
    #         first = False
    #     else:
    #         scoreArr = np.concatenate((scoreArr, np.around(fileScoreArr, decimals=5).astype(str)), axis=0)
    #         locationArr = np.concatenate((locationArr, fileLocationArr), axis=0)
    #         observationArr = np.concatenate((observationArr, fileObservationArr), axis=0)

    # print("Observation Calculation Time:", time.time() - tLoop)

    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    tLoop = time.time()
    first = True
    scoreTypes = (("U5", "int", "int") + tuple("int" for i in range(numStates)))
    scoreNames = (("chr", "binstart", "binend") + tuple("s{}".format(i) for i in range(numStates)))
    obsTypes = ("U5", "int", "int", "int", "float", "int", "float")
    obsNames = ("chr", "binstart", "binend", "maxloc", "maxval", "one", "totalscore")
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
            first = False
            scoreArrStructured = np.zeros(fileScoreArr.shape[0], dtype={"names":scoreNames, "formats":scoreTypes})
            scoreArrStructured["chr"] = fileLocationArr[:,0]
            scoreArrStructured["binstart"] = fileLocationArr[:,1]
            scoreArrStructured["binend"] = fileLocationArr[:,2]
            for i in range(numStates):
                scoreArrStructured["s{}".format(i)] = fileScoreArr[:,i]

            obsArrStructured = np.zeros(fileObservationArr.shape[0], dtype={"names":obsNames, "formats":obsTypes})
            obsArrStructured["chr"] = fileLocationArr[:,0]
            obsArrStructured["binstart"] = fileLocationArr[:,1]
            obsArrStructured["binend"] = fileLocationArr[:,2]
            obsArrStructured["maxloc"] = fileObservationArr[:,0]
            obsArrStructured["maxval"] = fileObservationArr[:,1]
            obsArrStructured["one"] = np.ones(fileObservationArr.shape[0], dtype=int)
            obsArrStructured["totalscore"] = fileObservationArr[:,2]
        else:
            fileScoreArrStructured = np.zeros(fileScoreArr.shape[0], dtype={"names":scoreNames, "formats":scoreTypes})
            fileScoreArrStructured["chr"] = fileLocationArr[:,0]
            fileScoreArrStructured["binstart"] = fileLocationArr[:,1]
            fileScoreArrStructured["binend"] = fileLocationArr[:,2]
            for i in range(numStates):
                fileScoreArrStructured["s{}".format(i)] = fileScoreArr[:,i]

            fileObsArrStructured = np.zeros(fileObservationArr.shape[0], dtype={"names":obsNames, "formats":obsTypes})
            fileObsArrStructured["chr"] = fileLocationArr[:,0]
            fileObsArrStructured["binstart"] = fileLocationArr[:,1]
            fileObsArrStructured["binend"] = fileLocationArr[:,2]
            fileObsArrStructured["maxloc"] = fileObservationArr[:,0]
            fileObsArrStructured["maxval"] = fileObservationArr[:,1]
            fileObsArrStructured["one"] = np.ones(fileObservationArr.shape[0], dtype=int)
            fileObsArrStructured["totalscore"] = fileObservationArr[:,2]

            scoreArrStructured = np.append(scoreArrStructured, fileScoreArrStructured)
            obsArrStructured = np.append(obsArrStructured, fileObsArrStructured)

    print("Observation Calculation Time:", time.time() - tLoop)

    # tConcat = time.time()
    # scoreConcatArr = np.concatenate((locationArr, scoreArr), axis=1)
    # obsConcatArr = np.concatenate((locationArr, observationArr), axis=1)

    # print("Concatenation Time:", time.time() - tConcat)


    # tScore = time.time()
    # np.savetxt(scoresTxtPath, scoreConcatArr, delimiter="\t")
    # print("Score SaveTxt Time:", time.time() - tScore)

    # tObs = time.time()
    # np.savetxt(observationsTxtPath, obsConcatArr, delimiter="\t")
    # print("Observation SaveTxt:", time.time() - tObs)

    tScoreStruct = time.time()
    np.savetxt(scoresTxtPath, scoreArrStructured, delimiter="\t")
    print("Score Structured SaveTxt Time:", time.time() - tScoreStruct)

    tObsStruct = time.time()
    np.savetxt(observationsTxtPath, obsArrStructured, delimiter="\t")
    print("Observation Structured SaveTxt Time:", time.time() - tObsStruct)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])