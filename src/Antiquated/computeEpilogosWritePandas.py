import os
import numpy as np
from pathlib import Path
import gzip
import sys
import time
import pandas as pd

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
    observationsTxtPath = outputDirPath / "observationsTOCSV_{}.txt.gz".format(fileTag)
    scoresTxtPath = outputDirPath / "scoresTOCSV_{}.txt.gz".format(fileTag)

    # Order matters to us when writing, so use sorted
    # Loop over all score files and write them all to scores and observations txt
    tLoop = time.time()
    first = True
    for file in sorted(outputDirPath.glob("temp_scores_{}_*.npy".format(fileTag))):
        combinedArr = np.load(file, allow_pickle=False)

        fileScoreArr = combinedArr[:, 3:].astype(float)
        fileLocationArr = combinedArr[:, 0:3]

        # Calculating observation values
        maxContributions = np.amax(fileScoreArr, axis=1)
        maxContributionsLocs = np.argmax(fileScoreArr, axis=1)
        totalScores = np.sum(fileScoreArr, axis=1)

        # Storing the per file arrays into the entire array
        if first:
            first = False
            observationDF = pd.concat([pd.DataFrame(fileLocationArr), pd.DataFrame({3: maxContributionsLocs, 4: maxContributions, 5: np.ones(len(maxContributions), dtype=int), 6: totalScores})], axis=1)
            scoreDF = pd.concat([pd.DataFrame(fileLocationArr), pd.DataFrame(fileScoreArr)], axis=1)
        else:
            observationDF = observationDF.append(pd.concat([pd.DataFrame(fileLocationArr), pd.DataFrame({3: maxContributionsLocs, 4: maxContributions, 5: np.ones(len(maxContributions), dtype=int), 6: totalScores})], axis=1))
            scoreDF = scoreDF.append(pd.concat([pd.DataFrame(fileLocationArr), pd.DataFrame(fileScoreArr)], axis=1))
    print("Observation Calculation Time:", time.time() - tLoop)


    tScore = time.time()
    scoreDF.to_csv(scoresTxtPath, sep="\t", float_format="%.5f", index=False, header=False)
    print("Score To_CSV Time:", time.time() - tScore)

    tObs = time.time()
    observationDF.to_csv(observationsTxtPath, sep="\t", float_format="%.5f", index=False, header=False)
    print("Observation To_CSV Time:", time.time() - tObs)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])