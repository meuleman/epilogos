import os
import numpy as np
from pathlib import Path
import gzip
import sys
import time

def main(file, numStates, outputDirectory, fileTag):
    tTotal = time.time()
    outputDirPath = Path(outputDirectory)
    filePath = Path(file)

    # chrName = filePath.name.split("_")[-1].split(".")[0]
    chrName = filePath.name.split(".")[0]

    writeScores(chrName, numStates, outputDirPath, fileTag)

    print("Total Time:", time.time() - tTotal)

# Helper to write the final scores to files
def writeScores(chrName, numStates, outputDirPath, fileTag):
    # Initializing writing process
    scoresTxtPath = outputDirPath / "scores_{}_{}.txt.gz".format(fileTag, chrName)
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    tRead = time.time()
    # Taking in the the score array
    filePath = outputDirPath / Path("temp_scores_{}_{}.npz".format(fileTag, chrName))
    npzFile = np.load(filePath)
    scoreArr = npzFile['scoreArr']
    locationArr = npzFile['locationArr']

    print("Read Time:", time.time() - tRead)

    # Create one string of all the scores to write out
    tCreate = time.time()
    scoresTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates-1)) + "{1[%d]:.5f}\n" % (numStates - 1)
    scoreStr = "".join(scoresTemplate.format(locationArr[i], scoreArr[i]) for i in range(scoreArr.shape[0]))
    print("String creation time:", time.time() - tCreate)

    # Write out the string
    tScore = time.time()
    scoresTxt.write(scoreStr)
    print("Score Write Time:", time.time() - tScore)
    
    scoresTxt.close()

    # Clean up temp file
    os.remove(filePath)

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4])