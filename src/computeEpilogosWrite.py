import os
import numpy as np
from pathlib import Path
import gzip
import sys
import time

def main(file, numStates, outputDirectory, fileTag, verbose):
    tTotal = time.time()
    outputDirPath = Path(outputDirectory)
    filePath = Path(file)

    filename = filePath.name.split(".")[0]

    # Initializing writing process
    scoresTxtPath = outputDirPath / "scores_{}_{}.txt.gz".format(fileTag, filename)
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    if verbose: print("\nReading in scores..."); tRead = time.time()
    # Taking in the the score array
    filePath = outputDirPath / Path("temp_scores_{}_{}.npz".format(fileTag, filename))
    npzFile = np.load(filePath)
    scoreArr = npzFile['scoreArr']
    locationArr = npzFile['locationArr']
    if verbose: print("    Time:", time.time() - tRead)

    # Create one string of all the scores to write out
    if verbose: print("Creating output string...")
    tCreate = time.time()
    scoresTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates-1)) + "{1[%d]:.5f}\n" % (numStates - 1)
    scoreStr = "".join(scoresTemplate.format(locationArr[i], scoreArr[i]) for i in range(scoreArr.shape[0]))
    if verbose: print("    Time:", time.time() - tCreate)

    # Write out the string
    if verbose: print("Writing scores...")
    tScore = time.time()
    scoresTxt.write(scoreStr)
    if verbose: print("    Time:", time.time() - tScore)
    
    scoresTxt.close()

    # Clean up temp file
    os.remove(filePath)

    print("Total Time:", time.time() - tTotal) if verbose else print("    {} Completed".format(filename))

# Helper for slurm to send boolean values
def strToBool(string):
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        raise ValueError("Invalid boolean string")

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], strToBool(sys.argv[5]))