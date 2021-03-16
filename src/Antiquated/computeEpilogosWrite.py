from os import remove
import numpy as np
from pathlib import Path
import gzip
from sys import argv
from time import time

def main(file, numStates, outputDirectory, fileTag, verbose):
    if verbose: tTotal = time()
    outputDirPath = Path(outputDirectory)
    filePath = Path(file)

    filename = filePath.name.split(".")[0]

    # Initializing writing process
    scoresTxtPath = outputDirPath / "scores_{}_{}.txt.gz".format(fileTag, filename)
    scoresTxt = gzip.open(scoresTxtPath, "wt")

    if verbose: print("\nReading in scores...", flush=True); tRead = time()
    # Taking in the the score array
    filePath = outputDirPath / Path("temp_scores_{}_{}.npz".format(fileTag, filename))
    npzFile = np.load(filePath)
    scoreArr = npzFile['scoreArr']
    locationArr = npzFile['locationArr']
    if verbose: print("    Time:", time() - tRead, flush=True)

    # Create one string of all the scores to write out
    if verbose: print("Creating output string...", flush=True)
    tCreate = time()
    scoresTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates-1)) + "{1[%d]:.5f}\n" % (numStates - 1)
    scoreStr = "".join(scoresTemplate.format(locationArr[i], scoreArr[i]) for i in range(scoreArr.shape[0]))
    if verbose: print("    Time:", time() - tCreate, flush=True)

    # Write out the string
    if verbose: print("Writing scores...", flush=True)
    tScore = time()
    scoresTxt.write(scoreStr)
    if verbose: print("    Time:", time() - tScore, flush=True)
    
    scoresTxt.close()

    # Clean up temp file
    remove(filePath)

    print("Total Time:", time() - tTotal, flush=True) if verbose else print("    {}\t\t[Done]".format(filename), flush=True)

# Helper for slurm to send boolean values
def strToBool(string):
    if string == 'True':
        return True
    elif string == 'False':
        return False
    else:
        raise ValueError("Invalid boolean string")

if __name__ == "__main__":
    main(argv[1], int(argv[2]), argv[3], argv[4], strToBool(argv[5]))