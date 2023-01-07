from sys import argv
import numpy as np
from os import remove
from pathlib import Path
from time import time
from epilogos.helpers import strToBool


def main(outputDirectory, storedExpInput, fileTag, verbose):
    """
    Combines all the temporary expected frequency numpy arrays into one expected frequency array

    Input:
    outputDirectory -- The output directory for epilogos
    storedExpInput  -- The path of the expected frequency array which is being created
    fileTag         -- A string which helps ensure outputed files are named similarly within an epilogos run
    verbose         -- Boolean which if True, causes much more detailed prints

    Output:
    A .npy file containing an array of the expected frequencies across all input files
    """
    if verbose: tTotal = time()

    outputDirPath = Path(outputDirectory)
    storedExpPath = Path(storedExpInput)

    # Loop over all the expected value arrays and add them up
    expFreqFileCount = 0
    expFreqArr = np.zeros((1, 1))
    for file in outputDirPath.glob("temp_exp_freq_{}_*.npy".format(fileTag)):
        if expFreqFileCount == 0:
            expFreqArr = np.load(file, allow_pickle=False)
        else:
            expFreqArr += np.load(file, allow_pickle=False)
        expFreqFileCount += 1

    # Clean up temp files
    for file in outputDirPath.glob("temp_exp_freq_*.npy"):
        remove(file)

    # normalize expected frequency array
    expFreqArr = (expFreqArr / np.sum(expFreqArr)).astype(np.float32)

    np.save(storedExpPath, expFreqArr, allow_pickle=False)

    print("Total Time:", time() - tTotal) if verbose else print("    [Done]")


if __name__ == "__main__":
    main(argv[1], argv[2], argv[3], strToBool(argv[4]))
