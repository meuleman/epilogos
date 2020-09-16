import sys
import numpy as np
import os
from pathlib import Path


def main(outputDirectory, fileTag, storeExp, storedExpInput):
    outputDirPath = Path(outputDirectory)
    storedExpPath = Path(storedExpInput)

    # Loop over all the expected value arrays and add them up (normalize for number of chromosomes)
    count = 0
    expFreqArr = np.zeros((1,1))
    for file in outputDirPath.glob("temp_exp_freq_*.npy"):
        if count == 0:
            expFreqArr = np.load(file, allow_pickle=False)
        else:
            expFreqArr += np.load(file, allow_pickle=False)
        count += 1
        # Delete file after we're done with it
        os.remove(file)
    expFreqArr /= count

    # If user desires, store away the expected frequency array
    if storeExp == "True":
        np.save(storedExpPath, expFreqArr, allow_pickle=False)

    # Save somewhere regardless for use in score calculation
    expFreqFilename = "temp_exp_freq_{}.npy".format(fileTag)
    expFreqPath = outputDirPath / expFreqFilename

    np.save(expFreqPath, expFreqArr, allow_pickle=False)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

