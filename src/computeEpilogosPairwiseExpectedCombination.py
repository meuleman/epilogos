import sys
import numpy as np
import os
from pathlib import Path
import time

def main(outputDirectory, storedExpInput):
    tTotal = time.time()

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
    
    # Clean up temp files
    for file in outputDirPath.glob("temp_exp_freq_*.npy"):
        os.remove(file)

    expFreqArr /= count

    # Save the expected frequency array
    np.save(storedExpPath, expFreqArr, allow_pickle=False)

    print("Total Time:", time.time() - tTotal)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

