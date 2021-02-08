import sys
import numpy as np
import os
from pathlib import Path
import time

def main(outputDirectory, storedExpInput, fileTag):
    tTotal = time.time()

    outputDirPath = Path(outputDirectory)
    storedExpPath = Path(storedExpInput)

    # Loop over all the expected value arrays and add them up
    expFreqFileCount = 0
    totalRows = 0
    expFreqArr = np.zeros((1,1))
    sums = []
    for file in outputDirPath.glob("temp_exp_freq_{}_*.npz".format(fileTag)):
        if expFreqFileCount == 0:
            npzFile = np.load(file)
            expFreqArr = npzFile["expFreqArr"]
            totalRows = npzFile["totalRows"][0]
            print("File: {}\t\tSum: {}\t\tRows:{}".format(file.name, np.sum(npzFile["expFreqArr"]) / npzFile["totalRows"][0], npzFile["totalRows"][0]))
            sums.append(np.sum(npzFile["expFreqArr"]) / npzFile["totalRows"][0])
        else:
            npzFile = np.load(file)
            expFreqArr += npzFile["expFreqArr"]
            totalRows += npzFile["totalRows"][0]
            print("File: {}\t\tSum: {}\t\tRows:{}".format(file.name, np.sum(npzFile["expFreqArr"]) / npzFile["totalRows"][0], npzFile["totalRows"][0]))
            sums.append(np.sum(npzFile["expFreqArr"]) / npzFile["totalRows"][0])
        expFreqFileCount += 1
    
    # Clean up temp files
    for file in outputDirPath.glob("temp_exp_freq_*.npz"):
        os.remove(file)

    # normalize for number of temp expected frequencies calculated
    expFreqArr /= totalRows

    print("Total Sum:", np.sum(expFreqArr))
    print("Sum of Sums:", sum(sums))

    # Save the expected frequency array
    np.save(storedExpPath, expFreqArr, allow_pickle=False)

    print("Total Time:", time.time() - tTotal)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])

