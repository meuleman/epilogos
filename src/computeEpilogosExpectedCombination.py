import sys
import numpy as np
import os
from pathlib import Path
import time

def main(outputDirectory, storedExpInput, saliency, fileTag):
    tTotal = time.time()

    outputDirPath = Path(outputDirectory)
    storedExpPath = Path(storedExpInput)

    # Loop over all the expected value arrays and add them up
    expFreqFileCount = 0
    numCols = 0
    expFreqArr = np.zeros((1,1))
    for file in outputDirPath.glob("temp_exp_freq_{}_*.npz".format(fileTag)):
        if expFreqFileCount == 0:
            npzFile = np.load(file)
            expFreqArr = npzFile["expFreqArr"]
            numCols = npzFile["numCols"][0]
        else:
            npzFile = np.load(file)
            expFreqArr += npzFile["expFreqArr"]
        expFreqFileCount += 1
    
    # Clean up temp files
    for file in outputDirPath.glob("temp_exp_freq_*.npz"):
        os.remove(file)

    # normalize expected frequency array

    expFreqArr = (expFreqArr / np.sum(expFreqArr)).astype(np.float32)

    # if saliency == 1:
    #     totalRows = np.sum(expFreqArr) / numCols
    #     print("TotalRows = ", totalRows)
    #     expFreqArr = expFreqArr.astype(np.float32) / (totalRows * numCols)
    # elif saliency == 2:
    #     combinations = ncr(numCols, 2) if sys.version_info < (3, 8) else math.comb(numCols, 2)
    #     totalRows = np.sum(expFreqArr) / combinations / 2
    #     print("Sum", np.sum(expFreqArr))
    #     print("combinations", combinations)
    #     print("TotalRows = ", totalRows)
    #     expFreqArr = expFreqArr.astype(np.float32) / (totalRows * combinations * 2)
    # elif saliency == 3:
    #     totalRows = np.sum(expFreqArr[0, 1])
    #     print("TotalRows = ", totalRows)
    #     expFreqArr = expFreqArr.astype(np.float32) / (totalRows * numCols * (numCols - 1))
    # else:
    #     raise ValueError("Inputted saliency value not supported")

    print("Total Sum:", np.sum(expFreqArr))

    # Save the expected frequency array
    np.save(storedExpPath, expFreqArr, allow_pickle=False)

    print("Total Time:", time.time() - tTotal)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])

