import sys
from pathlib import Path
import numpy as np
import pandas as pd
import time
import gzip

def main(file1, file2, outputDir1, outputDir2):
    file1Path = Path(file1)
    file2Path = Path(file2)
    outputDir1Path = Path(outputDir1)
    outputDir2Path = Path(outputDir2)
    filename1 = file1Path.name.split(".")[0] + ".txt.gz"
    filename2 = file2Path.name.split(".")[0] + ".txt.gz"

    # Read in the data
    print("\nReading data from file 1...")
    tRead1 = time.time()
    file1DF = pd.read_table(file1Path, header=None, sep="\t")
    print("    Time: ", time.time() - tRead1)

    print("\nReading data from file 2...")
    tRead2 = time.time()
    file2DF = pd.read_table(file2Path, header=None, sep="\t")
    print("    Time: ", time.time() - tRead2)


    # Converting to a np array for faster functions later
    print("Converting to numpy arrays...")
    tConvert = time.time()
    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=int)
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=int)
    locationArr = file1DF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    # Combining the arrays for per row shuffling
    combinedArr = np.concatenate((file1Arr, file2Arr), axis=1)

    # Row independent vectorized shuffling of the 2 arrays
    randomIndices = np.argpartition(np.random.rand(*combinedArr.shape), 1, axis=1)
    shuffledArr = np.take_along_axis(combinedArr, randomIndices, axis=1)
    file1OutputArr = np.concatenate((locationArr, shuffledArr[:,:file1Arr.shape[1]]), axis=1)
    file2OutputArr = np.concatenate((locationArr, shuffledArr[:,file1Arr.shape[1]:]), axis=1)

    print("input arrays sizes: ")
    print("    file1 =", file1Arr.shape)
    print("    file2 =", file2Arr.shape)
    print("    LocationArr =", locationArr.shape)
    print("Output arrays sizes: ")
    print("    file1 =", file1OutputArr.shape)
    print("    file2 =", file2OutputArr.shape)

    # Writing out the shuffled arrays
    np.savetxt(outputDir1Path / filename1, file1OutputArr, delimiter="\t", fmt="%s")
    np.savetxt(outputDir2Path / filename2, file2OutputArr, delimiter="\t", fmt="%s")

    # file1Txt = gzip.open(outputDir1Path / filename1, "wt")
    # file2Txt = gzip.open(outputDir2Path / filename2, "wt")

    # # Saving the files in the output directories
    # file1StrTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t" + "".join("{1[%d]:.5f}\t" % i for i in range(numStates-1)) + "{1[%d]:.5f}\n" % (numStates - 1)
    # file1Str = "".join(scoresTemplate.format(locationArr[i], scoreArr[i]) for i in range(scoreArr.shape[0]))
    
    # file1Txt.write(scoreStr)
    # file1Txt.close()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])