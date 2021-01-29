# import math
# import time
# import operator as op
# import numpy as np
# from functools import reduce
# import multiprocessing
# import ctypes as c
# import itertools
# import random
# from pathlib import Path
# import glob
# import pandas as pd
# import os
# import sys
# import subprocess
# from pathlib import PurePath
# import glob
# import gzip

import time
import gzip
from pathlib import Path
import pandas as pd
import numpy as np


def main():
    tTotal = time.time()
    # Taking in the the score array
    filePath = Path("~meuleman/public_html/scores_data_pyData_male_epilogos_matrix_chr1.txt.gz")
    # Read in the data
    print("\nReading data from file...")
    tRead = time.time()
    dataDF = pd.read_table(filePath, header=None, sep="\t")
    print("    Time: ", time.time() - tRead)

    # Converting to a np array for faster functions later
    print("Converting to numpy array...")
    tConvert = time.time()
    fileArr = dataDF.iloc[:,3:].to_numpy(dtype=np.float32)
    locationArr = dataDF.iloc[:,0:3].to_numpy(dtype=str)
    print("    Time: ", time.time() - tConvert)

    print(dataDF.iloc[:10])

    print(fileArr[:10])
    print(locationArr[:10])

    # Summing
    print("Summing up array...")
    tSum = time.time()
    fileArr = fileArr.sum(axis=1)
    print("    Time:", time.time() - tSum)

    print(fileArr[:10])

    # Create one string of all the scores to write out
    print("Creating Output String...")
    tCreate = time.time()
    scoresTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1:.5f}\n"
    scoreStr = "".join(scoresTemplate.format(locationArr[i], fileArr[i]) for i in range(fileArr.shape[0]))
    print("    Time:", time.time() - tCreate)

    print(scoreStr[:100])

    # Write out the string
    print("Writing Score...")
    tScore = time.time()
    scoresTxtPath = Path("/home/jquon/fortnightFridayContest/scores_test.txt.gz")
    scoresTxt = gzip.open(scoresTxtPath, "wt")
    scoresTxt.write(scoreStr)
    scoresTxt.close()
    print("    Time:", time.time() - tScore)
    
    print("Total Time:", time.time() - tTotal)

if __name__ == "__main__":
    main()