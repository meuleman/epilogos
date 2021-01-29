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

import timeit
import gzip
from pathlib import Path
import pandas as pd
import numpy as np


def main():
    timeArr = np.zeros((7, 10))
    for i in range(10):
        tTotal = timeit.default_timer()
        # Taking in the the score array
        filePath = Path("~meuleman/public_html/scores_data_pyData_male_epilogos_matrix_chr1.txt.gz")
        # Read in the data
        tRead = timeit.default_timer()
        # Challenge 1
        # dataDF = pd.read_table(filePath, header=None, sep="\t")

        # Challenge 2 500000 to 625,000
        dataDF = pd.read_table(filePath, skiprows=500000, nrows=625000-500000, header=None, sep="\t")
        timeArr[0][i] = timeit.default_timer() - tRead

        # Converting to a np array for faster functions later
        tConvert = timeit.default_timer()
        fileArr = dataDF.iloc[:,3:].to_numpy(dtype=np.float32)
        locationArr = dataDF.iloc[:,0:3].to_numpy(dtype=str)
        timeArr[1][i] = timeit.default_timer() - tConvert

        # Summing
        tSum = timeit.default_timer()
        fileArr = fileArr.sum(axis=1)
        timeArr[2][i] = timeit.default_timer() - tSum

        # Create one string of all the scores to write out
        tCreate = timeit.default_timer()
        scoresTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1:.5f}\n"
        scoreStr = "".join(scoresTemplate.format(locationArr[i], fileArr[i]) for i in range(fileArr.shape[0]))
        timeArr[3][i] = timeit.default_timer() - tCreate

        # Write out the string
        tScore = timeit.default_timer()
        # Challenge 1:
        # scoresTxtPath = Path("/home/jquon/fortnightFridayContest/scores_test.txt.gz")

        # Challenge 2:
        scoresTxtPath = Path("/home/jquon/fortnightFridayContest/scores_test_partial.txt.gz")
        scoresTxt = gzip.open(scoresTxtPath, "wt")
        scoresTxt.write(scoreStr)
        scoresTxt.close()
        timeArr[4][i] = timeit.default_timer() - tScore

        timeArr[5][i] = timeit.default_timer() - tTotal

        timeArr[6][i] = timeArr[:5][i].sum()

    print("\t\t" + "\t".join(i for i in range(10)))

    print("Reading Time:\t" + "\t".join(x for x in timeArr[0]))
    print("Convert Time:\t" + "\t".join(x for x in timeArr[1]))
    print("Summing Time:\t" + "\t".join(x for x in timeArr[2]))
    print("Out Str Time:\t" + "\t".join(x for x in timeArr[3]))
    print("Writing Time:\t" + "\t".join(x for x in timeArr[4]))
    print("Total   Time:\t" + "\t".join(x for x in timeArr[5]))
    print("Sum of Times:\t" + "\t".join(x for x in timeArr[6]))


if __name__ == "__main__":
    main()