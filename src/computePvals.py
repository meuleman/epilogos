from pathlib import Path
import numpy as np
import scipy.stats as st
import sys
import pandas as pd

def main(file1, file2, outputDir, numStates, a, b, loc, scale, index):
    file1Path = Path(file1)
    file2Path = Path(file2)
    outputPath = Path(outputDir)

    names = ["chr", "binStart", "binEnd"]
    for i in range(1, numStates + 1):
        names.append("s{}".format(i))
        
    chrOrder = []
    for i in range(1, 23):
        chrOrder.append("chr{}".format(i))
    chrOrder.append("chrX")
    chrOrder.append("chrY")

    file1DF = pd.read_table(file1Path, header=None, sep="\s+", names=names)
    file2DF = pd.read_table(file2Path, header=None, sep="\s+", names=names)

    file1DF["chr"] = pd.Categorical(file1DF["chr"], categories=chrOrder, ordered=True)
    file2DF["chr"] = pd.Categorical(file2DF["chr"], categories=chrOrder, ordered=True)

    file1DF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)
    file2DF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)

    file1Arr = file1DF.iloc[:,3:].to_numpy(dtype=float)
    file2Arr = file2DF.iloc[:,3:].to_numpy(dtype=float)
    locationArr = file1DF.iloc[:,0:3].to_numpy(dtype=str)


    file1SumArr = np.sum(file1Arr, axis=1)
    file2SumArr = np.sum(file2Arr, axis=1)

    distances = np.sum(np.square(file1Arr - file2Arr), axis=1)

    size = 1000000

    rv = st.betaprime.rvs(a,b,loc=loc,scale=scale, size=size)

    if index + 50000 < len(distances):
        subValues = distances[index:index+50000].reshape(50000, 1)
        pvals = np.concatenate((locationArr[index:index+50000], ((rv > subValues).sum(axis=1) / float(size)).reshape(50000, 1)), axis=1)
    else:
        subValues = distances[index:].reshape(len(distances) - index, 1)
        pvals = np.concatenate((locationArr[index:], ((rv > subValues).sum(axis=1) / float(size)).reshape(len(distances) - index, 1)), axis=1)

    saveName = outputPath / "pvalarr_{}.npy".format(index)
    np.save(saveName, pvals, allow_pickle=False)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]), int(sys.argv[9]))