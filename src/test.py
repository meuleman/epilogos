import sys
import numpy as np
from pathlib import Path
import pandas as pd
import scipy.stats as st
import warnings
from time import time
from multiprocessing import cpu_count, Pool
from contextlib import closing
from itertools import repeat
from epilogosHelpers import strToBool, getNumStates

def main(stateInfo, outputDir, numProcesses, verbose):
    tTotal = time()

    outputDirPath = Path(outputDir)
    np.random.seed(7032016)
    
    numStates = getNumStates(stateInfo)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    # Read in observation files
    if verbose: print("\nReading in observation files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, distanceArrReal, distanceArrNull, maxDiffArr, diffArr = readInData(outputDirPath, numProcesses, numStates)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    # Fitting a gennorm distribution to the distances
    if verbose: print("Fitting gennorm distribution to distances...", flush=True); tFit = time()
    else: print("    Fitting distances\t", end="", flush=True)
    fitDistances(outputDirPath, distanceArrReal, distanceArrNull, diffArr, numStates)
    if verbose: print("    Time:", time() - tFit, flush=True)
    else: print("\t[Done]", flush=True)


    if verbose: print("Total Time:", time() - tTotal, flush=True)


# Helper to read in the necessary data to fit and visualize pairwise results
def readInData(outputDirPath, numProcesses, numStates):
    # For keeping the data arrays organized correctly
    realNames = ["chr", "binStart", "binEnd"] + ["s{}".format(i) for i in range(1, numStates + 1)]

    # Data frame to dump inputed data into
    diffDF = pd.DataFrame(columns=realNames)

    # Multiprocess the reading
    with closing(Pool(numProcesses)) as pool:
        results = pool.starmap(readTableMulti, zip(outputDirPath.glob("pairwiseDelta_*.txt.gz"), \
            outputDirPath.glob("temp_nullDistances_*.npz"), repeat(realNames)))
    pool.join()

    # Concatenating all chunks to the real differences dataframe
    for diffDFChunk, _ in results:
        diffDF = pd.concat((diffDF, diffDFChunk), axis=0, ignore_index=True)

    # Figuring out chromosome order
    chromosomes = diffDF.loc[diffDF['binStart'] == 0]['chr'].values
    rawChrNamesInts = []
    rawChrNamesStrs = []
    for chromosome in chromosomes:
        try:
            rawChrNamesInts.append(int(chromosome.split("chr")[-1]))
        except ValueError:
            rawChrNamesStrs.append(chromosome.split("chr")[-1])
    rawChrNamesInts.sort()
    rawChrNamesStrs.sort()
    chrOrder = rawChrNamesInts + rawChrNamesStrs
    for i in range(len(chrOrder)):
        chrOrder[i] = "chr" + str(chrOrder[i])

    # Sorting the dataframes by chromosomal location
    diffDF["chr"] = pd.Categorical(diffDF["chr"], categories=chrOrder, ordered=True)
    diffDF.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)

    # Convert dataframes to np arrays for easier manipulation
    locationArr     = diffDF.iloc[:,0:3].to_numpy(dtype=str)
    diffArr         = diffDF.iloc[:,3:].to_numpy(dtype=float)

    # Creating array of null distances ordered by chromosome based on the read in chunks
    nullChunks = list(zip(*list(zip(*results))[1]))
    index = nullChunks[0].index(chrOrder[0])
    distanceArrNull = nullChunks[1][index]
    for chrName in chrOrder[1:]:
        index = nullChunks[0].index(chrName)
        distanceArrNull = np.concatenate((distanceArrNull, nullChunks[1][index]))

    # Cleaning up the temp files after we've read them
    # for file in outputDirPath.glob("temp_nullDistances_*.npz"):
    #     remove(file)

    # Calculate the distance array for the real data
    diffSign = np.sign(np.sum(diffArr, axis=1))
    distanceArrReal = np.sum(np.square(diffArr), axis=1) * diffSign

    # Calculate the maximum contributing state for each bin
    # In the case of a tie, the higher number state wins (e.g. last state wins if all states are 0)
    maxDiffArr = np.abs(np.argmax(np.abs(np.flip(diffArr, axis=1)), axis=1) - diffArr.shape[1]).astype(int)

    return locationArr, distanceArrReal, distanceArrNull, maxDiffArr, diffArr


def readTableMulti(realFile, nullFile, realNames):
    diffDFChunk = pd.read_table(Path(realFile), header=None, sep="\t", names=realNames)
    npzFile = np.load(Path(nullFile))

    return diffDFChunk, (npzFile['chrName'][0], npzFile['nullDistances'])


# Helper to fit the distances
def fitDistances(outputDirPath, distanceArrReal, distanceArrNull, diffArr, numStates):
    # Filtering out quiescent values (When there are exactly zero differences between both score arrays)
    idx = [i for i in range(len(distanceArrReal)) if round(distanceArrReal[i], 5) != 0 \
        or np.any(diffArr[i] != np.zeros((numStates)))]
    dataNull = pd.Series(distanceArrNull[idx])

    with open(outputDirPath / "fits.txt", "w") as f:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            params = st.gennorm.fit(dataNull)
            mle = st.gennorm.nnlf(params, pd.Series(dataNull))
            f.write("Gennorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.johnsonsu.fit(dataNull)
            mle = st.johnsonsu.nnlf(params, pd.Series(dataNull))
            f.write("johnsonsu:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.t.fit(dataNull)
            mle = st.t.nnlf(params, pd.Series(dataNull))
            f.write("t:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.norminvgauss.fit(dataNull)
            mle = st.norminvgauss.nnlf(params, pd.Series(dataNull))
            f.write("norminvgauss:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.cauchy.fit(dataNull)
            mle = st.cauchy.nnlf(params, pd.Series(dataNull))
            f.write("cauchy:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.laplace.fit(dataNull)
            mle = st.laplace.nnlf(params, pd.Series(dataNull))
            f.write("laplace:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.hypsecant.fit(dataNull)
            mle = st.hypsecant.nnlf(params, pd.Series(dataNull))
            f.write("hypsecant:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.genlogistic.fit(dataNull)
            mle = st.genlogistic.nnlf(params, pd.Series(dataNull))
            f.write("genlogistic:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.logistic.fit(dataNull)
            mle = st.logistic.nnlf(params, pd.Series(dataNull))
            f.write("logistic:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.lognorm.fit(dataNull)
            mle = st.lognorm.nnlf(params, pd.Series(dataNull))
            f.write("lognorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.powernorm.fit(dataNull)
            mle = st.powernorm.nnlf(params, pd.Series(dataNull))
            f.write("powernorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.norm.fit(dataNull)
            mle = st.norm.nnlf(params, pd.Series(dataNull))
            f.write("norm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.loglaplace.fit(dataNull)
            mle = st.loglaplace.nnlf(params, pd.Series(dataNull))
            f.write("loglaplace:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.gumbel_l.fit(dataNull)
            mle = st.gumbel_l.nnlf(params, pd.Series(dataNull))
            f.write("gumbel_l:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.gumbel_r.fit(dataNull)
            mle = st.gumbel_r.nnlf(params, pd.Series(dataNull))
            f.write("gumbel_r:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.exponnorm.fit(dataNull)
            mle = st.exponnorm.nnlf(params, pd.Series(dataNull))
            f.write("exponnorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.skewnorm.fit(dataNull)
            mle = st.skewnorm.nnlf(params, pd.Series(dataNull))
            f.write("skewnorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.foldnorm.fit(dataNull)
            mle = st.foldnorm.nnlf(params, pd.Series(dataNull))
            f.write("foldnorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.powerlognorm.fit(dataNull)
            mle = st.powerlognorm.nnlf(params, pd.Series(dataNull))
            f.write("powerlognorm:\tParams={}\t\tMLE={}".format(params, mle))

            params = st.gausshyper.fit(dataNull)
            mle = st.gausshyper.nnlf(params, pd.Series(dataNull))
            f.write("gausshyper:\tParams={}\t\tMLE={}".format(params, mle))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), strToBool(sys.argv[4]))