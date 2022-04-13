from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import scipy.stats as st
import warnings
from time import time
import gzip
import heapq
from multiprocessing import cpu_count, Pool
from contextlib import closing
from itertools import repeat
from os import remove
from statsmodels.stats.multitest import multipletests
from epilogos.helpers import strToBool, getStateNames, getStateColorsRGB, getNumStates


def main(group1Name, group2Name, stateInfo, outputDir, fileTag, numProcesses, pvalBool, diagnosticBool, numTrials,
         samplingSize, expFreqPath, exemplarWidth, verbose):
    """
    Takes in the scores for the 2 paired groups and finds the distance between them. Then fits a gennorm distribution to the
    distances between the null scores and uses this to calculate the pvalues of the distances. These pvalues are written out,
    used to generate txt file of most interesting loci, and used to create manhattan plots of the genome.

    Input:
    group1Name -- String of the name of the first epilogos group
    group2Name -- String of the name of the second epilogos group
    stateInfo -- State model tab seperated information file
    outputDir -- The output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    numProcesses -- The number of cores which to run on
    pvalBool -- Boolean which if True tells epilogos to generate pvalues for pairwise scores
    diagnosticBool -- Boolean which if True tells us to generate diagnostic plots of the gennorm fit on the null data and
                      comparisons between the null and real data
    numTrials -- The number of gennorm fits to do
    samplingSize -- The amount of null data to fit
    expFreqPath -- The location of the stored expected frequency array
    exemplarWidth -- 2*exemplarWidth+1 = size of exemplar regions
    verbose -- Boolean which if True, causes much more detailed prints
    """
    tTotal = time()

    outputDirPath = Path(outputDir)

    # Plotting setting
    plt.rcParams['agg.path.chunksize'] = 10000

    numStates = getNumStates(stateInfo)
    stateColorList = getStateColorsRGB(stateInfo)
    stateNameList = getStateNames(stateInfo)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    if pvalBool:
        # Fitting a gennorm distribution to the distances
        if verbose: print("\nFitting gennorm distribution to distances...", flush=True); tFit = time()
        else: print("    Fitting distances\t", end="", flush=True)
        params, distanceArrNull, nonQuiescentIdx = fitDistances(outputDirPath, numProcesses, numTrials, samplingSize)
        if verbose: print("    Time:", time() - tFit, flush=True)
        else: print("\t[Done]", flush=True)

    # Read in observation files
    if verbose: print("Reading in observation files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, distanceArrReal, maxDiffArr, chrDict = readInData(outputDirPath, numProcesses, numStates)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    if pvalBool:
        # Splitting the params up
        beta, loc, scale = params[:-2], params[-2], params[-1]

        # Creating Diagnostic Figures
        if diagnosticBool:
            if verbose: print("Creating diagnostic figures...", flush=True); tDiagnostic = time()
            else: print("    Diagnostic figures\t", end="", flush=True)
            createDiagnosticFigures(distanceArrReal, distanceArrNull, nonQuiescentIdx, beta, loc, scale, outputDirPath, fileTag)
            if verbose: print("    Time:", time() - tDiagnostic, flush=True)
            else: print("\t[Done]", flush=True)

        # Calculating PValues
        if verbose: print("Calculating P-Values...", flush=True); tPVal = time()
        else: print("    Calculating p-vals\t", end="", flush=True)
        pvals = calculatePVals(distanceArrReal, beta, loc, scale)
        if verbose: print("    Time:", time() - tPVal, flush=True)
        else: print("\t[Done]", flush=True)

        # Perform multiple hypothesis correction on pvals
        if verbose: print("Performing Benjamini-Hochberg procedure...", flush=True); tMH = time()
        else: print("    Benjamini-Hochberg procedure\t", end="", flush=True)
        mhPvals = multipletests(pvals, method="fdr_bh")[1]
        if verbose: print("    Time:", time() - tMH, flush=True)
        else: print("\t[Done]", flush=True)

    else:
        # Calculate the absolute value of zcores for creating the manhattan plots
        if verbose: print("Calculating Z-Scores...", flush=True); tZS = time()
        else: print("    Z-Scores\t", end="", flush=True)
        zScores = np.abs(st.zscore(distanceArrReal))
        if verbose:print("    Time:", time() - tZS, flush=True)
        else: print("\t[Done]", flush=True)

    # Create an output file which summarizes the results
    if verbose: print("Writing metrics file...", flush=True); tMetrics = time()
    else: print("    Writing metrics\t", end="", flush=True)
    if pvalBool:
        writeMetrics(locationArr, chrDict, maxDiffArr, stateNameList, distanceArrReal, outputDirPath, fileTag, pvalBool,
                     pvals=pvals, mhPvals=mhPvals)
    else:
        writeMetrics(locationArr, chrDict, maxDiffArr, stateNameList, distanceArrReal, outputDirPath, fileTag, pvalBool)
    if verbose: print("    Time:", time() - tMetrics, flush=True)
    else: print("\t[Done]", flush=True)

    if pvalBool:
        # Create txt file of significant loci with adjacent merged
        if verbose: print("Creating .txt file of top loci...", flush=True); tGreat = time()
        else: print("    Greatest hits txt\t", end="", flush=True)
        roiPath = outputDirPath / "greatestHits_{}.txt".format(fileTag)
        createTopScoresTxt(roiPath, locationArr, chrDict, distanceArrReal, maxDiffArr, stateNameList, pvals, False, mhPvals, exemplarWidth)
        if verbose: print("    Time:", time() - tGreat, flush=True)
        else: print("\t[Done]", flush=True)

        # Create txt file of significant loci
        if verbose: print("Creating .txt file of significant loci...", flush=True); tSig = time()
        else: print("    Significant loci txt\t", end="", flush=True)
        roiPath = outputDirPath / "significantLoci_{}.txt.gz".format(fileTag)
        createTopScoresTxt(roiPath, locationArr, chrDict, distanceArrReal, maxDiffArr, stateNameList, pvals, True, mhPvals, exemplarWidth)
        if verbose: print("    Time:", time() - tSig, flush=True)
        else: print("\t[Done]", flush=True)
    else:
        if verbose: print("Creating .txt file of top loci...", flush=True); tGreat = time()
        else: print("    Greatest hits txt\t", end="", flush=True)
        roiPath = outputDirPath / "greatestHits_{}.txt".format(fileTag)
        greatestHitsNoSignificance(roiPath, locationArr, chrDict, distanceArrReal, maxDiffArr, stateNameList, zScores, exemplarWidth)
        if verbose: print("    Time:", time() - tGreat, flush=True)
        else: print("\t[Done]", flush=True)

    # Create Chromosome Manhattan Plot
    if verbose: print("Creating Individual Chromosome Manhattan Plots", flush=True); tCManhattan = time()
    else: print("    Chromosome Manhattan\t", end="", flush=True)
    if pvalBool:
        createChromosomeManhattan(group1Name, group2Name, locationArr, chrDict, distanceArrReal, maxDiffArr,
                                  stateColorList, outputDirPath, fileTag, numProcesses, pvalBool, params=params,
                                  mhPvals=mhPvals)
    else:
        createChromosomeManhattan(group1Name, group2Name, locationArr, chrDict, distanceArrReal, maxDiffArr,
                                  stateColorList, outputDirPath, fileTag, numProcesses, pvalBool, zScores=zScores)
    if verbose: print("    Time:", time() - tCManhattan, flush=True)
    else: print("\t[Done]", flush=True)

    # Create Genome Manhattan Plot
    if verbose: print("Creating Genome-Wide Manhattan Plot", flush=True); tGManhattan = time()
    else: print("    Genome-wide Manhattan\t", end="", flush=True)
    if pvalBool:
        createGenomeManhattan(group1Name, group2Name, locationArr, chrDict, distanceArrReal, maxDiffArr,
                              stateColorList, outputDirPath, fileTag, pvalBool, beta=beta, loc=loc, scale=scale,
                              mhPvals=mhPvals)
    else:
        createGenomeManhattan(group1Name, group2Name, locationArr, chrDict, distanceArrReal, maxDiffArr,
                              stateColorList, outputDirPath, fileTag, pvalBool, zScores=zScores)
    if verbose: print("    Time:", time() - tGManhattan, flush=True)
    else: print("\t[Done]", flush=True)

    # Removing the expected frequency array
    remove(Path(expFreqPath))

    if verbose: print("Total Time:", time() - tTotal, flush=True)


def fitDistances(outputDirPath, numProcesses, numTrials, samplingSize):
    """
    Filters out quiescent bins and deploys the processes which fits the null distances. Then calculates the median fit based
    on the negative loglikelihood function

    Input:
    outputDirPath -- Path to the epilogos output directory (this contains the score files)
    numProcesses -- The number of cores to run on
    numTrials -- The number of fits to do
    samplingSize -- The amount of data to fit each time

    Output:
    (fitDF.iloc[medianIndex, 0], fitDF.iloc[medianIndex, 1], fitDF.iloc[medianIndex, 2]) -- Tuple with beta, loc, and scale
                                                                                            params of the median fit
    distanceArrNull -- Numpy array of the null distances
    nonQuiescentIdx -- indices of non-quiescent bins (specifically used for distance arrays)
    """
    # Filtering out quiescent values (When there are exactly zero differences between both score arrays)

    with closing(Pool(numProcesses)) as pool:
        results = pool.starmap(readNull, zip(outputDirPath.glob("temp_nullDistances_*.npz"),
                                             outputDirPath.glob("temp_quiescence_*.npz")))
    pool.join()

    # Figuring out chromosome order
    chromosomes = list(zip(*list(zip(*results))[0]))[0]
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

    # Creating array of null distances ordered by chromosome based on the read in chunks
    nullChunks = list(zip(*list(zip(*results))[0]))
    index = nullChunks[0].index(chrOrder[0])
    distanceArrNull = nullChunks[1][index]
    for chrName in chrOrder[1:]:
        index = nullChunks[0].index(chrName)
        distanceArrNull = np.concatenate((distanceArrNull, nullChunks[1][index]))

    # Creating quiescence array ordered by chromosome based on the read in chunks
    quiescenceChunks = list(zip(*list(zip(*results))[1]))
    index = quiescenceChunks[0].index(chrOrder[0])
    quiescenceArr = quiescenceChunks[1][index]
    for chrName in chrOrder[1:]:
        index = quiescenceChunks[0].index(chrName)
        quiescenceArr = np.concatenate((quiescenceArr, quiescenceChunks[1][index]))

    nonQuiescentIdx = np.where(quiescenceArr == False)[0]

    with closing(Pool(numProcesses)) as pool:
        results = pool.starmap(fitOnSample, zip(repeat(distanceArrNull[nonQuiescentIdx], numTrials),
                                                repeat(samplingSize, numTrials)))
    pool.join()

    # Creating dataframe of all params and nnlf so that we can figure out median
    index = [i for i in range(numTrials)]
    columns = ["beta", "loc", "scale", "nnlf"]
    fitDF = pd.DataFrame(index=index, columns=columns)
    for i in range(len(results)):
        fitDF.iloc[i, 0] = results[i][0][0]
        fitDF.iloc[i, 1] = results[i][0][1]
        fitDF.iloc[i, 2] = results[i][0][2]
        fitDF.iloc[i, 3] = results[i][1]
    fitDF.sort_values(by=["nnlf"], inplace=True)

    # return params, dataReal, dataNull
    medianIndex = int((numTrials-1)/2)
    return (fitDF.iloc[medianIndex, 0], fitDF.iloc[medianIndex, 1], fitDF.iloc[medianIndex, 2]),\
           distanceArrNull,\
           nonQuiescentIdx


def readNull(nullFile, quiescenceFile):
    """
    Reads in the null scores

    Input:
    nullFile -- The path to the file containing the null scores
    quiescenceFile -- The path to the file containing T/F regarding quiescence for each bin

    Output:
    (npzFile['chrName'][0], npzFile['nullDistances']) -- Tuple with the chromosome name and the null signed squared euclidean
                                                         distances
    (npzFileQuiescence['chrName'][0], npzFileQuiescence['quiescenceArr']) -- Tuple with chromosome name and the T/F value of
                                                                             quiescence for each bin
    """
    npzFileNull = np.load(Path(nullFile))
    npzFileQuiescence = np.load(Path(quiescenceFile))

    return (npzFileNull['chrName'][0], npzFileNull['nullDistances']), (npzFileQuiescence['chrName'][0],
                                                                       npzFileQuiescence['quiescenceArr'])


def fitOnSample(distanceArrNull, samplingSize):
    """
    Fits a sample of the null distances

    Input:
    distanceArrNull -- Numpy array containing the distances to sample for the fit
    samplingSize -- The size of the sample to fit

    Output:
    params -- The fit parameters obtained
    nnlf -- The negative loglikelihood function obtained by the fit
    """
    if len(distanceArrNull) <= samplingSize:
        sampleData = distanceArrNull
    else:
        np.random.seed()  # On linux, multiprocessing inherits the master seed and doesn't generate fully random numbers
        sampleData = pd.Series(np.random.choice(distanceArrNull, size=samplingSize, replace=False))

    # ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Fit the data
        params = st.gennorm.fit(sampleData)

        # Calculate SSE and MLE
        nnlf = st.gennorm.nnlf(params, pd.Series(distanceArrNull))

    return params, nnlf


def readInData(outputDirPath, numProcesses, numStates):
    """
    Reads all the epilogos score files in and combines them into a numpy array ordered by location

    Input:
    outputDirPath -- Path to the epilogos output directory (this contains the score files)
    numProcesses -- The number of cores used to read in score files
    numStates -- The number of states in the state model

    Output:
    locationArr -- Numpy array containing the genomic locations for all the scores
    distanceArrReal -- Numpy array containing binwise signed squared euclidean distances of the scores of the two groups
    maxDiffArr -- Numpy array containing the state which had the absolute distance in each bin
    chrDict -- Dictionary containing mappings between number values and chromosome names (helps locationArr use less memory)
    """
    # For keeping the data arrays organized correctly
    realNames = ["chr", "binStart", "binEnd"] + ["s{}".format(i) for i in range(1, numStates + 1)]

    # Data frame to dump inputed data into
    diffDF = pd.DataFrame(columns=realNames)

    # Multiprocess the reading
    with closing(Pool(numProcesses)) as pool:
        results = pool.starmap(readTableMulti, zip(outputDirPath.glob("pairwiseDelta_*.txt.gz"), repeat(realNames)))
    pool.join()

    # Concatenating all chunks to the real differences dataframe
    for diffDFChunk in results:
        diffDF = pd.concat((diffDF, diffDFChunk), axis=0, ignore_index=True)

    # Figuring out chromosome order
    chromosomes = diffDF['chr'].unique()
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

    # Creating a dictionary to make location array take less memory
    chrNumbers = [i for i in range(1, len(chrOrder) + 1)]
    chrDict    = dict(zip(chrNumbers, chrOrder))

    # Convert dataframes to np arrays for easier manipulation
    locationArr = diffDF.iloc[:, 0:3].replace({"chr": dict(zip(chrOrder, chrNumbers))}).to_numpy(dtype=np.int32)
    diffArr     = diffDF.iloc[:, 3:].to_numpy(dtype=np.float32)

    # Cleaning up the temp files after we've read them
    for file in outputDirPath.glob("temp_*.npz"):
        remove(file)

    # Calculate the signed squared euclidean distance array for the real data
    distanceArrReal = np.sum(np.square(diffArr), axis=1) * np.sign(np.sum(diffArr, axis=1))

    # Calculate the maximum contributing state for each bin
    # In the case of a tie, the higher number state wins (e.g. last state wins if all states are 0)
    maxDiffArr = np.abs(np.argmax(np.abs(np.flip(diffArr, axis=1)), axis=1) - diffArr.shape[1]).astype(np.int32)

    return locationArr, distanceArrReal, maxDiffArr, chrDict


def readTableMulti(realFile, realNames):
    """
    Reads in the real and null scores

    Input:
    realFile -- The path to the file containing the real scores

    Output:
    diffDFChunk -- Pandas dataframe containing the real scores
    """
    diffDFChunk = pd.read_table(Path(realFile), header=None, sep="\t", names=realNames)

    return diffDFChunk


def createDiagnosticFigures(distanceArrReal, distanceArrNull, nonQuiescentIdx, beta, loc, scale, outputDirPath, fileTag):
    """
    Generate diagnostic plots of the gennorm fit on the null data and comparisons between the null and real data

    Input:
    distanceArrReal -- Numpy array containing the real distances
    distanceArrNull -- Numpy array containing the null distances
    nonQuiescentIdx -- indices of non-quiescent bins (specifically used for distance arrays)
    beta -- gennorm fit parameter
    loc -- gennorm fit parameter
    scale -- gennorm fit parameter
    outputDirPath -- The path to the epilogos output directory
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    """
    diagnosticDirPath = outputDirPath / "diagnosticFigures_{}".format(fileTag)
    if not diagnosticDirPath.exists():
        diagnosticDirPath.mkdir(parents=True)

    dataReal = pd.Series(distanceArrReal[nonQuiescentIdx])
    dataNull = pd.Series(distanceArrNull[nonQuiescentIdx])

    # Fit on data (range=(min, max))
    y, x = np.histogram(dataNull, bins=400, range=(np.amin(distanceArrNull), np.amax(distanceArrNull)), density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0
    fig = plt.figure(figsize=(12, 8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta, loc, scale), legend=True)
    dataNull.plot(kind='hist', bins=400, range=(np.amin(distanceArrNull), np.amax(distanceArrNull)), density=True,
                  alpha=0.5, label='Null Data', legend=True, ax=ax)
    plt.title("Gennorm on null data (range=(min,max))")
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_minToMax.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)

    # Fit on data (range=(-1, 1))
    y, x = np.histogram(dataNull, bins=400, range=(-1, 1), density=True)
    x = (x + np.roll(x, -1))[:-1] / 2.0
    fig = plt.figure(figsize=(12, 8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta, loc, scale), legend=True)
    dataNull.plot(kind='hist', bins=400, range=(-1, 1), density=True, alpha=0.5, label='Null Data', legend=True, ax=ax)
    plt.title("Gennorm on null data (range=(-1,1))")
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_n1to1.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)

    # Real Data Histogram vs. Null Data Histogram (Range=(-1, 1))
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    dataReal.plot(kind='hist', bins=400, range=(-1, 1), density=True, alpha=0.5, label='Distances in Real Data', legend=True,
                  ax=ax)
    dataNull.plot(kind='hist', bins=400, range=(-1, 1), density=True, alpha=0.5, label='Distances in Null Data', legend=True,
                  ax=ax)
    plt.title("Real Data vs. Null Data (range=(-1, 1))")
    figPath = diagnosticDirPath / "real_vs_null_histogram_n1to1.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)

    # Real Data Histogram vs. Null Data Histogram (Range=(-max(abs), max(abs)))
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    rangeLim = np.amax(np.abs(dataReal))
    dataReal.plot(kind='hist', bins=400, range=(-rangeLim, rangeLim), density=True, alpha=0.5, label='Distances in Real Data',
                  legend=True, ax=ax)
    dataNull.plot(kind='hist', bins=400, range=(-rangeLim, rangeLim), density=True, alpha=0.5, label='Distances in Null Data',
                  legend=True, ax=ax)
    plt.title("Real Data vs. Null Data (range=(-max(abs), max(abs)))")
    figPath = diagnosticDirPath / "real_vs_null_histogram_minToMax.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)

    # Real vs Null distance scatter plot
    fig = plt.figure(figsize=(12, 12))
    plt.scatter(distanceArrReal, distanceArrNull, color='r', rasterized=True)
    plt.xlim(-rangeLim, rangeLim)
    plt.ylim(-rangeLim, rangeLim)
    plt.xlabel("Real Distances")
    plt.ylabel("Null Distances")
    plt.title("Real Distances vs Null Distances")
    figPath = diagnosticDirPath / "real_vs_null_scatter.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)

    # Box Plots Showcasing fit accuracy
    dist = st.gennorm.rvs(beta, loc=loc, scale=scale, size=dataNull.size)
    data = [dataNull, dist, dataReal]
    medianprops = dict(linewidth=2, color='black')
    boxprops = dict(linewidth=2, color='black')
    whiskerprops = dict(linewidth=2, color='black')
    capprops = dict(linewidth=2, color='black')
    fig = plt.figure(figsize=(12,8))
    bplot = plt.boxplot(data, patch_artist=True, medianprops=medianprops, boxprops=boxprops, whiskerprops=whiskerprops, capprops=capprops, showfliers=False)
    plt.xticks([1, 2, 3], ['Null', 'Fit', 'Real'])
    plt.xlabel("Data")
    plt.ylabel("Signed Squared Euclidean Distance")

    colors = ['#ff7f0e', '#bcbd22', '#d62728']
    for patch, color in zip(bplot['boxes'], colors):
        patch.set_facecolor(color)

    plt.title("Box Plots of Null and Real Data vs Fit")
    figPath = diagnosticDirPath / "null_vs_fit_vs_real_boxplots.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)


def calculatePVals(distanceArrReal, beta, loc, scale):
    """
    Calculates the pvalues of the real distances based on the fit of the null distances

    Input:
    distanceArrReal -- Numpy array containing real distances
    beta -- gennorm fit parameter
    loc -- gennorm fit parameter
    scale -- gennorm fit parameter

    Output:
    pvals -- Numpy array containing pvalues of the distances
    """
    pvalsBelowLoc = 2 * st.gennorm.cdf(distanceArrReal[np.where(distanceArrReal <= loc)[0]], beta, loc=loc, scale=scale)
    pvalsAboveLoc = 2 * (1 - st.gennorm.cdf(distanceArrReal[np.where(distanceArrReal > loc)[0]], beta, loc=loc, scale=scale))

    pvals = np.zeros(len(distanceArrReal))
    pvals[np.where(distanceArrReal <= loc)[0]] = pvalsBelowLoc
    pvals[np.where(distanceArrReal > loc)[0]]  = pvalsAboveLoc

    return pvals


def writeMetrics(locationArr, chrDict, maxDiffArr, nameArr, distanceArrReal, outputDirPath, fileTag, pvalBool,
                 pvals = [], mhPvals = []):
    """
    Writes metrics file to disk. Metrics file contains the following columns in order:
    chromosome, bin start, bin end, state with highest difference, signed squared euclidean distance, pvalue of distance

    Input:
    locationArr -- Numpy array containing the genomic locations of all the bins
    chrDict -- Dictionary containing mappings between number values and chromosome names (helps locationArr use less memory)
    maxDiffArr -- Numpy array containing the states which had the largest difference between the two groups in each bin
    nameArr -- Numpy array containing the names of all the states
    distanceArrReal -- Numpy array containing the real distances
    outputDirPath -- The path of the output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    pvalBool -- Boolean which if True tells epilogos to generate pvalues for pairwise scores
    pvals -- Numpy array containing the pvalues of all the distances
    mhPvals -- Multiple hypothesis corrected pvals using the Benjamini-Hochberg procedure
    """
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    metricsTxtPath = outputDirPath / "pairwiseMetrics_{}.txt.gz".format(fileTag)
    metricsTxt = gzip.open(metricsTxtPath, "wt")

    # Creating a string to write out the raw differences (faster than np.savetxt)
    if pvalBool:
        metricsTemplate = "{0}\t{1[0]}\t{1[1]}\t{2}\t{3:.5f}\t{4}\t{5:.5e}\t{6:.5e}\n"
        metricsStr = "".join(metricsTemplate.format(chrDict[locationArr[i, 0]], locationArr[i, 1:3],
                                                    nameArr[maxDiffArr[i] - 1], abs(distanceArrReal[i]),
                                                    findSign(distanceArrReal[i]), pvals[i], mhPvals[i])
                                                    for i in range(len(distanceArrReal)))
    else:
        metricsTemplate = "{0}\t{1[0]}\t{1[1]}\t{2}\t{3:.5f}\t{4}\n"
        metricsStr = "".join(metricsTemplate.format(chrDict[locationArr[i, 0]], locationArr[i, 1:3],
                                                    nameArr[maxDiffArr[i] - 1], abs(distanceArrReal[i]),
                                                    findSign(distanceArrReal[i]))
                                                    for i in range(len(distanceArrReal)))

    metricsTxt.write(metricsStr)
    metricsTxt.close()


def createTopScoresTxt(filePath, locationArr, chrDict, distanceArr, maxDiffArr, nameArr, pvals, onlySignificant, mhPvals, exemplarWidth):
    """
    Finds the either the 1000 largest distance bins and merges adjacent bins or finds all significant loci and does not merge.
    Then it outputs a txt containing these highest distances regions and some information about each (chromosome, bin start,
    bin end, state name, absolute value of distance, sign of score, pvalue of distance, stars repersenting significance)

    Input:
    filePath -- The path of the file to write to
    locationArr -- Numpy array containing the genomic locations of all the bins
    chrDict -- Dictionary containing mappings between number values and chromosome names (helps locationArr use less memory)
    distanceArr -- Numpy array containing the distance between the pairwise groups
    maxDiffArr -- Numpy array containing the states which had the largest difference in each bin
    nameArr -- Numpy array containing the names of all the states
    pvals -- Numpy array containing the pvalues of all the distances between the pairwise groups
    onlySignificant -- Boolean telling us whether to use significant loci or 1000 largest distance bins
    mhPvals -- Multiple hypothesis corrected pvals using the Benjamini-Hochberg procedure
    exemplarWidth -- 2*exemplarWidth+1 = size of exemplar regions
    """
    f = gzip.open(filePath, "wt") if onlySignificant else open(filePath, 'w')

    # Pick values above significance threshold
    if onlySignificant:
        indices = maxIndices = np.where(mhPvals <= 0.1)[0]
    else:
        # in the case of creating greatest hits indices should be a intersection of all significantLoci and regions generated by priorityqueue
        indices, maxIndices = meanAndMax(distanceArr, 100, exemplarWidth)
        maxIndices, _, comm2 = np.intersect1d(np.where(mhPvals <= 0.1)[0], maxIndices, assume_unique=True, return_indices=True)
        indices = np.array(indices)[comm2]

    if len(indices) == 0:
        f.close()
        return

    locations = pd.DataFrame(np.concatenate((locationArr[indices], distanceArr[maxIndices].reshape(len(maxIndices), 1),
                                             maxDiffArr[maxIndices].reshape(len(maxIndices), 1),
                                             pvals[maxIndices].reshape(len(maxIndices), 1),
                                             mhPvals[maxIndices].reshape(len(maxIndices), 1)), axis=1),
                             columns=["Chromosome", "Start", "End", "Score", "MaxDiffLoc", "Pval", "MhPval"])\
                  .astype({"Chromosome": np.int32, "Start": np.int32, "End": np.int32, "Score": np.float32,
                           "MaxDiffLoc": np.int32, "Pval": np.float32, "MhPval": np.float32})\
                  .replace({"Chromosome": chrDict})

    # Don't want to merge when creating significantLoci.txt
    if not onlySignificant:
        locations = locations.iloc[(-locations["Score"].abs()).argsort()]
        # Pad the bins
        locations["Start"] = locations["Start"].to_numpy(dtype=np.int32) - 200 * exemplarWidth
        locations["End"]   = locations["End"].to_numpy(dtype=np.int32) + 200 * exemplarWidth

    # Locations get 3 stars if they are significant at .01, 2 stars at .05, 1 star at .1, and a period if not significant
    stars = np.array(["***" if float(locations.iloc[i, 6]) <= 0.01 else
                      ("**" if float(locations.iloc[i, 6]) <= 0.05 else
                      ("*" if float(locations.iloc[i, 6]) <= 0.1 else "."))
                      for i in range(locations.shape[0])])

    # Write all the locations to the file for significantLoci.txt
    # Write only top 100 loci to file for greatestHits.txt
    outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3}\t{4:.5e}\t{5:.5e}\t{6}\n"
    outString = "".join(outTemplate.format(locations.iloc[i], nameArr[int(float(locations.iloc[i, 4])) - 1],
                                           abs(float(locations.iloc[i, 3])), findSign(float(locations.iloc[i, 3])),
                                           float(locations.iloc[i, 5]), float(locations.iloc[i, 6]), stars[i])
                        for i in range(locations.shape[0]))

    f.write(outString)
    f.close()


def findSign(x):
    """
    Returns a string containing the sign of the inputted number

    Input:
    x -- Any number

    Output:
    "+" or "-"
    """
    if (x >= 0):
        return "+"
    else:
        return "-"


def meanAndMax(scoreArr, k, w):
    """
    Takes the scores of all the bins in the genome and takes the top k regions of width w.
    If a region is within a higher scoring region it is removed from the list and a new region is added

    Input:
    scoreArr -- Numpy array containing the distance between the pairwise groups
    k -- the number of regions we want to end up with
    w -- the width (in number of bins) to search in each direction
    """
    # Generate the table used to
    n = len(scoreArr)
    hits = np.zeros(n, dtype=np.bool)
    priorityQueue = []

    # Have to make scores negative so they pop off in correct order
    for i, v in enumerate(np.abs(scoreArr)):
        heapq.heappush(priorityQueue, (-v, i))

    indices = []
    maxIndices = []
    while k > 0:
        try:
            (v, i) = heapq.heappop(priorityQueue)
            for j in findBestIndex(scoreArr, i, w):
                start = (j - w) if (j - w) > 0 else 0
                stop = (j + w + 1) if (j + w + 1) <= n else n
                if not np.any(hits[start:stop]):
                    hits[start:stop] = True
                    indices.append(j)
                    maxIndices.append(i)
                    k -= 1
                    break
        except IndexError:
            k = 0

    return indices, maxIndices


def findBestIndex(scoreArr, i, exemplarWidth):
    """
    Within a given range [i - exemplarWidth, i + exemplarWidth], find the index, j, in the score array with the highest mean
    score of the range [j - exemplarWidth, j + exemplarWidth + 1].
    ie. we want to find the index, j, which contains the highest mean score in its range while also containing i in said range.

    Input:
    scoreArr -- 1D numpy array containing the scores of each bin
    i -- The index in the score array around which we will find the optimal index
    exemplarWidth -- 2*exemplarWidth+1 = size of exemplar regions

    Output:
    Array of indices sorted in descending order of mean score
    """
    farthestLeft = i - exemplarWidth
    farthestRight = i + exemplarWidth
    means = []
    for j in range(farthestLeft, farthestRight + 1):
        means.append(-np.mean(np.abs(scoreArr[j - exemplarWidth:j + exemplarWidth + 1])))
    return np.arange(farthestLeft, farthestRight + 1)[np.argsort(means)]


def greatestHitsNoSignificance(filePath, locationArr, chrDict, distanceArr, maxDiffArr, nameArr, zScores, exemplarWidth):
    """
    Writes only the locations out to file which are specified by priorityQueueMerging. Used to write the results of
    priority queue merging

    Input:
    filePath -- The path of the file to write to
    locationArr -- Numpy array containing the genomic locations of all the bins
    chrDict -- Dictionary containing mappings between number values and chromosome names (helps locationArr use less memory)
    distanceArr -- Numpy array containing the distance between the pairwise groups
    maxDiffArr -- Numpy array containing the states which had the largest difference in each bin
    nameArr -- Numpy array containing the names of all the states
    zScores -- Absolute value of the ZScores of the distances
    exemplarWidth -- 2*exemplarWidth+1 = size of exemplar regions    """
    f = open(filePath, 'w')

    indices, maxIndices = meanAndMax(distanceArr, 100, exemplarWidth)

    locations = pd.DataFrame(np.concatenate((locationArr[indices], distanceArr[maxIndices].reshape(len(maxIndices), 1),
                                             maxDiffArr[maxIndices].reshape(len(maxIndices), 1),
                                             zScores[maxIndices].reshape(len(maxIndices), 1)), axis=1),
                             columns=["Chromosome", "Start", "End", "Score", "MaxDiffLoc", "ZScores"])\
                  .astype({"Chromosome": np.int32, "Start": np.int32, "End": np.int32, "Score": np.float32,
                           "MaxDiffLoc": np.int32, "ZScores": np.float32})\
                  .replace({"Chromosome": chrDict})

    # Pad the bins
    locations["Start"] = locations["Start"].to_numpy(dtype=np.int32) - 200 * exemplarWidth
    locations["End"]   = locations["End"].to_numpy(dtype=np.int32) + 200 * exemplarWidth

    # Locations get 3 stars if they are 3 sd away, 2 at 2 sd, 1 at 1 sd, and a period if withing 1 sd
    stars = np.array(["***" if float(locations.iloc[i, 5]) >= 3 else
                      ("**" if float(locations.iloc[i, 5]) >= 2 else
                      ("*" if float(locations.iloc[i, 5]) >= 1 else "."))
                      for i in range(locations.shape[0])])

    # Write all the locations to the file for significantLoci.txt
    # Write only top 100 loci to file for greatestHits.txt
    outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3}\t{4:.5f}\t{5}\n"
    outString = "".join(outTemplate.format(locations.iloc[i], nameArr[int(float(locations.iloc[i, 4])) - 1],
                                           abs(float(locations.iloc[i, 3])), findSign(float(locations.iloc[i, 3])),
                                           float(locations.iloc[i, 5]), stars[i])
                        for i in range(len(indices)))

    f.write(outString)
    f.close()

def createGenomeManhattan(group1Name, group2Name, locationArr, chrDict, distanceArrReal, maxDiffArr, stateColorList,
                          outputDirPath, fileTag, pvalBool, beta=0, loc=0, scale=0, mhPvals=[], zScores=[]):
    """
    Creates a manhattan plot based on the distances between the two groups for the entire genome

    Input:
    group1Name -- The name of the first epilogos group
    group2Name -- The name of the second epilogos group
    locationArr -- Numpy array containing the genomic locations of all the bins
    chrDict -- Dictionary containing mappings between number values and chromosome names (helps locationArr use less memory)
    distanceArrReal -- Numpy array containing the real distances
    maxDiffArr -- Numpy array containing the states which had the largest difference between the two groups in each bin
    stateColorList -- Numpy array containing the colors of each of the states in the state model
    outputDirPath -- The path of the output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    pvalBool -- A boolean telling us whether we are using pvals or not
    beta -- gennorm fit parameter
    loc -- gennorm fit parameter
    scale -- gennorm fit parameter
    mhPvals -- Multiple hypothesis corrected pvals using the Benjamini-Hochberg procedure
    zScores -- Absolute value of the Z-Scores of the scores in the case that we are not using pvals
    """
    manhattanDirPath = outputDirPath / "manhattanPlots_{}".format(fileTag)
    if not manhattanDirPath.exists():
        manhattanDirPath.mkdir(parents=True)

    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.set_facecolor("#FFFFFF")
    ax.set_axisbelow(True)
    ax.grid(True, axis='both', color='k', linewidth=.25, linestyle="-")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    plt.title("Differential epilogos between {} and {} biosamples".format(group1Name, group2Name))
    ax.set_ylabel("Distance")
    plt.xlabel("Chromosome")
    xticks = np.where(locationArr[:, 1] == 0)[0]
    plt.xticks(ticks=xticks, labels=list(map(lambda x: x.split("chr")[-1], [chrDict[x] for x in locationArr[:, 0][xticks]])))

    plt.margins(x=0)
    ylim = np.amax(np.abs(distanceArrReal)) * 1.1
    ax.set_ylim(-ylim, ylim)
    yticks, ytickLabels = pvalAxisScaling(ylim, beta, loc, scale) if pvalBool else zScoreAxisScaling(ylim, np.mean(distanceArrReal), np.std(distanceArrReal))

    ax.set_yticks(yticks)
    ax.set_yticklabels([str(np.abs(np.round(val, 1))) for val in yticks])

    axR = ax.twinx()
    axR.set_ylabel("P-Value" if pvalBool else "Z-Score")
    axR.spines["top"].set_visible(False)
    axR.spines["left"].set_visible(False)
    axR.spines["bottom"].set_visible(False)
    axR.set_yticks(yticks)
    axR.set_ylim(ax.get_ylim())
    axR.set_yticklabels(ytickLabels)

    ax.text(0.99, 0.99, group1Name, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=15)
    ax.text(0.99, 0.01, group2Name, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes,
            fontsize=15)

    locationOnGenome = np.arange(len(distanceArrReal))

    for i in range(len(xticks)):
        if i == len(xticks)-1:
            points = np.where((locationOnGenome >= xticks[i]) & (mhPvals > 0.1))[0] if pvalBool else np.where((locationOnGenome >= xticks[i]) & (zScores < 10))[0]
            color = "gray" if i % 2 == 0 else "black"
            plt.scatter(locationOnGenome[points], distanceArrReal[points],
                        s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color=color,
                        marker=".", alpha=0.1, edgecolors='none', rasterized=True)
        elif i % 2 == 0:
            points = np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1]) & (mhPvals > 0.1))[0] if pvalBool else np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1]) & (zScores < 10))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points],
                        s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray",
                        marker=".", alpha=0.1, edgecolors='none', rasterized=True)
        else:
            points = np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1]) & (mhPvals > 0.1))[0] if pvalBool else np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1]) & (zScores < 10))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points],
                        s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="black",
                        marker=".", alpha=0.1, edgecolors='none', rasterized=True)

    line1Indices = np.where(mhPvals <= 0.1)[0] if pvalBool else np.where(zScores >= 10)[0]
    line2Indices = np.where(mhPvals <= 0.05)[0] if pvalBool else np.where(zScores >= 20)[0]
    line3Indices = np.where(mhPvals <= 0.01)[0] if pvalBool else np.where(zScores >= 30)[0]

    colorArr = stateColorList[maxDiffArr[line1Indices].astype(int) - 1]
    opacityArr = np.array((np.abs(distanceArrReal[line1Indices]) /
                           np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[line1Indices]), 1)
    rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
    sizeArr = np.abs(distanceArrReal[line1Indices]) / np.amax(np.abs(distanceArrReal)) * 100

    plt.scatter(line1Indices, distanceArrReal[line1Indices], s=sizeArr, color=rgbaColorArr, marker=".",
                edgecolors='none', rasterized=True)

    if len(line3Indices) > 0:
        point1Line = np.min(np.abs(distanceArrReal[line1Indices]))
        point05Line = np.min(np.abs(distanceArrReal[line2Indices]))
        point01Line = np.min(np.abs(distanceArrReal[line3Indices]))
        plt.axhspan(point1Line, point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(point05Line, point01Line, facecolor='black', alpha=0.10)
        plt.axhspan(point01Line, ylim, facecolor='black', alpha=0.15)
        plt.axhspan(-point1Line, -point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(-point05Line, -point01Line, facecolor='black', alpha=0.10)
        plt.axhspan(-point01Line, -ylim, facecolor='black', alpha=0.15)
    elif len(line2Indices) > 0:
        point1Line = np.min(np.abs(distanceArrReal[line1Indices]))
        point05Line = np.min(np.abs(distanceArrReal[line2Indices]))
        plt.axhspan(point1Line, point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(point05Line, ylim, facecolor='black', alpha=0.10)
        plt.axhspan(-point1Line, -point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(-point05Line, -ylim, facecolor='black', alpha=0.10)
    elif len(line1Indices) > 0:
        point1Line = np.min(np.abs(distanceArrReal[line1Indices]))
        plt.axhspan(point1Line, ylim, facecolor='black', alpha=0.05)
        plt.axhspan(-point1Line, -ylim, facecolor='black', alpha=0.05)

    figPath = manhattanDirPath / "manhattan_plot_genome.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)


def _initChromosomeManhattan(group1Name_, group2Name_, locationArr_, distanceArrReal_, maxDiffArr_, params, mhPvals_, zScores_,
                             stateColorList_, manhattanDirPath_, pvalBool_):
    """
    Initializes global variables for multiprocessing in the single epilogos case

    Input:
    group1Name_ -- The name of the first epilogos group
    group2Name_ -- The name of the second epilogos group
    locationArr_ -- Numpy array containing the genomic locations of all the bins
    distanceArrReal_ -- Numpy array containing the real distances
    maxDiffArr_ -- Numpy array containing the states which had the largest difference between the two groups in each bin
    params -- gennorm fit parameters
    mhPvals_ -- Multiple hypothesis corrected pvals using the Benjamini-Hochberg procedure
    zScores_ -- Absolute values of the z-scores of the distances
    stateColorList_ -- Numpy array containing the colors of each of the states in the state model
    manhattanDirPath_ -- The path to directory to put manhattan plots
    pvalBool_ -- A boolean to tell us whether we will use pvals or zscores
    """
    global group1Name
    global group2Name
    global locationArr
    global distanceArrReal
    global maxDiffArr
    global beta
    global loc
    global scale
    global mhPvals
    global zScores
    global stateColorList
    global manhattanDirPath
    global pvalBool

    group1Name = group1Name_
    group2Name = group2Name_
    locationArr = locationArr_
    distanceArrReal = distanceArrReal_
    maxDiffArr = maxDiffArr_
    beta, loc, scale = params[:-2], params[-2], params[-1]
    mhPvals = mhPvals_
    zScores = zScores_
    stateColorList = stateColorList_
    manhattanDirPath = manhattanDirPath_
    pvalBool = pvalBool_


def createChromosomeManhattan(group1Name, group2Name, locationArr, chrDict, distanceArrReal, maxDiffArr, stateColorList,
                              outputDirPath, fileTag, numProcesses, pvalBool, params=[0,0,0], mhPvals=[], zScores=[]):
    """
    Creates a manhattan plot based on the distances between the two groups for the each chromosome

    Input:
    group1Name -- The name of the first epilogos group
    group2Name -- The name of the second epilogos group
    locationArr -- Numpy array containing the genomic locations of all the bins
    chrDict -- Dictionary containing mappings between number values and chromosome names (helps locationArr use less memory)
    distanceArrReal -- Numpy array containing the real distances
    maxDiffArr -- Numpy array containing the states which had the largest difference between the two groups in each bin
    stateColorList -- Numpy array containing the colors of each of the states in the state model
    outputDirPath -- The path of the output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    numProcesses -- The number of cores to run on
    pvalBool -- A boolean value to tell us whether we are going to use pvals or not
    params -- gennorm fit parameters
    mhPvals -- Multiple hypothesis corrected pvals using the Benjamini-Hochberg procedure
    zScores -- Absolute values of the Z-Scores of the distances
    """
    manhattanDirPath = outputDirPath / "manhattanPlots_{}".format(fileTag)
    if not manhattanDirPath.exists():
        manhattanDirPath.mkdir(parents=True)

    xticks = np.where(locationArr[:, 1] == 0)[0]
    startEnd = []
    for i in range(len(xticks)):
        if not i == len(xticks) - 1:
            startEnd.append((xticks[i], xticks[i+1]))
        else:
            startEnd.append((xticks[i], -1))

    chrOrder = list(map(lambda x: x.split("chr")[-1], [chrDict[x] for x in locationArr[:, 0][xticks]]))

    # Multiprocess the reading
    with closing(Pool(numProcesses, initializer=_initChromosomeManhattan, initargs=(group1Name, group2Name, locationArr,
                                                                                    distanceArrReal, maxDiffArr, params,
                                                                                    mhPvals, zScores, stateColorList,
                                                                                    manhattanDirPath, pvalBool))) as pool:
        pool.starmap(graphChromosomeManhattan, zip(chrOrder, startEnd))
    pool.join()


def graphChromosomeManhattan(chromosome, startEnd):
    """
    Generates the manhattan plots for each chromosome based on the distances between the two groups.
    Note there are global variables which are used across all chromosomes (see _init())

    Input:
    chromosome -- The chromosome which we are plotting
    startEnd -- The start and end indices for the chromosome on all the global numpy arrays

    Also see global variables from _initChromosomeManhattan above
    """
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(111)
    ax.set_facecolor("#FFFFFF")
    ax.set_axisbelow(True)
    ax.grid(True, axis='y', color='k', linewidth=.25, linestyle="-")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_ylabel("Distance")
    plt.xlabel("Location in Chromosome {} (Mb)".format(chromosome))
    plt.title("Differential epilogos between {} and {} donor biosamples (Chromosome {})".format(group1Name, group2Name,
                                                                                                chromosome))

    plt.margins(x=0)
    ylim = np.amax(np.abs(distanceArrReal)) * 1.1
    ax.set_ylim(-ylim, ylim)
    yticks, ytickLabels = pvalAxisScaling(ylim, beta, loc, scale) if pvalBool else zScoreAxisScaling(ylim, np.mean(distanceArrReal), np.std(distanceArrReal))

    ax.set_yticks(yticks)
    ax.set_yticklabels([str(np.abs(np.round(val, 1))) for val in yticks])

    axR = ax.twinx()
    axR.set_ylabel("P-Value" if pvalBool else "Z-Score")
    axR.spines["top"].set_visible(False)
    axR.spines["left"].set_visible(False)
    axR.spines["bottom"].set_visible(False)
    axR.set_yticks(yticks)
    axR.set_ylim(ax.get_ylim())
    axR.set_yticklabels(ytickLabels)

    ax.text(0.99, 0.99, group1Name, verticalalignment='top', horizontalalignment='right', transform=ax.transAxes, fontsize=15)
    ax.text(0.99, 0.01, group2Name, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes,
            fontsize=15)

    locationOnGenome = np.arange(len(distanceArrReal))

    if startEnd[1] == -1:
        realxticks = np.where((locationOnGenome >= startEnd[0]) & (locationArr[:, 1].astype(int) % 10000000 == 0))[0]
        plt.xticks(ticks=realxticks, labels=[str(int(int(locationArr[tick, 1])/1000000)) for tick in realxticks])

        points = np.where((locationOnGenome >= startEnd[0]) & (mhPvals > 0.1))[0] if pvalBool else np.where((locationOnGenome >= startEnd[0]) & (zScores < 10))[0]
        plt.scatter(locationOnGenome[points], distanceArrReal[points],
                    s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".",
                    alpha=0.1, edgecolors='none', rasterized=True)

        line1Indices = np.where((locationOnGenome >= startEnd[0]) & (mhPvals <= 0.1))[0] if pvalBool else np.where((locationOnGenome >= startEnd[0]) & (zScores >= 10))[0]
        line2Indices = np.where((locationOnGenome >= startEnd[0]) & (mhPvals <= 0.05))[0] if pvalBool else np.where((locationOnGenome >= startEnd[0]) & (zScores >= 20))[0]
        line3Indices = np.where((locationOnGenome >= startEnd[0]) & (mhPvals <= 0.01))[0] if pvalBool else np.where((locationOnGenome >= startEnd[0]) & (zScores >= 30))[0]

    else:
        realxticks = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome < startEnd[1]))
                              & (locationArr[:, 1].astype(int) % 10000000 == 0))[0]
        plt.xticks(ticks=realxticks, labels=[str(int(int(locationArr[tick, 1])/1000000)) for tick in realxticks])

        points = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome < startEnd[1])) & (mhPvals > 0.1))[0] if pvalBool else np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome < startEnd[1])) & (zScores < 10))[0]
        plt.scatter(locationOnGenome[points], distanceArrReal[points],
                    s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".",
                    alpha=0.1, edgecolors='none', rasterized=True)

        line1Indices = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome <= startEnd[1])) & (mhPvals <= 0.1))[0] if pvalBool else np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome <= startEnd[1])) & (zScores >= 10))[0]
        line2Indices = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome <= startEnd[1])) & (mhPvals <= 0.05))[0] if pvalBool else np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome <= startEnd[1])) & (zScores >= 20))[0]
        line3Indices = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome <= startEnd[1])) & (mhPvals <= 0.01))[0] if pvalBool else np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome <= startEnd[1])) & (zScores >= 30))[0]

    colorArr = stateColorList[maxDiffArr[line1Indices].astype(int) - 1]
    opacityArr = np.array((np.abs(distanceArrReal[line1Indices]) /
                           np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[line1Indices]), 1)
    rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
    sizeArr = np.abs(distanceArrReal[line1Indices]) / np.amax(np.abs(distanceArrReal)) * 100

    plt.scatter(line1Indices, distanceArrReal[line1Indices], s=sizeArr, color=rgbaColorArr, marker=".",
                edgecolors='none', rasterized=True)

    if len(line3Indices) > 0:
        point1Line = np.min(np.abs(distanceArrReal[line1Indices]))
        point05Line = np.min(np.abs(distanceArrReal[line2Indices]))
        point01Line = np.min(np.abs(distanceArrReal[line3Indices]))
        plt.axhspan(point1Line, point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(point05Line, point01Line, facecolor='black', alpha=0.10)
        plt.axhspan(point01Line, ylim, facecolor='black', alpha=0.15)
        plt.axhspan(-point1Line, -point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(-point05Line, -point01Line, facecolor='black', alpha=0.10)
        plt.axhspan(-point01Line, -ylim, facecolor='black', alpha=0.15)
    elif len(line2Indices) > 0:
        point1Line = np.min(np.abs(distanceArrReal[line1Indices]))
        point05Line = np.min(np.abs(distanceArrReal[line2Indices]))
        plt.axhspan(point1Line, point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(point05Line, ylim, facecolor='black', alpha=0.10)
        plt.axhspan(-point1Line, -point05Line, facecolor='black', alpha=0.05)
        plt.axhspan(-point05Line, -ylim, facecolor='black', alpha=0.10)
    elif len(line1Indices) > 0:
        point1Line = np.min(np.abs(distanceArrReal[line1Indices]))
        plt.axhspan(point1Line, ylim, facecolor='black', alpha=0.05)
        plt.axhspan(-point1Line, -ylim, facecolor='black', alpha=0.05)

    figPath = manhattanDirPath / "manhattan_plot_chr{}.pdf".format(chromosome)
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    fig.clear()
    plt.close(fig)


def pvalAxisScaling(ylim, beta, loc, scale):
    """
    Generates list containing proper tick marks for the manhattan plots

    Input:
    ylim -- The y limit of the manhattan plot
    beta -- gennorm fit parameter
    loc -- gennorm fit parameter
    scale -- gennorm fit parameter
    """
    yticks = []
    ytickLabels = ["$10^{%d}$" % i for i in range(-16, -3)] + ["$1$"] + ["$10^{-%d}$" % i for i in range(4, 17)]

    for i in range(-16, -3):
        yticks.append(-st.gennorm.isf(10**i/2, beta, loc=loc, scale=scale))
        yticks.append(st.gennorm.isf(10**i/2, beta, loc=loc, scale=scale))
    yticks.append(0)
    yticks.sort()

    yticksFinal = []
    ytickLabelsFinal = []

    for i in range(len(yticks)):
        if yticks[i] >= -ylim and yticks[i] <= ylim:
            yticksFinal.append(float(yticks[i]))
            ytickLabelsFinal.append(ytickLabels[i])

    return (yticksFinal, ytickLabelsFinal)

def zScoreAxisScaling(ylim, mean, stanDev):
    """
    Generates list containing proper tick marks for the manhattan plots according to z-score

    Input:
    ylim -- The y limit of the manhattan plot
    mean -- the mean of the scores
    stanDev -- the standard deviation of the scores
    """

    yticks = []
    ytickLabels = ["{0:.1f}".format(i) for i in np.linspace(-ylim, ylim, 11)]

    for i in np.linspace(-ylim, ylim, 11):
        yticks.append(i * stanDev + mean)

    yticksFinal = []
    ytickLabelsFinal = []

    for i in range(len(yticks)):
        if yticks[i] >= -ylim and yticks[i] <= ylim:
            yticksFinal.append(float(yticks[i]))
            ytickLabelsFinal.append(ytickLabels[i])

    return (yticksFinal, ytickLabelsFinal)


if __name__ == "__main__":
    main(argv[1], argv[2], argv[3], argv[4], argv[5], int(argv[6]), strToBool(argv[7]), strToBool(argv[8]), int(argv[9]),
         int(argv[10]), argv[11], int(argv[12]), strToBool(argv[13]))
