from sys import argv
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import scipy.stats as st
import warnings
from time import time
import gzip
from multiprocessing import cpu_count, Pool
from contextlib import closing
from itertools import repeat
from os import remove
from helpers import strToBool, getStateNames, getStateColorsRGB, getNumStates
import pyranges as pr


def main(group1Name, group2Name, stateInfo, outputDir, fileTag, numProcesses, diagnosticBool, numTrials, samplingSize,
         expFreqPath, verbose):
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
    diagnosticBool -- Boolean which if True tells us to generate diagnostic plots of the gennorm fit on the null data and 
                      comparisons between the null and real data
    numTrials -- The number of gennorm fits to do
    samplingSize -- The amount of null data to fit
    verbose -- Boolean which if True, causes much more detailed prints
    """
    tTotal = time()

    outputDirPath = Path(outputDir)
    np.random.seed(7032016)

    # Plotting setting
    plt.rcParams['agg.path.chunksize'] = 10000
    
    numStates = getNumStates(stateInfo)
    stateColorList = getStateColorsRGB(stateInfo)
    stateNameList = getStateNames(stateInfo)

    # If user doesn't want to choose number of cores use as many as available
    if numProcesses == 0:
        numProcesses = cpu_count()

    # Read in observation files
    if verbose: print("\nReading in observation files...", flush=True); tRead = time()
    else: print("    Reading in files\t", end="", flush=True)
    locationArr, distanceArrReal, distanceArrNull, maxDiffArr, diffArr, quiescenceArr = readInData(outputDirPath, numProcesses, numStates)
    if verbose: print("    Time:", time() - tRead, flush=True)
    else: print("\t[Done]", flush=True)

    # Fitting a gennorm distribution to the distances
    if verbose: print("Fitting gennorm distribution to distances...", flush=True); tFit = time()
    else: print("    Fitting distances\t", end="", flush=True)
    params, dataReal, dataNull = fitDistances(distanceArrReal, distanceArrNull, quiescenceArr, diffArr, numStates, numProcesses,
        numTrials, samplingSize)
    if verbose: print("    Time:", time() - tFit, flush=True)
    else: print("\t[Done]", flush=True)

    # Splitting the params up
    beta, loc, scale = params[:-2], params[-2], params[-1]

    # Creating Diagnostic Figures
    if diagnosticBool:
        if verbose: print("Creating diagnostic figures...", flush=True); tDiagnostic = time()
        else: print("    Diagnostic figures\t", end="", flush=True)
        createDiagnosticFigures(dataReal, dataNull, distanceArrReal, distanceArrNull, beta, loc, scale, outputDirPath, fileTag)
        if verbose: print("    Time:", time() - tDiagnostic, flush=True)
        else: print("\t[Done]", flush=True)

    # Calculating PValues
    if verbose: print("Calculating P-Values...", flush=True); tPVal = time()
    else: print("    Calculating p-vals\t", end="", flush=True)
    pvals = calculatePVals(distanceArrReal, beta, loc, scale)
    if verbose: print("    Time:", time() - tPVal, flush=True)
    else: print("\t[Done]", flush=True)

    # Create an output file which summarizes the results
    if verbose: print("Writing metrics file...", flush=True); tMetrics = time()
    else: print("    Writing metrics\t", end="", flush=True)
    writeMetrics(locationArr, maxDiffArr, distanceArrReal, pvals, outputDirPath, fileTag)
    if verbose: print("    Time:", time() - tMetrics, flush=True)
    else: print("\t[Done]", flush=True)

    # Determine Significance Threshold (based on n*)
    genomeAutoCorrelation = 0.987
    nStar = len(distanceArrReal) * ((1 - genomeAutoCorrelation) / (1 + genomeAutoCorrelation))
    significanceThreshold = .1 / nStar

    # Create txt file of top 1000 loci with adjacent merged
    if verbose: print("Creating .txt file of top loci...", flush=True); t1000 = time()
    else: print("    Greatest hits txt\t", end="", flush=True)
    roiPath = outputDirPath / "greatestHits_{}.txt".format(fileTag)
    createTopScoresTxt(roiPath, locationArr, distanceArrReal, maxDiffArr, stateNameList, pvals, nStar, False)
    if verbose: print("    Time:", time() - t1000, flush=True)
    else: print("\t[Done]", flush=True)

    # Create txt file of significant loci
    if verbose: print("Creating .txt file of significant loci...", flush=True); tSig = time()
    else: print("    Significant loci txt\t", end="", flush=True)
    roiPath = outputDirPath / "signficantLoci_{}.txt".format(fileTag)
    createTopScoresTxt(roiPath, locationArr, distanceArrReal, maxDiffArr, stateNameList, pvals, nStar, True)
    if verbose: print("    Time:", time() - tSig, flush=True)
    else: print("\t[Done]", flush=True)

    # Create Genome Manhattan Plot
    if verbose: print("Creating Genome-Wide Manhattan Plot", flush=True); tGManhattan = time()
    else: print("    Genome-wide Manhattan\t", end="", flush=True)
    createGenomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArr, beta, loc, scale,
        significanceThreshold, pvals, stateColorList, outputDirPath, fileTag)
    if verbose: print("    Time:", time() - tGManhattan, flush=True)
    else: print("\t[Done]", flush=True)
    
    # Create Chromosome Manhattan Plot
    if verbose: print("Creating Individual Chromosome Manhattan Plots", flush=True); tCManhattan = time()
    else: print("    Chromosome Manhattan\t", end="", flush=True)
    createChromosomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArr, params,
        significanceThreshold, pvals, stateColorList, outputDirPath, fileTag, numProcesses)
    if verbose: print("    Time:", time() - tCManhattan, flush=True)
    else: print("\t[Done]", flush=True)

    # Removing the expected frequency array
    remove(Path(expFreqPath))

    if verbose: print("Total Time:", time() - tTotal, flush=True)


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
    distanceArrNull -- Numpy array containing binwise signed squared euclidean distances of the null scores of the two groups
    maxDiffArr -- Numpy array containing the state which had the absolute distance in each bin
    diffArr -- Numpy array containing the raw binwise differences between the two groups
    """
    # For keeping the data arrays organized correctly
    realNames = ["chr", "binStart", "binEnd"] + ["s{}".format(i) for i in range(1, numStates + 1)]

    # Data frame to dump inputed data into
    diffDF = pd.DataFrame(columns=realNames)

    # Multiprocess the reading
    with closing(Pool(numProcesses)) as pool:
        results = pool.starmap(readTableMulti, zip(outputDirPath.glob("pairwiseDelta_*.txt.gz"),
            outputDirPath.glob("temp_nullDistances_*.npz"), outputDirPath.glob("temp_quiescence_*.npz"), repeat(realNames)))
    pool.join()

    # Concatenating all chunks to the real differences dataframe
    for diffDFChunk, _, _ in results:
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

    # Creating quiescence array ordered by chromosome based on the read in chunks
    quiescenceChunks = list(zip(*list(zip(*results))[2]))
    index = quiescenceChunks[0].index(chrOrder[0])
    quiescenceArr = quiescenceChunks[1][index]
    for chrName in chrOrder[1:]:
        index = quiescenceChunks[0].index(chrName)
        quiescenceArr = np.concatenate((quiescenceArr, quiescenceChunks[1][index]))

    # # Cleaning up the temp files after we've read them
    # for file in outputDirPath.glob("temp_*.npz"):
    #     remove(file)

    # Calculate the distance array for the real data
    diffSign = np.sign(np.sum(diffArr, axis=1))
    distanceArrReal = np.sum(np.square(diffArr), axis=1) * diffSign

    # Calculate the maximum contributing state for each bin
    # In the case of a tie, the higher number state wins (e.g. last state wins if all states are 0)
    maxDiffArr = np.abs(np.argmax(np.abs(np.flip(diffArr, axis=1)), axis=1) - diffArr.shape[1]).astype(int)

    return locationArr, distanceArrReal, distanceArrNull, maxDiffArr, diffArr, quiescenceArr


def readTableMulti(realFile, nullFile, quiescenceFile, realNames):
    """
    Reads in the real and null scores

    Input:
    realFile -- The path to the file containing the real scores
    nullFile -- The path to the file containing the null scores
    realNames -- The names for the columns of the real dataframe

    Output:
    diffDFChunk -- Pandas dataframe containing the real scores
    (npzFile['chrName'][0], npzFile['nullDistances']) -- Tuple with the chromosome name and the null signed squared euclidean 
                                                         distances
    """
    diffDFChunk = pd.read_table(Path(realFile), header=None, sep="\t", names=realNames)
    npzFileNull = np.load(Path(nullFile))
    npzFileQuiescence = np.load(Path(quiescenceFile))

    return diffDFChunk, (npzFileNull['chrName'][0], npzFileNull['nullDistances']), (npzFileQuiescence['chrName'][0], 
                                                                                    npzFileQuiescence['quiescenceArr'])


def fitDistances(distanceArrReal, distanceArrNull, quiescenceArr, diffArr, numStates, numProcesses, numTrials, samplingSize):
    """
    Filters out quiescent bins and deploys the processes which fits the null distances. Then calculates the median fit based
    on the negative loglikelihood function

    Input:
    distanceArrReal -- Numpy array containing the distances between the real scores
    distanceArrNull -- Numpy array containing the distances between the null scores
    quiescenceArr -- Numpy array containing booleans informing us which bins to filter out
    diffArr -- Numpy array containing the raw differences between the real scores
    numStates -- The number of states in the state model
    numProcesses -- The number of cores to run on
    numTrials -- The number of fits to do
    samplingSize -- The amount of data to fit each time

    Output:
    (fitDF.iloc[medianIndex, 0], fitDF.iloc[medianIndex, 1], fitDF.iloc[medianIndex, 2]) -- Tuple with beta, loc, and scale
                                                                                            params of the median fit
    dataReal -- Pandas series containing the real distances filtered for quiescence
    dataNull -- Pandas series containing the null distances filtered for quiescence
    """
    # Filtering out quiescent values (When there are exactly zero differences between both score arrays)
    idx = np.where(quiescenceArr == False)[0]
    dataReal = pd.Series(distanceArrReal[idx])
    dataNull = pd.Series(distanceArrNull[idx])

    with closing(Pool(numProcesses)) as pool:
        results = pool.starmap(fitOnSample, zip(repeat(distanceArrNull[idx], numTrials), repeat(samplingSize, numTrials)))
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
    return (fitDF.iloc[medianIndex, 0], fitDF.iloc[medianIndex, 1], fitDF.iloc[medianIndex, 2]), dataReal, dataNull


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
        np.random.seed() # On linux, multiprocessing inherits the master seed and doesn't generate fully random numbers
        sampleData = pd.Series(np.random.choice(distanceArrNull, size=samplingSize, replace=False))

    # ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Fit the data
        params = st.gennorm.fit(sampleData)

        # Calculate SSE and MLE
        nnlf = st.gennorm.nnlf(params, pd.Series(distanceArrNull))

    return params, nnlf


def createDiagnosticFigures(dataReal, dataNull, distanceArrReal, distanceArrNull, beta, loc, scale, outputDirPath, fileTag):
    """
    Generate diagnostic plots of the gennorm fit on the null data and comparisons between the null and real data

    Input:
    dataReal -- Pandas series containing the real distances filtered for quiescence
    dataNull -- Pandas series containing the null distances filtered for quiescence
    distanceArrReal -- Numpy array containing the real distances
    distanceArrNull -- Numpy array containing the null distances
    beta -- gennorm fit parameter
    loc -- gennorm fit parameter
    scale -- gennorm fit parameter
    outputDirPath -- The path to the epilogos output directory
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    """
    diagnosticDirPath = outputDirPath / "diagnosticFigures_{}".format(fileTag)
    if not diagnosticDirPath.exists():
        diagnosticDirPath.mkdir(parents=True)

    # Real Data Histogram vs. Null Data Histogram (Range=(-1, 1))
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    dataReal.plot(kind='hist', bins=200, range=(-1, 1), density=True, alpha=0.5, label='Non-Random Distances', legend=True,
        ax=ax, rasterized=True)
    dataNull.plot(kind='hist', bins=200, range=(-1, 1), density=True, alpha=0.5, label='Random Distances', legend=True,
        ax=ax, rasterized=True)
    plt.title("Real Data vs. Null Data (range=(-1, 1))")
    figPath = diagnosticDirPath / "real_vs_null_histogram_n1to1.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    plt.close(fig)

    # Real Data Histogram vs. Null Data Histogram (Range=(-max(abs), max(abs)))
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    rangeLim = np.amax(np.abs(dataReal))
    dataReal.plot(kind='hist', bins=200, range=(-rangeLim, rangeLim), density=True, alpha=0.5, label='Non-Random Distances',
        legend=True, ax=ax, rasterized=True)
    dataNull.plot(kind='hist', bins=200, range=(-rangeLim, rangeLim), density=True, alpha=0.5, label='Random Distances',
        legend=True, ax=ax, rasterized=True)
    plt.title("Real Data vs. Null Data (range=(-max(abs), max(abs)))")
    figPath = diagnosticDirPath / "real_vs_null_histogram_minToMax.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    plt.close(fig)

    # Real vs Null distance scatter plot
    fig = plt.figure(figsize=(12,12))
    plt.scatter(distanceArrReal, distanceArrNull, color='r', rasterized=True)
    plt.xlim(-rangeLim, rangeLim)
    plt.ylim(-rangeLim, rangeLim)
    plt.xlabel("Real Distances")
    plt.ylabel("Null Distances")
    plt.title("Real Distances vs Null Distances")
    figPath = diagnosticDirPath / "real_vs_null_scatter.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    plt.close(fig)

    # Fit on data (range=(min, max))
    y, x = np.histogram(dataNull, bins=20000, range=(np.amin(distanceArrNull), np.amax(distanceArrNull)), density=True,
        rasterized=True);
    x = (x + np.roll(x, -1))[:-1] / 2.0
    fig = plt.figure(figsize=(12,8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta,loc,scale), legend=True,
        rasterized=True)
    dataNull.plot(kind='hist', bins=20000, range=(np.amin(distanceArrNull), np.amax(distanceArrNull)), density=True,
        alpha=0.5, label='Data', legend=True, ax=ax, rasterized=True)
    plt.title("Gennorm on data (range=(min,max))")
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_minToMax.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    plt.close(fig)

    # Fit on data (range=(-1, 1))
    y, x = np.histogram(dataNull, bins=20000, range=(-1, 1), density=True, rasterized=True);
    x = (x + np.roll(x, -1))[:-1] / 2.0
    fig = plt.figure(figsize=(12,8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta,loc,scale), legend=True,
        rasterized=True)
    dataNull.plot(kind='hist', bins=20000, range=(-1, 1), density=True, alpha=0.5, label='Data', legend=True, ax=ax,
        rasterized=True)
    plt.title("Gennorm on data (range=(-1,1))")
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_n1to1.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    plt.close(fig)

    # Fit on data (range=(-0.1, 0.1))
    y, x = np.histogram(dataNull, bins=20000, range=(-1, 1), density=True, rasterized=True);
    x = (x + np.roll(x, -1))[:-1] / 2.0
    fig = plt.figure(figsize=(12,8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta,loc,scale), legend=True,
        rasterized=True)
    dataNull.plot(kind='hist', bins=20000, range=(-1, 1), density=True, alpha=0.5, label='Data', legend=True, ax=ax,
        rasterized=True)
    plt.title("Gennorm on data (range=(-0.1,0.1))")
    plt.xlim(-.1, .1)
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_0p1to0p1.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
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


def createGenomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArr, beta, loc, scale,
    significanceThreshold, pvals, stateColorList, outputDirPath, fileTag):
    """
    Creates a manhattan plot based on the distances between the two groups for the entire genome

    Input:
    group1Name -- The name of the first epilogos group
    group2Name -- The name of the second epilogos group
    locationArr -- Numpy array containing the genomic locations of all the bins
    distanceArrReal -- Numpy array containing the real distances
    maxDiffArr -- Numpy array containing the states which had the largest difference between the two groups in each bin
    beta -- gennorm fit parameter
    loc -- gennorm fit parameter
    scale -- gennorm fit parameter
    significanceThreshold -- The value below which pvalues are considered significant
    pvals -- Numpy array containing the pvalues of all the distances
    stateColorList -- Numpy array containing the colors of each of the states in the state model
    outputDirPath -- The path of the output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    """
    manhattanDirPath = outputDirPath / "manhattanPlots_{}".format(fileTag)
    if not manhattanDirPath.exists():
        manhattanDirPath.mkdir(parents=True)

    logSignificanceThreshold = -np.log10(significanceThreshold)

    fig = plt.figure(figsize=(16,9))
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
    xticks = np.where(locationArr[:, 1] == "0")[0]
    plt.xticks(ticks=xticks, labels=list(map(lambda x: x.split("chr")[-1], list(locationArr[:, 0][xticks]))))

    plt.margins(x=0)
    ylim = np.amax(np.abs(distanceArrReal)) * 1.1
    ax.set_ylim(-ylim, ylim)
    yticks, ytickLabels = pvalAxisScaling(ylim, beta, loc, scale)

    ax.set_yticks(yticks)
    ax.set_yticklabels([str(np.abs(np.round(val, 1))) for val in yticks])

    axR = ax.twinx()
    axR.set_ylabel("P-Value")
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
    pvalsGraph = -np.log10(pvals.astype(float)) * np.sign(distanceArrReal)

    for i in range(len(xticks)):
        if i == len(xticks)-1:
            points = np.where((locationOnGenome >= xticks[i]) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points],
                s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".",
                    alpha=0.1, edgecolors='none', rasterized=True)
        elif i%2 == 0:
            points = np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1])
                & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points],
                s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".",
                    alpha=0.1, edgecolors='none', rasterized=True)
        else:
            points = np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1])
                & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points],
                s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="black", marker=".",
                    alpha=0.1, edgecolors='none', rasterized=True)
            
    opaqueSigIndices = np.where(np.abs(pvalsGraph) >= logSignificanceThreshold)[0]

    colorArr=stateColorList[maxDiffArr[opaqueSigIndices].astype(int) - 1]
    opacityArr=np.array((np.abs(distanceArrReal[opaqueSigIndices]) /
        np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[opaqueSigIndices]), 1)
    rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
    sizeArr = np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)) * 100

    plt.scatter(opaqueSigIndices, distanceArrReal[opaqueSigIndices], s=sizeArr, color=rgbaColorArr, marker=".",
        edgecolors='none', rasterized=True)
    ax.axhline(st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
    ax.axhline(-st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")

    figPath = manhattanDirPath / "manhattan_plot_genome.pdf"
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
    plt.close(fig)


def _init(group1Name_, group2Name_, locationArr_, distanceArrReal_, maxDiffArr_, params, significanceThreshold_, pvalsGraph_,
    stateColorList_, manhattanDirPath_):
    """
    Initializes global variables for multiprocessing in the single epilogos case

    Input:
    group1Name_ -- The name of the first epilogos group
    group2Name_ -- The name of the second epilogos group
    locationArr_ -- Numpy array containing the genomic locations of all the bins
    distanceArrReal_ -- Numpy array containing the real distances
    maxDiffArr_ -- Numpy array containing the states which had the largest difference between the two groups in each bin
    params -- gennorm fit parameters
    significanceThreshold_ -- The value below which pvalues are considered significant
    pvalsGraph_ -- Numpy array containing the signed (according to distance) -log10 pvalues of all the distances
    stateColorList_ -- Numpy array containing the colors of each of the states in the state model
    manhattanDirPath_ -- The path to directory to put manhattan plots
    """
    global group1Name
    global group2Name
    global locationArr
    global distanceArrReal
    global maxDiffArr
    global beta
    global loc
    global scale
    global significanceThreshold
    global pvalsGraph
    global stateColorList
    global manhattanDirPath

    group1Name = group1Name_
    group2Name = group2Name_
    locationArr = locationArr_
    distanceArrReal = distanceArrReal_
    maxDiffArr = maxDiffArr_
    beta, loc, scale = params[:-2], params[-2], params[-1]
    significanceThreshold = significanceThreshold_
    pvalsGraph = pvalsGraph_
    stateColorList = stateColorList_
    manhattanDirPath = manhattanDirPath_


def createChromosomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArr, params,
    significanceThreshold, pvals, stateColorList, outputDirPath, fileTag, numProcesses):
    """
    Creates a manhattan plot based on the distances between the two groups for the each chromosome

    Input:
    group1Name -- The name of the first epilogos group
    group2Name -- The name of the second epilogos group
    locationArr -- Numpy array containing the genomic locations of all the bins
    distanceArrReal -- Numpy array containing the real distances
    maxDiffArr -- Numpy array containing the states which had the largest difference between the two groups in each bin
    params -- gennorm fit parameters
    significanceThreshold -- The value below which pvalues are considered significant
    pvals -- Numpy array containing the pvalues of all the distances
    stateColorList -- Numpy array containing the colors of each of the states in the state model
    outputDirPath -- The path of the output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    numProcesses -- The number of cores to run on
    """
    manhattanDirPath = outputDirPath / "manhattanPlots_{}".format(fileTag)
    if not manhattanDirPath.exists():
        manhattanDirPath.mkdir(parents=True)

    pvalsGraph = -np.log10(pvals.astype(float)) * np.sign(distanceArrReal)

    xticks = np.where(locationArr[:, 1] == "0")[0]
    startEnd = []
    for i in range(len(xticks)):
        if not i == len(xticks) - 1:
            startEnd.append((xticks[i], xticks[i+1]))
        else:
            startEnd.append((xticks[i], -1))

    chrOrder = list(map(lambda x: x.split("chr")[-1], list(locationArr[:, 0][xticks])))

    # Multiprocess the reading
    with closing(Pool(numProcesses, initializer=_init, initargs=(group1Name, group2Name, locationArr, distanceArrReal,
        maxDiffArr, params, significanceThreshold, pvalsGraph, stateColorList, manhattanDirPath))) as pool:
        pool.starmap(graphChromosomeManhattan, zip(chrOrder, startEnd))
    pool.join()


def graphChromosomeManhattan(chromosome, startEnd):
    """
    Generates the manhattan plots for each chromosome based on the distances between the two groups.
    Note there are global variables which are used across all chromosomes (see _init())

    Input:
    chromosome -- The chromosome which we are plotting
    startEnd -- The start and end indices for the chromosome on all the global numpy arrays
    """
    logSignificanceThreshold = -np.log10(significanceThreshold)    

    fig = plt.figure(figsize=(16,9))
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
    yticks, ytickLabels = pvalAxisScaling(ylim, beta, loc, scale)

    ax.set_yticks(yticks)
    ax.set_yticklabels([str(np.abs(np.round(val, 1))) for val in yticks])

    axR = ax.twinx()
    axR.set_ylabel("P-Value")
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
        realxticks = np.where((locationOnGenome >= startEnd[0]) & (locationArr[:, 1].astype(int)%10000000 == 0))[0]
        plt.xticks(ticks=realxticks, labels=[str(int(int(locationArr[tick, 1])/1000000)) for tick in realxticks])

        points = np.where((locationOnGenome >= startEnd[0]) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
        plt.scatter(locationOnGenome[points], distanceArrReal[points],
            s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".",
                alpha=0.1, edgecolors='none', rasterized=True)

        opaqueSigIndices = np.where((locationOnGenome >= startEnd[0]) & (np.abs(pvalsGraph) >= logSignificanceThreshold))[0]
    else:
        realxticks = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome < startEnd[1]))
            & (locationArr[:, 1].astype(int)%10000000 == 0))[0]
        plt.xticks(ticks=realxticks, labels=[str(int(int(locationArr[tick, 1])/1000000)) for tick in realxticks])

        points = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome < startEnd[1]))
            & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
        plt.scatter(locationOnGenome[points], distanceArrReal[points],
            s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".",
                alpha=0.1, edgecolors='none', rasterized=True)

        opaqueSigIndices = np.where(((locationOnGenome >= startEnd[0]) & (locationOnGenome < startEnd[1]))
            & (np.abs(pvalsGraph) >= logSignificanceThreshold))[0]

    colorArr=stateColorList[maxDiffArr[opaqueSigIndices].astype(int) - 1]
    opacityArr=np.array((np.abs(distanceArrReal[opaqueSigIndices]) /
        np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[opaqueSigIndices]), 1)
    rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
    sizeArr = np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)) * 100

    plt.scatter(opaqueSigIndices, distanceArrReal[opaqueSigIndices], s=sizeArr, color=rgbaColorArr, marker=".",
        edgecolors='none', rasterized=True)
    ax.axhline(st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
    ax.axhline(-st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")

    figPath = manhattanDirPath / "manhattan_plot_chr{}.pdf".format(chromosome)
    fig.savefig(figPath, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
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
    

def writeMetrics(locationArr, maxDiffArr, distanceArrReal, pvals, outputDirPath, fileTag):
    """
    Writes metrics file to disk. Metrics file contains the following columns in order:
    chromosome, bin start, bin end, state with highest difference, signed squared euclidean distance, pvalue of distance

    Input:
    locationArr -- Numpy array containing the genomic locations of all the bins
    maxDiffArr -- Numpy array containing the states which had the largest difference between the two groups in each bin
    distanceArrReal -- Numpy array containing the real distances
    pvals -- Numpy array containing the pvalues of all the distances
    outputDirPath -- The path of the output directory for epilogos
    fileTag -- A string which helps ensure outputed files are named similarly within an epilogos run
    """
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)

    metricsTxtPath = outputDirPath / "pairwiseMetrics_{}.txt.gz".format(fileTag)
    metricsTxt = gzip.open(metricsTxtPath, "wt")

    # Creating a string to write out the raw differences (faster than np.savetxt)
    metricsTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3:.5e}\n"
    metricsStr = "".join(metricsTemplate.format(locationArr[i], maxDiffArr[i], distanceArrReal[i], pvals[i])
        for i in range(len(distanceArrReal)))

    metricsTxt.write(metricsStr)
    metricsTxt.close()


def createTopScoresTxt(filePath, locationArr, distanceArr, maxDiffArr, nameArr, pvals, nStar, onlySignificant):
    """
    Finds the either the 1000 largest distance bins and merges adjacent bins or finds all significant loci and does not merge.
    Then it outputs a txt containing these highest distances regions and some information about each (chromosome, bin start, 
    bin end, state name, absolute value of distance, sign of score, pvalue of distance, stars repersenting significance)

    Input:
    filePath -- The path of the file to write to
    locationArr -- Numpy array containing the genomic locations of all the bins
    distanceArr -- Numpy array containing the distance between the pairwise groups
    maxScoreArr -- Numpy array containing the states which had the largest difference in each bin
    nameArr -- Numpy array containing the names of all the states
    pvals -- Numpy array containing the pvalues of all the distances between the pairwise groups
    nStar -- Number of bins scored adjusted for autocorrelation
    onlySignificant -- Boolean telling us whether to use significant loci or 1000 largest distance bins
    """
    
    significantAt1 = .1 / nStar
    significantAt05 = .05 / nStar
    significantAt01 = .01 / nStar

    with open(filePath, 'w') as f:
        # Pick values above significance threshold and then sort
        indices = np.where(pvals <= significantAt1)[0][(-np.abs(distanceArr[np.where(pvals <= significantAt1)[0]])).argsort()]
        # Make sure that there are at least 1000 values if creating greatestHits.txt
        if not onlySignificant and len(indices) < 1000:
            indices = (-np.abs(distanceArr)).argsort()[:1000]

        locations = pd.DataFrame(np.concatenate((locationArr[indices], distanceArr[indices].reshape(len(indices), 1),
            maxDiffArr[indices].reshape(len(indices), 1), pvals[indices].reshape(len(indices), 1)), axis=1), 
            columns=["chr", "binStart", "binEnd", "distance", "maxDiffLoc", "pval"])\
                .astype({"chr": str, "binStart": np.int32, "binEnd": np.int32, "distance": np.float32, "maxDiffLoc": np.int32, "pval": np.float32})

        # print("premerge lenght", locations.shape[0])

        # locations = pd.DataFrame(np.concatenate((locationArr[indices], distanceArr[indices].reshape(len(indices), 1),
        #     maxDiffArr[indices].reshape(len(indices), 1), pvals[indices].reshape(len(indices), 1)), axis=1), 
        #     columns=["Chromosome", "Start", "End", "distance", "maxDiffLoc", "pval"])\
        #         .astype({"Chromosome": str, "Start": np.int32, "End": np.int32, "distance": np.float32, "maxDiffLoc": np.int32, "pval": np.float32})

        # Figuring out chromosome order
        chromosomes = locations['chr'].unique()
        # chromosomes = locations['Chromosome'].unique()
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
        locations["chr"] = pd.Categorical(locations["chr"], categories=chrOrder, ordered=True)
        locations.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)

        # locations["chr"] = pd.Categorical(locations["Chromosome"], categories=chrOrder, ordered=True)
        # locations.sort_values(by=["Chromosome", "Start", "End"], inplace=True)


        # Iterate until all is merged, but only for the general case
        if not onlySignificant:
            locations = mergeAdjacent(locations)

        # print("postmerge lenght", locations.shape[0])


        # if not onlySignificant:
        #     tMerge = time()
        #     locations = pr.PyRanges(locations)
        #     locations = mergeAdjacent(locations)
        #     print("TIME TO MERGE:", time() - tMerge)

        locations = locations.iloc[-locations.iloc[:, 3].abs().argsort()]

        print("end locations shape", locations.shape)
        print("end locations", locations.head())

        # locations.sort_values(by=["distance", "chr", "binStart", "binEnd"], inplace=True, ascending=False)

        # Locations get 3 stars if they are significant at .01, 2 stars at .05, 1 star at .1, and a period if not significant
        stars = np.array(["***" if float(locations.iloc[i, 5]) <= significantAt01 else
            ("**" if float(locations.iloc[i, 5]) <= significantAt05 else
                ("*" if float(locations.iloc[i, 5]) <= significantAt1 else "."))
                    for i in range(locations.shape[0])]).reshape(locations.shape[0], 1)

        print("stars shape", stars.shape)


        # Write all the locations to the file for significantLoci.txt
        # Write only top 100 loci to file for greatestHits.txt
        outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2:.5f}\t{3}\t{4:.5e}\t{5}\n"
        outString = "".join(outTemplate.format(locations.iloc[i], nameArr[int(float(locations.iloc[i, 4])) - 1],
            abs(float(locations.iloc[i, 3])), findSign(float(locations.iloc[i, 3])), float(locations.iloc[i, 5]), stars[i, 0])
                for i in range(locations.shape[0] if onlySignificant else min((locations.shape[0], 100))))
        f.write(outString)


def mergeAdjacent(originalLocations):
    """
    Takes a pandas dataframe sorted by genomic location and merges all adjacent loci

    Input:
    locationArr -- pandas dataframe containing genomic loci in the first 3 columns

    Output:
    dataframe with merged loci
    """
    i = 0
    mergedLocations = []
    while i < len(originalLocations) - 1:
        j = 1
        while i + j < len(originalLocations) and originalLocations.iloc[i, 1] == originalLocations.iloc[i+j, 1] - 200 * j \
            and originalLocations.iloc[i, 0] == originalLocations.iloc[i+j, 0]:
            j += 1
        maxDistIndex = i + originalLocations.iloc[i:i+j, 3].argmax()
        mergedLocations.append([originalLocations.iloc[i, 0], originalLocations.iloc[i, 1], originalLocations.iloc[i+j-1, 2]] 
                               + list(originalLocations.iloc[maxDistIndex, 3:]))
        i += j
    return pd.DataFrame(mergedLocations)

# def mergeAdjacent(originalLocations):
#     merged_data = originalLocations.merge()
#     join_merged_to_all = merged_data.join(originalLocations)

#     # res is a pyranges object
#     res = join_merged_to_all.apply(max_scoring_element)
#     # res.df is a pandas dataframe
#     return res.df

# def max_scoring_element(df):
#     return df \
#         .sort_values('distance', ascending=False) \
#         .drop_duplicates(['Chromosome', 'Start', 'End', 'maxDiffLoc', 'pval'], keep='first') \
#         .sort_index() \
#         .reset_index(drop=True)



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


if __name__ == "__main__":
    main(argv[1], argv[2], argv[3], argv[4], argv[5], int(argv[6]), strToBool(argv[7]), int(argv[8]), int(argv[9]),
         argv[10], strToBool(argv[11]))
