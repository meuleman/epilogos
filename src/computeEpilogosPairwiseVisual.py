import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
import math
import scipy.stats as st
import statsmodels as sm
import warnings
import time

def main(group1Name, group2Name, numStates, outputDir):
    tTotal = time.time()

    outputDirPath = Path(outputDir)

    # Plotting setting
    plt.rcParams['agg.path.chunksize'] = 10000
    
    if numStates == 18:
        stateColorList = np.array([(1.0, 0.0, 0.0), (1.0, 0.2706, 0.0), (1.0, 0.2706, 0.0), (1.0, 0.2706, 0.0), (0.0, 0.502, 0.0), (0.0, 0.3922, 0.0), (0.7608, 0.8824, 0.0196), (0.7608, 0.8824, 0.0196), (1.0, 0.7647, 0.302), (1.0, 0.7647, 0.302), (1.0, 1.0, 0.0), (0.4, 0.8039, 0.6667), (0.5412, 0.5686, 0.8157), (0.8039, 0.3608, 0.3608), (0.7412, 0.7176, 0.4196), (0.502, 0.502, 0.502), (0.7529, 0.7529, 0.7529), (1.0, 1.0, 1.0)])
        stateNameList = np.array(["TssA", "TssFlnk", "TssFlnkU", "TssFlnkD", "Tx", "TxWk", "EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "ZNF/Rpts", "Het", "TssBiv", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"])
    elif numStates == 15:
        stateColorList = np.array([(1.0, 0.0, 0.0), (1.0, 0.27059, 0.0), (0.19608, 0.80392, 0.19608), (0.0, 0.50196, 0.0), (0.0, 0.39216, 0.0), (0.76078, 0.88235, 0.01961), (1.0, 1.0, 0.0), (0.4, 0.80392, 0.66667), (0.54118, 0.56863, 0.81569), (0.80392, 0.36078, 0.36078), (0.91373, 0.58824, 0.47843), (0.74118, 0.71765, 0.41961), (0.50196, 0.50196, 0.50196), (0.75294, 0.75294, 0.75294), (1.0, 1.0, 1.0)])
        stateNameList = np.array(["TssA", "TssAFlnk", "TxFlnk", "Tx", "TxWk", "EnhG", "Enh", "ZNF/Rpts", "Het", "TssBiv", "BivFlnk", "EnhBiv", "ReprPC", "ReprPCWk", "Quies"])
    else:
        print("State model not supported for plotting")
        return

    # Read in observation files
    print("\nReading in observation files...")
    tRead = time.time()
    locationArr, distanceArrReal, distanceArrNull, maxDiffArrReal, maxDiffArrNull = readInData(outputDirPath)
    print("    Time:", time.time() - tRead)

    # Fitting a gennorm distribution to the distances
    print("Fitting gennorm distribution to distances...")
    tFit = time.time()
    params, dataReal, dataNull = fitDistances(distanceArrReal, distanceArrNull)
    print("    Time:", time.time() - tFit)

    # Splitting the params up
    beta, loc, scale = params[:-2], params[-2], params[-1]

    # # Creating Diagnostic Figures
    # print("Creating diagnostic figures...")
    # tDiagnostic = time.time()
    # createDiagnosticFigures(dataReal, dataNull, distanceArrReal, distanceArrNull, beta, loc, scale, outputDirPath)
    # print("    Time:", time.time() - tDiagnostic)

    # Calculating PValues
    print("Calculating P-Values...")
    tPVal = time.time()
    pvals = calculatePVals(distanceArrReal, beta, loc, scale)
    print("    Time:", time.time() - tPVal)

    # Determine Significance Threshold (based on n*)
    genomeAutoCorrelation = 0.987
    nStar = len(distanceArrReal) * ((1 - genomeAutoCorrelation) / (1 + genomeAutoCorrelation))
    significanceThreshold = .1 / nStar

    # Create Genome Manhattan Plot
    print("Creating Genome-Wide Manhattan Plot")
    tGManhattan = time.time()
    createGenomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArrReal, beta, loc, scale, significanceThreshold, pvals, stateColorList, outputDirPath)
    print("    Time:", time.time() - tGManhattan)
    
    # Create Chromosome Manhattan Plot
    print("Creating Individual Chromosome Manhattan Plots")
    tCManhattan = time.time()
    createChromosomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArrReal, beta, loc, scale, significanceThreshold, pvals, stateColorList, outputDirPath)
    print("    Time:", time.time() - tCManhattan)

    # Create roiURL
    roiPath = Path("/home/jquon/public_html/roiUrl/largestDistanceLoci_{}_{}.bed".format(group1Name, group2Name))
    sendRoiUrl(roiPath, locationArr, distanceArrReal, maxDiffArrReal, stateNameList)
    print("https://epilogos.altius.org/?application=viewer&sampleSet=vC&mode=single&genome=hg19&model=18&complexity=KL&group=all&roiPaddingAbsolute=10000&roiURL=https%3A%2F%2Fencode%3Acollabor8%40resources.altius.org%2F~jquon%2FroiUrl/{}".format(roiPath.name))

    print("Total Time:", time.time() - tTotal)

# Helper to read in the necessary data to fit and visualize pairwise results
def readInData(outputDirPath):
    # For keeping the data arrays organized correctly
    names = ["chr", "binStart", "binEnd", "maxDiffState", "distance"]
    chrOrder = []
    for i in range(1, 23):
        chrOrder.append("chr{}".format(i))
    chrOrder.append("chrX")

    # Data frames to dump inputed data into
    observationDFReal = pd.DataFrame(columns=names)
    observationDFNull = pd.DataFrame(columns=names)
    
    # Take in all the real distances
    for file in outputDirPath.glob("pairwiseObservations_*.txt.gz"):
        observationDFChunk = pd.read_table(Path(file), header=None, sep="\s+", names=names)
        observationDFReal = pd.concat((observationDFReal, observationDFChunk), axis=0, ignore_index=True)

    # Take in all the null distances
    for file in outputDirPath.glob("pairwiseObservationsNull_*.txt.gz"):
        observationDFChunk = pd.read_table(Path(file), header=None, sep="\s+", names=names)
        observationDFNull = pd.concat((observationDFNull, observationDFChunk), axis=0, ignore_index=True)

    # Sorting the dataframes by chromosomal location
    observationDFReal["chr"] = pd.Categorical(observationDFReal["chr"], categories=chrOrder, ordered=True)
    observationDFNull["chr"] = pd.Categorical(observationDFNull["chr"], categories=chrOrder, ordered=True)
    observationDFReal.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)
    observationDFNull.sort_values(by=["chr", "binStart", "binEnd"], inplace=True)

    # Split the locations, distances, and maximum difference state
    locationArr     = observationDFReal.iloc[:,0:3].to_numpy(dtype=str)
    distanceArrReal = observationDFReal.iloc[:,4].to_numpy(dtype=float).flatten()
    distanceArrNull = observationDFNull.iloc[:,4].to_numpy(dtype=float).flatten()
    maxDiffArrReal  = observationDFReal.iloc[:,3].to_numpy(dtype=float).flatten()
    maxDiffArrNull  = observationDFNull.iloc[:,3].to_numpy(dtype=float).flatten()

    return locationArr, distanceArrReal, distanceArrNull, maxDiffArrReal, maxDiffArrNull

# Helper to fit the distances
def fitDistances(distanceArrReal, distanceArrNull):
    # Filtering out quiescent values (We know that the first entry from chromosome 1 is quiescent so we use that as base)
    quiescentVal = round(distanceArrReal[0], 5)
    idx = [i for i in range(len(distanceArrReal)) if round(distanceArrReal[i]) != quiescentVal]
    dataReal = pd.Series(distanceArrReal[idx])
    dataNull = pd.Series(distanceArrNull[idx])

    # y, x = np.histogram(data.values, bins=100, range=(np.amin(data), np.amax(data)), density=True)
    # x = (x + np.roll(x, -1))[:-1] / 2.0

    # ignore warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # Fit the data
        params = st.gennorm.fit(dataNull)

        # # Separate parts of parameters
        # distArgs = params[:-2]
        # loc = params[-2]
        # scale = params[-1]
        # # Calculate SSE and MLE
        # pdf = st.gennorm.pdf(x, loc=loc, scale=scale, *distArgs)
        # sse = np.sum(np.power(y - pdf, 2.0))
        # mle = st.gennorm.nnlf(params, data)

    return params, dataReal, dataNull


# Helper for creating and saving diagnostic figures
def createDiagnosticFigures(dataReal, dataNull, distanceArrReal, distanceArrNull, beta, loc, scale, outputDirPath):
    diagnosticDirPath = outputDirPath / "diagnosticFigures"
    if not diagnosticDirPath.exists():
        diagnosticDirPath.mkdir(parents=True)
    
    # Real Data Histogram vs. Null Data Histogram (Range=(-1, 1))
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    dataReal.plot(kind='hist', bins=200, range=(-1, 1), density=True, alpha=0.5, label='Non-Random Distances', legend=True, ax=ax)
    dataNull.plot(kind='hist', bins=200, range=(-1, 1), density=True, alpha=0.5, label='Random Distances', legend=True, ax=ax)
    plt.title("Real Data vs. Null Data (range=(-1, 1))")
    figPath = diagnosticDirPath / "real_vs_null_histogram_n1to1.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    # Real Data Histogram vs. Null Data Histogram (Range=(-max(abs), max(abs)))
    fig = plt.figure(figsize=(16,9))
    ax = fig.add_subplot(111)
    rangeLim = np.amax(np.abs(dataReal))
    dataReal.plot(kind='hist', bins=200, range=(-rangeLim, rangeLim), density=True, alpha=0.5, label='Non-Random Distances', legend=True, ax=ax)
    dataNull.plot(kind='hist', bins=200, range=(-rangeLim, rangeLim), density=True, alpha=0.5, label='Random Distances', legend=True, ax=ax)
    plt.title("Real Data vs. Null Data (range=(-max(abs), max(abs)))")
    figPath = diagnosticDirPath / "real_vs_null_histogram_minToMax.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    # Real vs Null distance scatter plot
    fig = plt.figure(figsize=(12,12))
    print("distanceArrReal Length:", distanceArrReal.shape)
    print("distanceArrNull Length:", distanceArrNull.shape)
    plt.scatter(distanceArrReal, distanceArrNull, color='r')
    plt.xlim(-rangeLim, rangeLim)
    plt.ylim(-rangeLim, rangeLim)
    plt.xlabel("Real Distances")
    plt.ylabel("Null Distances")
    plt.title("Real Distances vs Null Distances")
    figPath = diagnosticDirPath / "real_vs_null_scatter.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    # Fit on data (range=(min, max))
    y, x = np.histogram(dataNull, bins=20000, range=(np.amin(distanceArrNull), np.amax(distanceArrNull)), density=True);
    x = (x + np.roll(x, -1))[:-1] / 2.0
    plt.figure(figsize=(12,8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta,loc,scale), legend=True)
    dataNull.plot(kind='hist', bins=20000, range=(np.amin(distanceArrNull), np.amax(distanceArrNull)), density=True, alpha=0.5, label='Data', legend=True, ax=ax)
    plt.title("Gennorm on data (range=(min,max))")
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_minToMax.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    # Fit on data (range=(-1, 1))
    y, x = np.histogram(dataNull, bins=20000, range=(-1, 1), density=True);
    x = (x + np.roll(x, -1))[:-1] / 2.0
    plt.figure(figsize=(12,8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta,loc,scale), legend=True)
    dataNull.plot(kind='hist', bins=20000, range=(-1, 1), density=True, alpha=0.5, label='Data', legend=True, ax=ax)
    plt.title("Gennorm on data (range=(-1,1))")
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_n1to1.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    # Fit on data (range=(-0.1, 0.1))
    y, x = np.histogram(dataNull, bins=20000, range=(-1, 1), density=True);
    x = (x + np.roll(x, -1))[:-1] / 2.0
    plt.figure(figsize=(12,8))
    pdf = st.gennorm.pdf(x, beta, loc=loc, scale=scale)
    ax = pd.Series(pdf, x).plot(label='gennorm(beta={}, loc={}, scale={})'.format(beta,loc,scale), legend=True)
    dataNull.plot(kind='hist', bins=20000, range=(-1, 1), density=True, alpha=0.5, label='Data', legend=True, ax=ax)
    plt.title("Gennorm on data (range=(-0.1,0.1))")
    plt.xlim(-.1, .1)
    plt.xlabel("Signed Squared Euclidean Distance")
    figPath = diagnosticDirPath / "gennorm_on_data_0p1to0p1.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)


# Helper to find the P-Values of all the distances
def calculatePVals(distanceArrReal, beta, loc, scale):
    pvalsBelowLoc = 2 * st.gennorm.cdf(distanceArrReal[np.where(distanceArrReal <= loc)[0]], beta, loc=loc, scale=scale)
    pvalsAboveLoc = 2 * (1 - st.gennorm.cdf(distanceArrReal[np.where(distanceArrReal > loc)[0]], beta, loc=loc, scale=scale))
    
    pvals = np.zeros(len(distanceArrReal))
    pvals[np.where(distanceArrReal <= loc)[0]] = pvalsBelowLoc
    pvals[np.where(distanceArrReal > loc)[0]]  = pvalsAboveLoc

    return pvals


# Helper to create a genome-wide manhattan plot
def createGenomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArrReal, beta, loc, scale, significanceThreshold, pvals, stateColorList, outputDirPath):
    manhattanDirPath = outputDirPath / "manhattanPlots"
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
    plt.xticks(ticks=xticks, labels=[i for i in range(1, 23)] + ["X"])

    plt.margins(x=0)
    ylim = np.amax(np.abs(distanceArrReal)) * 1.1
    ax.set_ylim(-ylim, ylim)
    yticks, ytickLabels = pvalAxisScaling(ylim, beta, loc, scale)

    ax.set_yticks(yticks)
    print(yticks)
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
    ax.text(0.99, 0.01, group2Name, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=15)

    locationOnGenome = np.arange(len(distanceArrReal))
    pvalsGraph = -np.log10(pvals.astype(float)) * np.sign(distanceArrReal)

    for i in range(len(xticks)):
        if i == len(xticks)-1:
            points = np.where((locationOnGenome >= xticks[i]) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points], s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".", alpha=0.1, edgecolors='none')
        elif i%2 == 0:
            points = np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1]) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points], s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".", alpha=0.1, edgecolors='none')
        else:
            points = np.where((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1]) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points], s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="black", marker=".", alpha=0.1, edgecolors='none')
            
    opaqueSigIndices = np.where(np.abs(pvalsGraph) >= logSignificanceThreshold)[0]

    colorArr=np.array(stateColorList)[maxDiffArrReal[opaqueSigIndices].astype(int) - 1]
    opacityArr=np.array((np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[opaqueSigIndices]), 1)
    rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
    sizeArr = np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)) * 100

    plt.scatter(opaqueSigIndices, distanceArrReal[opaqueSigIndices], s=sizeArr, color=rgbaColorArr, marker=".", edgecolors='none')
    ax.axhline(st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
    ax.axhline(-st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")

    figPath = manhattanDirPath / "manhattan_plot_genome.png"
    fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)


# Helper for generating individual chromosome manhattan plots
def createChromosomeManhattan(group1Name, group2Name, locationArr, distanceArrReal, maxDiffArrReal, beta, loc, scale, significanceThreshold, pvals, stateColorList, outputDirPath):
    manhattanDirPath = outputDirPath / "manhattanPlots"
    if not manhattanDirPath.exists():
        manhattanDirPath.mkdir(parents=True)

    logSignificanceThreshold = -np.log10(significanceThreshold)    
    pvalsGraph = -np.log10(pvals.astype(float)) * np.sign(distanceArrReal)
    xticks = np.where(locationArr[:, 1] == "0")[0]

    for i in range(len(xticks)):
        if i == len(xticks)-1:
            fig = plt.figure(figsize=(16,9))
            ax = fig.add_subplot(111)
            ax.set_facecolor("#FFFFFF")
            ax.set_axisbelow(True)
            ax.grid(True, axis='y', color='k', linewidth=.25, linestyle="-")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.set_ylabel("Distance")
            plt.xlabel("Location in Chromosome X (Mb)")
            plt.title("Differential epilogos between {} and {} donor biosamples (Chromosome X)".format(group1Name, group2Name))
            
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
            ax.text(0.99, 0.01, group2Name, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=15)

            locationOnGenome = np.arange(len(distanceArrReal))

            realxticks = np.where((locationOnGenome >= xticks[i]) & (locationArr[:, 1].astype(int)%10000000 == 0))[0]
            plt.xticks(ticks=realxticks, labels=[str(int(int(locationArr[tick, 1])/1000000)) for tick in realxticks])

            points = np.where((locationOnGenome >= xticks[i]) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points], s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".", alpha=0.1, edgecolors='none')

            opaqueSigIndices = np.where((locationOnGenome >= xticks[i]) & (np.abs(pvalsGraph) >= logSignificanceThreshold))[0]

            colorArr=np.array(stateColorList)[maxDiffArrReal[opaqueSigIndices].astype(int) - 1]
            opacityArr=np.array((np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[opaqueSigIndices]), 1)
            rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
            sizeArr = np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)) * 100

            plt.scatter(opaqueSigIndices, distanceArrReal[opaqueSigIndices], s=sizeArr, color=rgbaColorArr, marker=".", edgecolors='none')
            ax.axhline(st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
            ax.axhline(-st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
            
            figPath = manhattanDirPath / "manhattan_plot_chrX.png"
            fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)
        
        else:
            fig = plt.figure(figsize=(16,9))
            ax = fig.add_subplot(111)
            ax.set_facecolor("#FFFFFF")
            ax.set_axisbelow(True)
            ax.grid(True, axis='y', color='k', linewidth=.25, linestyle="-")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["bottom"].set_visible(False)
            ax.set_ylabel("Distance")
            plt.xlabel("Location in Chromosome {} (Mb)".format(i+1))
            plt.title("Differential epilogos between {} and {} donor biosamples (Chromosome {})".format(group1Name, group2Name, i+1))
            
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
            ax.text(0.99, 0.01, group2Name, verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=15)

            locationOnGenome = np.arange(len(distanceArrReal))
            
            realxticks = np.where(((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1])) & (locationArr[:, 1].astype(int)%10000000 == 0))[0]
            plt.xticks(ticks=realxticks, labels=[str(int(int(locationArr[tick, 1])/1000000)) for tick in realxticks])

            points = np.where(((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1])) & (np.abs(pvalsGraph) < logSignificanceThreshold))[0]
            plt.scatter(locationOnGenome[points], distanceArrReal[points], s=(np.abs(distanceArrReal[points]) / np.amax(np.abs(distanceArrReal)) * 100), color="gray", marker=".", alpha=0.1, edgecolors='none')

            opaqueSigIndices = np.where(((locationOnGenome >= xticks[i]) & (locationOnGenome < xticks[i+1])) & (np.abs(pvalsGraph) >= logSignificanceThreshold))[0]

            colorArr=np.array(stateColorList)[maxDiffArrReal[opaqueSigIndices].astype(int) - 1]
            opacityArr=np.array((np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)))).reshape(len(distanceArrReal[opaqueSigIndices]), 1)
            rgbaColorArr = np.concatenate((colorArr, opacityArr), axis=1)
            sizeArr = np.abs(distanceArrReal[opaqueSigIndices]) / np.amax(np.abs(distanceArrReal)) * 100

            plt.scatter(opaqueSigIndices, distanceArrReal[opaqueSigIndices], s=sizeArr, color=rgbaColorArr, marker=".", edgecolors='none')
            ax.axhline(st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
            ax.axhline(-st.gennorm.isf(significanceThreshold/2, beta, loc=loc, scale=scale), linewidth=.25, linestyle="-")
            
            figPath = manhattanDirPath / "manhattan_plot_chr{}.png".format(i+1)
            fig.savefig(figPath, bbox_inches='tight', dpi=300, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)


# Helper function for generating proper tick marks on the manhattan plots
def pvalAxisScaling(ylim, beta, loc, scale):
    yticks = []
    ytickLabels = ["$10^{-16}$", "$10^{-15}$", "$10^{-14}$", "$10^{-13}$", "$10^{-12}$", "$10^{-11}$", "$10^{-10}$", "$10^{-9}$", "$10^{-8}$", "$10^{-7}$", "$10^{-6}$", "$10^{-5}$", "$10^{-4}$", "$1$", "$10^{-4}$", "$10^{-5}$", "$10^{-6}$", "$10^{-7}$", "$10^{-8}$", "$10^{-9}$", "$10^{-10}$", "$10^{-11}$", "$10^{-12}$", "$10^{-13}$", "$10^{-14}$", "$10^{-15}$", "$10^{-16}$"]
    
    yticks.append(-st.gennorm.isf(10**-16/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-15/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-14/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-13/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-12/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-11/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-10/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-9/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-8/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-7/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-6/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-5/2, beta, loc=loc, scale=scale))
    yticks.append(-st.gennorm.isf(10**-4/2, beta, loc=loc, scale=scale))
    yticks.append(0)
    yticks.append(st.gennorm.isf(10**-4/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-5/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-6/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-7/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-8/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-9/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-10/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-11/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-12/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-13/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-14/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-15/2, beta, loc=loc, scale=scale))
    yticks.append(st.gennorm.isf(10**-16/2, beta, loc=loc, scale=scale))
    
    yticksFinal = []
    ytickLabelsFinal = []
    
    for i in range(len(yticks)):
        if yticks[i] >= -ylim and yticks[i] <= ylim:
            yticksFinal.append(yticks[i])
            ytickLabelsFinal.append(ytickLabels[i])
            
    print(yticks)
    print(yticksFinal)
    
    return (yticksFinal, ytickLabelsFinal)
    

# Helper function to create a roiURL bed file of the top 1000 loci (Adjacent loci are merged)
def sendRoiUrl(filePath, locationArr, distanceArr, maxDiffArr, nameArr):
    with open(filePath, 'w') as f:
        # Sort the significant values
        sortedIndices = (-np.abs(distanceArr)).argsort()[:1000]
        locations = np.concatenate((locationArr[sortedIndices], distanceArr[sortedIndices].reshape(len(sortedIndices), 1), maxDiffArr[sortedIndices].reshape(len(sortedIndices), 1)), axis=1)
    
        # Iterate until all is merged
        while(hasAdjacent(locations)):
            locations = mergeAdjacent(locations)
            
        print(len(locations))
        # Write all the locations to the file
        outTemplate = "{0[0]}\t{0[1]}\t{0[2]}\t{1}\t{2}\t{3}\n"
        outString = "".join(outTemplate.format(locations[i], nameArr[int(float(locations[i, 4])) - 1], abs(float(locations[i, 3])), findSign(float(locations[i, 3]))) for i in range(locations.shape[0]))
        f.write(outString)


# Helper function for determining when to stop merging roiURL loci
def hasAdjacent(locationArr):
    # Checks each location against every other location
    for i in range(locationArr.shape[0]):
        for j in range(locationArr.shape[0]):
            # If the chromosomes are the same and they are adjacent return True
            # Also check if the distance is in the same direction
            if locationArr[i, 0] == locationArr[j, 0] and (int(locationArr[i, 2]) - int(locationArr[j, 1]) == 0 or int(locationArr[j, 2]) - int(locationArr[i, 1]) == 0) and np.sign(float(locationArr[i, 3])) == np.sign(float(locationArr[j, 3])):
                return True
    # If we have gotten through everything and not found adjacent locations, return false
    return False


# Helper function for merging adjacent loci in the roiURL bed file
def mergeAdjacent(locationArr):
    for i in range(locationArr.shape[0]):
        for j in range(locationArr.shape[0]):
            # If the chromosomes are the same, they are adjacent, and their distances are in the same direction merge and delete the originals
            if locationArr[i, 0] == locationArr[j, 0] and int(locationArr[i, 2]) - int(locationArr[j, 1]) == 0 and np.sign(float(locationArr[i, 3])) == np.sign(float(locationArr[j, 3])):
                mergedLocation = np.array([locationArr[i, 0], locationArr[i, 1], locationArr[j, 2], locationArr[i, 3], locationArr[i, 4]]).reshape(1, 5)
                locationArr = np.delete(locationArr, [i, j], axis=0)
                locationArr = np.insert(locationArr, i, mergedLocation, axis=0)
                return locationArr
            elif locationArr[i, 0] == locationArr[j, 0] and int(locationArr[j, 2]) - int(locationArr[i, 1]) == 0 and np.sign(float(locationArr[i, 3])) == np.sign(float(locationArr[j, 3])):
                mergedLocation = np.array([locationArr[i, 0], locationArr[j, 1], locationArr[i, 2], locationArr[i, 3], locationArr[i, 4]]).reshape(1, 5)
                locationArr = np.delete(locationArr, [i, j], axis=0)
                locationArr = np.insert(locationArr, i, mergedLocation, axis=0)
                return locationArr
    return locationArr


# Helper function for sendRoiUrl
def findSign(x):
    if (x >= 0):
        return "+"
    else:
        return "-"


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])