"""
Similarity Search For Epilogos

Written by Jacob Quon, Nalu Tripician, Alex Reynolds, and Wouter Meuleman
"""
import numpy as np
import pandas as pd
import sys
from time import time
from pathlib import Path
import os
import tempfile
import json
import csv
import pysam
import click
from multiprocessing import cpu_count, Pool, RawArray
from contextlib import closing
from epilogos.helpers import generateRegionArr, maxMean, sharedToNumpy, splitRows
from sklearn.metrics.pairwise import euclidean_distances


# Custom help format for click to separate the build and query commands
class CustomClickHelpFormat(click.Command):
    def format_help(self, ctx, formatter):
        click.echo("Usage: simsearch [OPTIONS]\n\
\n\
Options:\n\
To Build Similarity Search Data:\n\
  -b, --build                     If true builds the similarity search files\n\
                                  needed to query regions\n\
  -s, --scores TEXT               Path to scores file to be used in\n\
                                  similarity search\n\
  -o, --output-directory TEXT     Path to desired similarity search output\n\
                                  directory\n\
  -w, --window-bp INTEGER         Window size (in BP) on which to perform\n\
                                  similarity search  [default: 25000]\n\
  -j, --num-jobs INTEGER          Number of cores to be used in simsearch\n\
                                  calculation. If set to 0, uses all cores.\n\
                                  [default: 0=All cores]\n\
  -n, --num-matches INTEGER       Number of matches to be found by simsearch\n\
                                  for each query region [default: 100]\n\
  -f, --filter-state INTEGER      If the max signal within a region is from the\n\
                                  filter state, it is removed from the region\n\
                                  list. The purpose of this is to be used to \n\
                                  filter out 'quiescent regions'. If set to 0,\n\
                                  filtering is not done. [default: last state]\
\n\
\n\
To Query Similarity Search Data:\n\
  -q, --query TEXT                Query region formatted as chr:start-end or\n\
                                  path to tab-separated bed file containing\n\
                                  query regions\n\
  -m, --matches-file TEXT\n\
                                  Path to previously built\n\
                                  simsearch.bed.gz file to be queried\n\
                                  for matches\n\
  -o, --output-directory TEXT     Path to desired similarity search output\n\
                                  directory")

@click.command(context_settings=dict(help_option_names=['-h', '--help']), cls=CustomClickHelpFormat)
@click.option("-b", "--build", "buildBool", is_flag=True,
              help="If true builds the similarity search files needed to query regions")
@click.option("-s", "--scores", "scoresPath", type=str, help="Path to scores file to be used in similarity search")
@click.option("-o", "--output-directory", "outputDir", required=True, type=str,
              help="Path to desired similarity search output directory")
@click.option("-w", "--window-bp", "windowBP", type=int, default=-1, show_default=True,
              help="Window size (in BP) on which to perform similarity search")
@click.option("-c", "--num-cores", "nCores", type=int, default=0, show_default=True,
              help="Number of cores to be used in simsearch calculation. If set to 0, uses all cores. [default: 0=All cores]")
@click.option("-n", "--num-matches", "nDesiredMatches", type=int, default=100, show_default=True,
              help="Number of matches to be found by simsearch for each query region [default: 100]")
@click.option("-f", "--filter-state", "filterState", type=int, default=-1,
              help="If the max signal within a region is from the filter state, it is removed from the region list. " +
                   "The purpose of this is to be used to filter out 'quiescent regions'." +
                   "If set to 0, filtering is not done. [default: last state]")
@click.option("-q", "--query", "query", type=str, default="",
              help="Query region formatted as chr:start-end or path to tab-separated bed file containing query regions")
@click.option("-m", "--matches-file", "simSearchPath", type=str,
              help="Path to previously built simsearch.bed.gz file to be queried for matches")
def main(buildBool, scoresPath, outputDir, windowBP, nCores, nDesiredMatches, filterState, query, simSearchPath):
    print("""\n
                 d8b                                                           888
                 Y8P                                                           888
                                                                               888
        .d8888b  888 88888b.d88b.   .d8888b   .d88b.   8888b.  888d888 .d8888b 88888b.
        88K      888 888 "888 "88b  88K      d8P  Y8b     "88b 888P"  d88P"    888 "88b
        "Y8888b. 888 888  888  888  "Y8888b. 88888888 .d888888 888    888      888  888
             X88 888 888  888  888       X88 Y8b.     888  888 888    Y88b.    888  888
         88888P' 888 888  888  888   88888P'  "Y8888  "Y888888 888     "Y8888P 888  888
    """, flush=True)

    if not buildBool and query == "":
        raise ValueError("Either -b or -q flag must be used to run simsearch")
    elif buildBool and query != "":
        raise ValueError("Both -b and -q flags cannot be used at the same time")

    outputDir = Path(outputDir)
    if not outputDir.exists():
        outputDir.mkdir(parents=True)
    if not outputDir.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDir)))

    if buildBool:
        if nCores == 0: nCores = cpu_count()
        buildSimSearch(scoresPath, outputDir, windowBP, nCores, nDesiredMatches, filterState)
    else:
        querySimSearch(query, simSearchPath, outputDir)


def buildSimSearch(scoresPath, outputDir, windowBP, nCores, nDesiredMatches, filterState):
    """
    Builds the similarity search files for a set of scores

    Input:
    scoresPath        -- Path to scores file to be used in similarity search
    outputDir         -- pathlib Path to the output directory for similarity search files
    windowBP          -- Size of the similarity search window in bp
    nCores             -- Number of jobs to be used in nearest neighbor algorithm
    nDesiredMatches -- Number of matches to be found by nearest neighbor algorithm (first match is always the query region)

    Output:
    simsearch_cube.npz     -- NPZ file containing all the reduced regions and their coords (scores=scores, coords=coords)
    simsearch_knn.npz      -- NPZ file containing the top nDesiredMatches matches for each region
                              (arr=coords, idx=indices, dist=distances)
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredMatches matches for each of regions
    """

    print("\n\n\n        Reading in data...", flush=True); readTime = time()

    # Bins are assumed to be 200bp or 20bp
    if determineBinSize(scoresPath) == 200:
        if windowBP == -1: windowBP = 25000
        windowBins = int(windowBP / 200)
        blockSize = determineBlockSize200(windowBP)
    elif determineBinSize(scoresPath) == 20:
        if windowBP == -1: windowBP = 2500
        windowBins = int(windowBP / 20)
        blockSize = determineBlockSize20(windowBP)
    else:
        raise ValueError("Similarity Search is only compatible with bins of size 200bp or 20bp")

    stateScores, inputArr = readScores(scoresPath)

    # The maximum number of regions chosen should depend on the window size
    # We want to allow for full coverage of the genome if possible (maxRegions is chosen accordingly)
    maxRegions = int(stateScores.shape[0] // windowBins)

    print("            Time:", format(time() - readTime,'.0f'), "seconds\n", flush=True)
    print("        Finding regions of size {}kb...".format(windowBP // 1000), flush=True); cubeTime = time()

    # Filter-regions package to perform maxmean algorithm & pull out top X regions
    rois, _ = maxMean(inputArr, windowBins, maxRegions)

    # Seperates just the coordinates
    roiCoords = rois.iloc[:,:3].reset_index(drop=True, inplace=True)
    # Generate slices of reduced data for the cube
    roiCube = rois["OriginalIdx"].apply(lambda x: makeSlice(stateScores, x, windowBins, blockSize))
    roiCube = np.stack(roiCube.to_numpy())

    # Dropping overlapping chromosomes
    chrOverlapIndices = np.where(roiCoords["Start"] >= roiCoords["End"])[0]
    roiCoords.drop(labels=chrOverlapIndices, inplace=True)
    roiCoords.reset_index(drop=True, inplace=True)
    roiCube = np.delete(roiCube, chrOverlapIndices, axis=0)

    # Dropping regions where the quiescent state is the maximally contributing state
    if filterState != 0:
        filterState = roiCube.shape[2] - 1 if filterState == -1 else filterState - 1
        quiescentIndices = np.where(np.argmax(np.max(roiCube, axis=1), axis=1) == filterState)[0]
        roiCoords.drop(labels=quiescentIndices, inplace=True)
        roiCoords.reset_index(drop=True, inplace=True)
        roiCube = np.delete(roiCube, quiescentIndices, axis=0)

    # Save the cube
    np.savez_compressed(file=outputDir / 'simsearch_cube', scores=roiCube, coords=roiCoords.values)

    print("            Time:", format(time() - cubeTime,'.0f'), "seconds\n", flush=True)
    print("        Finding simsearch matches...", flush=True); knnTime = time()

    # for each region find the 101 nearest neighbors to be used as matches
    findSimilarRegions(outputDir, pd.DataFrame(inputArr[:,:3], columns=["Chromosome", "Start", "End"]), stateScores, roiCoords, roiCube, windowBins, blockSize, nCores, nDesiredMatches)

    print("            Time:", format(time() - knnTime,'.0f'), "seconds\n", flush=True)


def determineBinSize(scoresPath):
    """
    Determines the size of an individual bin using the first row of the scores file

    Input:
    scoresPath -- Path to scores file to be used in similarity search

    Output:
    int(row1.iloc[0, 2] - row1.iloc[0, 1]) -- size of an individual bin
    """
    row1 = pd.read_table(scoresPath, sep="\t", header=None, usecols=[0,1,2], nrows=1)
    return int(row1.iloc[0, 2] - row1.iloc[0, 1])


def querySimSearch(query, simSearchPath, outputDir):
    """
    Queries the simsearch.bed.gz file a set of input regions and outputs each regions matches to an individual output file

    Input:
    query         -- A 2d numpy array containing the coordinates of user queried regions
    simSearchPath -- The path to the simsearch.bed.gz file to be queried for matches
    outputDir     -- The directory to output each regions matches
    """
    print("\n\n\n        Reading in data...", flush=True); readTime = time()

    # Parse query to generate array
    queryArr = generateRegionArr(query)

    # Read in matches file
    matchesDF = pd.read_table(Path(simSearchPath), sep="\t", header=None)

    print("            Time:", format(time() - readTime,'.0f'), "seconds\n", flush=True)
    print("        Querying regions...", flush=True)

    # For each query output top 100 matches to a bed file
    for chr, start, end in queryArr:
        index = np.where((matchesDF.iloc[:, 0] == chr) & (matchesDF.iloc[:, 1] >= start) & (matchesDF.iloc[:, 2] <= end))[0]
        if index.size > 0:
            index = index[0]

            # Find source coords for matches
            regionChr, regionStart, regionEnd = matchesDF.iloc[index, :3]

            # write to bed file
            outfile = outputDir / "similarity_search_region_{}_{}_{}_recs.bed".format(regionChr, regionStart, regionEnd)
            with open(outfile, "w+") as f:
                matchesStr = ""

                recs = matchesDF.iloc[index,3][2:-2] # trim off brackets
                recs = recs.split('", "')[1:] # split matches
                for i in range(len(recs)):
                    # Append bed formatted coords to file
                    matchesStr += "{0[0]}\t{0[1]}\t{0[2]}\n".format(recs[i].split(":"))

                f.write(matchesStr)

            # Print out information to user
            print("            Found region {}:{}-{} within user query {}:{}-{}".format(regionChr, regionStart, regionEnd,
                                                                                        chr, start, end))
            print("                See {} for matches\n".format(outfile), flush=True)
        else:
            print("            ValueError: Could not find region in given query range: {}:{}-{}\n".format(chr, start, end))


def determineBlockSize20(windowBP):
    """
    Similarity search reduces scores over search windows to mimic what the human eye might focus on.
    This is meant to provide regions which contain similar characteristics/trends
    This function determines the reduction factor for a input file which has bins of 20bp

    Input:
    windowBP - The size of the search window in BP

    Output:
    blockSize - The reduction factor for the similarity search
    """

    # We only support some set window sizes, these each have a reduction factor to
    # create an overall data size of 25 points per window
    if windowBP == 500:
        blockSize = 1
    elif windowBP == 1000:
        blockSize = 2
    elif windowBP == 2500:
        blockSize = 5
    elif windowBP == 5000:
        blockSize = 10
    elif windowBP == 7500:
        blockSize = 15
    elif windowBP == 10000:
        blockSize = 20
    else:
        raise ValueError("Error: window size must be either 500, 1000, 2500, 5000, 7500, or 10000 (in bp)")

    return blockSize


def determineBlockSize200(windowBP):
    """
    Similarity search reduces scores over search windows to mimic what the human eye might focus on.
    This is meant to provide regions which contain similar characteristics/trends
    This function determines the reduction factor for a input file which has bins of 200bp

    Input:
    windowBP - The size of the search window in BP

    Output:
    blockSize - The reduction factor for the similarity search
    """

    # We only support some set window sizes, these each have a reduction factor to
    # create an overall data size of 25 points per window
    if windowBP == 5000:
        blockSize = 1
    elif windowBP == 10000:
        blockSize = 2
    elif windowBP == 25000:
        blockSize = 5
    elif windowBP == 50000:
        blockSize = 10
    elif windowBP == 75000:
        blockSize = 15
    elif windowBP == 100000:
        blockSize = 20
    else:
        raise ValueError("Error: window size must be either 5000, 10000, 25000, 50000, 75000, or 100000 (in bp)")

    return blockSize


def readScores(scoresPath):
    """
    Reads in the scores file provided by the user and computes a per bin sum

    Input:
    scoresPath - The path to the scores file provided by the user

    Output:
    scores - Just the scores from the scores file
    inputArr - A numpy array to be read in by the filter-regions package (contains coords and summed scores)
    """

    # Read in the scores file for regions of interest generation
    scores = pd.read_table(scoresPath, sep="\t", header=None)

    # Sum up the scores for all states for maxMean algorithm
    inputDF = scores.iloc[:, :3]
    inputDF[3] = scores.iloc[:,3:].sum(axis=1)
    inputDF.columns = ["Chromosome", "Start", "End", "Score"]
    inputArr = inputDF.to_numpy()

    # Pull out just the scores for later window reduction
    scores = scores.iloc[:,3:]

    return scores, inputArr


def makeSlice(genome, idx, windowBins, blockSize):
    """
    Creates a slice of a reduced scores track used for construction of the scores "cube"
    This will take a location and use it as a center point for each "slice" (region) of the cube

    Input:
    genome     -- The scores across the whole genome
    idx        -- The centerpoint of a similarity search window
    windowBins -- The number of bins to have in each window
    blockSize  -- The reduction factor for the similiarty search

    Output:
    genomeWindow.loc[reducedIndices[0].values].to_numpy() -- A numpy array representing the reduced scores of a window
    """
    genomeWindow = genome.iloc[idx - windowBins // 2:idx + windowBins // 2 + 1] if windowBins % 2\
                   else genome.iloc[idx - windowBins // 2:idx + windowBins // 2]
    sumsOverWindow = genomeWindow.sum(axis=1).to_frame()
    reducedIndices = sumsOverWindow.groupby(np.arange(len(sumsOverWindow), dtype=np.float32)//blockSize,).idxmax()
    return genomeWindow.loc[reducedIndices[0].values].to_numpy()


def _initMultiprocessing(genomeCoords_, reducedGenome_, roiCoords_, roiCube_, sharedArr_, windowBins_, blockSize_, nDesiredMatches_):
    """
    Initializes global variables for multiprocessing the convolution calculation

    Input:
    genomeCoords_   -- chr, start, and end coordinates for the non-reduced genome
    reducedGenome_  -- The reduced state scores
    roiCoords_      -- chr, start, and end coordinates for the maxmean chosen regions of interest
    roiCube_        -- Reduced scores for the regions of interest
    sharedArr_      -- A tuple containing relevant information about the shared array (stores simsearch results)
    windowBins_     -- The number of bins to have in each window
    blockSize_      -- The reduction factor for the similiarty search
    nDesiredMatches_ -- The number of simsearch results to output for each region of interest
    """
    global genomeCoords
    global reducedGenome
    global roiCoords
    global roiCube
    global sharedArr
    global windowBins
    global blockSize
    global nDesiredMatches

    genomeCoords = genomeCoords_
    reducedGenome = reducedGenome_
    roiCoords = roiCoords_
    roiCube = roiCube_
    sharedArr = sharedArr_
    windowBins = windowBins_
    blockSize = blockSize_
    nDesiredMatches = nDesiredMatches_


def runEuclideanDistance(rowsToCalc):
    """
    Runs the euclidean distance function for the rows in roiCube present in rowlist

    Input:
    rowsToCalc -- A tuple containing the first row (inclusive) and last row (exclusive) to investigate
    To see rest of input see _initMultiprocessing description

    Output:
    Writes indices of the similarity search results to the shared array for each of the rows in rowList
    """

    similarRegionArr = sharedToNumpy(*sharedArr)

    # Run convolvultion on all regions
    for row in range(*rowsToCalc):
        """
        euclidean_distances(reducedGenome, roiCube[row], squared=True)
            Calculates the euclidean distances all of the single bin vectors in the the region and single bin positions in the
            genome. This returns a len(reducedGenome)x(windowBins // blockSize) size array. The sum of each diagonal represents
            the squared euclidean distance of that location in the genome.
            i.e. assume a minimal case where the distance array is 30x4,
            sum(distarr[5,0], distarr[6,1], distarr[7,2], distarr[8,3])
            would represent the euclidean distance between the query and the window created by indexing reducedGenome[5,6,7,8]

        np.add(*np.broadcast_arrays(np.arange(25), np.arange(distanceSize).reshape(distanceSize, 1)))
            Calculates the primary indices required to create a 2d array which summed across its secondary axis gives the
            euclidean distance for each position

            np.broadcast_arrays(np.arange(25), np.arange(distanceSize).reshape(distanceSize, 1))
                Creates two arrays of size len(reducedGenome)-24x(windowBins // blockSize) and returns them in a list
                The first of these arrays is simply len(reducedGenome)-24 rows containing the ints [0,24] inclusive
                The second  is len(reducedGenome)-24 rows with each row containing its primary index repeated 25 times

            np.add(*np.broadcast_arrays())
                Takes the two generated arrays and adds them together so that each row counts up 25 ints from its start
                i.e. [[0,1,2,3, ... 24],
                      [1,2,3,4, ... 25]
                      [2,3,4,5, ... 26],
                      ...]
                This represents the primary indices necessary to easily sum each diagonal

            The reason we use this format as opposed to the inbuilt python broadcasting of
                np.arange(25) + np.arange(len(reducedGenome_25kb) - 24).reshape(len(reducedGenome_25kb) - 24, 1)
            is because this method is slightly faster which matters when we are running it so many times across the genome

        np.broadcast_to(np.arange(25), (distanceSize, 25))
            Calculates the secondary indices required to create a 2d array which summed across its secondary axis gives the
            euclidean distance for each position. It generates a len(reducedGenome)-24x(windowBins // blockSize) size array
            where each row simply contains the ints [0,24] inclusive

        np.sum(euclidean_distances()[], axis=1)
            Takes the results of the euclidean distance calculation and indexing and sums them so that we obtain the proper
            euclidean distances
        """
        distanceSize = len(reducedGenome) - 24
        euclideanDistances = np.sum(euclidean_distances(reducedGenome, roiCube[row], squared=True)[np.add(*np.broadcast_arrays(np.arange(25), np.arange(distanceSize).reshape(distanceSize, 1))), np.broadcast_to(np.arange(25), (distanceSize, 25))], axis=1)

        # Only want to recommend non-overlapping regions
        overlap_arr = np.zeros(len(reducedGenome))

        # Do not want to recommend region itself
        regionStart = np.where((genomeCoords["Chromosome"] == roiCoords.iloc[row, 0]) & (genomeCoords["Start"] == roiCoords.iloc[row, 1]))[0][0] // blockSize # Genome coords is the coords for the whole unreduced genoem
        overlap_arr[regionStart:regionStart + windowBins // blockSize] = 1

        # Loop to find nDesiredMatches similar regions
        numMatches = 0
        for hitIndex in np.argsort(euclideanDistances):
            # Don't want overlapping regions
            if np.any(overlap_arr[hitIndex:hitIndex + windowBins // blockSize]):
                continue
            else:
                similarRegionArr[row, numMatches] = hitIndex  # Store the indices of these similar regions (can convert later)
                overlap_arr[hitIndex:hitIndex + windowBins // blockSize] = 1
                numMatches += 1
                if numMatches >= nDesiredMatches:
                    break


def findSimilarRegions(outputDir, genomeCoords, stateScores, roiCoords, roiCube, windowBins, blockSize, nCores, nDesiredMatches):
    """
    Runs a convolution funnction on normalized data to determine the most similar regions for all ROIs across the genome

    Input:
    outputDir      -- The output directory for the results
    genomeCoords   -- chr, start, and end coordinates for the non-reduced genome
    stateScores    -- The scores for each state in each bin of the genom
    roiCoords      -- Pandas Dataframe of the coordinates for the mean max generated windows
    roiCube        -- The reduced scores over the mean max generated windows
    windowBins     -- The number of bins to have in each window
    blockSize      -- The reduction factor for the similiarty search
    nCores         -- The number of cores to use for the multiprocessing of the convolution calculation
    nDesiredMatches -- Number of similar regions to be returned for each roi

    Output:
    simsearch_knn.npz      -- for each of regions, the top nDesiredMatches matches (arr=coords, idx=indices, dist=distances)
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredMatches matches for each of regions
    """

    # Reduce the genome for the appropriate region size
    # Give the same blockIndex for every sequential group of blockSize
    sumsOverGenome = stateScores.sum(axis=1).to_frame(name="scores")
    sumsOverGenome["blockIndices"] = np.arange(len(sumsOverGenome), dtype=np.int32) // blockSize
    # For each blockIndex group, output only the index of the bin with the highest overall score. use these indices to get reduced scores
    reducedIndices = np.array(sumsOverGenome.sort_values('scores').drop_duplicates(['blockIndices'], keep='last').sort_values('blockIndices').index)
    reducedGenome = stateScores.iloc[reducedIndices].to_numpy()

    # Normalize the genome using (x-mean)/std for each superbin
    # reducedGenomeNormalized  = ((reducedGenome.transpose() - np.mean(reducedGenome, axis=1)) / np.std(reducedGenome, axis=1)).transpose()

    # Split into multiprocessing
    # Shared arrays across the multiple processes for storing the recommendations
    # We avoid race conditions by writing to separate parts of the array in each process
    numRegions = roiCube.shape[0]
    sharedArr = RawArray(np.ctypeslib.as_ctypes_type(np.int32), numRegions * nDesiredMatches)

    #########################################################################################################################################################################################
    ##                                                                                                                                                                                     ##
    ## If there are memory issues might have to do something where we write cube to disk earlier and then read in only relevevant parts (hopefully garbage collection will deal with cube) ##
    ##                                                                                                                                                                                     ##
    #########################################################################################################################################################################################

    rowList = splitRows(numRegions, nCores)

    # Start the processes
    with closing(Pool(nCores, initializer=_initMultiprocessing,
                      initargs=(genomeCoords, reducedGenome, roiCoords, roiCube,
                                (sharedArr, numRegions, nDesiredMatches), windowBins, blockSize, nDesiredMatches))) as pool:
        pool.map(runEuclideanDistance, rowList)
    pool.join()

    # Generate coords for the reduced genome
    firstIndices = np.array(sorted(sumsOverGenome.drop_duplicates(['blockIndices'], keep='first').index))
    lastIndices = np.array(sorted(sumsOverGenome.drop_duplicates(['blockIndices'], keep='last').index))
    reducedGenomeCoords = pd.concat((genomeCoords.iloc[firstIndices, :2].reset_index(drop=True), genomeCoords.iloc[lastIndices, 2].reset_index(drop=True)), axis=1, names=["Chromsome", "Start", "End"])

    # take the shared array and convert all of the indices into locations
    resultsChrAndStart = reducedGenomeCoords.iloc[sharedToNumpy(sharedArr, numRegions, nDesiredMatches).reshape(numRegions * nDesiredMatches), :2].values
    resultsEnd = reducedGenomeCoords.iloc[sharedToNumpy(sharedArr, numRegions, nDesiredMatches).reshape(numRegions * nDesiredMatches) + windowBins // blockSize - 1, 2].values.reshape(numRegions * nDesiredMatches, 1)

    searchResults = np.concatenate((resultsChrAndStart, resultsEnd), axis=1).reshape(numRegions, nDesiredMatches, 3)

    # Store the lcoation of the query at the beginning of array
    searchResults = np.concatenate((roiCoords.values, searchResults), axis=1)

    i = 0
    final = ['' for i in range(numRegions)]
    for row in searchResults:
        query_v = roiCoords.iloc[i,].values
        query = '{}:{}:{}'.format(query_v[0], query_v[1], query_v[2])
        recs = []
        query_added = False
        for chrom, start, end in row:
            hit = '{}:{}:{}'.format(chrom, start, end)
            if query != hit:
                recs.append(hit)
            else:
                recs = [query] + recs
                query_added = True
        # if the query region was not added to the beginning of the final list,
        # we push it at the front here and remove a trailing element
        if not query_added:
            recs = [query] + recs[:-1]
        assert(len(recs) == nDesiredMatches)
        final[i] = json.dumps(recs)
        i += 1

    tbx = pd.DataFrame(roiCoords)
    tbx[3] = final
    tbx.columns = ['chrom', 'start', 'stop', 'matches']
    tbx = tbx.sort_values(by=['chrom','start'])

    simsearch_fn = os.path.join(outputDir, "simsearch.bed.gz")
    simsearch_idx_fn = os.path.join(outputDir, "simsearch.bed.gz.tbi")

    if os.path.exists(simsearch_fn): os.remove(simsearch_fn)
    if os.path.exists(simsearch_idx_fn): os.remove(simsearch_idx_fn)

    # save as tabix
    temp_csv_fn = None
    with tempfile.NamedTemporaryFile(mode='w+b', delete=False, dir=os.path.dirname(os.path.realpath(__file__))) as temp_fh:
        tbx.to_csv(temp_fh, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
        temp_fh.seek(0)
        temp_csv_fn = temp_fh.name

    pysam.tabix_compress(temp_csv_fn, simsearch_fn, force=True)
    if not os.path.exists(simsearch_fn) or os.stat(simsearch_fn).st_size == 0:
        raise Exception("Error: Could not create bgzip archive [{}]".format(simsearch_fn))

    pysam.tabix_index(simsearch_fn, force=True, zerobased=True, preset="bed")
    if not os.path.exists(simsearch_idx_fn) or os.stat(simsearch_idx_fn).st_size == 0:
        raise Exception("Error: Could not create index of bgzip archive [{}]".format(simsearch_idx_fn))

    if os.path.exists(temp_csv_fn): os.remove(temp_csv_fn)


if __name__ == "__main__":
    main()