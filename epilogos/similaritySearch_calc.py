import numpy as np
import pandas as pd
import sys
from time import time
from pathlib import Path
from multiprocessing import cpu_count, Pool, RawArray
from contextlib import closing
from epilogos.helpers import sharedToNumpy, splitRows
from sklearn.metrics.pairwise import euclidean_distances

def main(outputDir, windowBins, blockSize, nCores, nDesiredMatches, nJobs, processTag):
    print("Calculating search results...", flush=True); t = time()
    genomeCoords = pd.DataFrame(np.load(outputDir / "genome_stats.npz", allow_pickle=True)["coords"], columns=["Chromosome", "Start", "End"])

    simsearchCube = np.load(outputDir / "simsearch_cube.npz", allow_pickle=True)
    roiCube = simsearchCube["scores"]
    roiCoords = pd.DataFrame(simsearchCube["coords"], columns=["Chromosome", "Start", "End"])

    if nCores == 0: nCores = cpu_count()

    rowsToCalc = splitRows(roiCube.shape[0], nJobs)[processTag]

    euclideanDistanceMulti(outputDir, genomeCoords, roiCoords, roiCube, windowBins, blockSize, nCores, nDesiredMatches, rowsToCalc, processTag)

    print("    Time:", format(time() - t,'.0f'), "seconds\n", flush=True)

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


def runEuclideanDistance(rowsToCalcMulti):
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
    for row in range(*rowsToCalcMulti):
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


def euclideanDistanceMulti(outputDir, genomeCoords, roiCoords, roiCube, windowBins, blockSize, nCores, nDesiredMatches, rowsToCalc, processTag):
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

    reducedGenome = np.load(outputDir / "reduced_genome.npy", allow_pickle=True)

    # Split into multiprocessing
    # Shared arrays across the multiple processes for storing the recommendations
    # We avoid race conditions by writing to separate parts of the array in each process
    nRegions = rowsToCalc[1] - rowsToCalc[0]
    sharedArr = RawArray(np.ctypeslib.as_ctypes_type(np.int32), nRegions * nDesiredMatches)

    rowList = splitRows(nRegions, nCores)

    # Start the processes
    with closing(Pool(nCores, initializer=_initMultiprocessing,
                      initargs=(genomeCoords, reducedGenome, roiCoords.iloc[rowsToCalc[0]:rowsToCalc[1]], roiCube[rowsToCalc[0]:rowsToCalc[1]],
                                (sharedArr, nRegions, nDesiredMatches), windowBins, blockSize, nDesiredMatches))) as pool:
        pool.map(runEuclideanDistance, rowList)
    pool.join()

    np.save(outputDir / "simsearch_regions_{}.npy".format(processTag), sharedToNumpy(sharedArr, nRegions, nDesiredMatches), allow_pickle=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]))