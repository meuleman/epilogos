import numpy as np
import pandas as pd
import sys
from time import time
from pathlib import Path
from multiprocessing import cpu_count, Pool, RawArray
from contextlib import closing
from epilogos.helpers import sharedToNumpy, splitRows
from sklearn.metrics.pairwise import euclidean_distances
import scipy.stats as st


def main(outputDir, windowBins, blockSize, nCores, nDesiredMatches, nJobs, processTag):
    print("Calculating search results...", flush=True); t = time()

    genomeCoords = pd.DataFrame(np.load(outputDir / "genome_stats.npz", allow_pickle=True)["coords"],
                                columns=["Chromosome", "Start", "End"])

    simsearchCube = np.load(outputDir / "simsearch_cube.npz", allow_pickle=True)
    roiCube = simsearchCube["scores"]
    roiCoords = pd.DataFrame(simsearchCube["coords"], columns=["Chromosome", "Start", "End"])

    if nCores == 0: nCores = cpu_count()

    rowsToCalc = splitRows(roiCube.shape[0], nJobs)[processTag]

    euclideanDistanceMulti(outputDir, genomeCoords, roiCoords, roiCube, windowBins, blockSize, nCores, nDesiredMatches,
                           rowsToCalc, processTag)

    print("    Time:", format(time() - t, '.0f'), "seconds\n", flush=True)


def _initEuclideanDistance(genomeCoords_, reducedGenome_, roiCoords_, roiCube_, sharedSimilarRegionArr_, windowBins_,
                           blockSize_, nDesiredMatches_):
    """
    Initializes global variables for multiprocessing the convolution calculation

    Input:
    genomeCoords_           -- chr, start, and end coordinates for the non-reduced genome
    reducedGenome_          -- The reduced state scores
    roiCoords_              -- chr, start, and end coordinates for the maxmean chosen regions of interest
    roiCube_                -- Reduced scores for the regions of interest
    sharedSimilarRegionArr_ -- A tuple containing relevant information about the shared array (stores simsearch results)
    windowBins_             -- The number of bins to have in each window
    blockSize_              -- The reduction factor for the similiarty search
    nDesiredMatches_        -- The max number of simsearch results to output for each region of interest
    """
    global genomeCoords
    global reducedGenome
    global roiCoords
    global roiCube
    global sharedSimilarRegionArr
    global windowBins
    global blockSize
    global nDesiredMatches

    genomeCoords = genomeCoords_
    reducedGenome = reducedGenome_
    roiCoords = roiCoords_
    roiCube = roiCube_
    sharedSimilarRegionArr = sharedSimilarRegionArr_
    windowBins = windowBins_
    blockSize = blockSize_
    nDesiredMatches = nDesiredMatches_


def runEuclideanDistance(rowsToCalcMulti):
    """
    Runs the euclidean distance function for the rows in roiCube present in rowlist

    Input:
    rowsToCalcMulti -- A tuple containing the first row (inclusive) and last row (exclusive) to investigate
    To see rest of input see _initMultiprocessing description

    Output:
    Writes indices of the similarity search results to the shared array for each of the rows in rowList
    """

    if type(sharedSimilarRegionArr) is tuple:
        similarRegionArr = sharedToNumpy(*sharedSimilarRegionArr)
    else:
        similarRegionArr = sharedSimilarRegionArr

    # Precalculate values for euclidean distance indexing
    nSuperBins = windowBins // blockSize
    distanceSize = len(reducedGenome) - (nSuperBins - 1)
    diagonalIndices0 = np.add(*np.broadcast_arrays(np.arange(nSuperBins),
                                                   np.arange(distanceSize).reshape(distanceSize, 1)))
    diagonalIndices1 = np.broadcast_to(np.arange(nSuperBins), (distanceSize, nSuperBins))

    # Run convolvultion on all regions
    for row in range(*rowsToCalcMulti):

        euclideanDistances = np.sum(euclidean_distances(reducedGenome, roiCube[row], squared=True)[diagonalIndices0,
                                                                                                   diagonalIndices1],
                                    axis=1)

        halfMode = st.mode(euclideanDistances, keepdims=False)[0] / 2

        # Only want to recommend non-overlapping regions
        overlapArr = np.zeros(len(reducedGenome))

        # Do not want to recommend region itself
        regionStart = np.where((genomeCoords["Chromosome"] == roiCoords.iloc[row, 0])
                               & (genomeCoords["Start"] == roiCoords.iloc[row, 1]))[0][0] // blockSize
        overlapArr[regionStart:regionStart + nSuperBins] = 1

        # Loop to find nDesiredMatches similar regions
        numMatches = 0
        for hitIndex in np.argsort(euclideanDistances):
            # Don't want overlapping regions
            if np.any(overlapArr[hitIndex:hitIndex + nSuperBins]):
                continue
            elif euclideanDistances[hitIndex] > halfMode:
                similarRegionArr[row, numMatches:] = -1
                break
            else:
                # Store the indices of these similar regions (can convert later)
                similarRegionArr[row, numMatches] = hitIndex
                overlapArr[hitIndex:hitIndex + nSuperBins] = 1
                numMatches += 1
                if numMatches >= nDesiredMatches:
                    break


def euclideanDistanceMulti(outputDir, genomeCoords, roiCoords, roiCube, windowBins, blockSize, nCores, nDesiredMatches,
                           rowsToCalc, processTag):
    """
    Runs a convolution funnction on normalized data to determine the most similar regions for all ROIs across the genome

    Input:
    outputDir       -- The output directory for the results
    genomeCoords    -- chr, start, and end coordinates for the non-reduced genome
    roiCoords       -- Pandas Dataframe of the coordinates for the mean max generated windows
    roiCube         -- The reduced scores over the mean max generated windows
    windowBins      -- The number of bins to have in each window
    blockSize       -- The reduction factor for the similiarty search
    nCores          -- The number of cores to use for the multiprocessing of the convolution calculation
    nDesiredMatches -- Max Number of similar regions to be returned for each roi
    rowsToCalc      -- The rows within roiCube over which the running slurm job should calculate simsearch values
    processTag      -- An integer which indicates which slurm job we are running (and thus which rows to calculate)

    Output:
    simsearch_knn.npz      -- for each of regions, the top nDesiredMatches matches
                              (arr=coords, idx=indices, dist=distances)
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredMatches matches for each of regions
    """

    reducedGenome = np.load(outputDir / "reduced_genome.npy", allow_pickle=True)

    nRegions = rowsToCalc[1] - rowsToCalc[0]

    if nCores > 1:

        # Split into multiprocessing
        # Shared arrays across the multiple processes for storing the recommendations
        # We avoid race conditions by writing to separate parts of the array in each process
        sharedSimilarRegionArr = RawArray(np.ctypeslib.as_ctypes_type(np.int32), nRegions * nDesiredMatches)

        rowList = splitRows(nRegions, nCores)

        # Start the processes
        with closing(Pool(nCores, initializer=_initEuclideanDistance,
                          initargs=(genomeCoords, reducedGenome, roiCoords.iloc[rowsToCalc[0]:rowsToCalc[1]],
                                    roiCube[rowsToCalc[0]:rowsToCalc[1]],
                                    (sharedSimilarRegionArr, nRegions, nDesiredMatches), windowBins, blockSize,
                                    nDesiredMatches))) as pool:
            pool.map(runEuclideanDistance, rowList)
        pool.join()

        np.save(outputDir / "simsearch_indices_{}.npy".format(processTag),
                sharedToNumpy(sharedSimilarRegionArr, nRegions, nDesiredMatches), allow_pickle=True)
    else:
        similarRegionArr = np.zeros((nRegions, nDesiredMatches), dtype=np.int32)

        _initEuclideanDistance(genomeCoords, reducedGenome, roiCoords.iloc[rowsToCalc[0]:rowsToCalc[1]],
                               roiCube[rowsToCalc[0]:rowsToCalc[1]], similarRegionArr, windowBins, blockSize,
                               nDesiredMatches)

        runEuclideanDistance((0, nRegions))

        np.save(outputDir / "simsearch_indices_{}.npy".format(processTag), similarRegionArr, allow_pickle=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]),
         int(sys.argv[7]))
