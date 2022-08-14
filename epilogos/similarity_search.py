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
import filter_regions as fr
from sklearn.neighbors import NearestNeighbors
import tempfile
import json
import csv
import pysam
import click

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-s", "--scores-path", "epilogosScoresPath", type=str, help="Path to epilogos scores file to be used in similarity search")
@click.option("-o", "--output-directory", "outputDir", type=str, help="Path to desired similarity search output directory")
@click.option("-w", "--window-KB", "windowKB", type=int, default=25, show_default=True, help="Window size (in KB) on which to perform similarity search")
@click.option("-j", "--number-of-jobs", "nJobs", type=int, default=8, show_default=True, help="Number of jobs to be used in nearest neighbor algorithm")
@click.option("-n", "--number-of-neighbors", "nDesiredNeighbors", type=int, default=101, show_default=True, help="Number of neighbors to be found by nearest neighbor algorithm (note that first neighbor is always the query region)")
def main(epilogosScoresPath, outputDir, windowKB, nJobs, nDesiredNeighbors):
    outputDir = Path(outputDir)

    # Bins are assumed to be 200bp, thus there are 5 bins per KB
    windowBins = windowKB * 5

    blockSize = determineBlockSize(windowKB)

    try:
        os.makedirs(outputDir)
    except OSError:
        print("Directory exists, or cannot be created")

    cubeTime = time()

    epilogosScores, inputArr = readScores(epilogosScoresPath)

    # The maximum number of regions chosen should depend on the window size
    # We want to allow for full coverage of the genome if possible (maxRegions is chosen accordingly)
    maxRegions = epilogosScores.shape[0] // windowBins

    # Filter-regions package to perform maxmean algorithm & pull out top X regions
    f = fr.Filter(method='maxmean', input=inputArr, input_type='bedgraph', aggregation_method='max', window_bins=windowBins, max_elements=maxRegions, preserve_cols=True, quiet=False)
    f.read()
    f.filter()
    exemplars = f.output_df

    # Build cube which represents reduced epilogos scores for the top X regions
    exemplars = exemplars.sort_values(by=["RollingMax", "RollingMean", "Score"], ascending=False).reset_index(drop=True)
    # Seperates just the coordinates
    coords = exemplars.iloc[:,:3]
    # Generate slices of reduced data for the cube
    cube = exemplars["OriginalIdx"].apply(lambda x: makeSlice(epilogosScores, x, windowBins, blockSize))
    cube = np.stack(cube.to_numpy())

    # Save the cube
    np.savez_compressed(file=outputDir / 'epilogos_cube', scores=cube, coords=coords.values)

    print("Cube Time:", format(time() - cubeTime,'.0f'), "seconds\n")

    # for each region find the 101 nearest neighbors to be used as recommendations
    knnTime = time()
    knn(outputDir, cube, pd.DataFrame(coords.values), nJobs, nDesiredNeighbors)
    print("KNN time:", format(time() - knnTime,'.0f'), "seconds\n")


def determineBlockSize(windowKB):
    """
    The similarity search reduces epilogos scores over search windows to mimic what the human eye might focus on.
    This is meant to provide regions which contain similar characteristics/trends

    Input:
    windowKB - The size of the search window in KB

    Output:
    blockSize - The reduction factor for the similarity search
    """

    # We only support some set window sizes, these each have a reduction factor to
    # create an overall data size of 25 points per window
    if windowKB == 5:
        blockSize = 1
    elif windowKB == 10:
        blockSize = 2
    elif windowKB == 25:
        blockSize = 5
    elif windowKB == 50:
        blockSize = 10
    elif windowKB == 75:
        blockSize = 15
    elif windowKB == 100:
        blockSize = 20
    else:
        raise ValueError("Error: window size must be either 5, 10, 25, 50, 75, or 100 (in kb)")

    return blockSize


def readScores(epilogosScoresPath):
    """
    Reads in the epilogos scores file provided by the user

    Input:
    epilogosScoresPath - The path to the epilogos scores file provided by the user

    Output:
    epilogosScores - Just the scores from the epilogos scores file
    inputArr - A numpy array to be read in by the filter-regions package (contains coords and summed scores)
    """

    # Read in the scores file for exemplar generation
    epilogosScores = pd.read_table(epilogosScoresPath, sep="\t", header=None)

    # Sum up the scores for all states for maxMean algorithm
    inputDF = epilogosScores.iloc[:, :3]
    inputDF[3] = epilogosScores.iloc[:,3:].sum(axis=1)
    inputDF.columns = ["Chromosome", "Start", "End", "Score"]
    inputArr = inputDF.to_numpy()

    # Pull out just the scores for later window reduction
    epilogosScores = epilogosScores.iloc[:,3:]

    return epilogosScores, inputArr


def makeSlice(genome, idx, windowBins, blockSize):
    """
    Creates a slice of a reduced epilogos track used for construction of the epilogos "cube"
    This will take a location and use it as a center point for each "slice" (region) of the cube

    Input:
    genome - The epilogos scores across the whole genome
    idx - The centerpoint of a similarity search window
    windowBins - The number of bins to have in each window
    blockSize - The reduction factor for the similiarty search

    Output:
    genomeWindow.loc[reducedIndices[0].values].to_numpy() - A numpy array representing the reduced epilogos scores of a window
    """
    genomeWindow = genome.iloc[idx - windowBins // 2:idx + windowBins // 2 + 1] if windowBins % 2 else genome.iloc[idx - windowBins // 2:idx + windowBins // 2]
    sumsOverWindow = genomeWindow.sum(axis=1).to_frame()
    reducedIndices = sumsOverWindow.groupby(np.arange(len(sumsOverWindow), dtype=np.float32)//blockSize,).idxmax()
    return genomeWindow.loc[reducedIndices[0].values].to_numpy()


def knn(outputDir, cube, locs, nJobs, nDesiredNeighbors):
    """
    Finds KNearestNeighbors that will be used as recommendations

    Input:
    outputDir - The output directory for the results
    cube - The reduced epilogos scores over the mean max generated windows
    locs - Pandas Dataframe of the coordinates for the mean max generated windows
    nJobs - The number of jobs to use for the kNearestNeighbors algorithm
    nDesiredNeighbors - Number of neighbors to be found by kNearestNeighbors (note the 1st neighbor is the query region itself)

    Output:
    epilogos_knn.npz -- for each of regions, the top nDesiredNeighbors nearest neighbor regions
    recommendations.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredNeighbors nearest neighbors for each of regions
    """
    np.set_printoptions(threshold=sys.maxsize)

    samples, length, states = cube.shape

    # flattens the cube for NearestNeighbors algorithm
    cube = cube.reshape(samples, length * states)

    neighbors = NearestNeighbors(n_neighbors=nDesiredNeighbors, n_jobs=nJobs).fit(cube)

    # retuns array of distances and indices of the nearest neighbors
    dist, idx = neighbors.kneighbors(cube)

    # number of regions, number of nearest neighbors (should equal n_desired_neighbors)
    nRegions, nNearestNeighbors = idx.shape

    # retrieves the locations for each recommendation from the original dataframe
    res = locs.iloc[idx.reshape(nRegions * nNearestNeighbors), :3].values

    # reshapes results to be 3D array of k (101) [chr, start,end] coordinates
    res = res.reshape(nRegions, nNearestNeighbors, 3)

    # saved in order: Highest scoring -> Lowest scoring
    np.savez_compressed(outputDir / "epilogos_knn.npz", arr=res, idx=idx, dist=dist)

    i = 0
    final = ['' for i in range(nRegions)]
    for row in res:
        query_v = locs.iloc[i,].values
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
        # if the query region was not added to the beginning of the final list, we push it at the front here and remove a trailing element
        if not query_added:
            recs = [query] + recs[:-1]
        assert(len(recs) == nDesiredNeighbors)
        final[i] = json.dumps(recs)
        i += 1

    tbx = pd.DataFrame(locs)
    tbx[3] = final
    tbx.columns = ['chrom', 'start', 'stop', 'recommendations']
    tbx = tbx.sort_values(by=['chrom','start'])

    recommendations_fn = os.path.join(outputDir, "recommendations.bed.gz")
    recommendations_idx_fn = os.path.join(outputDir, "recommendations.bed.gz.tbi")

    if os.path.exists(recommendations_fn): os.remove(recommendations_fn)
    if os.path.exists(recommendations_idx_fn): os.remove(recommendations_idx_fn)

    # save as tabix
    temp_csv_fn = None
    with tempfile.NamedTemporaryFile(mode='w+b', delete=False, dir=os.path.dirname(os.path.realpath(__file__))) as temp_fh:
        tbx.to_csv(temp_fh, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)
        temp_fh.seek(0)
        temp_csv_fn = temp_fh.name

    pysam.tabix_compress(temp_csv_fn, recommendations_fn, force=True)
    if not os.path.exists(recommendations_fn) or os.stat(recommendations_fn).st_size == 0:
        raise Exception("Error: Could not create bgzip archive [{}]".format(recommendations_fn))

    pysam.tabix_index(recommendations_fn, force=True, zerobased=True, preset="bed")
    if not os.path.exists(recommendations_idx_fn) or os.stat(recommendations_idx_fn).st_size == 0:
        raise Exception("Error: Could not create index of bgzip archive [{}]".format(recommendations_idx_fn))

    if os.path.exists(temp_csv_fn): os.remove(temp_csv_fn)


if __name__ == "__main__":
    main()