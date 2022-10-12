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
from sklearn.neighbors import NearestNeighbors
import tempfile
import json
import csv
import pysam
import click
from epilogos.helpers import generateRegionArr, maxMean

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
  -w, --window-KB INTEGER         Window size (in KB) on which to perform\n\
                                  similarity search  [default: 25]\n\
  -j, --num-jobs INTEGER          Number of jobs to be used in nearest\n\
                                  neighbor algorithm  [default: 8]\n\
  -n, --num-matches INTEGER       Number of neighbors to be found by nearest\n\
                                  neighbor algorithm (note that first neighbor\n\
                                  is always the query region)  [default: 101]\n\
\n\
\n\
To Query Similarity Search Data:\n\
  -q, --query TEXT                Query region formatted as chr:start-end or\n\
                                  path to bed file containing query regions\n\
  -m, --matches-file TEXT\n\
                                  Path to previously built\n\
                                  simsearch.bed.gz file to be queried\n\
                                  for matches")

@click.command(context_settings=dict(help_option_names=['-h', '--help']), cls=CustomClickHelpFormat)
@click.option("-b", "--build", "buildBool", is_flag=True,
              help="If true builds the similarity search files needed to query regions")
@click.option("-s", "--scores", "scoresPath", type=str, help="Path to scores file to be used in similarity search")
@click.option("-o", "--output-directory", "outputDir", type=str, help="Path to desired similarity search output directory")
@click.option("-w", "--window-bp", "windowBP", type=int, default=25000, show_default=True,
              help="Window size (in BP) on which to perform similarity search")
@click.option("-j", "--num-jobs", "nJobs", type=int, default=8, show_default=True,
              help="Number of jobs to be used in nearest neighbor algorithm")
@click.option("-n", "--num-matches", "nDesiredNeighbors", type=int, default=101, show_default=True,
              help="Number of matches to be found by nearest neighbor algorithm"+
                   "(note that first match is always the query region)")
@click.option("-q", "--query", "query", type=str,
              help="Query region formatted as chr:start-end or path to bed file containing query regions")
@click.option("-m", "--matches-file", "simSearchPath", type=str,
              help="Path to previously built simsearch.bed.gz file to be queried for matches")
def main(buildBool, scoresPath, outputDir, windowBP, nJobs, nDesiredNeighbors, query, simSearchPath):
    outputDir = Path(outputDir)
    if not outputDir.exists():
        outputDir.mkdir(parents=True)
    if not outputDir.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDir)))

    if buildBool:
        buildSimSearch(scoresPath, outputDir, windowBP, nJobs, nDesiredNeighbors)
    else:
        querySimSearch(query, simSearchPath, outputDir)


def buildSimSearch(scoresPath, outputDir, windowBP, nJobs, nDesiredNeighbors):
    """
    Builds the similarity search files for a set of scores

    Input:
    scoresPath        -- Path to scores file to be used in similarity search
    outputDir         -- pathlib Path to the output directory for similarity search files
    windowBP          -- Size of the similarity search window in bp
    nJobs             -- Number of jobs to be used in nearest neighbor algorithm
    nDesiredNeighbors -- Number of matches to be found by nearest neighbor algorithm (first match is always the query region)

    Output:
    simsearch_cube.npz     -- NPZ file containing all the reduced regions and their coords (scores=scores, coords=coords)
    simsearch_knn.npz      -- NPZ file containing the top nDesiredNeighbors matches for each region
                              (arr=coords, idx=indices, dist=distances)
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredNeighbors matches for each of regions
    """
    # Bins are assumed to be 200bp or 20bp
    if determineBinSize(scoresPath) == 200:
        windowBins = windowBP / 200
        blockSize = determineBlockSize200(windowBP)
    elif determineBinSize(scoresPath) == 20:
        windowBins = windowBP / 20
        blockSize = determineBlockSize20(windowBP)
    else:
        raise ValueError("Similarity Search is only compatible with bins of size 200bp or 20bp")

    cubeTime = time()

    scores, inputArr = readScores(scoresPath)

    # The maximum number of regions chosen should depend on the window size
    # We want to allow for full coverage of the genome if possible (maxRegions is chosen accordingly)
    maxRegions = scores.shape[0] // windowBins

    # Filter-regions package to perform maxmean algorithm & pull out top X regions
    rois, _ = maxMean(inputArr, windowBins, maxRegions)

    # Seperates just the coordinates
    coords = rois.iloc[:,:3]
    # Generate slices of reduced data for the cube
    cube = rois["OriginalIdx"].apply(lambda x: makeSlice(scores, x, windowBins, blockSize))
    cube = np.stack(cube.to_numpy())

    # Save the cube
    np.savez_compressed(file=outputDir / 'simsearch_cube', scores=cube, coords=coords.values)

    print("Cube Time:", format(time() - cubeTime,'.0f'), "seconds\n", flush=True)

    # for each region find the 101 nearest neighbors to be used as matches
    knnTime = time()
    knn(outputDir, cube, pd.DataFrame(coords.values), nJobs, nDesiredNeighbors)
    print("KNN time:", format(time() - knnTime,'.0f'), "seconds\n", flush=True)


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
    # Parse query to generate array
    queryArr = generateRegionArr(query)

    # Read in matches file
    matchesDF = pd.read_table(Path(simSearchPath), sep="\t", header=None)

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
            print("Found region {}:{}-{} within user query {}:{}-{}".format(regionChr, regionStart, regionEnd,
                                                                            chr, start, end))
            print("\tSee {} for matches\n".format(outfile), flush=True)
        else:
            raise ValueError("ERROR: Could not find region in given range")


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


def knn(outputDir, cube, locs, nJobs, nDesiredNeighbors):
    """
    Finds KNearestNeighbors that will be used as matches

    Input:
    outputDir         -- The output directory for the results
    cube              -- The reduced scores over the mean max generated windows
    locs              -- Pandas Dataframe of the coordinates for the mean max generated windows
    nJobs             -- The number of jobs to use for the kNearestNeighbors algorithm
    nDesiredNeighbors -- Number of neighbors to be found by kNearestNeighbors (note 1st neighbor is the query region itself)

    Output:
    simsearch_knn.npz      -- for each of regions, the top nDesiredNeighbors matches (arr=coords, idx=indices, dist=distances)
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredNeighbors matches for each of regions
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

    # retrieves the locations for each match from the original dataframe
    res = locs.iloc[idx.reshape(nRegions * nNearestNeighbors), :3].values

    # reshapes results to be 3D array of k (101) [chr, start,end] coordinates
    res = res.reshape(nRegions, nNearestNeighbors, 3)

    # saved in order: Highest scoring -> Lowest scoring
    np.savez_compressed(outputDir / "simsearch_knn.npz", arr=res, idx=idx, dist=dist)

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
        # if the query region was not added to the beginning of the final list,
        # we push it at the front here and remove a trailing element
        if not query_added:
            recs = [query] + recs[:-1]
        assert(len(recs) == nDesiredNeighbors)
        final[i] = json.dumps(recs)
        i += 1

    tbx = pd.DataFrame(locs)
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