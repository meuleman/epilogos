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
from epilogos.helpers import splitRows


def main(outputDir, windowBins, blockSize, nJobs, nDesiredMatches):
    print("Reducing genome coordinates...", flush=True); reductionTime = time(); t = time()
    reducedGenomeCoords = reduceGenomeCoords(outputDir, blockSize)

    print("    Time:", format(time() - reductionTime, '.0f'), "seconds\n", flush=True)
    print("Reading in search results...", flush=True); readTime = time()

    simsearch_npz = np.load(outputDir / "simsearch_cube.npz", allow_pickle=True)
    nRegions = simsearch_npz["scores"].shape[0]
    roiCoords = pd.DataFrame(simsearch_npz["coords"], columns=["Chromosome", "Start", "End"])

    simsearchArr = readSimsearchIndices(outputDir, nRegions, nDesiredMatches, nJobs)

    searchResults = convertIndicesToCoords(simsearchArr, reducedGenomeCoords, roiCoords, windowBins, blockSize,
                                           nRegions, nDesiredMatches)

    print("    Time:", format(time() - readTime, '.0f'), "seconds\n", flush=True)
    print("Writing search results...", flush=True); writeTime = time()

    writeResults(outputDir, searchResults, simsearchArr, roiCoords, nRegions)

    print("    Time:", format(time() - writeTime, '.0f'), "seconds\n", flush=True)
    print("Cleaning up temp files...", flush=True); cleanTime = time()

    cleanUpFiles(outputDir, simsearchArr)

    print("    Time:", format(time() - cleanTime, '.0f'), "seconds\n", flush=True)
    print("Total time:", format(time() - t, '.0f'), "seconds\n", flush=True)


def reduceGenomeCoords(inputDir, blockSize):
    """
    Reduces the entire genome by a factor of blockSize

    Input:
    inputDir  -- The directory in which to read the unreduced genome from
    blockSize -- The reduction factor

    Output:
    The coordinates of the reduced genome
    """
    genome_stats = np.load(inputDir / "genome_stats.npz", allow_pickle=True)
    stateScores = genome_stats["scores"]
    genomeCoords = pd.DataFrame(genome_stats["coords"], columns=["Chromosome", "Start", "End"])
    sumsOverGenome = pd.DataFrame(stateScores).sum(axis=1).to_frame(name="scores")
    sumsOverGenome["blockIndices"] = np.arange(len(sumsOverGenome), dtype=np.int32) // blockSize

    firstIndices = np.array(sorted(sumsOverGenome.drop_duplicates(['blockIndices'], keep='first').index))
    lastIndices = np.array(sorted(sumsOverGenome.drop_duplicates(['blockIndices'], keep='last').index))
    return pd.concat((genomeCoords.iloc[firstIndices, :2].reset_index(drop=True),
                      genomeCoords.iloc[lastIndices, 2].reset_index(drop=True)),
                     axis=1, names=["Chromsome", "Start", "End"])


def readSimsearchIndices(inputDir, nRegions, nDesiredMatches, nJobs):
    """
    Reads in the similarity search results in index form (indices of the reduced_genome.npy file)

    Input:
    inputDir        -- The directory from which to read the simsearch results
    nRegions        -- The number of regions which results were calculated for
    nDesiredMatches -- The number of matches for each region
    nJobs           -- The number of slurm jobs used to calculate results

    Output:
    simsearchArr -- The combined results of simsearch indices across all slurm jobs
    """
    simsearchArr = np.zeros((nRegions, nDesiredMatches), dtype=np.int32)

    rowList = splitRows(nRegions, nJobs)
    for file in inputDir.glob("simsearch_indices_*.npy"):
        i = int(file.stem.split("_")[-1])
        simsearchArr[rowList[i][0]:rowList[i][1]] = np.load(file, allow_pickle=True)

    return simsearchArr


def convertIndicesToCoords(simsearchArr, reducedGenomeCoords, roiCoords, windowBins, blockSize, nRegions,
                           nDesiredMatches):
    """
    Takes the similarity search results in index form and converts them to genomic coordinates

    Input:
    simsearchArr        -- The similarity search results in index form
    reducedGenomeCoords -- The coordinates of the reduced genome
    roiCoords           -- The coordinates of the regions for which results were calculated
    windowBins          -- The number of unreduced bins to have in each window
    blockSize           -- The reduction factor for the similarity search
    nRegions            -- The number of regions which results were calculated for
    nDesiredMatches     -- The maximum number of matches for each region

    Output:
    The similarity search results in coordinate form with the coordinates of the query region at the start of each row
    """
    # take the shared array and convert all of the indices into locations
    resultsChrAndStart = reducedGenomeCoords.iloc[simsearchArr.flatten(), :2].values\
        .reshape(nRegions * nDesiredMatches, 2)
    resultsEnd = reducedGenomeCoords.iloc[simsearchArr.flatten() + windowBins // blockSize - 1, 2].values\
        .reshape(nRegions * nDesiredMatches, 1)

    searchResults = np.concatenate((resultsChrAndStart, resultsEnd), axis=1).reshape(nRegions, nDesiredMatches, 3)

    # Store the lcoation of the query at the beginning of array
    return np.concatenate((roiCoords.values.reshape((nRegions, 1, 3)), searchResults), axis=1)


def writeResults(outputDir, searchResults, simsearchArr, roiCoords, nRegions):
    """
    Written by Alex Reynolds - https://github.com/alexpreynolds

    Writes the similarity search results to gzipped bed and tabix files

    Input:
    outputDir     -- The output directory for the results
    searchResults -- The similarity search results in coordinate form (including the coordinates of the query region)
    simsearchArr  -- The similarity search results in index form
    roiCoords     -- The coordinates of the regions for which results were calculated
    nRegions      -- The number of regions which results were calculated for

    Output:
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredMatches matches for each of regions
    """
    i = 0
    final = ['' for i in range(nRegions)]
    for resultsRow, arrRow in zip(searchResults,
                                  np.concatenate((np.ones((len(simsearchArr), 1)), simsearchArr), axis=1)):
        recs = []
        for chrom, start, end in resultsRow[np.where(arrRow != -1)[0]]:
            hit = '{}:{}:{}'.format(chrom, start, end)
            recs.append(hit)
        final[i] = json.dumps(recs)
        i += 1

    tbx = pd.DataFrame(roiCoords)
    tbx[3] = final
    tbx.columns = ['chrom', 'start', 'stop', 'matches']
    tbx = tbx.sort_values(by=['chrom', 'start'])

    simsearch_fn = os.path.join(outputDir, "simsearch.bed.gz")
    simsearch_idx_fn = os.path.join(outputDir, "simsearch.bed.gz.tbi")

    if os.path.exists(simsearch_fn): os.remove(simsearch_fn)
    if os.path.exists(simsearch_idx_fn): os.remove(simsearch_idx_fn)

    # save as tabix
    temp_csv_fn = None
    with tempfile.NamedTemporaryFile(mode='w+b', delete=False, dir=os.path.dirname(os.path.realpath(__file__)))\
            as temp_fh:
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


def cleanUpFiles(outputDir, simsearchArr):
    """
    Removes and combines some of the files used for similarity search calculatio

    Input:
    outputDir    -- The output directory for the results
    simsearchArr -- The similarity search results in index form
    """
    os.remove(outputDir / "genome_stats.npz")
    for file in outputDir.glob("simsearch_indices_*.npy"):
        os.remove(file)
    np.save(outputDir / "simsearch_indices.npy", simsearchArr, allow_pickle=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))
