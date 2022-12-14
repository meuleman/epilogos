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

    print("    Time:", format(time() - reductionTime,'.0f'), "seconds\n", flush=True)
    print("Reading in search results...", flush=True); readTime = time()

    simsearch_npz = np.load(outputDir / "simsearch_cube.npz", allow_pickle=True)
    nRegions = simsearch_npz["scores"].shape[0]
    roiCoords = pd.DataFrame(simsearch_npz["coords"], columns=["Chromosome", "Start", "End"])

    simsearch_arr = readSimsearchIndices(outputDir, nRegions, nDesiredMatches, nJobs)

    searchResults = convertIndicesToCoords(simsearch_arr, reducedGenomeCoords, roiCoords, windowBins, blockSize, nRegions, nDesiredMatches)

    print("    Time:", format(time() - readTime,'.0f'), "seconds\n", flush=True)
    print("Writing search results...", flush=True); writeTime = time()

    writeResults(outputDir, searchResults, roiCoords, nRegions, nDesiredMatches)

    print("    Time:", format(time() - writeTime,'.0f'), "seconds\n", flush=True)
    print("Cleaning up temp files...", flush=True); cleanTime = time()

    cleanUpTempFiles(outputDir)

    print("    Time:", format(time() - cleanTime,'.0f'), "seconds\n", flush=True)
    print("Total time:", format(time() - t,'.0f'), "seconds\n", flush=True)


def reduceGenomeCoords(outputDir, blockSize):
    genome_stats = np.load(outputDir / "genome_stats.npz", allow_pickle=True)
    stateScores = genome_stats["scores"]
    genomeCoords = pd.DataFrame(genome_stats["coords"], columns=["Chromosome", "Start", "End"])
    sumsOverGenome = pd.DataFrame(stateScores).sum(axis=1).to_frame(name="scores")
    sumsOverGenome["blockIndices"] = np.arange(len(sumsOverGenome), dtype=np.int32) // blockSize

    firstIndices = np.array(sorted(sumsOverGenome.drop_duplicates(['blockIndices'], keep='first').index))
    lastIndices = np.array(sorted(sumsOverGenome.drop_duplicates(['blockIndices'], keep='last').index))
    return pd.concat((genomeCoords.iloc[firstIndices, :2].reset_index(drop=True), genomeCoords.iloc[lastIndices, 2].reset_index(drop=True)), axis=1, names=["Chromsome", "Start", "End"])


def readSimsearchIndices(outputDir, nRegions, nDesiredMatches, nJobs):
    simsearch_arr = np.zeros((nRegions, nDesiredMatches), dtype=np.int32)

    rowList = splitRows(nRegions, nJobs)
    for file in outputDir.glob("simsearch_regions_*.npy"):
        i = int(file.stem.split("_")[-1])
        simsearch_arr[rowList[i][0]:rowList[i][1]] = np.load(file, allow_pickle=True)

    return simsearch_arr


def convertIndicesToCoords(simsearch_arr, reducedGenomeCoords, roiCoords, windowBins, blockSize, nRegions, nDesiredMatches):
    # take the shared array and convert all of the indices into locations
    resultsChrAndStart = reducedGenomeCoords.iloc[simsearch_arr.flatten(), :2].values.reshape(nRegions * nDesiredMatches, 2)
    resultsEnd = reducedGenomeCoords.iloc[simsearch_arr.flatten() + windowBins // blockSize - 1, 2].values.reshape(nRegions * nDesiredMatches, 1)

    searchResults = np.concatenate((resultsChrAndStart, resultsEnd), axis=1).reshape(nRegions, nDesiredMatches, 3)

    # Store the lcoation of the query at the beginning of array
    return np.concatenate((roiCoords.values.reshape((1,nRegions,3)), searchResults), axis=1)


def writeResults(outputDir, searchResults, roiCoords, nRegions, nDesiredMatches):
    i = 0
    final = ['' for i in range(nRegions)]
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


def cleanUpTempFiles(outputDir):
    os.remove(outputDir / "genome_coords.npy")
    os.remove(outputDir / "reduced_genome.npy")
    for file in outputDir.glob("simsearch_regions_*.npy"):
        os.remove(file)


if __name__ == "__main__":
    main(Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]))