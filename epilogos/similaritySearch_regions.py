import numpy as np
import pandas as pd
import sys
from time import time
from pathlib import Path
from epilogos.helpers import maxMean


def main(outputDir, scoresPath, windowBins, blockSize, windowBP, filterState):
    print("Reading in data...", flush=True); readTime = time(); t = time()
    stateScores, inputArr = readScores(scoresPath)

    np.savez_compressed(outputDir / "genome_stats", scores=stateScores, coords=inputArr[:,:3])

    # The maximum number of regions chosen should depend on the window size
    # We want to allow for full coverage of the genome if possible (maxRegions is chosen accordingly)
    maxRegions = int(stateScores.shape[0] // windowBins)

    print("    Time:", format(time() - readTime,'.0f'), "seconds\n", flush=True)
    print("Finding regions of size {}kb...".format(windowBP // 1000), flush=True); maxmeanTime = time()

    # Filter-regions package to perform maxmean algorithm & pull out top X regions
    rois, _ = maxMean(inputArr, windowBins, maxRegions)

    # Make sure rois are sorted properly for indexing later
    rois.reset_index(drop=True, inplace=True)
    roiCoords = rois.iloc[:,:3]

    print("    Time:", format(time() - maxmeanTime,'.0f'), "seconds\n", flush=True)
    print("Reducing region scores by factor of {}...".format(blockSize), flush=True); reductionTime = time()

    # Generate slices of reduced data for the cube
    roiCube = np.stack(rois["OriginalIdx"].apply(lambda x: makeSlice(stateScores, x, windowBins, blockSize)).to_numpy())

    print("    Time:", format(time() - reductionTime,'.0f'), "seconds\n", flush=True)
    print("Filtering out uninteresting regions...", flush=True); filterTime = ()

    roiCoords, roiCube = removeRegions(roiCoords, roiCube, filterState)

    np.savez_compressed(file=outputDir / 'simsearch_cube', scores=roiCube, coords=roiCoords.values)

    print("    Time:", format(time() - filterTime,'.0f'), "seconds\n", flush=True)
    print("Reducing genome scores by factor of {}...".format(blockSize), flush=True); genomeTime = time()

    reduceGenome(outputDir, stateScores, blockSize)

    print("    Time:", format(time() - genomeTime,'.0f'), "seconds\n", flush=True)
    print("Total time:", format(time() - t,'.0f'), "seconds\n", flush=True)


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


def removeRegions(roiCoords, roiCube, filterState):
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

    return roiCoords, roiCube

def reduceGenome(outputDir, stateScores, blockSize):
    # Reduce the genome for the appropriate region size
    # Give the same blockIndex for every sequential group of blockSize
    sumsOverGenome = stateScores.sum(axis=1).to_frame(name="scores")
    sumsOverGenome["blockIndices"] = np.arange(len(sumsOverGenome), dtype=np.int32) // blockSize
    # For each blockIndex group, output only the index of the bin with the highest overall score. use these indices to get reduced scores
    reducedIndices = np.array(sumsOverGenome.sort_values('scores').drop_duplicates(['blockIndices'], keep='last').sort_values('blockIndices').index)
    reducedGenome = stateScores.iloc[reducedIndices].to_numpy()
    np.save(outputDir / "reduced_genome.npy", reducedGenome, allow_pickle=True)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))


