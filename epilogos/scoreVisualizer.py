import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import numpy.lib.recfunctions as nlr
import click
from epilogos.helpers import getStateNames, getStateColorsRGB
from epilogos.similarity_search import generateQueryArr
import re

@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-r", "--regions", "regions", type=str, help="Region formatted as chr:start-end or path to bed file containing regions to visualize")
@click.option("-s", "--scores-file", "epilogosScoresPath", type=str, help="Path to epilogos scores file to be used for region visualization")
@click.option("-j", "--state-info", "metadataPath", type=str, help="Path to state metadata file to be used for region coloring")
@click.option("-o", "--output-directory", "outputDir", type=str, help="Path to desired output directory")
def main(regions, epilogosScoresPath, metadataPath, outputDir):
    # Read in regions
    regionArr = generateQueryArr(regions)

    # Read in epilogos scores
    epilogosScores = pd.read_table(Path(epilogosScoresPath), sep="\t", header=None)

    # Determine state names
    stateNames = getStateNames(metadataPath)

    # Determine state colors
    stateColors = getStateColorsRGB(metadataPath)

    # Draw each of the desired regions
    for chr, start, end in regionArr:
        # Process query scores for graphing
        regionScoresSorted, regionColorsSorted = processEpilogosScoresForDrawing(chr, start, end, epilogosScores, stateColors)

        # Draw the query region
        drawEpilogosScores(chr, start, end, regionScoresSorted, regionColorsSorted, stateNames, stateColors, Path(outputDir) / "epilogos_region_{}_{}_{}.pdf".format(chr, start, end))


# Takes in a region, genome wide epilogos scores, and state colors and outputs a numpy array containing scores for the region
# with each bin independently sorted by its scores. It also generates a corresponding state color array for the score array
def processEpilogosScoresForDrawing(chr, start, end, epilogosScores, stateColors):
    # Find the epilogos scores for the region
    scoresStartIndex = np.where((epilogosScores.iloc[:,0] == chr) & (epilogosScores.iloc[:,1] == start))[0][0]
    scoresEndIndex = np.where((epilogosScores.iloc[:,0] == chr) & (epilogosScores.iloc[:,1] == end))[0][0]
    regionScores = epilogosScores.iloc[scoresStartIndex:scoresEndIndex, 3:].to_numpy(dtype=np.float64).T

    # Generate state color array for the region
    state_colors_2d = nlr.unstructured_to_structured(np.swapaxes(np.array([stateColors for i in range(regionScores.shape[1])]), 0, 1)).astype('O')

    # Sort each bin by the scores
    sortedScoresIndices = np.argsort(regionScores, axis=0)
    regionScoresSorted = np.take_along_axis(regionScores, sortedScoresIndices, axis=0)
    regionColorsSorted = np.take_along_axis(state_colors_2d, sortedScoresIndices, axis=0)

    return regionScoresSorted, regionColorsSorted


# Takes in a region, epilogos scores for that region (each bin is sorted by scores),
#  state colors for the region (sorted same order as scores), and ylims and draws the epilogos graph for the region
def drawEpilogosScores(chr, start, end, scores, colors, stateNames, stateColors, file):
    # create the bar chart
    fig, axs = plt.subplots(1, 1, figsize=(24,5))

    # Determine ylim
    ymax, ymin = ylim(scores)

    # Formatting the graph
    axs.set_ylim(ymin=(ymin * 1.1), ymax=(ymax * 1.1))
    axs.set_xticks([0, scores.shape[1]])
    axs.set_xticklabels([start, end])
    axs.set_xlabel(chr)
    axs.set_facecolor("black")
    axs.set_title(file.name)

    x = np.arange(scores.shape[1])

    # Graph the positive values
    # Find the column with the most positives and split the score array based on that
    max_positives = max(np.count_nonzero(scores > 0, axis=0))
    positive_scores = scores[-max_positives:]
    positive_colors = colors[-max_positives:]

    # Set any negative scores to 0 so they don't mess up coloring
    positive_scores[positive_scores < 0] = 0

    # Graph all the positives with the base at 0
    bottoms = np.zeros(scores.shape[1])
    for i in range(max_positives):
        axs.bar(x, positive_scores[i, :], bottom=bottoms, color=positive_colors[i, :], align='edge', width=1)
        bottoms += positive_scores[i, :].flatten()

    # Graph the negative values
    # Find the column with the most negatives and split the score array based on that
    max_negatives = max(np.count_nonzero(scores < 0, axis=0))
    negative_scores = scores[:max_negatives]
    negative_colors = colors[:max_negatives]

    # Set any positive scores to 0 so they don't mess up coloring
    negative_scores[negative_scores > 0] = 0

    # Graph all the negatives with the base at 0
    # We graph backwards so that the bottoms works out properly
    bottoms = np.zeros(scores.shape[1])
    for i in range(max_negatives - 1, -1, -1):
        axs.bar(x, negative_scores[i, :], bottom=bottoms, color=negative_colors[i, :], align='edge', width=1)
        bottoms += negative_scores[i, :].flatten()

    # Draw the legend
    custom_lines = [Line2D([0], [0], color=stateColors[i], label=stateNames[i], lw=4) for i in range(len(stateColors))]
    axs.legend(handles=custom_lines, loc='upper left', fontsize=6.25)

    # Draw a zero line
    axs.axhline(0, color='purple', linewidth=0.5)

    fig.savefig(file, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)


# Determines what the bounds should be for all the graphs so that they are all graphed on the same scale
def ylim(regionScores):
    # Have to treat positive and negative scores separately because of the way they stack in opposite directions
    positiveScores = regionScores.copy()
    positiveScores[positiveScores < 0] = 0
    max_score = max(np.sum(positiveScores, axis=0))

    negativeScores = regionScores.copy()
    negativeScores[negativeScores > 0] = 0
    min_score = min(np.sum(negativeScores, axis=0))

    return max_score, min_score



