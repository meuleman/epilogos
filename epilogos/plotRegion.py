import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import numpy.lib.recfunctions as nlr
import click
from epilogos.helpers import getStateNames, getStateColorsRGB, generateRegionArr


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-r", "--regions", "regions", type=str,
              help="Region formatted as chr:start-end or path to bed file containing regions to visualize")
@click.option("-s", "--scores-file", "epilogosScoresPath", type=str,
              help="Path to epilogos scores file to be used for region visualization")
@click.option("-j", "--state-info", "metadataPath", type=str,
              help="Path to state metadata file to be used for region coloring")
@click.option("-o", "--output-directory", "outputDir", type=str, help="Path to desired output directory")
@click.option("-y", "--individual-ylims", "individualYlims", is_flag=True,
              help="If true each region is plotted on its own y-axis")
@click.option("-f", "--file-format", "fileFormat", type=str, default="pdf",
              help="File format for the output images [default: pdf]")
def main(regions, epilogosScoresPath, metadataPath, outputDir, individualYlims, fileFormat):
    print("""\n
                 888          888                              d8b
                 888          888                              Y8P
                 888          888
        88888b.  888  .d88b.  888888  888d888 .d88b.   .d88b.  888  .d88b.  88888b.
        888 "88b 888 d88""88b 888     888P"  d8P  Y8b d88P"88b 888 d88""88b 888 "88b
        888  888 888 888  888 888     888    88888888 888  888 888 888  888 888  888
        888 d88P 888 Y88..88P Y88b.   888    Y8b.     Y88b 888 888 Y88..88P 888  888
        88888P"  888  "Y88P"   "Y888  888     "Y8888   "Y88888 888  "Y88P"  888  888
        888                                                888
        888                                           Y8b d88P
        888                                            "Y88P"
    """, flush=True)

    print("\n\n\nReading in data...", flush=True)
    # Read in regions
    regionArr = generateRegionArr(regions)

    # Read in epilogos scores
    epilogosScores = pd.read_table(Path(epilogosScoresPath), sep="\t", header=None)

    # Determine state names
    stateNames = getStateNames(metadataPath)

    # Determine state colors
    stateColors = getStateColorsRGB(metadataPath)

    print("Processing region scores...", flush=True)
    # Process query scores for graphing
    processedScoresList = []
    processedColorsList = []
    for chr, start, end in regionArr:
        processedRegion = processEpilogosScoresForDrawing(chr, start, end, epilogosScores, stateColors)
        processedScoresList.append(processedRegion[0])
        processedColorsList.append(processedRegion[1])

    # Precalculate ylims for drawing if using same axes
    ymin, ymax = ylim(processedScoresList) if individualYlims else (np.nan, np.nan)

    # Strip period off file format
    if fileFormat.startswith('.'): fileFormat = fileFormat[1:]

    # Draw the query regions
    "Drawing regions..."
    for i, (chr, start, end) in enumerate(regionArr):
        print("\tRegion {}...".format(i+1))
        drawEpilogosScores(chr, start, end, processedScoresList[i], processedColorsList[i], ymin, ymax, stateNames,\
                           stateColors, Path(outputDir) / "epilogos_region_{}_{}_{}.{}".format(chr, start, end, fileFormat))


def processEpilogosScoresForDrawing(chr, start, end, epilogosScores, stateColors):
    """
    Takes in a region, genome wide epilogos scores, and state colors and outputs a numpy array containing scores for the
    region with each bin independently sorted by its scores.
    It also generates a corresponding state color array for the score array

    Input:
    chr            -- The region's chromosome
    start          -- The region's start position on the chromosome (in bp)
    end            -- The region's end position on the chromosome (in bp)
    epilogosScores -- Numpy array containing epilogos scores across the whole genome
    stateColors    -- Numpy array containing the rgba colors for each state

    Output:
    regionScoresSorted -- Numpy array containing scores across the region with each bin individually sorted by scores
    regionColorsSorted -- Numpy array containing state colors across the region with each bin individually sorted by scores
    """
    # Find the epilogos scores for the region
    scoresStartIndex = np.where((epilogosScores.iloc[:,0] == chr) & (epilogosScores.iloc[:,1] == start))[0][0]
    scoresEndIndex = np.where((epilogosScores.iloc[:,0] == chr) & (epilogosScores.iloc[:,1] == end))[0][0]
    regionScores = epilogosScores.iloc[scoresStartIndex:scoresEndIndex, 3:].to_numpy(dtype=np.float64).T

    # Generate state color array for the region
    state_colors_2d = nlr.unstructured_to_structured(\
                          np.swapaxes(np.array([stateColors for i in range(regionScores.shape[1])]), 0, 1)).astype('O')

    # Sort each bin by the scores
    sortedScoresIndices = np.argsort(regionScores, axis=0)
    regionScoresSorted  = np.take_along_axis(regionScores, sortedScoresIndices, axis=0)
    regionColorsSorted  = np.take_along_axis(state_colors_2d, sortedScoresIndices, axis=0)

    return regionScoresSorted, regionColorsSorted


def drawEpilogosScores(chr, start, end, scores, colors, ymin, ymax, stateNames, stateColors, file):
    """
    Draws the epilogos scores for a region. State scores are drawn with lowest scores at the bottom and highest at the top.
    This sorting is done individually for each bin

    Input:
    chr         -- The region's chromosome
    start       -- The region's start position on the chromosome (in bp)
    end         -- The region's end position on the chromosome (in bp)
    scores      -- Numpy array containing scores across the region with each bin individually sorted by scores
    colors      -- Numpy array containing state colors across the region with each bin individually sorted by scores
    ymin        -- The minimum lower bound for the y-axis
    ymax        -- The minimum upper bound for the y-axis
    stateNames  -- Numpy array containing names of states
    stateColors -- Numpy array containg colors of states
    file        -- Name for the output file

    Output:
    A image containing the epilogos scores drawn
    """
    # create the bar chart
    fig, axs = plt.subplots(1, 1, figsize=(24,5))

    # Calculate individual ylims if generic ylims haven't been calculated
    if np.isnan(ymin) and np.isnan(ymax):
        ymin, ymax = ylim([scores])

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


def ylim(regionScoresList):
    """
    Determines what the bounds should be for all the graphs so that they are all graphed on the same scale

    Input:
    regionScoresList -- List of Numpy arrays containg epilogos scores over queried regions

    Output:
    minScore -- The minimum lower bound for the y axis
    maxScore -- The minimum upper bound for the y axis
    """

    minScore = np.finfo(np.float32).max
    maxScore = np.finfo(np.float32).min

    # Compare the max/min ylims needed for each region
    # Evaluate negative and positive scores separately due to drawing stracks in opposite directions
    for regionScores in regionScoresList:
        negativeScores = regionScores.copy()
        negativeScores[negativeScores > 0] = 0
        minCandidate = min(np.sum(negativeScores, axis=0))
        minScore = minCandidate if minCandidate < minScore else minScore

        positiveScores = regionScores.copy()
        positiveScores[positiveScores < 0] = 0
        maxCandidate = max(np.sum(positiveScores, axis=0))
        maxScore = maxCandidate if maxCandidate > maxScore else maxScore

    return minScore, maxScore

if __name__ == "__main__":
    main()