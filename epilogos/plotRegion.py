import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import numpy.lib.recfunctions as nlr
from time import time
import click
from epilogos.helpers import getStateNames, getStateColorsRGB, generateRegionArr


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-r", "--regions", "regions", type=str, required=True,
              help="Region formatted as chr:start-end or path to bed file containing regions to visualize")
@click.option("-s", "--scores-file", "epilogosScoresPath", type=str, default="",
              help="Path to epilogos scores file to be used for region visualization")
@click.option("-a", "--scores-a", "scoresPathGroupA", type=str, default="",
              help="Path to epilogos scores file to be used for region visualization")
@click.option("-b", "--scores-b", "scoresPathGroupB", type=str, default="",
              help="Path to epilogos scores file to be used for region visualization")
@click.option("-c", "--scores-diff", "scoresDiffPath", type=str, default="",
              help="Path to epilogos scores file to be used for region visualization")
@click.option("-j", "--state-info", "metadataPath", type=str, required=True,
              help="Path to state metadata file to be used for region coloring")
@click.option("-o", "--output-directory", "outputDir", type=str, required=True, help="Path to desired output directory")
@click.option("-y", "--individual-ylims", "individualYlims", is_flag=True,
              help="If true each region is plotted on its own y-axis")
@click.option("-f", "--file-format", "fileFormat", type=str, default="pdf",
              help="File format for the output images [default: pdf]")
def main(regions, epilogosScoresPath, scoresPathGroupA, scoresPathGroupB, scoresDiffPath, metadataPath, outputDir,
         individualYlims, fileFormat):
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

    outputDir = Path(outputDir)
    if not outputDir.exists():
        outputDir.mkdir(parents=True)
    if not outputDir.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDir)))

    print("\n\n\n        Reading in data...", flush=True); readTime = time()
    # Read in regions
    regionArr = generateRegionArr(regions)

    # Determine state names
    stateNames = getStateNames(metadataPath)

    # Determine state colors
    stateColors = getStateColorsRGB(metadataPath)

    if epilogosScoresPath:
        plotOneTrackRegion(epilogosScoresPath, regionArr, outputDir, stateColors, stateNames, individualYlims,
                           fileFormat, readTime)
    elif scoresPathGroupA and scoresPathGroupB:
        plotMultiTrackRegion(scoresPathGroupA, scoresPathGroupB, scoresDiffPath, regionArr, outputDir, stateColors,
                             stateNames, individualYlims, fileFormat, readTime)
    else:
        raise ValueError("Missing scores file input")


def plotOneTrackRegion(epilogosScoresPath, regionArr, outputDir, stateColors, stateNames, individualYlims, fileFormat,
                       readTime):
    """
    Wrapper function for plotting single epilogos regions

    Input:
    epilogosScoresPath -- Path to the single epilogos scores file
    regionArr          -- 2D numpy array containing coordinates (chr, start, end) of all the regions to plot
    outputDir          -- Path to the directory to output drawings to
    stateColors        -- Numpy array of colors for each state in rgb format
    stateNames         -- Numpy array of names for each state
    individualYlims    -- Boolean determining whether plots will share the same y-axis or not
    fileFormat         -- String of file format for the drawings
    readTime           -- time.time() value used for outputing time information to users

    Output:
    matplotlib graphs of single epilogos scores for each of the user input regions to the output directory
    """
    # Read in epilogos scores
    epilogosScores = pd.read_table(Path(epilogosScoresPath), sep="\t", header=None)

    print("            Time:", format(time() - readTime, '.0f'), "seconds\n", flush=True)
    print("        Processing region scores...", flush=True); processTime = time()

    # Process query scores for graphing
    processedScoresList, processedColorsList = [], []
    for chr, start, end in regionArr:
        processedRegion = processEpilogosScoresForDrawing(chr, start, end, epilogosScores, stateColors)
        processedScoresList.append(processedRegion[0])
        processedColorsList.append(processedRegion[1])

    # Precalculate ylims for drawing if using same axes
    ymin, ymax = (np.nan, np.nan) if individualYlims else ylim(processedScoresList)

    print("            Time:", format(time() - processTime, '.0f'), "seconds\n", flush=True)
    print("        Drawing regions...", flush=True)

    # Strip period off file format
    if fileFormat.startswith('.'): fileFormat = fileFormat[1:]

    # Draw the query regions
    for i, (chr, start, end) in enumerate(regionArr):
        print("            Region {}...".format(i + 1), flush=True); regionTime = time()
        drawOneTrackEpilogosScores(chr, start, end, processedScoresList[i], processedColorsList[i], ymin, ymax,
                                   stateNames, stateColors,
                                   outputDir / "epilogos_region_{}_{}_{}.{}".format(chr, start, end, fileFormat))
        print("                Time: ", format(time() - regionTime, '.0f'), "seconds\n", flush=True)


def plotMultiTrackRegion(scoresPathGroupA, scoresPathGroupB, scoresDiffPath, regionArr, outputDir, stateColors,
                         stateNames, individualYlims, fileFormat, readTime):
    """
    Wrapper function for plotting single epilogos regions

    Input:
    epilogosPathGroupA -- Path to the group A epilogos scores file
    epilogosPathGroupB -- Path to the group B epilogos scores file
    epilogosDiffPath   -- Path to file containing the differences between the two single epilogos scores
                          if empty string, this difference is calculated from epilogosPathGroupA and epilogosPathGroupB
    regionArr          -- 2D numpy array containing coordinates (chr, start, end) of all the regions to plot
    outputDir          -- Path to the directory to output drawings to
    stateColors        -- Numpy array of colors for each state in rgb format
    stateNames         -- Numpy array of names for each state
    individualYlims    -- Boolean determining whether plots will share the same y-axis or not
    fileFormat         -- String of file format for the drawings
    readTime           -- time.time() value used for outputing time information to users

    Output:
    matplotlib graphs of single epilogos scores for each of the user input regions to the output directory
    """
    # Read in epilogos scores
    epilogosScoresGroupA = pd.read_table(Path(scoresPathGroupA), sep="\t", header=None)
    epilogosScoresGroupB = pd.read_table(Path(scoresPathGroupB), sep="\t", header=None)
    epilogosScoresDiff = pd.read_table(Path(scoresDiffPath), sep="\t", header=None) if scoresDiffPath \
        else pd.concat((epilogosScoresGroupA.iloc[:, :3],
                        epilogosScoresGroupA.iloc[:, 3:] - epilogosScoresGroupB.iloc[:, 3:]), axis=1)

    print("            Time:", format(time() - readTime, '.0f'), "seconds\n", flush=True)
    print("        Processing region scores...", flush=True); processTime = time()

    # Process query scores for graphing
    processedScoresListA, processedColorsListA = [], []
    processedScoresListB, processedColorsListB = [], []
    processedScoresListDiff, processedColorsListDiff = [], []
    for chr, start, end in regionArr:
        processedRegionGroupA = processEpilogosScoresForDrawing(chr, start, end, epilogosScoresGroupA, stateColors)
        processedRegionGroupB = processEpilogosScoresForDrawing(chr, start, end, epilogosScoresGroupB, stateColors)
        processedRegionDiff   = processEpilogosScoresForDrawing(chr, start, end, epilogosScoresDiff, stateColors)

        processedScoresListA.append(processedRegionGroupA[0])
        processedScoresListB.append(processedRegionGroupB[0])
        processedScoresListDiff.append(processedRegionDiff[0])

        processedColorsListA.append(processedRegionGroupA[1])
        processedColorsListB.append(processedRegionGroupB[1])
        processedColorsListDiff.append(processedRegionDiff[1])

    # Precalculate ylims for drawing if using same axes
    ymin, ymax = ylim(processedScoresListA + processedScoresListB + processedScoresListDiff) if individualYlims \
        else (np.nan, np.nan)

    print("            Time:", format(time() - processTime, '.0f'), "seconds\n", flush=True)
    print("        Drawing regions...", flush=True)

    # Strip period off file format
    if fileFormat.startswith('.'): fileFormat = fileFormat[1:]

    # Draw the query regions
    for i, (chr, start, end) in enumerate(regionArr):
        print("            Region {}...".format(i + 1), flush=True); regionTime = time()
        drawMultiTrackEpilogosScores(chr, start, end, processedScoresListA[i], processedColorsListA[i],
                                     processedScoresListB[i], processedColorsListB[i], processedScoresListDiff[i],
                                     processedColorsListDiff[i], ymin, ymax, stateNames, stateColors,
                                     outputDir / "epilogos_region_{}_{}_{}.{}".format(chr, start, end, fileFormat))
        print("                Time: ", format(time() - regionTime, '.0f'), "seconds\n", flush=True)


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
    regionColorsSorted -- Numpy array containing state colors across the region with each bin individually sorted by
                          scores
    """
    # Find the epilogos scores for the region
    scoresStartIndex = np.where((epilogosScores.iloc[:, 0] == chr) & (epilogosScores.iloc[:, 1] == start))[0][0]
    scoresEndIndex = np.where((epilogosScores.iloc[:, 0] == chr) & (epilogosScores.iloc[:, 1] == end))[0][0]
    regionScores = epilogosScores.iloc[scoresStartIndex:scoresEndIndex, 3:].to_numpy(dtype=np.float64).T

    # Generate state color array for the region
    state_colors_2d = nlr.unstructured_to_structured(
        np.swapaxes(np.array([stateColors for i in range(regionScores.shape[1])]), 0, 1)).astype('O')

    # Sort each bin by the scores
    sortedScoresIndices = np.argsort(regionScores, axis=0)
    regionScoresSorted  = np.take_along_axis(regionScores, sortedScoresIndices, axis=0)
    regionColorsSorted  = np.take_along_axis(state_colors_2d, sortedScoresIndices, axis=0)

    return regionScoresSorted, regionColorsSorted


def drawOneTrackEpilogosScores(chr, start, end, scores, colors, ymin, ymax, stateNames, stateColors, file):
    """
    Draws the single group epilogos scores for a region.
    State scores are drawn with lowest scores at the bottom and highest at the top.
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
    fig, axs = plt.subplots(1, 1, figsize=(24, 5))

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

    plotPosNeg(scores, colors, axs)

    # Draw the legend
    custom_lines = [Line2D([0], [0], color=stateColors[i], label=stateNames[i], lw=4) for i in range(len(stateColors))]
    axs.legend(handles=custom_lines, loc='upper left', fontsize=6.25)

    # Draw a zero line
    axs.axhline(0, color='purple', linewidth=0.5)

    fig.savefig(file, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    plt.close()


def drawMultiTrackEpilogosScores(chr, start, end, scoresA, colorsA, scoresB, colorsB, scoresDiff, colorsDiff, ymin,
                                 ymax, stateNames, stateColors, file):
    """
    Draws the epilogos scores for a region. State scores are drawn with lowest scores at the bottom and highest at the
    top. This sorting is done individually for each bin

    Input:
    chr         -- The region's chromosome
    start       -- The region's start position on the chromosome (in bp)
    end         -- The region's end position on the chromosome (in bp)
    scoresA     -- Numpy array containing group A scores across the region with each bin individually sorted by scores
    colorsA     -- Numpy array containing group A state colors across the region with each bin individually sorted by
                   scores
    scoresB     -- Numpy array containing group B scores across the region with each bin individually sorted by scores
    colorsB     -- Numpy array containing group B state colors across the region with each bin individually sorted by
                   scores
    scoresDiff  -- Numpy array containing group A & B scores diff across the region with each bin individually sorted by
                   scores
    colorsDiff  -- Numpy array containing score diff state colors across the region with each bin individually sorted by
                   scores
    ymin        -- The minimum lower bound for the y-axis
    ymax        -- The minimum upper bound for the y-axis
    stateNames  -- Numpy array containing names of states
    stateColors -- Numpy array containg colors of states
    file        -- Name for the output file

    Output:
    A image containing the epilogos scores drawn for group A, group B, and their difference stacked upon one another
    """
    # create the bar chart
    fig, axs = plt.subplots(3, 1, figsize=(24, 15))

    # Calculate individual ylims if generic ylims haven't been calculated
    if np.isnan(ymin) and np.isnan(ymax):
        ymin, ymax = ylim([scoresA, scoresB, scoresDiff])

    # Formatting the graph
    for ax in axs:
        ax.set_ylim(ymin=(ymin * 1.1), ymax=(ymax * 1.1))
        ax.set_facecolor("black")
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.xaxis.set_ticks_position('none')

    axs[0].text(0.99, 0.99, "Group A", verticalalignment='top', horizontalalignment='right', transform=axs[0].transAxes,
                color='w', fontsize=15)
    axs[1].text(0.99, 0.99, "Group B", verticalalignment='top', horizontalalignment='right', transform=axs[1].transAxes,
                color='w', fontsize=15)
    axs[2].text(0.99, 0.99, "Group A vs. Group B", verticalalignment='top', horizontalalignment='right',
                transform=axs[2].transAxes, color='w', fontsize=15)

    axs[2].set_xticks([0, scoresA.shape[1] / 2, scoresA.shape[1]])
    axs[2].set_xticklabels([start, chr, end])

    axs[0].set_title(file.name)

    plt.subplots_adjust(hspace=0.01)

    plotPosNeg(scoresA, colorsA, axs[0])
    plotPosNeg(scoresB, colorsB, axs[1])
    plotPosNeg(scoresDiff, colorsDiff, axs[2])

    # Draw the legend
    custom_lines = [Line2D([0], [0], color=stateColors[i], label=stateNames[i], lw=4) for i in range(len(stateColors))]
    axs[0].legend(handles=custom_lines, loc='upper left', fontsize=6.25)

    # Draw a zero line
    for ax in axs:
        ax.axhline(0, color='purple', linewidth=0.5)

    fig.savefig(file, bbox_inches='tight', dpi=400, facecolor="#FFFFFF", edgecolor="#FFFFFF", transparent=False)

    plt.close()


def plotPosNeg(scores, colors, ax):
    """
    Helper function to plot the positive & negative epilogos scores
    (must be done separately as we don't want to stack them)

    Input:
    scores -- Numpy array containing scores across the region with each bin individually sorted by scores
    colors -- Numpy array containing state colors across the region with each bin individually sorted by scores
    ax     -- matplotlib axis to plot scores on
    """
    x = np.arange(scores.shape[1])

    # Graph the positive values
    # Find the column with the most positives and split the score array based on that
    maxPositives = max(np.count_nonzero(scores > 0, axis=0))
    positiveScores = scores[-maxPositives:]
    positiveColors = colors[-maxPositives:]

    # Set any negative scores to 0 so they don't mess up coloring
    positiveScores[positiveScores < 0] = 0

    # Graph all the positives with the base at 0
    bottoms = np.zeros(scores.shape[1])
    for i in range(maxPositives):
        ax.bar(x, positiveScores[i, :], bottom=bottoms, color=positiveColors[i, :], align='edge', width=1)
        bottoms += positiveScores[i, :].flatten()

    # Graph the negative values
    # Find the column with the most negatives and split the score array based on that
    maxNegatives = max(np.count_nonzero(scores < 0, axis=0))
    negativeScores = scores[:maxNegatives]
    negativeColors = colors[:maxNegatives]

    # Set any positive scores to 0 so they don't mess up coloring
    negativeScores[negativeScores > 0] = 0

    # Graph all the negatives with the base at 0
    # We graph backwards so that the bottoms works out properly
    bottoms = np.zeros(scores.shape[1])
    for i in range(maxNegatives - 1, -1, -1):
        ax.bar(x, negativeScores[i, :], bottom=bottoms, color=negativeColors[i, :], align='edge', width=1)
        bottoms += negativeScores[i, :].flatten()


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
