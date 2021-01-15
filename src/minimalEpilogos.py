import numpy as np
from pathlib import Path
from pathlib import PurePath
import click
import computeEpilogosExpected
import computeEpilogosExpectedCombination
import computeEpilogosScores
import computeEpilogosWrite


@click.command()
@click.option("-f", "--file-directory", "fileDirectory", type=str, required=True, help="Path to directory that contains files to read from (All files in this directory will be read in)")
@click.option("-s", "--state-model", "numStates", type=int, required=True, help="Number of states in chromatin state model")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-l", "--saliency-level", "saliency", type=int, default=1, show_default=True, help="Desired saliency level (1, 2, or 3)")
@click.option("-m", "--mode-of-operation", "modeOfOperation", type=click.Choice(["bg", "s", "both"]), default="both", show_default=True, help="bg for background, s for scores, both for both")
@click.option("-b", "--background-directory", "expFreqDir", type=str, default="null", help="Path to where the background frequency array is read from (-m s) or written to (-m bg, -m both) [default: output-directory]") # default output directory
def main(fileDirectory, numStates, saliency, outputDirectory, modeOfOperation, expFreqDir):
    dataFilePath = Path(fileDirectory)
    outputDirPath = Path(outputDirectory)

    fileTag = "_".join(str(dataFilePath).split("/")[-5:])

    print()
    print("Input Directory =", dataFilePath)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)
    print("Mode of Operation =", modeOfOperation)
    print("Background Directory =", expFreqDir)

    if expFreqDir == "null":
        expFreqDir = outputDirectory

    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath

    if not PurePath(dataFilePath).is_absolute():
        dataFilePath = Path.cwd() / dataFilePath

    if saliency != 1 and saliency != 2 and saliency != 3:
        print()
        print("ERROR: Ensure that saliency metric is either 1, 2, or 3")
        print()
        return

    # Check that paths are valid before doing anything
    if not dataFilePath.exists() or not dataFilePath.is_dir():
        print()
        print("ERROR: Given file path does not exist or is not a directory")
        print()
        return

    if not list(dataFilePath.glob("*")):
        print()
        print("ERROR: Ensure that file directory is not empty")
        print()
        return

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)
    
    if not outputDirPath.is_dir():
        print()
        print("ERROR: Output directory is not a directory")
        print()
        return

    # For slurm output and error later
    (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)

    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to path of the input file directory
    expFreqFilename = "exp_freq_{}.npy".format(fileTag)
    storedExpPath = Path(expFreqDir) / expFreqFilename
    print()
    print("Background Frequency Array Location:", storedExpPath)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    # Check if user wants to calculate it
    if modeOfOperation == "s":
        try:
            expFreqArr = np.load(storedExpPath, allow_pickle=False)
        except IOError:
            print("ERROR: Could not load stored expected value array.\n\tPlease check that the directory is correct or that the file exists")
            return
    else:
        print()
        print("Calculating Per Datafile Background Frequency Arrays...")
        for file in dataFilePath.glob("*"):
            if not file.is_dir():
                computeEpilogosExpected.main(file, numStates, saliency, outputDirPath, fileTag)

        print()
        print("Combining Per Datafile Background Frequency Arrays....")
        computeEpilogosExpectedCombination.main(outputDirPath, fileTag, storedExpPath)

    if modeOfOperation == "s" or modeOfOperation == "both":
        print()
        print("Calculating Per Datafile Scores...")
        for file in dataFilePath.glob("*"):
            if not file.is_dir():
                computeEpilogosScores.main(file, numStates, saliency, outputDirPath, storedExpPath, fileTag)

        print()
        print("Writing to Score Files....")
        for file in dataFilePath.glob("*"):
            if not file.is_dir():
                filename = file.name.split(".")[0]
                computeEpilogosWrite.main(filename, numStates, outputDirPath, fileTag)

    

if __name__ == "__main__":
    main()