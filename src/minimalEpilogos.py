from multiprocessing import cpu_count
from pathlib import Path
from pathlib import PurePath
import click
import computeEpilogosExpectedMaster
import computeEpilogosExpectedCombination
import computeEpilogosScoresMaster
import computeEpilogosGreatestHits
import errno

print("""\n
                  d8b 888                                     
                  Y8P 888                                     
                      888                                     
 .d88b.  88888b.  888 888  .d88b.   .d88b.   .d88b.  .d8888b  
d8P  Y8b 888 "88b 888 888 d88""88b d88P"88b d88""88b 88K      
88888888 888  888 888 888 888  888 888  888 888  888 "Y8888b. 
Y8b.     888 d88P 888 888 Y88..88P Y88b 888 Y88..88P      X88 
 "Y8888  88888P"  888 888  "Y88P"   "Y88888  "Y88P"   88888P' 
         888                            888                   
         888                       Y8b d88P                   
         888                        "Y88P"                    
""")

@click.command()
@click.option("-i", "--input-directory", "inputDirectory", type=str, required=True, multiple=True, help="Path to directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, multiple=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-n", "--state-model", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency-level", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1, 2, or 3)")
@click.option("-m", "--mode-of-operation", "modeOfOperation", type=click.Choice(["bg", "s", "both"]), default=["both"], show_default=True, multiple=True, help="bg for background, s for scores, both for both")
@click.option("-b", "--background-directory", "expFreqDir", type=str, default=["null"], multiple=True, help="Path to where the background frequency array is read from (-m s) or written to (-m bg, -m both) [default: output-directory]")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-v", "--verbose", "verbose", is_flag=True, multiple=True, help="If flag is enabled, program will print verbose outputs on current processes")
def main(inputDirectory, outputDirectory, numStates, saliency, modeOfOperation, expFreqDir, numProcesses, verbose):
    """
    This script computes scores for chromatin states across the genome.
    """

    # Handling case if user inputs flag multiples times
    if len(inputDirectory) > 1:
        raise ValueError("Too many [-f, --file-directory] arguments provided")
    elif len(outputDirectory) > 1:
        raise ValueError("Too many [-o, --output-directory] arguments provided")
    elif len(numStates) > 1:
        raise ValueError("Too many [-s, --state-model] arguments provided")
    elif len(saliency) > 1:
        raise ValueError("Too many [-l, --saliency-level] arguments provided")
    elif len(modeOfOperation) > 1:
        raise ValueError("Too many [-m, --mode-of-operation] arguments provided")
    elif len(expFreqDir) > 1:
        raise ValueError("Too many [-b, --background-directory] arguments provided")
    elif len(numProcesses) > 1:
        raise ValueError("Too many [-c, --num-cores] arguments provided")
    elif len(verbose) > 1:
        raise ValueError("Too many [-v, --verbose] arguments provided")
    inputDirectory = inputDirectory[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    modeOfOperation = modeOfOperation[0]
    expFreqDir = expFreqDir[0]
    numProcesses = numProcesses[0]
    if verbose:
        verbose = verbose[0]
        
    inputDirPath = Path(inputDirectory)
    outputDirPath = Path(outputDirectory)

    print()
    print("Input Directory =", inputDirPath)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)
    print("Mode of Operation =", modeOfOperation)
    if expFreqDir == "null":
        print("Background Directory =", outputDirPath)
    else:
        print("Background Directory =", expFreqDir)
    if numProcesses == 0:
        print("Number of Cores = All available")
    else:
        print("Number of Cores =", numProcesses)

    # If user does not specificy a directory to look for expected frequencies default the output directory
    if expFreqDir == "null":
        expFreqDir = outputDirectory

    # Making paths absolute
    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath
    if not PurePath(inputDirPath).is_absolute():
        inputDirPath = Path.cwd() / inputDirPath

    # For making sure all files are consistently named
    fileTag = "{}_s{}".format(inputDirPath.name, saliency)

    if saliency != 1 and saliency != 2 and saliency != 3:
        raise ValueError("Please ensure that saliency metric is either 1, 2, or 3")

    # Check that paths are valid before doing anything
    if not inputDirPath.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath)))
    if not inputDirPath.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath)))
    if not list(inputDirPath.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty:", str(inputDirPath))

    # If the output directory does not exist yet, make it for the user 
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)
    
    if not outputDirPath.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDirPath)))

    if numProcesses < 0:
        raise ValueError("Number of cores must be positive or zero (0 means use all cores)")
    elif numProcesses == 0:
        numProcesses = cpu_count()

    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to path of the input file directory
    storedExpPath = Path(expFreqDir) / "exp_freq_{}.npy".format(fileTag)

    # Only calculate the expected frequencies if user asks for it, otherwise just load from where the user said
    if modeOfOperation != "s":
        print("\nSTEP 1: Per data file background frequency calculation", flush=True)
        for file in inputDirPath.glob("*"):
            if file.name.split(".")[1] == "genome":
                continue
            if not file.is_dir():
                computeEpilogosExpectedMaster.main(file, "null", numStates, saliency, outputDirPath, fileTag, numProcesses, verbose)

        print("\nSTEP 2: Background frequency combination", flush=True)
        computeEpilogosExpectedCombination.main(outputDirPath, storedExpPath, fileTag, verbose)

    if modeOfOperation == "s" or modeOfOperation == "both":
        print("\nSTEP 3: Score calculation", flush=True)
        for file in inputDirPath.glob("*"):
            if file.name.split(".")[1] == "genome":
                continue
            if not file.is_dir():
                computeEpilogosScoresMaster.main(file, "null", numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses, verbose)
                computeEpilogosGreatestHits.main(outputDirPath, numStates, fileTag, numProcesses)

if __name__ == "__main__":
    main()