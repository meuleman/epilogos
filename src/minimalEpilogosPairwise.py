from multiprocessing import cpu_count
from pathlib import Path
import click
from pathlib import PurePath
import errno
import computeEpilogosExpectedMaster
import computeEpilogosExpectedCombination
import computeEpilogosScoresMaster
import computeEpilogosPairwiseVisual

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
@click.option("-a", "--directory-one", "inputDirectory1", type=str, required=True, multiple=True, help="Path to first directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-b", "--directory-two", "inputDirectory2", type=str, required=True, multiple=True, help="Path to second directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, multiple=True, help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-n", "--number-of-states", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1 or 2)")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-v", "--verbose", "verbose", is_flag=True, multiple=True, help="If flag is enabled, program will print verbose outputs on current processes")
@click.option("-d", "--diagnostic-figures", "diagnosticBool", is_flag=True, multiple=True, help="If flag is enabled, program will output some diagnostic figures of the gennorm fit on the null data and comparisons between the null and real data")
@click.option("-t", "--num-trials", "numTrials", type=int, default=[101], show_default=True, multiple=True, help="The number of times subsamples of the scores are fit")
@click.option("-z", "--sampling-size", "samplingSize", type=int, default=[100000], show_default=True, multiple=True, help="The size of the subsamples on which the scores are fit")
def main(inputDirectory1, inputDirectory2, outputDirectory, numStates, saliency, numProcesses, verbose, diagnosticBool, numTrials, samplingSize):
    """
    This script computes and visualizes differences between chromatin states across epigenomes
    """
    # Handling case if user inputs flag multiples times
    if len(inputDirectory1) > 1:
        raise ValueError("Too many [-i, --directory-one] arguments provided")
    elif len(inputDirectory2) > 1:
        raise ValueError("Too many [-n, --directory-two] arguments provided")
    elif len(outputDirectory) > 1:
        raise ValueError("Too many [-o, --output-directory] arguments provided")
    elif len(numStates) > 1:
        raise ValueError("Too many [-s, --state-model] arguments provided")
    elif len(saliency) > 1:
        raise ValueError("Too many [-l, --saliency-level] arguments provided")
    elif len(numProcesses) > 1:
        raise ValueError("Too many [-c, --num-cores] arguments provided")
    elif len(verbose) > 1:
        raise ValueError("Too many [-x, --exit-when-complete] arguments provided")
    elif len(diagnosticBool) > 1:
        raise ValueError("Too many [-d, --diagnostic-figures] arguments provided")
    elif len(numTrials) > 1:
        raise ValueError("Too many [-t, --num-trials] arguments provided")
    elif len(samplingSize) > 1:
        raise ValueError("Too many [-z, --sampling-size] arguments provided")
    inputDirectory1 = inputDirectory1[0]
    inputDirectory2 = inputDirectory2[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    numProcesses = numProcesses[0]
    if diagnosticBool:
        diagnosticBool = True
    else:
        diagnosticBool = False
    numTrials = numTrials[0]
    samplingSize = samplingSize[0]
    if verbose:
        verbose = verbose[0]

    inputDirPath1 = Path(inputDirectory1)
    inputDirPath2 = Path(inputDirectory2)
    outputDirPath = Path(outputDirectory)

    print()
    print("File Directory 1 =", inputDirPath1)
    print("File Directory 2 =", inputDirPath2)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)
    if numProcesses == 0:
        print("Number of Cores = All available")
    else:
        print("Number of Cores =", numProcesses)

    if saliency != 1 and saliency != 2:
        raise ValueError("Saliency Metric Invalid: {} Please ensure that saliency metric is either 1 or 2 (Saliency of 3 is unsupported for pairwise comparison".format(saliency))

    # Making paths absolute
    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath
    if not PurePath(inputDirPath1).is_absolute():
        inputDirPath1 = Path.cwd() / inputDirPath1
    if not PurePath(inputDirPath2).is_absolute():
        inputDirPath2 = Path.cwd() / inputDirPath2

    # For making sure all files are consistently named
    fileTag = "{}_{}_s{}".format(inputDirPath1.name, inputDirPath2.name, saliency)

    # Check that paths are valid before doing anything
    if not inputDirPath1.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath1)))
    if not inputDirPath1.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath1)))
    if not list(inputDirPath1.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty:", str(inputDirPath1))
    if not inputDirPath2.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath2)))
    if not inputDirPath2.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath2)))
    if not list(inputDirPath2.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty:", str(inputDirPath2))

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
    storedExpPath = outputDirPath / "exp_freq_{}.npy".format(fileTag)
 
    print("\nSTEP 1: Per data file background frequency calculation", flush=True)
    for file1 in inputDirPath1.glob("*"):
        # Skip over ".genome" files
        if file1.name.split(".")[-1] == "genome":
            continue
        # Find matching file in other directory
        if not list(inputDirPath2.glob(file1.name)):
            raise FileNotFoundError("File not found: {} Please ensure corresponding files within input directories directories 1 and 2 have the same name".format(str(inputDirPath2 / file1.name)))
        else:
            file2 = next(inputDirPath2.glob(file1.name))
        if not file1.is_dir() and not file2.is_dir():
            computeEpilogosExpectedMaster.main(file1, file2, numStates, saliency, outputDirPath, fileTag, numProcesses, verbose)


    print("\nSTEP 2: Background frequency combination", flush=True)
    computeEpilogosExpectedCombination.main(outputDirPath, storedExpPath, fileTag, verbose)

    # Calculate the observed frequencies and scores
    print("\nSTEP 3: Score calculation", flush=True)
    for file1 in inputDirPath1.glob("*"):
        # Skip over ".genome" files
        if file1.name.split(".")[-1] == "genome":
            continue

        # Find matching file in other directory
        if not list(inputDirPath2.glob(file1.name)):
            raise FileNotFoundError("File not found: {} Please ensure corresponding files within input directories directories 1 and 2 have the same name".format(str(inputDirPath2 / file1.name)))
        else:
            file2 = next(inputDirPath2.glob(file1.name))

        if not file1.is_dir() and not file2.is_dir():
            computeEpilogosScoresMaster.main(file1, file2, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses, verbose)
            
    # Fitting, calculating p-values, and visualizing pairiwse differences
    print("\nSTEP 4: Generating p-values and figures", flush=True)
    computeEpilogosPairwiseVisual.main(inputDirPath1.name, inputDirPath2.name, numStates, outputDirPath, fileTag, numProcesses, diagnosticBool, numTrials, samplingSize)


if __name__ == "__main__":
    main()