import sys
from pathlib import Path
from time import sleep
from os import remove
import subprocess
from pathlib import PurePath
import errno
import click

from epilogos.expected import main as expected
from epilogos.expectedCombination import main as expectedCombination
from epilogos.scores import main as scores
from epilogos.roiSingle import main as roiSingle
from epilogos.roiAndVisualPairwise import main as roiPairwise
from epilogos.helpers import getNumStates


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
@click.option("-m", "--mode", "mode", type=click.Choice(["single", "paired"]), default="single", show_default=True,
              help="single for single group epilogos and paired for 2 group epilogos")
@click.option("-l", "--local", "commandLineBool", is_flag=True,
              help="If enabled, Epilogos will run locally in your terminal rather than on a SLURM cluster")
@click.option("-i", "--input-directory", "inputDirectory", type=str,
              help="Path to directory that contains files to read from (ALL files in this directory will be read in)")
@click.option("-a", "--directory-one", "inputDirectory1", type=str,
              help="Path to first directory that contains files to read from "
                   + "(ALL files in this directory will be read in)")
@click.option("-b", "--directory-two", "inputDirectory2", type=str,
              help="Path to second directory that contains files to read from "
                   + "(ALL files in this directory will be read in)")
@click.option("-o", "--output-directory", "outputDirectory", type=str,
              help="Output Directory (CANNOT be the same as input directory)\n")
@click.option("-j", "--state-info", "stateInfo", type=str, help="State model info file")
@click.option("-s", "--saliency", "saliency", type=int, default=1, show_default=True,
              help="Desired saliency level (1, 2, or 3)")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=1,
              help="The number of cores to run on (0=Uses all cores) [default: 1]")
@click.option("-x", "--exit", "exitBool", is_flag=True,
              help="If flag is enabled program will exit upon submission of all SLURM processes rather than "
                   + "completion of all processes")
@click.option("-d", "--diagnostic-figures", "diagnosticBool", is_flag=True,
              help="If flag is enabled, program will output some diagnostic figures of the gennorm fit on the null "
                   + "data and comparisons between the null and real data")
@click.option("-t", "--num-trials", "numTrials", type=int, default=101, show_default=True,
              help="The number of times subsamples of the scores are fit when using a null distribution")
@click.option("-z", "--sampling-size", "samplingSize", type=int, default=100000, show_default=True,
              help="The size of the subsamples on which the scores are fit when using a null distribution")
@click.option("-q", "--quiescent-state", "quiescentState", type=int, default=-1,
              help="If a bin contains only states of this value, it is treated as quiescent and not factored into "
                   + "fitting. If set to 0, filtering is not done. [default: last state]")
@click.option("-g", "--group-size", "groupSize", type=int, default=-1, show_default=True,
              help="In pairwise epilogos controls the sizes of the shuffled arrays. "
                   + "Default is sizes of the input groups")
@click.option("-v", "--version", "version", is_flag=True,
              help="If flag is enabled epilogos will print the installed version number and exit")
@click.option("-p", "--partition", "partition", type=str,
              help="Request a specific partition for the SLURM resource allocation. If not specified, uses the "
                   + "default partition as designated by the system administrator")
@click.option("-n", "--null-distribution", "pvalBool", is_flag=True,
              help="If flag is enabled, epilogos will calculate p-values for pairwise scores")
@click.option("-w", "--roi-width", "roiWidth", type=int, default=0,
              help="The number of bins in a region of interest [default: 50(single)/125(paired)]")
@click.option("-f", "--file-tag", "fileTag", type=str, default="null",
              help="Tag to be appended in output filenames [default: input-directory_saliency]")
@click.option("--exp-freq-mem", "expFreqMem", type=int, default=20000,
              help="Memory (in MB) for the expected frequency calcuation jobs [default: 20000]")
@click.option("--exp-comb-mem", "expCombMem", type=int, default=8000,
              help="Memory (in MB) for the expected frequency combination job [default: 8000]")
@click.option("--score-mem", "scoreMem", type=int, default=40000,
              help="Memory (in MB) for the expected frequency calcuation jobs [default: 40000]")
@click.option("--roi-mem", "roiMem", type=int, default=-1,
              help="Memory (in MB) for the expected frequency calcuation jobs "
                   + "[default: 20000 (single) / 100000 (paired)]")
def main(mode, commandLineBool, inputDirectory, inputDirectory1, inputDirectory2, outputDirectory, stateInfo, saliency,
         numProcesses, exitBool, diagnosticBool, numTrials, samplingSize, quiescentState, groupSize, version, partition,
         pvalBool, roiWidth, fileTag, expFreqMem, expCombMem, scoreMem, roiMem):
    """
    Information-theoretic navigation of multi-tissue functional genomic annotations

    Written by Jacob Quon and Wouter Meuleman
    """

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
    """, flush=True)

    if len(sys.argv) == 1:
        print("Run 'epilogos -h' for help")
        sys.exit()

    if version:
        from epilogos import __version__
        print("Version:", __version__)
        sys.exit()

    # Make sure all flags are submitted as expected
    checkFlags(mode, commandLineBool, inputDirectory, inputDirectory1, inputDirectory2, outputDirectory, stateInfo,
               exitBool, diagnosticBool, partition, pvalBool)

    verbose = False if commandLineBool else True
    numStates = getNumStates(stateInfo)
    # Quiescent value is user input - 1 because states are read in to be -1 from their values
    quiescentState = numStates - 1 if quiescentState == -1 else quiescentState - 1

    # setting roiWidth from default of 0
    if roiWidth == 0:
        roiWidth = 50 if mode == "single" else 125

    # Get paths from arguments and turn them into absolute paths
    if mode == "single":
        inputDirPath = Path(inputDirectory)
        inputDirPath2 = ""
    else:
        inputDirPath = Path(inputDirectory1)
        inputDirPath2 = Path(inputDirectory2)
        if not PurePath(inputDirPath2).is_absolute():
            inputDirPath2 = Path.cwd() / inputDirPath2
    if not PurePath(inputDirPath).is_absolute():
        inputDirPath = Path.cwd() / inputDirPath
    outputDirPath = Path(outputDirectory)
    if not PurePath(outputDirPath).is_absolute():
        outputDirPath = Path.cwd() / outputDirPath

    # Make sure argments are valid
    checkArguments(mode, saliency, inputDirPath, inputDirPath2, outputDirPath, numProcesses, numStates, numTrials,
                   samplingSize, quiescentState, groupSize, roiWidth)

    # Informing user of their inputs
    print()
    if mode == "single":
        print("Input Directory =", inputDirPath)
    else:
        print("Input Directory 1 =", inputDirPath)
        print("Input Directory 2 =", inputDirPath2)
    print("State Model =", numStates)
    print("Saliency level =", saliency)
    print("Output Directory =", outputDirPath)
    if numProcesses == 0:
        print("Number of Cores = All available", flush=True)
    else:
        print("Number of Cores =", numProcesses, flush=True)
    if mode == "paired" and quiescentState == -1:  # If quiescentState was user input as 0 it is now -1
        print("Quiescent State = No quiescent filtering")
    elif mode == "paired":
        print("Quiescent State =", quiescentState + 1)

    # For making sure all files are consistently named
    if fileTag == "null":
        if mode == "single":
            fileTag = "{}_s{}".format(inputDirPath.name, saliency)
        else:
            fileTag = "{}_{}_s{}".format(inputDirPath.name, inputDirPath2.name, saliency)

    # Path for storing/retrieving the expected frequency array
    storedExpPath = outputDirPath / "exp_freq_{}.npy".format(fileTag)

    if not commandLineBool:
        # Variable for the sbatch submission in case we are using slurm
        expFreqMem, expCombMem, scoreMem, roiMem = determineMemories(numProcesses, expFreqMem, expCombMem, scoreMem,
                                                                     roiMem, mode)

        if partition:
            partition = "--partition=" + partition
        else:
            partition = ""

        # Creating directories for slurm output and error logs
        (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
        (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)
        print("\nSlurm .out log files are located at: {}".format(outputDirPath / ".out/"))
        print("Slurm .err log files are located at: {}".format(outputDirPath / ".err/"), flush=True)

        # Finding the location of the .py files that must be run
        if PurePath(__file__).is_absolute():
            pythonFilesDir = Path(__file__).parents[1] / "epilogos/"
        else:
            pythonFilesDir = (Path.cwd() / Path(__file__)).parents[1] / "epilogos/"
            print("Path generated from current working directory. May cause errors")

    # Calculate the expected frequency for each file in the input directory
    expJobIDArr = []
    print("\nSTEP 1: Per data file background frequency calculation", flush=True)
    for file in inputDirPath.glob("*"):
        if mode == "single":
            if commandLineBool:
                expected(file, "null", numStates, saliency, outputDirPath, fileTag, numProcesses, verbose)
            else:
                computeExpectedPy = pythonFilesDir / "expected.py"
                pythonCommand = "python {} {} null {} {} {} {} {} {}".format(computeExpectedPy, file, numStates,
                                                                             saliency, outputDirPath, fileTag,
                                                                             numProcesses, verbose)
                expJobIDArr.append(submitSlurmJob("_" + file.name.split(".")[0], "exp_calc", fileTag, outputDirPath,
                                                  pythonCommand, saliency, partition, expFreqMem, ""))
        else:
            # Find matching file in other directory
            if not list(inputDirPath2.glob(file.name)):
                raise FileNotFoundError("File not found: {}".format(str(inputDirPath2 / file.name))
                                        + "Please ensure corresponding files within input directories "
                                        + "directories 1 and 2 have the same name")
            else:
                file2 = next(inputDirPath2.glob(file.name))

            if commandLineBool:
                expected(file, file2, numStates, saliency, outputDirPath, fileTag, numProcesses, verbose)
            else:
                computeExpectedPy = pythonFilesDir / "expected.py"
                pythonCommand = "python {} {} {} {} {} {} {} {} {}".format(computeExpectedPy, file, file2, numStates,
                                                                           saliency, outputDirPath, fileTag,
                                                                           numProcesses, verbose)
                expJobIDArr.append(submitSlurmJob("_" + file.name.split(".")[0], "exp_calc", fileTag, outputDirPath,
                                                  pythonCommand, saliency, partition, expFreqMem, ""))

    if not commandLineBool:
        # Create a string for slurm dependency to work and to print more nicely
        expJobIDStr = str(expJobIDArr).strip('[]').replace(" ", "")
        print("    JobIDs:", expJobIDStr, flush=True)

    # Combining all the different chromosome expected frequency arrays into one
    print("\nSTEP 2: Background frequency combination", flush=True)
    if commandLineBool:
        expectedCombination(outputDirPath, storedExpPath, fileTag, verbose)
    else:
        computeExpectedCombinationPy = pythonFilesDir / "expectedCombination.py"
        pythonCommand = "python {} {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, storedExpPath,
                                                       fileTag, verbose)
        combinationJobID = submitSlurmJob("", "exp_comb", fileTag, outputDirPath, pythonCommand, saliency, partition,
                                          expCombMem, "--dependency=afterok:{}".format(expJobIDStr))
        print("    JobID:", combinationJobID, flush=True)

    scoreJobIDArr = []
    # Calculate the observed frequencies and scores
    print("\nSTEP 3: Score calculation", flush=True)
    for file in inputDirPath.glob("*"):
        if mode == "single":
            if commandLineBool:
                scores(file, "null", numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses,
                       quiescentState, groupSize, verbose)
            else:
                computeScorePy = pythonFilesDir / "scores.py"
                pythonCommand = "python {} {} null {} {} {} {} {} {} {} {} {}".format(computeScorePy, file, numStates,
                                                                                      saliency, outputDirPath,
                                                                                      storedExpPath, fileTag,
                                                                                      numProcesses, quiescentState,
                                                                                      groupSize, verbose)
                scoreJobIDArr.append(submitSlurmJob("_" + file.name.split(".")[0], "score", fileTag, outputDirPath,
                                                    pythonCommand, saliency, partition, scoreMem,
                                                    "--dependency=afterok:{}".format(combinationJobID)))
        else:
            # Find matching file in other directory
            if not list(inputDirPath2.glob(file.name)):
                raise FileNotFoundError("File not found: {}".format(str(inputDirPath2 / file.name))
                                        + "Please ensure corresponding files within input directories "
                                        + "directories 1 and 2 have the same name")
            else:
                file2 = next(inputDirPath2.glob(file.name))

            if commandLineBool:
                scores(file, file2, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses,
                       quiescentState, groupSize, verbose)
            else:
                computeScorePy = pythonFilesDir / "scores.py"
                pythonCommand = "python {} {} {} {} {} {} {} {} {} {} {} {}".format(computeScorePy, file, file2,
                                                                                    numStates, saliency, outputDirPath,
                                                                                    storedExpPath, fileTag,
                                                                                    numProcesses, quiescentState,
                                                                                    groupSize, verbose)
                scoreJobIDArr.append(submitSlurmJob("_" + file.name.split(".")[0], "score", fileTag, outputDirPath,
                                                    pythonCommand, saliency, partition, scoreMem,
                                                    "--dependency=afterok:{}".format(combinationJobID)))

    if not commandLineBool:
        # Create a string for slurm dependency to work
        scoreJobIDStr = str(scoreJobIDArr).strip('[]').replace(" ", "")
        print("    JobIDs:", scoreJobIDStr, flush=True)

    if mode == "single":
        # Create a regions of interest text file
        print("\nSTEP 4: Finding regions of interest", flush=True)
        if commandLineBool:
            roiSingle(outputDirPath, stateInfo, fileTag, storedExpPath, roiWidth, verbose)
        else:
            computeROISinglePy = pythonFilesDir / "roiSingle.py"
            pythonCommand = "python {} {} {} {} {} {} {}".format(computeROISinglePy, outputDirPath, stateInfo, fileTag,
                                                                 storedExpPath, roiWidth, verbose)
            summaryJobID = submitSlurmJob("", "roi", fileTag, outputDirPath, pythonCommand, saliency, partition,
                                          roiMem, "--dependency=afterok:{}".format(scoreJobIDStr))
            print("    JobID:", summaryJobID, flush=True)
    else:
        # Fitting, calculating p-values, and visualizing pairiwse differences
        print("\nSTEP 4: Generating p-values, rois, & figures", flush=True)
        if commandLineBool:
            roiPairwise(inputDirPath.name, inputDirPath2.name, stateInfo, outputDirPath, fileTag, numProcesses,
                        pvalBool, diagnosticBool, numTrials, samplingSize, storedExpPath, roiWidth, verbose)
        else:
            computeROIPairwisePy = pythonFilesDir / "roiAndVisualPairwise.py"
            pythonCommand = "python {} {} {} {} {} {} {} {} {} {} {} {} {} {}".format(computeROIPairwisePy,
                                                                                      inputDirPath.name,
                                                                                      inputDirPath2.name, stateInfo,
                                                                                      outputDirPath, fileTag,
                                                                                      numProcesses, pvalBool,
                                                                                      diagnosticBool, numTrials,
                                                                                      samplingSize, storedExpPath,
                                                                                      roiWidth, verbose)
            summaryJobID = submitSlurmJob("", "visual", fileTag, outputDirPath, pythonCommand, saliency, partition,
                                          roiMem, "--dependency=afterok:{}".format(scoreJobIDStr))
            print("    JobID:", summaryJobID, flush=True)

    if not commandLineBool:
        allJobIDs = "{},{},{},{}".format(expJobIDStr, combinationJobID, scoreJobIDStr, summaryJobID)
        print("\nAll JobIDs:\n    ", allJobIDs, flush=True)

    # If the user wants to exit upon job completion rather than submission
    # If a job fails, it cancels all other jobs
    if not commandLineBool and not exitBool:
        checkExit(mode, allJobIDs, expJobIDArr, scoreJobIDArr, outputDirPath, saliency)


def checkFlags(mode, commandLineBool, inputDirectory, inputDirectory1, inputDirectory2, outputDirectory, stateInfo,
               exitBool, diagnosticBool, partition, pvalBool):
    """
    Checks all the input flags are makes sure that there are not duplicates, required flags are present, and
    incompatible flags are not present together

    Input:
    See click options at top of script
    """
    # Checking that all required flags are inputted
    if mode == "single" and not inputDirectory:
        print("ERROR: [-i, --input-directory] required in 'single' group mode")
        sys.exit()
    elif mode == "paired" and not inputDirectory1:
        print("ERROR: [-a, --directory-one] required in 'paired' group mode")
        sys.exit()
    elif mode == "paired" and not inputDirectory2:
        print("ERROR: [-b, --directory-two] required in 'paired' group mode")
        sys.exit()
    elif not outputDirectory:
        print("ERROR: [-o, --output-directory] required")
        sys.exit()
    elif not stateInfo:
        print("ERROR: [-n, --state-info] required")
        sys.exit()

    # Checking that incompatible flags are not input together
    if mode == "single" and inputDirectory1:
        print("ERROR: [-m, --mode] 'single' not compatible with [-a, --directory-one] option")
        sys.exit()
    elif mode == "single" and inputDirectory2:
        print("ERROR: [-m, --mode] 'single' not compatible with [-b, --directory-two] option")
        sys.exit()
    elif mode == "single" and diagnosticBool:
        print("ERROR: [-m, --mode] 'single' not compatible with [-d, --diagnostic-figures] flag")
        sys.exit()
    elif mode == "single" and pvalBool:
        print("ERROR: [-m, --mode] 'single' not compatible with [-n, --null-distribution] flag")
        sys.exit()
    elif mode == "paired" and inputDirectory:
        print("ERROR: [-m, --mode] 'paired' not compatible with [-i, --input-directory] option")
        sys.exit()
    elif commandLineBool and exitBool:
        print("ERROR: [-l, --cli] flag not compatible with [-x, --exit] flag")
        sys.exit()
    elif commandLineBool and partition:
        print("Error: [-l, --cli] flag not compatible with [-p, --partition] option")
        sys.exit()


def checkArguments(mode, saliency, inputDirPath, inputDirPath2, outputDirPath, numProcesses, numStates, numTrials,
                   samplingSize, quiescentState, groupSize, roiWidth):
    """
    Checks whether user submitted arguments have valid values

    Input:
    mode           -- 'single' or 'paired'; tells us which version of epilogos we are running
    saliency       -- The saliency metric input by the user
    inputDirPath   -- The path to the only input directory in single epilogos and the first input directory in paired
                      epilogos
    inputDirPath2  -- The path to the second input directory in the paired epilogos case
    outputDirPath  -- The path to the output directory for epilogos
    numProcesses   -- The number of cores the user would like to use
    numStates      -- The number of states in the state model
    numTrials      -- The number of gennorm fits to do in pairwise mode
    samplingSize   -- The amount of null data to fit in pairwise mode
    quiescentState -- The state used to filter out quiescent bins
    groupSize      -- The size of the null (shuffled) score arrays, -1 means inputed sizes
    roiWidth       -- size of regions of interest in bins
    """
    # Check validity of saliency
    if mode == "single" and saliency != 1 and saliency != 2 and saliency != 3:
        raise ValueError("Saliency Metric Invalid: {}".format(saliency)
                         + "Please ensure that saliency metric is either 1, 2, or 3")
    elif mode == "paired" and saliency != 1 and saliency != 2:
        raise ValueError("Saliency Metric Invalid: {}".format(saliency)
                         + "Please ensure that saliency metric is either 1 or 2 "
                         + "(Saliency of 3 is unsupported for pairwise comparison")

    # Check validity of directories
    if not inputDirPath.exists():
        raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath)))
    if not inputDirPath.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath)))
    if not list(inputDirPath.glob("*")):
        raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty", str(inputDirPath))
    if mode == "paired":
        if not inputDirPath2.exists():
            raise FileNotFoundError("Given path does not exist: {}".format(str(inputDirPath2)))
        if not inputDirPath2.is_dir():
            raise NotADirectoryError("Given path is not a directory: {}".format(str(inputDirPath2)))
        if not list(inputDirPath2.glob("*")):
            raise OSError(errno.ENOTEMPTY, "Ensure given directory is not empty", str(inputDirPath2))
    if not outputDirPath.exists():
        outputDirPath.mkdir(parents=True)
    if not outputDirPath.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDirPath)))

    if numProcesses < 0:
        print("ERROR: Number of cores must be positive or zero (0 means use all cores)")
        sys.exit()

    if numTrials <= 0:
        print("ERROR: Number of trials must be greater than zero")
        sys.exit()

    if samplingSize <= 0:
        print("ERROR: Sampling size must be greater than zero")
        sys.exit()

    if quiescentState < -1:
        print("ERROR: Quiescent state value must be positive or zero (0 means do not filter)")
        sys.exit()
    elif quiescentState >= numStates:
        print("ERROR: Quiescent state value must be a state provided in the state model")
        sys.exit()

    if groupSize < -1:
        print("ERROR: Group size value must be positive or -1 (-1 means use inputted group sizes)")
        sys.exit()

    if roiWidth < 0:
        print("ERROR: Group size value must be greater than 0")
        sys.exit()


def submitSlurmJob(filename, jobPrefix, fileTag, outputDirPath, pythonCommand, saliency, partition, memory, dependency):
    """
    Submits a epilogos job to a SLURM cluster

    Input:
    filename      -- The name of the file epilogos is using
    jobPrefix     -- The step of epilogos we are on ('exp_calc', 'exp_comb', 'score', 'visual', 'roi')
    fileTag       -- String that helps us ensure jobs are named consistenly across one epilogos run
    outputDirPath -- The output directory for epilogos will be used for error and output logs
    pythonCommand -- The python command to run on the SLURM job
    saliency      -- The saliency metric being used in this epilogos run
    paritition    -- The slurm paritition to use if user elected to choose a specific one
    memory        -- The amount of memory to be allocated to the slurm job
    dependency    -- Which jobs the SLURM job must wait for before starting

    Output:
    The job number
    """
    jobName = "{}_{}{}".format(jobPrefix, fileTag, filename)
    jobOutPath = outputDirPath / (".out/" + jobName + ".out")
    jobErrPath = outputDirPath / (".err/" + jobName + ".err")

    # Creating the out and err files for the batch job
    if jobOutPath.exists():
        remove(jobOutPath)
    if jobErrPath.exists():
        remove(jobErrPath)
    try:
        jobOutPath.touch()
        jobErrPath.touch()
    except FileExistsError as err:
        # This error should never occur because we are deleting the files first
        print(err)
        return

    # Create a string for the slurm command
    if saliency == 1:
        slurmCommand = "sbatch {} --job-name={}.job --output={} --error={} {} {} --wrap='{}'"\
            .format(dependency, jobName, jobOutPath, jobErrPath, partition, memory, pythonCommand)
    elif saliency == 2:
        slurmCommand = "sbatch {} --job-name={}.job --output={} --error={} {} {} --wrap='{}'"\
            .format(dependency, jobName, jobOutPath, jobErrPath, partition, memory, pythonCommand)
    elif saliency == 3:
        slurmCommand = "sbatch {} --job-name={}.job --output={} --error={} {} {} --wrap='{}'"\
            .format(dependency, jobName, jobOutPath, jobErrPath, partition, memory, pythonCommand)

    sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    if not sp.stdout.startswith("Submitted batch"):
        raise ChildProcessError("SlurmError: sbatch not submitted correctly")

    return int(sp.stdout.split()[-1])


def checkExit(mode, allJobIDs, expJobIDArr, scoreJobIDArr, outputDirPath, saliency):
    """
    Prevents exiting of the main script until there is an error in one of the SLURM jobs or all SLURM jobs have finished

    Input:
    mode          -- 'single' or 'paired' depending on which version of epilogos we are running
    allJobIDs     -- String of all the SLURM job IDs
    expJobIDArr   -- List of the SLURM job IDs for the expected frequency step
    scoreJobIDArr -- List of the SLURM job IDs for the score calculation step
    outputDirPath -- Epilogos' output directory
    saliency      -- The saliency metric epilogos is using
    """

    jobCheckStr = "sacct --format=JobID%18,JobName%50,State%10 --jobs {}".format(allJobIDs)

    completedJobs = []
    calculationStep = ""

    # Every x seconds check what jobs are done and print accordingly
    while True:
        # Check the jobs which are done
        sp = subprocess.run(jobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
        spLines = sp.stdout.split("\n")

        # Printing separate step headers
        if len(completedJobs) == 0 and calculationStep == "":
            print("\n Step 1: Per data file background frequency calculation\n{}\n{}\n{}"
                  .format("-" * 80, spLines[0], spLines[1]), flush=True)
            calculationStep = "exp_calc"
        elif len(completedJobs) >= len(expJobIDArr) and calculationStep == "exp_calc":
            print("\n Step 2: Background frequency combination\n{}\n{}\n{}"
                  .format("-" * 80, spLines[0], spLines[1]), flush=True)
            calculationStep = "exp_comb"
        elif len(completedJobs) >= (len(expJobIDArr) + 1) and calculationStep == "exp_comb":
            print("\n Step 3: Score calculation\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]), flush=True)
            calculationStep = "score"
        elif mode == "single" and len(completedJobs) >= (len(expJobIDArr) + 1 + len(scoreJobIDArr)) \
                and calculationStep == "score":
            print("\n Step 4: Finding regions of interest\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]),
                  flush=True)
            calculationStep = "roi"
        elif mode == "paired" and len(completedJobs) >= (len(expJobIDArr) + 1 + len(scoreJobIDArr)) \
                and calculationStep == "score":
            print("\n Step 4: Generating p-values and figures\n{}\n{}\n{}"
                  .format("-" * 80, spLines[0], spLines[1]), flush=True)
            calculationStep = "visual"

        # Print out jobs when they are completed
        for line in spLines[2:]:
            if "COMPLETED" in line and "allocation" not in line:
                jobID = line.split()[0]
                # Don't want to print if we have already printed
                # Also ensure we're only printing jobs from correct step
                if jobID not in completedJobs and ".batch" not in jobID and calculationStep in line:
                    completedJobs.append(jobID)
                    print(line, flush=True)

        # Check if there was an error, if so cancel everything and exit the program
        if "FAILED" in sp.stdout or "CANCELLED" in sp.stdout:
            print("\nERROR RUNNING JOBS: CANCELLING ALL REMAINING JOBS\n")
            print("Please check error logs in: {}/.err/\n".format(outputDirPath))
            subprocess.run("scancel {}".format(allJobIDs), shell=True)
            break

        # If final job is done, exit the program
        # Checks are if the 3rd line is not empty, if there are no more running or pending values,
        # and if an "allocation" job is not in the output
        if spLines[2] and not ("RUNNING" in sp.stdout or "PENDING" in sp.stdout) and "allocation" not in sp.stdout:
            print("\nAll jobs finished successfully. Please find output in: {}".format(outputDirPath))
            print("\nPlease find output and error logs in {} and /.err respectively\n".format(outputDirPath / ".out"))
            break

        # Don't want to spam cluster with commands so sleep
        # Sleep less time for saliency one since it is much faster
        if saliency == 1:
            sleep(2)
        else:
            sleep(10)


def determineMemories(numProcesses, expFreqMem, expCombMem, scoreMem, roiMem, mode):
    """
    Generates memory strings for slurm commands

    Input:
    numProcesses -- The number of cores the user would like to use
    expFreqMem   -- The memory (in MB) to be used for the expected frequency calculations
    expCombMem   -- The memory (in MB) to be used for the expected frequency calculations
    scoreMem     -- The memory (in MB) to be used for the expected frequency calculations
    roiMem       -- The memory (in MB) to be used for the expected frequency calculations
    mode         -- 'single' or 'paired'; tells us which version of epilogos we are running

    Output:
    expFreqMem -- A string for slurm commands containing information on ntasks and memory (for exp freq calculation)
    expCombMem -- A string for slurm commands containing information on ntasks and memory (for exp comb calculation)
    scoreMem   -- A string for slurm commands containing information on ntasks and memory (for score calculation)
    roiMem     -- A string for slurm commands containing information on ntasks and memory (for roi calculation)
    """

    # 0 cores is defined to take all resources on node
    if numProcesses == 0:
        expFreqMem, expCombMem, scoreMem, roiMem = "--exclusive --mem=0"
    else:
        expFreqMem = "--ntasks={} --mem={}".format(numProcesses, expFreqMem)
        expCombMem = "--ntasks={} --mem={}".format(numProcesses, expCombMem)
        scoreMem   = "--ntasks={} --mem={}".format(numProcesses, scoreMem)

        # Default roi memory is different for single/paired mode
        if roiMem == -1:
            if mode == "single":
                roiMem = "--ntasks={} --mem=20000".format(numProcesses)
            else:
                roiMem = "--ntasks={} --mem=100000".format(numProcesses)
        else:
            roiMem = "--ntasks={} --mem={}".format(numProcesses, roiMem)

    return expFreqMem, expCombMem, scoreMem, roiMem


if __name__ == "__main__":
    main()
