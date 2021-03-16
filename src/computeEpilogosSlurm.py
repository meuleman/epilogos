from pathlib import Path
from time import sleep
import click
from os import remove
import subprocess
from pathlib import PurePath
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
@click.option("-n", "--number-of-states", "numStates", type=int, required=True, multiple=True, help="Number of states in chromatin state model")
@click.option("-s", "--saliency", "saliency", type=int, default=[1], show_default=True, multiple=True, help="Desired saliency level (1, 2, or 3)")
@click.option("-m", "--mode-of-operation", "modeOfOperation", type=click.Choice(["bg", "s", "both"]), default=["both"], show_default=True, multiple=True, help="bg for background, s for scores, both for both")
@click.option("-b", "--background-directory", "expFreqDir", type=str, default=["null"], multiple=True, help="Path to where the background frequency array is read from (-m s) or written to (-m bg, -m both) [default: output-directory]")
@click.option("-c", "--num-cores", "numProcesses", type=int, default=[0], multiple=True, help="The number of cores to run on [default: 0 = Uses all cores]")
@click.option("-x", "--exit", "exitBool", is_flag=True, multiple=True, help="If flag is enabled program will exit upon submission of all SLURM processes rather than completion of all processes")
def main(inputDirectory, outputDirectory, numStates, saliency, modeOfOperation, expFreqDir, numProcesses, exitBool):
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
    elif len(exitBool) > 1:
        raise ValueError("Too many [-x, --exit-when-complete] arguments provided")
    inputDirectory = inputDirectory[0]
    outputDirectory = outputDirectory[0]
    numStates = numStates[0]
    saliency = saliency[0]
    modeOfOperation = modeOfOperation[0]
    expFreqDir = expFreqDir[0]
    numProcesses = numProcesses[0]
    verbose=True

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
        numTasks = "--exclusive"
    else:
        numTasks = "--ntasks={}".format(numProcesses)

    # For slurm output and error later
    (outputDirPath / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / ".err/").mkdir(parents=True, exist_ok=True)
    print("\nSlurm .out log files are located at: {}".format(outputDirPath / ".out/"))
    print("Slurm .err log files are located at: {}".format(outputDirPath / ".err/"))


    # Path for storing/retrieving the expected frequency array
    # Expected frequency arrays are stored according to path of the input file directory
    storedExpPath = Path(expFreqDir) / "exp_freq_{}.npy".format(fileTag)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[1] / "src/"
    else:
        pythonFilesDir = (Path.cwd() / Path(__file__)).parents[1] / "src/"
        print("Path generate from current working directory. May cause errors")

    # Only calculate the expected frequencies if user asks for it, otherwise just load from where the user said
    expJobIDArr = []   
    if modeOfOperation != "s":   
        print("\nSTEP 1: Per data file background frequency calculation")
        for file in inputDirPath.glob("*"):
            # Skip over ".genome" files
            if file.name.split(".")[-1] == "genome":
                continue
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "exp_calc_{}_{}".format(fileTag, filename)
                jobOutPath = outputDirPath / (".out/" + jobName + ".out")
                jobErrPath = outputDirPath / (".err/" + jobName + ".err")

                # Creating the out and err files for the batch job
                if jobOutPath.exists():
                    remove(jobOutPath)
                if jobErrPath.exists():
                    remove(jobErrPath)
                try:
                    jout = open(jobOutPath, 'x')
                    jout.close()
                    jerr = open(jobErrPath, 'x')
                    jerr.close()
                except FileExistsError as err:
                    # This error should never occur because we are deleting the files first
                    print(err)
                    return

                computeExpectedPy = pythonFilesDir / "computeEpilogosExpectedMaster.py"

                pythonCommand = "python {} {} null {} {} {} {} {} {}".format(computeExpectedPy, file, numStates, saliency, outputDirPath, fileTag, numProcesses, verbose)

                if saliency == 1:
                    slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                elif saliency == 2:
                    slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                elif saliency == 3:
                    slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    raise ChildProcessError("SlurmError: sbatch not submitted correctly")
                
                expJobIDArr.append(int(sp.stdout.split()[-1]))

        # Combining all the different chromosome expected frequency arrays into one
        # create a string for slurm dependency to work
        expJobIDStr = str(expJobIDArr).strip('[]').replace(" ", "")

        print("    JobIDs:", expJobIDStr)

        print("\nSTEP 2: Background frequency combination")

        jobName = "exp_comb_{}".format(fileTag)
        jobOutPath = outputDirPath / (".out/" + jobName + ".out")
        jobErrPath = outputDirPath / (".err/" + jobName + ".err")

        # Creating the out and err files for the batch job
        if jobOutPath.exists():
            remove(jobOutPath)
        if jobErrPath.exists():
            remove(jobErrPath)
        try:
            jout = open(jobOutPath, 'x')
            jout.close()
            jerr = open(jobErrPath, 'x')
            jerr.close()
        except FileExistsError as err:
            # This error should never occur because we are deleting the files first
            print(err)
            return

        computeExpectedCombinationPy = pythonFilesDir / "computeEpilogosExpectedCombination.py"

        pythonCommand = "python {} {} {} {} {}".format(computeExpectedCombinationPy, outputDirPath, storedExpPath, fileTag, verbose)

        if saliency == 1:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 2:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 3:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(expJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            raise ChildProcessError("SlurmError: sbatch not submitted correctly")
        
        combinationJobID = int(sp.stdout.split()[-1])

        print("    JobID:", combinationJobID)

    scoreJobIDArr = []
    if modeOfOperation == "s" or modeOfOperation == "both":
        # Calculate the observed frequencies and scores
        print("\nSTEP 3: Score calculation")
        for file in inputDirPath.glob("*"):
            # Skip over ".genome" files
            if file.name.split(".")[-1] == "genome":
                continue
            if not file.is_dir():
                filename = file.name.split(".")[0]
                jobName = "score_{}_{}".format(fileTag, filename)
                jobOutPath = outputDirPath / (".out/" + jobName + ".out")
                jobErrPath = outputDirPath / (".err/" + jobName + ".err")

                # Creating the out and err files for the batch job
                if jobOutPath.exists():
                    remove(jobOutPath)
                if jobErrPath.exists():
                    remove(jobErrPath)
                try:
                    jout = open(jobOutPath, 'x')
                    jout.close()
                    jerr = open(jobErrPath, 'x')
                    jerr.close()
                except FileExistsError as err:
                    # This error should never occur because we are deleting the files first
                    print(err)
                    return
                
                computeScorePy = pythonFilesDir / "computeEpilogosScoresMaster.py"

                pythonCommand = "python {} {} null {} {} {} {} {} {} {}".format(computeScorePy, file, numStates, saliency, outputDirPath, storedExpPath, fileTag, numProcesses, verbose)

                if modeOfOperation == "s":
                    if saliency == 1:
                        slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 2:
                        slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 3:
                        slurmCommand = "sbatch --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                else:
                    if saliency == 1:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 2:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)
                    elif saliency == 3:
                        slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} {} --mem=0 --wrap='{}'".format(combinationJobID, jobName, jobOutPath, jobErrPath, numTasks, pythonCommand)

                sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

                if not sp.stdout.startswith("Submitted batch"):
                    raise ChildProcessError("SlurmError: sbatch not submitted correctly")
                
                scoreJobIDArr.append(int(sp.stdout.split()[-1]))

        # Create a string for slurm dependency to work
        scoreJobIDStr = str(scoreJobIDArr).strip('[]').replace(" ", "")
        
        print("    JobIDs:", scoreJobIDStr)

        # # Write scores out to gzipped text files
        # print("\nSTEP 4: Writing score files")
        # writeJobIDArr = []
        # for file in inputDirPath.glob("*"):
        #     # Skip over ".genome" files
        #     if file.name.split(".")[-1] == "genome":
        #         continue
        #     if not file.is_dir():
        #         filename = file.name.split(".")[0]
        #         jobName = "write_{}_{}".format(fileTag, filename)
        #         jobOutPath = outputDirPath / (".out/" + jobName + ".out")
        #         jobErrPath = outputDirPath / (".err/" + jobName + ".err")

        #         # Creating the out and err files for the batch job
        #         if jobOutPath.exists():
        #             remove(jobOutPath)
        #         if jobErrPath.exists():
        #             remove(jobErrPath)
        #         try:
        #             jout = open(jobOutPath, 'x')
        #             jout.close()
        #             jerr = open(jobErrPath, 'x')
        #             jerr.close()
        #         except FileExistsError as err:
        #             # This error should never occur because we are deleting the files first
        #             print(err)
        #             return

        #         computeExpectedWritePy = pythonFilesDir / "computeEpilogosWrite.py"

        #         pythonCommand = "python {} {} {} {} {} {}".format(computeExpectedWritePy, file, numStates, outputDirPath, fileTag, verbose)

        #         slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=64000 --wrap='{}'".format(scoreJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

        #         sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        #         if not sp.stdout.startswith("Submitted batch"):
        #             raise ChildProcessError("SlurmError: sbatch not submitted correctly")

        #         writeJobIDArr.append(int(sp.stdout.split()[-1]))

        # writeJobIDStr = str(writeJobIDArr).strip('[]').replace(" ", "")
        # print("    JobIDs:", writeJobIDStr)
        # if modeOfOperation == "both":
        #     allJobIDs = "{},{},{},{}".format(expJobIDStr, combinationJobID, scoreJobIDStr, writeJobIDStr)
        #     print("\nAll JobIDs:\n    ", allJobIDs)
        # elif modeOfOperation == "s":
        #     allJobIDs = "{},{}".format(scoreJobIDStr, writeJobIDStr)
        #     print("\nAll JobIDs:\n    ", allJobIDs)
        # elif modeOfOperation == "bg":
        #     allJobIDs = "{},{}".format(expJobIDStr, combinationJobID)
        #     print("\nAll JobIDs:\n    ", allJobIDs)


        # Create a greatest hits text file
        print("\nSTEP 4: Finding greatest hits")

        jobName = "hits_{}".format(fileTag)
        jobOutPath = outputDirPath / (".out/" + jobName + ".out")
        jobErrPath = outputDirPath / (".err/" + jobName + ".err")

        # Creating the out and err files for the batch job
        if jobOutPath.exists():
            remove(jobOutPath)
        if jobErrPath.exists():
            remove(jobErrPath)
        try:
            jout = open(jobOutPath, 'x')
            jout.close()
            jerr = open(jobErrPath, 'x')
            jerr.close()
        except FileExistsError as err:
            # This error should never occur because we are deleting the files first
            print(err)
            return

        computeGreatestHitsPy = pythonFilesDir / "computeEpilogosGreatestHits.py"
        pythonCommand = "python {} {} {} {} {} {}".format(computeGreatestHitsPy, outputDirPath, numStates, fileTag, numProcesses, verbose)

        if saliency == 1:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(scoreJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 2:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(scoreJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)
        elif saliency == 3:
            slurmCommand = "sbatch --dependency=afterok:{} --job-name={}.job --output={} --partition=queue1 --error={} --ntasks=1 --mem-per-cpu=8000 --wrap='{}'".format(scoreJobIDStr, jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            raise ChildProcessError("SlurmError: sbatch not submitted correctly")
        
        greatestHitsJobID = int(sp.stdout.split()[-1])

        print("    JobID:", greatestHitsJobID)

        if modeOfOperation == "both":
            allJobIDs = "{},{},{},{}".format(expJobIDStr, combinationJobID, scoreJobIDStr, greatestHitsJobID)
            print("\nAll JobIDs:\n    ", allJobIDs)
        elif modeOfOperation == "s":
            allJobIDs = "{},{}".format(scoreJobIDStr, greatestHitsJobID)
            print("\nAll JobIDs:\n    ", allJobIDs)
        elif modeOfOperation == "bg":
            allJobIDs = "{},{}".format(expJobIDStr, combinationJobID)
            print("\nAll JobIDs:\n    ", allJobIDs)

    # If the user wants to exit upon job completion rather than submission
    # If a job fails, it cancels all other jobs
    if not exitBool:
        jobCheckStr = "sacct --format=JobID%18,JobName%50,State%10 --jobs {}".format(allJobIDs)

        completedJobs = []
        calculationStep = 0

        # Every ten seconds check what jobs are done and print accordingly
        while True:
            # Check the jobs which are done
            sp = subprocess.run(jobCheckStr, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)
            spLines = sp.stdout.split("\n")

            # Printing separate step headers
            if len(completedJobs) == 0 and calculationStep == 0 and (modeOfOperation == "bg" or modeOfOperation == "both"):
                print("\n Step 1: Per data file background frequency calculation\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]))
                calculationStep += 1
            elif len(completedJobs) >= len(expJobIDArr) and calculationStep == 1 and (modeOfOperation == "bg" or modeOfOperation == "both"):
                print("\n Step 2: Background frequency combination\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]))
                calculationStep += 1
            elif modeOfOperation == "s":
                if len(completedJobs) == 0 and calculationStep == 0:
                    print("\n Step 3: Score calculation\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]))
                    calculationStep += 1
                elif len(completedJobs) >= len(scoreJobIDArr) and calculationStep == 1:
                    print("\n Step 4: Finding greatest hits\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]))
                    calculationStep += 1
            elif modeOfOperation == "both":
                if len(completedJobs) >= (len(expJobIDArr) + 1) and calculationStep == 2:
                    print("\n Step 3: Score calculation\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]))
                    calculationStep += 1
                elif len(completedJobs) >= (len(expJobIDArr) + 1 + len(scoreJobIDArr)) and calculationStep == 3:
                    print("\n Step 4: Finding greatest hits\n{}\n{}\n{}".format("-" * 80, spLines[0], spLines[1]))
                    calculationStep += 1

            # Print out jobs when they are completed
            for line in spLines[2:]:
                if "COMPLETED" in line and "allocation" not in line:
                    jobID = line.split()[0]
                    # Don't want to print if we have already printed
                    if jobID not in completedJobs and ".batch" not in jobID:
                        completedJobs.append(jobID)
                        print(line, flush=True)

            # Check if there was an error, if so cancel everything and exit the program
            if "FAILED" in sp.stdout or "CANCELLED" in sp.stdout:
                print("\nERROR RUNNING JOBS: CANCELLING ALL REMAINING JOBS\n")
                print("Please check error logs in: {}/.err/\n".format(outputDirPath))
                subprocess.run("scancel {}".format(allJobIDs), shell=True)
                break

            # If final job is done, exit the program
            # Checks are if the 3rd line is not empty, if there are no more running or pending values and if an "allocation" job is not in the output
            if spLines[2] and not ("RUNNING" in sp.stdout or "PENDING" in sp.stdout) and "allocation" not in sp.stdout:
                print("\nAll jobs finished successfully. Please find output in: {}".format(outputDirPath))
                print("\nPlease find output and error logs in {} and {} respectively\n".format(outputDirPath / ".out", "/.err"))
                break
            
            if saliency == 1:
                sleep(2)
            else:
                sleep(10)
                
if __name__ == "__main__":
    main()
