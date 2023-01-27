"""
Similarity Search For Epilogos

Written by Jacob Quon, Nalu Tripician, Alex Reynolds, and Wouter Meuleman
"""
import numpy as np
import pandas as pd
from time import time
from pathlib import Path, PurePath
import os
import click
from epilogos.helpers import generateRegionArr
import subprocess


# Custom help format for click to separate the build and query commands
class CustomClickHelpFormat(click.Command):
    def format_help(self, ctx, formatter):
        click.echo("Usage: simsearch [OPTIONS]\n\
\n\
Options:\n\
To Build Similarity Search Data:\n\
  -b, --build                     If true builds the similarity search files\n\
                                  needed to query regions\n\
  -s, --scores TEXT               Path to scores file to be used in\n\
                                  similarity search\n\
  -o, --output-directory TEXT     Path to desired similarity search output\n\
                                  directory\n\
  -w, --window-bp INTEGER         Window size (in BP) on which to perform\n\
                                  similarity search [default: 25000]\n\
  -j, --num-jobs INTEGER          Number of slurm jobs to be used in\n\
                                  parallelization of sim search calculation\n\
                                  [default: 10]\n\
  -c, --num-cores INTEGER         Number of cores to be used in simsearch\n\
                                  calculation. If set to 0, uses all cores.\n\
                                  [default: 0=All cores]\n\
  -n, --num-matches INTEGER       Number of matches to be found by simsearch\n\
                                  for each query region [default: 100]\n\
  -f, --filter-state INTEGER      If the max signal within a region is from the\n\
                                  filter state, the region is removed from the\n\
                                  region list. The purpose of this is to be\n\
                                  used to filter out 'quiescent regions'. If\n\
                                  set to 0, filtering is not done.\n\
                                  [default: last state]\n\
  --filter-score FLOAT            If the max signal within a region is less\n\
                                  than the filter score, the region is removed\n\
                                  from the region list. The purpose of this is\n\
                                  to reduce computational complexity by\n\
                                  eliminating potentially uninteresting\n\
                                  regions. By default no filtering is done on\n\
                                  score. [default: -1==no filtering]\n\
  -p, --partition TEXT            Request a specific partition for the SLURM\n\
                                  resource allocation. If not specified, uses\n\
                                  the default partition as designated by the\n\
                                  system administrator\n\
  --mm-mem TEXT                   Memory (in MB) for the simsearch max mean\n\
                                  region selection job [default: 10000]\n\
  --calc-mem INTEGER              Memory (in MB) for the simsearch calcuation\n\
                                  jobs [default: 50000]\n\
  --write-mem INTEGER             Memory (in MB) for the simsearch write job\n\
                                  [default: 5000]\
\n\
\n\
To Query Similarity Search Data:\n\
  -q, --query TEXT                Query region formatted as chr:start-end or\n\
                                  path to tab-separated bed file containing\n\
                                  query regions\n\
  -m, --matches-file TEXT         Path to previously built\n\
                                  simsearch.bed.gz file to be queried\n\
                                  for matches\n\
  -o, --output-directory TEXT     Path to desired similarity search output\n\
                                  directory")


@click.command(context_settings=dict(help_option_names=['-h', '--help']), cls=CustomClickHelpFormat)
@click.option("-b", "--build", "buildBool", is_flag=True,
              help="If true builds the similarity search files needed to query regions")
@click.option("-s", "--scores", "scoresPath", type=str, help="Path to scores file to be used in similarity search")
@click.option("-o", "--output-directory", "outputDir", required=True, type=str,
              help="Path to desired similarity search output directory")
@click.option("-w", "--window-bp", "windowBP", type=int, default=-1, show_default=True,
              help="Window size (in BP) on which to perform similarity search [default: 25000]")
@click.option("-j", "--num-jobs", "nJobs", type=int, default=10, show_default=True,
              help="Number of slurm jobs to be used in parallelization of simsearch calculation. [default: 10]")
@click.option("-c", "--num-cores", "nCores", type=int, default=1, show_default=True,
              help="Number of cores to be used in simsearch calculation. If set to 0, uses all cores. "
                    + "[default: 1]")
@click.option("-n", "--num-matches", "nDesiredMatches", type=int, default=100, show_default=True,
              help="Number of matches to be found by simsearch for each query region [default: 100]")
@click.option("-f", "--filter-state", "filterState", type=int, default=-1,
              help="If the max signal within a region is from the filter state, the region is removed from the region "
                   + "list. The purpose of this is to be used to filter out 'quiescent regions'."
                   + "If set to 0, filtering is not done. [default: last state]")
@click.option("--filter-score", "filterScore", type=float, default=-1,
              help="If the max signal within a region is less than the filter score, the region is removed from the "
                   + "region list. The purpose of this is to reduce computational complexity by eliminating "
                   + "potentially uninteresting regions. By default no filtering is done on score. "
                   + "[default: -1 == no filtering]")
@click.option("-p", "--partition", "partition", type=str,
              help="Request a specific partition for the SLURM resource allocation. If not specified, uses the default "
                   + "partition as designated by the system administrator")
@click.option("-t", "--tag", "jobTag", type=str, default="", help="Tag to add to slurm job names [default: None]")
@click.option("--mm-mem", "mmMem", type=str, default=10000,
              help="Memory (in MB) for the simsearch max mean region selection job [default: 10000]")
@click.option("--calc-mem", "calcMem", type=int, default=50000,
              help="Memory (in MB) for the simsearch calcuation jobs [default: 50000]")
@click.option("--write-mem", "writeMem", type=int, default=5000,
              help="Memory (in MB) for the simsearch write job [default: 5000]")
@click.option("-q", "--query", "query", type=str, default="",
              help="Query region formatted as chr:start-end or path to tab-separated bed file containing query regions")
@click.option("-m", "--matches-file", "simSearchPath", type=str,
              help="Path to previously built simsearch.bed.gz file to be queried for matches")
def main(buildBool, scoresPath, outputDir, windowBP, nJobs, nCores, nDesiredMatches, filterState, filterScore,
         partition, jobTag, mmMem, calcMem, writeMem, query, simSearchPath):
    print("""\n
                 d8b                                                           888
                 Y8P                                                           888
                                                                               888
        .d8888b  888 88888b.d88b.   .d8888b   .d88b.   8888b.  888d888 .d8888b 88888b.
        88K      888 888 "888 "88b  88K      d8P  Y8b     "88b 888P"  d88P"    888 "88b
        "Y8888b. 888 888  888  888  "Y8888b. 88888888 .d888888 888    888      888  888
             X88 888 888  888  888       X88 Y8b.     888  888 888    Y88b.    888  888
         88888P' 888 888  888  888   88888P'  "Y8888  "Y888888 888     "Y8888P 888  888
    """, flush=True)

    if not buildBool and query == "":
        raise ValueError("Either -b or -q flag must be used to run simsearch")
    elif buildBool and query != "":
        raise ValueError("Both -b and -q flags cannot be used at the same time")

    outputDir = Path(outputDir)
    if not outputDir.exists():
        outputDir.mkdir(parents=True)
    if not outputDir.is_dir():
        raise NotADirectoryError("Given path is not a directory: {}".format(str(outputDir)))

    if buildBool:
        buildSimSearch(scoresPath, outputDir, windowBP, nJobs, nCores, nDesiredMatches, filterState, filterScore,
                       partition, jobTag, mmMem, calcMem, writeMem)
    else:
        querySimSearch(query, simSearchPath, outputDir)


def buildSimSearch(scoresPath, outputDir, windowBP, nJobs, nCores, nDesiredMatches, filterState, filterScore,
                   partition, jobTag, mmMem, calcMem, writeMem):
    """
    Deploys the slurm jobs to build the similarity search files for a set of scores

    Input:
    scoresPath      -- Path to scores file to be used in similarity search
    outputDir       -- pathlib Path to the output directory for similarity search files
    windowBP        -- Size of the similarity search window in bp
    nJobs           -- Number of slurm jobs to be use in parallelization of simsearch euclidean distance calculations
    nCores          -- Number of cores to be used for multiprocessing calculation of simsearch euclidean distance
    nDesiredMatches -- Number of matches to be found by simsearch algorithm
    filterState     -- If the max signal of a region is in this state, the region is removed
    filterScore     -- If the max signal of a region is less than filterScore, the region is removed
    partition       -- Slurm partition to use
    jobTag          -- Optional naming addition for slurm jobs
    mmMem           -- Memory in MB to be allocated to the maxmean region selection job
    calcMem         -- Memory in MB to be allocated to the simsearch euclidean distance calculations
    writeMem        -- Memory in MB to be allocated to writing the results of the simsearch calculation


    Output:
    simsearch_cube.npz     -- NPZ file containing all the reduced regions and their coords
                              (scores=scores, coords=coords)
    reduced_genome.npy     -- NPY file containing the reduced scores for the genome (used in simsearch eucl dist calc)
    simsearch_indices.npy  -- NPY file containing the simsearch results as indices of the reduced_genome.npy array
                              for each of the regions in simsearch_cube.npz
    simsearch.bed.gz(.tbi) -- Tabix-formatted file with top nDesiredMatches matches for each of regions
    """

    print("\n\n\n        Building Similarity Search Results...")

    # Bins are assumed to be 200bp or 20bp
    if determineBinSize(scoresPath) == 200:
        if windowBP == -1: windowBP = 25000
        windowBins = int(windowBP / 200)
        blockSize = determineBlockSize200(windowBP)
    elif determineBinSize(scoresPath) == 20:
        if windowBP == -1: windowBP = 2500
        windowBins = int(windowBP / 20)
        blockSize = determineBlockSize20(windowBP)
    else:
        raise ValueError("Similarity Search is only compatible with bins of size 200bp or 20bp")

    pythonFilesDir, partition, jobTag, mmMemCLO, calcMemCLO, writeMemCLO = setUpSlurm(outputDir, partition, jobTag,
                                                                                      nCores, mmMem, calcMem, writeMem)

    print("\n        STEP 1: Salient Region Selection", flush=True)
    submitMaxMeanCommand = "python {} {} {} {} {} {} {} {}".format(pythonFilesDir / "similaritySearch_max_mean.py",
                                                                   outputDir, scoresPath, windowBins, blockSize,
                                                                   windowBP, filterState, filterScore)
    mmJob = submitSlurmJob("simsearch_max_mean{}".format(jobTag), submitMaxMeanCommand, outputDir, partition, mmMemCLO,
                           "")
    print("            JobID:", mmJob, flush=True)

    print("\n        STEP 2: Similarity Search Calculation", flush=True)
    calcJobs = []
    for i in range(nJobs):
        submitCalcCommand = "python {} {} {} {} {} {} {} {}".format(pythonFilesDir / "similaritySearch_calc.py",
                                                                    outputDir, windowBins, blockSize, nCores,
                                                                    nDesiredMatches, nJobs, i)
        calcJobs.append(submitSlurmJob("simsearch_calc{}_{}".format(jobTag, i), submitCalcCommand, outputDir, partition,
                                       calcMemCLO, "--dependency=afterok:{}".format(mmJob)))
    calcJobsStr = str(calcJobs).strip('[]').replace(" ", "")
    print("            JobID:", calcJobsStr, flush=True)

    print("\n        STEP 3: Writing results", flush=True)
    submitWriteCommand = "python {} {} {} {} {} {}".format(pythonFilesDir / "similaritySearch_write.py", outputDir,
                                                           windowBins, blockSize, nJobs, nDesiredMatches)
    writeJob = submitSlurmJob("simsearch_write{}".format(jobTag), submitWriteCommand, outputDir, partition, writeMemCLO,
                              "--dependency=afterok:{}".format(calcJobsStr))
    print("            JobID:", writeJob, flush=True)

    allJobIDs = "{},{},{}".format(mmJob, calcJobsStr, writeJob)
    print("\nAll JobIDs:\n    ", allJobIDs, flush=True)


def determineBinSize(scoresPath):
    """
    Determines the size of an individual bin using the first row of the scores file

    Input:
    scoresPath -- Path to scores file to be used in similarity search

    Output:
    int(row1.iloc[0, 2] - row1.iloc[0, 1]) -- size of an individual bin
    """
    row1 = pd.read_table(scoresPath, sep="\t", header=None, usecols=[0, 1, 2], nrows=1)
    return int(row1.iloc[0, 2] - row1.iloc[0, 1])


def querySimSearch(query, simSearchPath, outputDir):
    """
    Queries the simsearch.bed.gz file a set of input regions and outputs each regions matches to an individual output
    file

    Input:
    query         -- A 2d numpy array containing the coordinates of user queried regions
    simSearchPath -- The path to the simsearch.bed.gz file to be queried for matches
    outputDir     -- The directory to output each regions matches
    """
    print("\n\n\n        Reading in data...", flush=True); readTime = time()

    # Parse query to generate array
    queryArr = generateRegionArr(query)

    # Read in matches file
    matchesDF = pd.read_table(Path(simSearchPath), sep="\t", header=None)

    print("            Time:", format(time() - readTime, '.0f'), "seconds\n", flush=True)
    print("        Querying regions...", flush=True)

    # For each query output top 100 matches to a bed file
    for chr, start, end in queryArr:
        index = np.where((matchesDF.iloc[:, 0] == chr) & (matchesDF.iloc[:, 1] >= start)
                         & (matchesDF.iloc[:, 2] <= end))[0]
        if index.size > 0:
            index = index[0]

            # Find source coords for matches
            regionChr, regionStart, regionEnd = matchesDF.iloc[index, :3]

            # write to bed file
            outfile = outputDir / "similarity_search_region_{}_{}_{}_recs.bed".format(regionChr, regionStart, regionEnd)
            with open(outfile, "w+") as f:
                matchesStr = ""

                recs = matchesDF.iloc[index, 3][2:-2]  # trim off brackets
                recs = recs.split('", "')[1:]  # split matches
                for i in range(len(recs)):
                    # Append bed formatted coords to file
                    matchesStr += "{0[0]}\t{0[1]}\t{0[2]}\n".format(recs[i].split(":"))

                f.write(matchesStr)

            # Print out information to user
            print("            Found region {}:{}-{} within user query {}:{}-{}".format(regionChr, regionStart,
                                                                                        regionEnd, chr, start, end))
            print("                See {} for matches\n".format(outfile), flush=True)
        else:
            print("            ValueError: Could not find region in given query range: {}:{}-{}\n".format(chr, start,
                                                                                                          end))


def determineBlockSize20(windowBP):
    """
    Similarity search reduces scores over search windows to mimic what the human eye might focus on.
    This is meant to provide regions which contain similar characteristics/trends
    This function determines the reduction factor for a input file which has bins of 20bp

    Input:
    windowBP - The size of the search window in BP

    Output:
    blockSize - The reduction factor for the similarity search
    """

    # We only support some set window sizes, these each have a reduction factor to
    # create an overall data size of 25 points per window
    if windowBP == 500:
        blockSize = 1
    elif windowBP == 1000:
        blockSize = 2
    elif windowBP == 2500:
        blockSize = 5
    elif windowBP == 5000:
        blockSize = 10
    elif windowBP == 7500:
        blockSize = 15
    elif windowBP == 10000:
        blockSize = 20
    else:
        raise ValueError("Error: window size must be either 500, 1000, 2500, 5000, 7500, or 10000 (in bp)")

    return blockSize


def determineBlockSize200(windowBP):
    """
    Similarity search reduces scores over search windows to mimic what the human eye might focus on.
    This is meant to provide regions which contain similar characteristics/trends
    This function determines the reduction factor for a input file which has bins of 200bp

    Input:
    windowBP - The size of the search window in BP

    Output:
    blockSize - The reduction factor for the similarity search
    """

    # We only support some set window sizes, these each have a reduction factor to
    # create an overall data size of 25 points per window
    if windowBP == 5000:
        blockSize = 1
    elif windowBP == 10000:
        blockSize = 2
    elif windowBP == 25000:
        blockSize = 5
    elif windowBP == 50000:
        blockSize = 10
    elif windowBP == 75000:
        blockSize = 15
    elif windowBP == 100000:
        blockSize = 20
    else:
        raise ValueError("Error: window size must be either 5000, 10000, 25000, 50000, 75000, or 100000 (in bp)")

    return blockSize


def setUpSlurm(outputDir, partition, jobTag, nCores, mmMem, calcMem, writeMem):
    """
    Generates strings which are used as flags in the slurm job deployment

    Input:
    outputDir -- pathlib Path to the output directory for similarity search files
    partition -- Slurm partition to use
    jobTag    -- Optional naming addition for slurm jobs
    nCores    -- Number of cores to be used for multiprocessing calculation of simsearch euclidean distance
    mmMem     -- Memory in MB to be allocated to the maxmean region selection job
    calcMem   -- Memory in MB to be allocated to the simsearch euclidean distance calculations
    writeMem  -- Memory in MB to be allocated to writing the results of the simsearch calculation

    Output:
    pythonFilesDir -- Path to the directory in which the python files to run are located
    partition      -- Slurm partition to use reformatted with slurm flag
    mmMemCLO       -- Memory to be allocated to the maxmean region selection job reformatted with slurm flag
    calcMemCLO     -- Memory to be allocated to the euclidean distance jobs reformatted with slurm flag
    writeMemCLO    -- Memory to be allocated to the writing job reformatted with slurm flag
    """
    (outputDir / ".out/").mkdir(parents=True, exist_ok=True)
    (outputDir / ".err/").mkdir(parents=True, exist_ok=True)
    print("\nSlurm .out log files are located at: {}".format(outputDir / ".out/"))
    print("Slurm .err log files are located at: {}".format(outputDir / ".err/"), flush=True)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[1] / "epilogos/"
    else:
        pythonFilesDir = (Path.cwd() / Path(__file__)).parents[1] / "epilogos/"
        print("Path generated from current working directory. May cause errors")

    partition = "--partition=" + partition if partition else ""

    jobTag = "_" + jobTag if jobTag else jobTag

    mmMemCLO = "--exclusive --mem=0" if nCores == 0 else "--ntasks={} --mem={}".format(nCores, mmMem)
    calcMemCLO = "--exclusive --mem=0" if nCores == 0 else "--ntasks={} --mem={}".format(nCores, calcMem)
    writeMemCLO = "--exclusive --mem=0" if nCores == 0 else "--ntasks={} --mem={}".format(nCores, writeMem)

    return pythonFilesDir, partition, jobTag, mmMemCLO, calcMemCLO, writeMemCLO


def submitSlurmJob(jobName, pythonCommand, outputDir, partition, memory, dependency):
    """
    Submits a slurm job to the cluster using the subprocess module

    Input:
    jobName       -- Name of the job to run
    pythonCommand -- Command to be submitted using sbatch
    outputDir     -- pathlib Path to the output directory for similarity search files
    partition     -- Slurm partition to use as string formatted with slurm flag
    memory        -- Memory to be allocated to the job as a string formatted with slurm flag
    dependency    -- Job dependencies if any

    Output:
    int(sp.stdout.split()[-1]) -- Slurm jobID number
    """
    jobOutPath = outputDir / (".out/" + jobName + ".out")
    jobErrPath = outputDir / (".err/" + jobName + ".err")

    # Creating the out and err files for the batch job
    if jobOutPath.exists():
        os.remove(jobOutPath)
    if jobErrPath.exists():
        os.remove(jobErrPath)
    try:
        jobOutPath.touch()
        jobErrPath.touch()
    except FileExistsError as err:
        # This error should never occur because we are deleting the files first
        print(err)
        return

    command = "sbatch {} --job-name={}.job --output={} --error={} {} {} --wrap='{}'".format(dependency, jobName,
                                                                                            jobOutPath, jobErrPath,
                                                                                            partition, memory,
                                                                                            pythonCommand)
    sp = subprocess.run(command, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

    if not sp.stdout.startswith("Submitted batch"):
        raise ChildProcessError("SlurmError: sbatch not submitted correctly")

    return int(sp.stdout.split()[-1])


if __name__ == "__main__":
    main()
