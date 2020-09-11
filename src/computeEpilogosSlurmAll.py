from pathlib import Path
import click
import subprocess

@click.command()
@click.option("-f", "--file", "combinationsFile", type=str, required=True, help="Path to file which contains all the combination of data (tab seperated) which you want to be calculated by epilogos\n(e.g. ROADMAP\thg19\t15\tFemale_donors\tS1)")
@click.option("-d", "--parent-directory", "parentDirectory", type=str, required=True, help="Path to the parent directory whose subdirectories contain all the data you want to be calculated by epilogos\n(e.g. PARENT_DIRECTORY/ROADMAP/hg19/15/Female_donors/S1/)")
@click.option("-o", "--output-directory", "outputDirectory", type=str, required=True, help="Path to the directory in which you want the outputed epilogos calculations to go")
def main(combinationsFile, parentDirectory, outputDirectory):
    outputDirPath = Path(outputDirectory)

    # Finding the location of the .py files that must be run
    if Path(__file__).is_absolute:
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    with open(combinationsFile, "r") as cf:
        line = cf.readline()
        while line:
            if not line.strip() == "":
                fileDirectory = Path(parentDirectory) / line.replace("\t", "/")
                jobName = line.replace("\t", "_")
                jobOutPath = outputDirPath / (".out/" + jobName + ".out")
                jobErrPath = outputDirPath / (".err/" + jobName + ".err")

                # Creating the out and err files for the batch job
                try:
                    jout = open(jobOutPath, "x")
                    jout.close()
                    jerr = open(jobErrPath, "x")
                    jerr.close()
                except FileExistsError:
                    print("ERROR: sbatch '.out' or '.err' file already exists")

                pythonFile = pythonFilesDir / "computeEpilogosSlurm.py"

                # Input for the python file is fileDirectory, state model, saliency metric, outputDirectory
                pythonCommand = "python {} -f {} -m {} -s {} -o {}".format(pythonFile, fileDirectory, line.split("\t")[2], line.split("\t")[4][1], outputDirectory)

                slurmCommand = "sbatch --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)

                print(pythonCommand)
                print(slurmCommand)

                subprocess.run(slurmCommand, shell=True)

if __name__ == "__main__":
    main()