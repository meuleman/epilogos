import sys
from pathlib import Path
from pathlib import PurePath
import os
import subprocess
import scipy.stats as st

def main(file1, file2, observationFile, filterBool, outputDir, binEnd):
    if filterBool == "ERROR: INVALID BOOL SUBMITTED":
        print("ERROR: INVALID BOOL SUBMITTED")
        return

    print("Submitting Slurm Jobs....")

    distributions = [st.betaprime, st.halfgennorm, st.pareto, st.lomax, st.genpareto, st.gamma, 
                    st.genexpon, st.expon, st.mielke, st.exponweib, st.loglaplace, st.chi, st.chi2,
                    st.nakagami, st.burr, st.ncx2, st.pearson3]

    outputDirPath = Path(outputDir)
    jobIDArr = []

    (outputDirPath / "out/").mkdir(parents=True, exist_ok=True)
    (outputDirPath / "err/").mkdir(parents=True, exist_ok=True)

    # Finding the location of the .py files that must be run
    if PurePath(__file__).is_absolute():
        pythonFilesDir = Path(__file__).parents[0]
    else:
        pythonFilesDir = Path.cwd() / Path(__file__).parents[0]

    for distNum in range(len(distributions)):
        jobName = "{}_{}{}".format(distributions[distNum].name, file1.split("/")[-2], file2.split("/")[-2])
        jobOutPath = outputDirPath / ("out/" + jobName + ".out")
        jobErrPath = outputDirPath / ("err/" + jobName + ".err")

        # Creating the out and err files for the batch job
        if jobOutPath.exists():
            os.remove(jobOutPath)
        if jobErrPath.exists():
            os.remove(jobErrPath)
        try:
            jout = open(jobOutPath, 'x')
            jout.close()
            jerr = open(jobErrPath, 'x')
            jerr.close()
        except FileExistsError:
            # This error should never occur because we are deleting the files first
            print("ERROR: sbatch '.out' or '.err' file already exists")

        fitDistributionPy = pythonFilesDir / "fitDistribution.py"

        pythonCommand = "python {} {} {} {} {} {} {} {}".format(fitDistributionPy, file1, file2, observationFile, filterBool, distNum, binEnd, outputDir)

        slurmCommand = "sbatch --job-name={}.job --output={} --error={} --nodes=1 --ntasks=1 --mem-per-cpu=16000 --wrap='{}'".format(jobName, jobOutPath, jobErrPath, pythonCommand)

        sp = subprocess.run(slurmCommand, shell=True, check=True, universal_newlines=True, stdout=subprocess.PIPE)

        if not sp.stdout.startswith("Submitted batch"):
            print("ERROR: sbatch not submitted correctly")

        jobIDArr.append(int(sp.stdout.split()[-1]))

    jobIDStr = str(jobIDArr).strip('[]').replace(" ", "")
    print("    JobIDs:", jobIDStr)

def strToBool(string):
    if string in ["True", "true", "T", "t", "y", "Y", "yes", "Yes"]:
        return True
    elif string in ["False", "false", "F", "f", "y", "Y", "yes", "Yes"]:
        return False
    else:
        return "ERROR: INVALID BOOL SUBMITTED"

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], strToBool(sys.argv[4]), sys.argv[5], sys.argv[6])