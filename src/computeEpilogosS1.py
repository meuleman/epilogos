import gzip
import numpy as np
import sys
from pathlib import Path
import math

def main(filename, numstates, saliency, outputDirectory, columnSpecification = "All"):
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    #Determine the number of columns in the data
    # with open(dataFilePath, "r") as f:
    #     nColsData = len(f.readline().split())

    nColsData = len(np.loadtxt(fname = dataFilePath, dtype="str", max_rows= 1))

    # Read in the data
    dataArr = np.loadtxt(fname = dataFilePath, dtype="str", usecols=columnSpecificationAsTuple(columnSpecification, nColsData))
    num_rows, num_cols = dataArr.shape

    # Calculate the expected frequencies of each state
    uniqueStates, stateCounts = np.unique(dataArr[:,3:], return_counts=True)
    expFreqArr = np.zeros(numstates)
    for i in range(len(uniqueStates)):
        expFreqArr[int(uniqueStates[i]) - 1] = stateCounts[i] / dataArr[:,3:].size

    print(expFreqArr)

    # Calculate the observed frequencies of each state for each row
    obsFreqArr = np.zeros((num_rows, numstates))
    for row in range(num_rows):
        uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
        for i in range(len(uniqueStates)):
            obsFreqArr[row, int(uniqueStates[i]) - 1] = stateCounts[i] / num_cols

    print(obsFreqArr)
        
    # Calculate the KL-divergence scores for each state in each location
    scoreArr = np.zeros((num_rows, numstates))
    for i in range(num_rows):
        for j in range(numstates):
            scoreArr[i, j] = logBase2(obsFreqArr[i, j], expFreqArr[j])


    # Writing the data to the files
    observationsTxtPath = outputDirPath / "observations.txt"
    scoresTxtPath = outputDirPath / "scores.txt"

    observationsTxt = open(observationsTxtPath, "w")
    scoresTxt = open(scoresTxtPath, "w")

    # Write each row in both observations and scores
    for i in range(num_rows):
        # Write in the coordinates
        locationData = dataArr[i, :3]
        for location in locationData:
            observationsTxt.write("%s " % location)
            scoresTxt.write("%s " % location)
        
        # Write to observations
        maxContribution = np.amax(scoreArr[i])
        maxContributionLoc = np.argmax(scoreArr[i]) + 1
        totalScore = np.sum(scoreArr[i])

        observationsTxt.write("%d " % maxContributionLoc)
        observationsTxt.write("%.5f " % maxContribution)
        observationsTxt.write("1 ")
        observationsTxt.write("%.5f \n" % totalScore)

        # Write to scores
        for j in range(numstates):
            scoresTxt.write("%.5f " % scoreArr[i, j])
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()



# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def logBase2(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to create a sequence from the inputed column specification
def columnSpecificationAsTuple(columnSpecification, ncols):
    columnListExpanded = []
    if columnSpecification == "All":
        # Include all columns except the first 3
        columnListExpanded = list(range(0, ncols))
    else:
        columnList = columnSpecification.split(',')
        for str in columnList:
            # Expand the hyphenated parts into a full list (e.g. "33-36" -> [33, 34, 35, 36])
            # Add 3 because of 3 not data columns in file
            if "-" in str:
                strList = str.split("-")
                columnListExpanded.append(list(range(int(strList[0]) + 3, int(strList[1] + 1 + 3))))
            else:
                columnListExpanded.append(int(str) + 3)
    return tuple(columnListExpanded)

if __name__ == "__main__":
    # Checking that the arguments are all correct
    if len(sys.argv) - 1 < 4:
        # Argument info if wrong number
        print("You must provide at least 4 arguments:\n")
        print("   1. Filename to read from\n")
        print("   2. Number of states in chromatin state model\n")
        print("   3. Saliency metric (1-3)\n")
        print("   4. Output directory\n")
        print("   5. (Optional) Specification of columns to read (all by default)\n")
    elif int(sys.argv[2]) != 15:
        print("We currently only offer support for a 15-state chromatin state model")
    elif int(sys.argv[3]) != 1:
        print("We currently only offer support for a saliency of one")
    elif len(sys.argv) - 1 == 4:
        # if the last argument is not used use the default
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])
    else:
        main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5])