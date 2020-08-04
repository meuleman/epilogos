import gzip
import numpy as np
import sys
from pathlib import Path
import math
import pandas as pd
import time

def main(filename, numStates, saliency, outputDirectory, columnSpecification = "All"):
    dataFilePath = Path(filename)
    outputDirPath = Path(outputDirectory)

    #Determine the number of columns in the data
    nColsData = len(np.loadtxt(fname = dataFilePath, dtype="str", max_rows= 1)) - 3

    # Read in the data
    t0 = time.time()
    dataDF = pd.read_table(dataFilePath, header = None, sep="\t", usecols=columnSpecificationAsTuple(columnSpecification, nColsData))
    numRows, numCols = dataDF.shape
    print("Time to read in: ", time.time() - t0)

    t1 = time.time()
    dataArr = dataDF.to_numpy()
    print("Time to convert to numpy: ", time.time() - t1)
    print("Total reading time: ", time.time() - t0)
    print(dataArr)


    # whitespaceTimes = []
    # regexWhitespaceTimes = []
    # tabTimes = []

    # for i in range(5):
    #     t0 = time.time()
    #     # Read in the data
    #     dataDF = pd.read_table(dataFilePath, header = None, sep=r"\s+", usecols=columnSpecificationAsTuple(columnSpecification, nColsData))
    #       numRows, numCols = dataDF.shape

    #     t1 = time.time() - t0
    #     regexWhitespaceTimes.append(t1)

    # for i in range(5):
    #     t0 = time.time()
    #     # Read in the data
    #     dataDF = pd.read_table(dataFilePath, header = None, delim_whitespace=True, usecols=columnSpecificationAsTuple(columnSpecification, nColsData))
    #       numRows, numCols = dataDF.shape

    #     t1 = time.time() - t0
    #     whitespaceTimes.append(t1)

    # for i in range(5):
    #     t0 = time.time()
    #     # Read in the data
    #     dataDF = pd.read_table(dataFilePath, header = None, sep="\t", usecols=columnSpecificationAsTuple(columnSpecification, nColsData))
    #       numRows, numCols = dataDF.shape

    #     t1 = time.time() - t0
    #     tabTimes.append(t1)

    # print(regexWhitespaceTimes, "  Average: ", sum(regexWhitespaceTimes) / 5)
    # print(whitespaceTimes, "  Average: ", sum(whitespaceTimes) / 5)
    # print(tabTimes, "  Average: ", sum(tabTimes) / 5)


    # Calculate the expected frequencies of each state
    t0 = time.time()
    stateIndices = list(range(1, 16))
    expFreqSeries = pd.Series(np.zeros(15), index=stateIndices)
    # Raw counts
    for i in range(3, numCols):
        stateCounts = dataDF[i].value_counts()
        for state, count in stateCounts.items():
            expFreqSeries.loc[state] += count
    # Computing the frequencies
    dfSize = expFreqSeries.sum()
    for state, count in expFreqSeries.items():
        expFreqSeries.loc[state] = count / dfSize

    print("Time to calculate exp freq: ", time.time() - t0)
    print(expFreqSeries)

    # Calculate the observed frequencies and final scores in one loop
    t0 = time.time()
    expFreqArr = expFreqSeries.to_numpy()
    scoreArr = np.zeros((numRows, numStates))
    for row in range(numRows):
        uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
        for i in range(len(uniqueStates)):
            # Function input is obsFreq and expFreq
            scoreArr[row, int(uniqueStates[i]) - 1] = klScore(stateCounts[i] / nColsData, expFreqArr[int(uniqueStates[i]) - 1])
    print("time for score values (one loop): ", time.time() - t0)

    # Writing the data to the files
    observationsTxtPath = outputDirPath / "observations.txt"
    scoresTxtPath = outputDirPath / "scores.txt"

    observationsTxt = open(observationsTxtPath, "w")
    scoresTxt = open(scoresTxtPath, "w")

    # Write each row in both observations and scores
    for i in range(numRows):
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
        for j in range(numStates):
            scoresTxt.write("%.5f " % scoreArr[i, j])
        scoresTxt.write("\n")

    observationsTxt.close()
    scoresTxt.close()



    # # Calculate the observed frequencies and final scores in one loop
    # t0 = time.time()
    # scoreArr2 = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
    #     for i in range(len(uniqueStates)):
    #         # Function input is obsFreq and expFreq
    #         scoreArr2[row, int(uniqueStates[i]) - 1] = klScore(stateCounts[i] / nColsData, expFreqSeries.values[int(uniqueStates[i]) - 1])
    # print("time for score values (one loop): ", time.time() - t0)

    # # Calculate the observed frequencies and final scores in two loops
    # t0 = time.time()
    # obsFreqArr2 = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
    #     for i in range(len(uniqueStates)):
    #         obsFreqArr2[row, int(uniqueStates[i]) - 1] = stateCounts[i] / nColsData
    # print("Time for score values (two loops | loop one): ", time.time() -t0)
    # t1 = time.time()
    # # Calculate the KL-divergence scores for each state in each location
    # scoreArr3 = np.zeros((numRows, numStates))
    # for i in range(numRows):
    #     for j in range(numStates):
    #         scoreArr3[i, j] = klScore(obsFreqArr2[i, j], expFreqArr[j])
    # print("Time for score values (two loops | loop two): ", time.time() - t1)
    # print("Time for score values (two loops): ", time.time() - t0)


    # Calculate the observed frequencies and final scores in two loops
    # t0 = time.time()
    # obsFreqArr = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
    #     for i in range(len(uniqueStates)):
    #         obsFreqArr[row, int(uniqueStates[i]) - 1] = stateCounts[i] / nColsData
    # print("Time for score values (two loops | loop one): ", time.time() -t0)
    # t1 = time.time()
    # # Calculate the KL-divergence scores for each state in each location
    # scoreArr = np.zeros((numRows, numStates))
    # for i in range(numRows):
    #     for j in range(numStates):
    #         scoreArr[i, j] = klScore(obsFreqArr[i, j], expFreqSeries.iloc[j])
    # print("Time for score values (two loops | loop two): ", time.time() - t1)
    # print("Time for score values (two loops): ", time.time() - t0)

    # # Calculate the observed frequencies and final scores in one loop
    # t0 = time.time()
    # scoreArr2 = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
    #     for i in range(len(uniqueStates)):
    #         # Function input is obsFreq and expFreq
    #         scoreArr2[row, int(uniqueStates[i]) - 1] = klScore(stateCounts[i] / nColsData, expFreqSeries.iloc[int(uniqueStates[i]) - 1])
    # print("time for score values (one loop): ", time.time() - t0)

    # print((scoreArr == scoreArr2).all())


    # Calculate the observed frequences of each state for each row
    # obsFreqDF = pd.DataFrame(index=list(range numRows)), columns=stateIndices)
    # dataDFTranspose = dataDF.transpose()
    # for row in range  numRows):
    #     stateCounts = dataDFTranspose.iloc[3:, row].value_counts()
    #     for state, count in stateCounts.items():
    #         obsFreqDF.loc[row, state] = count / nColsData
    # obsFreqDF.fillna(0)

    # t0 = time.time()
    # for index, row in dataDF.iterrows():
    #     if index >= 100:
    #         break
    #     stateCounts = row.value_counts()
    #     for state, count in stateCounts.items():
    #         obsFreqDF.loc[index, state] = count / nColsData
    # print("OBS value time: ", numRows/100*(time.time()-t0))

    # print(obsFreqDF.head())

    # t100 = time.time()
    # dataDF.apply(lambda row: calculateObsFreq(row[3:]), axis=1, result_type="expand")
    # print("Apply time: ", time.time() - t100)
    # print(dataDF.head())

    # t0 = time.time()
    # stateCounts = dataDF.iloc[1, 3:].value_counts()
    # print("Obs Value Counts Time: ",  numRows*(time.time() - t0))
    # t1 = time.time()
    # for state, count in stateCounts.items():
    #     t2 = time.time()
    #     obsFreqDF.loc[1, state] = count / nColsData
    #     print("    Freq Calc Time: ", numRows*(time.time() - t2))
    # print("Obs Freq loop time: ", numRows*(time.time() - t1))
    # print("Obs Freq Total Time: ", numRows*(time.time() - t0))
    
    # for row in range(numRows):
    #     t0 = time.time()
    #     stateCounts = dataDF.iloc[row, 3:].value_counts()
    #     print("ValueCounts Time: ",time.time() - t0)
    #     t0 = time.time()
    #     for state, count in stateCounts.items():
    #         t1 = time.time()
    #         obsFreqDF.loc[row, state] = count / nColsData
    #         print("    Freq Calc Time: ", time.time() - t0)
    #     print("Freq loop time: ", time.time() - t0)
    
    # print(obsFreqDF.head())
        



    # # Calculate the observed frequencies of each state for each row
    # obsFreqArr = np.zeros((numRows, numStates))
    # for row in range(numRows):
    #     uniqueStates, stateCounts = np.unique(dataArr[row,3:], return_counts=True)
    #     for i in range(len(uniqueStates)):
    #         obsFreqArr[row, int(uniqueStates[i]) - 1] = stateCounts[i] / numCols

    # print(obsFreqArr)
        
    # # Calculate the KL-divergence scores for each state in each location
    # scoreArr = np.zeros((numRows, numStates))
    # for i in range(numRows):
    #     for j in range(numStates):
    #         scoreArr[i, j] = logBase2(obsFreqArr[i, j], expFreqArr[j])


# def calculateObsFreq(row):
#     obsFreqSeries = pd.Series(np.zeros(15), index=list(range(1, 16)))
#     # Computing the frequencies
#     stateCounts = row.value_counts()
#     for state, count in stateCounts.items():
#         obsFreqSeries.loc[state] = count / row.size
#     return obsFreqSeries

# Helper to calculate KL-score (used because math.log2 errors out if obsFreq = 0)
def klScore(obs, exp):
    if obs == 0.0:
        return 0.0
    else:
        return obs * math.log2(obs / exp)

# Helper to create a sequence from the inputed column specification
def columnSpecificationAsTuple(columnSpecification, ncols):
    columnListExpanded = [0, 1, 2]
    if columnSpecification == "All":
        # Include all columns except the first 3
        columnListExpanded = list(range(0, ncols + 3))
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