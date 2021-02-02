import sys
from pathlib import Path
import random
import pandas as pd

def main(filename, randomSampling=False, equalSize=-1, column="Group", type1="", type2=""):
    if type(randomSampling) != type(True):
        return

    # Taking in Data
    filePath = Path(filename)
    dataDF = pd.read_table(filePath, header=0, sep="\t", keep_default_na=False)
    dataDF.rename(columns=str.lower, inplace=True)
    column = column.lower()
    type1 = type1.lower()
    type2 = type2.lower()

    if column == "":
        # If a column was not given, give back a list of all the columns
        print("COLUMNS:")
        print("-------------------------")
        for col in dataDF.columns: 
            print(col)
        return
    else:
        # Getting the column numbers for each of the group types
        typeDictionary = {}
        idDictionary = {}
        for i in range(len(dataDF[column])):
            if type(dataDF[column][i]) == type("HI"):
                typeVal = dataDF[column][i].lower()
            else:
                typeVal = dataDF[column][i]
            if typeVal not in typeDictionary:
                typeDictionary[typeVal] = [i + 4]
            else:
                typeDictionary[typeVal].append(i + 4)

            if typeVal not in idDictionary:
                idDictionary[typeVal] = [dataDF.iloc[i, 0]]
            else:
                idDictionary[typeVal].append(dataDF.iloc[i, 0])

        for k in typeDictionary:
            print("{} (Size={}): \t\t{}".format(k, len(typeDictionary[k]),typeDictionary[k]))

    if type1 != "":
        print()
        print()

        if type2 == "":
            otherList = []
            for k in typeDictionary:
                if k != type1:
                    otherList.extend(typeDictionary[k])
            if randomSampling:
                shuffledList = typeDictionary[type1] + otherList
                random.shuffle(shuffledList)
                if equalSize == 0:
                    allOtherLen = 0
                    for k in typeDictionary:
                        if k != type1:
                            allOtherLen += len(typeDictionary[k])
                    size = min((len(typeDictionary[type1]), allOtherLen))
                    print("{}\t1,2,3,{}".format(type1, str(shuffledList[:size]).strip('[]').replace(" ", "")))
                    print()
                    print("All Others\t1,2,3,{}".format(str(shuffledList[size:2*size]).strip('[]').replace(" ", "")))
                elif equalSize > 0:
                    print("{}\t1,2,3,{}".format(type1, str(shuffledList[:equalSize]).strip('[]').replace(" ", "")))
                    print()
                    print("All Others\t1,2,3,{}".format(str(shuffledList[equalSize:2*equalSize]).strip('[]').replace(" ", "")))
                else:
                    print("{}\t1,2,3,{}".format(type1, str(shuffledList[:len(typeDictionary[type1])]).strip('[]').replace(" ", "")))
                    print()
                    print("All Others\t1,2,3,{}".format(str(shuffledList[len(typeDictionary[type1]):]).strip('[]').replace(" ", "")))
            else:
                print("{}\t1,2,3,{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
                print()
                print("All Others\t1,2,3,{}".format(str(otherList).strip('[]').replace(" ", "")))
        else:
            if randomSampling:
                shuffledList = typeDictionary[type1] + typeDictionary[type2]
                random.shuffle(shuffledList)
                if equalSize == 0:
                    size = min((len(typeDictionary[type1]), len(typeDictionary[type2])))
                    print("{}\t1,2,3,{}".format(type1, str(shuffledList[:size]).strip('[]').replace(" ", "")))
                    print()
                    print("{}\t1,2,3,{}".format(type2, str(shuffledList[size:2*size]).strip('[]').replace(" ", "")))
                elif equalSize > 0:
                    print("{}\t1,2,3,{}".format(type1, str(shuffledList[:equalSize]).strip('[]').replace(" ", "")))
                    print()
                    print("{}\t1,2,3,{}".format(type2, str(shuffledList[equalSize:2*equalSize]).strip('[]').replace(" ", "")))
                else:
                    print("{}\t1,2,3,{}".format(type1, str(shuffledList[:len(typeDictionary[type1])]).strip('[]').replace(" ", "")))
                    print()
                    print("{}\t1,2,3,{}".format(type2, str(shuffledList[len(typeDictionary[type1]):]).strip('[]').replace(" ", "")))
            else:
                print("{}\t1,2,3,{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
                print()
                print("{}\t{}".format(type1, str(idDictionary[type1]).strip('[]').replace(" ", "")))
                print()
                print("{}\t1,2,3,{}".format(type2, str(typeDictionary[type2]).strip('[]').replace(" ", "")))
                print()
                print("{}\t{}".format(type2, str(idDictionary[type2]).strip('[]').replace(" ", "")))



# Helper to convert taken in string into a boolean 
def strToBool(string):
    if string in ["true", "t", "T", "y", "Y", "True", "yes", "Yes"]:
        return True
    elif string in ["false", "f", "F", "n", "N", "False", "no", "No"]:
        return False
    else:
        print("ERROR: STRTOBOOL FAILED")
        return "ERROR: STRTOBOOL FAILED"


if __name__ == "__main__":
    if len(sys.argv) <= 2:
        main(sys.argv[1], column="")
    elif len(sys.argv) <= 5:
        main(sys.argv[1], column=str(sys.argv[-1]))
    elif len(sys.argv) <= 6:
        main(sys.argv[1], randomSampling=strToBool(sys.argv[2]), equalSize=int(sys.argv[3]), column=sys.argv[4], type1=sys.argv[5])
    else:
        main(sys.argv[1], randomSampling=strToBool(sys.argv[2]), equalSize=int(sys.argv[3]), column=sys.argv[4], type1=sys.argv[5], type2=sys.argv[6])