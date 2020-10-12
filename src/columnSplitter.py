import sys
from pathlib import Path
import random

def main(filename, randomSampling=False, type1="", type2=""):
    if type(randomSampling) != type(True):
        return

    file = Path(filename)
    with open(file, 'r') as f:
        lines = f.readlines()
        genomeNumber = 4
        typeDictionary = {}
        for line in lines:
            genomeType = line.split("|")[-1].strip()
            if genomeType not in typeDictionary:
                typeDictionary[genomeType] = [genomeNumber]
            else:
                typeDictionary[genomeType].append(genomeNumber)
            genomeNumber += 1

    for k in typeDictionary:
        print("{}, {}\t\t{}".format(k, len(typeDictionary[k]),typeDictionary[k]))

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
                print("{}\t1,2,3,{}".format(type1, str(shuffledList[:len(typeDictionary[type1])]).strip('[]').replace(" ", "")))
                print()
                print("{}\t1,2,3,{}".format(type2, str(shuffledList[len(typeDictionary[type1]):]).strip('[]').replace(" ", "")))
            else:
                print("{}\t1,2,3,{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
                print()
                print("{}\t1,2,3,{}".format(type2, str(typeDictionary[type2]).strip('[]').replace(" ", "")))

def strToBool(string):
    if string in ["true", "t", "T", "y", "Y", "True", "yes", "Yes"]:
        return True
    elif string in ["false", "f", "F", "n", "N", "False", "no", "No"]:
        return False
    else:
        print("ERROR: STRTOBOOL FAILED")
        return "ERROR: STRTOBOOL FAILED"


if __name__ == "__main__":
    if len(sys.argv) <= 3:
        main(sys.argv[1])
    elif len(sys.argv) <= 4:
        main(sys.argv[1], strToBool(sys.argv[2]), sys.argv[3])
    else:
        main(sys.argv[1], strToBool(sys.argv[2]), sys.argv[3], sys.argv[4])