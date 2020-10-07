import sys
from pathlib import Path

def main(filename, type1="", type2=""):
    file = Path(filename)
    with open(file, 'r') as f:
        lines = f.readlines()
        genomeNumber = 4
        typeDictionary = {}
        for line in lines:
            genomeType = line.strip().split()[-1]
            if genomeType not in typeDictionary:
                typeDictionary[genomeType] = [1, 2, 3, genomeNumber]
            else:
                typeDictionary[genomeType].append(genomeNumber)
            genomeNumber += 1

    for k in typeDictionary:
        print("{}, {}\t\t{}".format(k, len(typeDictionary[k]),typeDictionary[k]))

    if type1 != "" and type2 != "":
        print()
        print()

        if type2 == "allOthers":
            otherList = []
            for k in typeDictionary:
                if k != type1:
                    otherList.extend(typeDictionary[k])
            print("{}\t{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
            print()
            print("{}\t{}".format(type2, str(otherList).strip('[]').replace(" ", "")))
        else:
            print("{}\t{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
            print()
            print("{}\t{}".format(type2, str(typeDictionary[type2]).strip('[]').replace(" ", "")))

if __name__ == "__main__":
    if len(sys.argv) < 4:
        main(sys.argv[1])
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])