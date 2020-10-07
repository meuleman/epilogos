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
            print("{}\t1,2,3,{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
            print()
            print("All Others\t1,2,3,{}".format(str(otherList).strip('[]').replace(" ", "")))
        else:
            print("{}\t1,2,3,{}".format(type1, str(typeDictionary[type1]).strip('[]').replace(" ", "")))
            print()
            print("{}\t1,2,3,{}".format(type2, str(typeDictionary[type2]).strip('[]').replace(" ", "")))

if __name__ == "__main__":
    if len(sys.argv) < 3:
        main(sys.argv[1])
    elif len(sys.argv) <= 3:
        main(sys.argv[1], sys.argv[2])
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3])