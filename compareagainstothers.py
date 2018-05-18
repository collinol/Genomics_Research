import re
import glob


# This script finds all the patient files of the same cancer type
# where the genes exist in the patient-in-question's gene-fusions, but not
# within the fusions of a same-cancer patient.

# Builds a dictionary of
# key = one patient
# value = (other patient file, missing-from-fusions gene)

# Needed as a reference for the RNA data in the individual patient to
# determine the difference in RNA scores

def UniqueGenes(patientDirectory, cancerType):
    MasterGenes = {}
    for name in glob.glob("./cancer_type_groups/" + cancerType + "/*"):
        fileheader = name.split("groups/" + cancerType + "/")[1]
        filename = open(name, 'r')
        patientDict = {}
        GeneA = 0
        GeneB = 0
        countGenes = 0
        start = 0
        exonsToSort = []
        for i, line in enumerate(filename):
            if line.split('\t')[0] == 'geneA':
                prevLine = i
            if i == prevLine + 1:
                GeneA = line.split('\t')[0]
                GeneB = line.split('\t')[1].split('\n')[0]
            else:
                line = line.split('\t')

                if len(line) > 2:
                    if line[1] == GeneA:
                        exonsToSort = []
                        exonsToSort.append((line[0], line[1], line[2], int(line[3].split('\n')[0])))

                    if line[1] == GeneB:
                        exonsToSort = []
                        exonsToSort.append((line[0], line[1], line[2], int(line[3].split('\n')[0])))
                if "chr" in line[0]:
                    exonsToSort.append(line[1].split('\n')[0])
                    countGenes += 1
                    start = i + 1
            if i >= start:
                try:
                    int(line[0])
                    exonsToSort.append((int(line[0]), int(line[1]), int(line[2].split('\n')[0])))
                except Exception as e:
                    pass
            if line[0] == '\n':
                patientDict[countGenes] = sorted(exonsToSort, key=lambda x: x[0])

        MasterGenes[fileheader] = (
        [int(re.sub('[^0-9]', '', v[-1][0])) for v in patientDict.values()], [v[-1][1] for v in patientDict.values()])
    CompareTo = {}
    for currentKey, currentValue in MasterGenes.iteritems():
        compare = []
        for checkKey, checkValue in MasterGenes.iteritems():

            if currentKey == checkKey: continue
            for i, gene in enumerate(currentValue[0]):
                if gene not in checkValue[0]:
                    if (checkKey, gene) not in compare:
                        compare.append((checkKey, gene, currentValue[1][i]))

        CompareTo[currentKey] = compare

    return CompareTo
