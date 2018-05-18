import glob
import os


# uses gene annotation file to find all
# exons and exon locations within the genes listed on
# patient fusions

def geneIdRef(filename):
    GeneFile = filename
    Genes = {}
    previousChrom = ""
    count = 0
    keyname = ""
    genename = ''
    exList = []
    excount = 0
    print "running gene file"
    for line in GeneFile:

        if line[0] == 'c':
            gene = line.split('\t')[8].split(';')[0].split('"')[1].split('.')[0]
            nametype = line.split('\t')[1]
            chrom = line.split('\t')[0]
            strandtype = line.split('\t')[2]

            if strandtype == "gene":
                if keyname != "":
                    if len(exList) > 0:
                        Genes[keyname] = (genename, previousChrom, excount, exList)
                    genename = line.split('\t')[8].split(';')[4].split('"')[1]

                keyname = gene
                exList = []
                excount = 0
            else:

                # I'm not sure what the nametypes are used for
                # Theres HAVANA and ENSEMBLE in this file, some of the exons repeat
                # within each type and some don't. This choice is arbitrary for now

                if nametype == "HAVANA":

                    transcript_type = line.split('\t')[8].split(';')[5].split('"')[1]
                    transcript_status = line.split('\t')[8].split(';')[6].split('"')[1]

                    if strandtype == "exon" and transcript_type == "protein_coding" and transcript_status == "KNOWN":
                        exonnumber = line.split('\t')[8].split(';')[8].split(' ')[2]

                        exList.append((exonnumber, line.split('\t')[3], line.split('\t')[4]))
                        excount += 1

            count += 1
            previousChrom = chrom

    return Genes


# These next two functions aren't needed for building results,
# but I wrote it to examine what kind of exon results I was getting from
# the file.
'''
def getExonFreq(filename):
    print "running exons"
    ExonsWithCount = {}
    previousExon = ""
    for line in ExonFile:
        exon = line.split('\t')[0]
        main = exon.split('.')[0]
        if main != previousExon:
            ExonsWithCount[main] = 1
        else:
            ExonsWithCount[main] += 1
        previousExon = main
    return ExonsWithCount


def checkExonNumbers(Genedata):
    Genes = Genedata
    totalgenes = 0
    totalgenes_with_exonerror = 0
    totalExons_with_error = 0
    total_exons_on_errorgene = 0
    percents = []
    Glist = []

    huh = 0
    with open("Gene_exon_errors.tsv", "w") as Geneswitherror:
        for k, v in Genes.iteritems():
            totalgenes += 1
            exonChecker = {}
            flag = 0
            exonError_specific = 0
            Wrong = []
            for exon in sorted(v[3], key=lambda x: x[0]):

                if exon[0] not in exonChecker.keys():
                    exonChecker[exon[0]] = (int(exon[1]), int(exon[2]))
                    if v[0] not in Wrong: Wrong.append(v[0])
                    Wrong.append(exon)
                else:
                    if int(exon[1]) != int(exonChecker[exon[0]][0]) or int(exon[2]) != int(exonChecker[exon[0]][1]):
                        if v[0] not in Wrong: Wrong.append(v[0])
                        Wrong.append(exon)
                        totalExons_with_error += 1
                        exonError_specific += 1
                        flag = 1
            for x in Wrong:
                if type(x) != tuple:
                    Geneswitherror.write(x + '\n')
                else:
                    Geneswitherror.write(x[0] + '\t' + x[1] + '\t' + x[2] + '\n')
            if flag == 1:
                Glist.append(v[0])

'''


def RNADonorId(ExonAnnotationFile, patient_directory):
    total = 0
    Genes = geneIdRef(ExonAnnotationFile)
    # create new patient_fusion files, this time with the gene fusions
    # containing the exons annotated along the gene
    os.system("mkdir patient_fusions_with_exons")
    for name in glob.glob("./" + patient_directory + "/*"):
        filename = open(name, "r")
        name = name.split("data/")[1]
        print name

        with open("./patient_fusions_with_exons/" + name, 'w') as newFile:

            filename.next()
            filename.next()
            for line in filename:
                total += 1
                geneA = line.split('\t')[0]
                geneB = line.split('\t')[2]
                breakpoint1 = int(line.split('\t')[5])
                breakpoint2 = int(line.split('\t')[8])

                matched = 0
                newFile.write("geneA" + '\t' + "geneB" + '\n')
                newFile.write(geneA + '\t' + geneB + '\n')
                for k in Genes.keys():
                    v = Genes[k]
                    if v[0] == geneA:
                        matched += 1
                        newFile.write(
                            k + '\t' + v[0] + '\t' + str(line.split('\t')[1]) + '\t' + str(breakpoint1) + '\n')
                        newFile.write(v[1] + '\t' + line.split('\t')[6] + '\n')
                        for j in v[3]:
                            newFile.write(j[0].split('"')[0] + '\t' + j[1] + '\t' + j[2] + '\n')
                        newFile.write('\n')

                    elif v[0] == geneB:
                        matched += 1
                        newFile.write(k + '\t' + v[0] + '\t' + line.split('\t')[3] + '\t' + str(breakpoint2) + '\n')
                        newFile.write(v[1] + '\t' + line.split('\t')[9] + '\n')
                        for j in v[3]:
                            newFile.write(j[0].split('"')[0] + '\t' + j[1] + '\t' + j[2] + '\n')
                        newFile.write('\n')
