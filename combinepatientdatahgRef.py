import glob
import csv
import time
import os

'''
This script converts original .bedpe patient files that contain chromosome
start and stop positions and orientation, into patient files consisting only
of gene fusions. They'll be placed in your current working directory under a new directory
called referenced_patient_data.

on ~2000 files, it takes roughly 10-12 minutes

column headers:

Gene A name, Gene A orientation, Gene B name, Gene A orientation, Gene B orientation, chrom1, chrom 1 position, chrom 1 orientation, chrom2, chrom 2 position, chrome 2 orientation 

'''


def createList(file):
    #	global reference_genome_list
    reference_genome_list = []
    for line in file:
        reference_genome_list.append(line.split('\t'))
    return reference_genome_list


# Converts chromosome strings to integer values
# eg. "chrX" -> 23
#	  "chrY" -> 24
def stripChr(chrname):
    try:
        return int(chrname.split("chr")[1])
    except ValueError:
        fix = (chrname.split("chr")[1]).split("_")[0]
        try:
            return int(fix)
        except ValueError:
            pass
        if fix != "Un":
            if fix == "X":
                return 23
            if fix == "Y":
                return 24
        else:
            return 25


'''
chromosome listing isn't sorted, but they're grouped. 
all 1s are together, followed by all 12s, etc...
Save the start position for each first occurrence of chromosomes, 
so we can just search from that position when we know what we're looking for
'''


def chromosomeLineIndexReference(List):
    # global Line_index_reference
    Line_index_reference = {}  # store start and end line numbers for the span of that chrome
    for i, line in enumerate(List):
        chrome = stripChr(line[2])
        if chrome not in Line_index_reference.keys() and len(Line_index_reference) == 0:
            Line_index_reference[chrome] = (i, '')
        elif chrome not in Line_index_reference.keys():
            if Line_index_reference[stripChr(List[i - 1][2])][1] == '':
                Line_index_reference.update(
                    {stripChr(List[i - 1][2]): (Line_index_reference[stripChr(List[i - 1][2])][0], i)})
            else:
                # hard coded
                Line_index_reference.update(
                    {stripChr(List[i - 8][2]): (Line_index_reference[stripChr(List[i - 8][2])][0], i)})
            Line_index_reference[chrome] = (i, '')

    return Line_index_reference


# scan through reference file from the chromsome number index
def run_from_position(reference_genome_list, chrom_start, current_chrom, start, stop, file):
    current_chrom = int(current_chrom)
    file.seek(start)

    if stop == '': stop = len(reference_genome_list)
    for i in reference_genome_list[start:stop]:

        chr_on_gene = stripChr(i[2])
        tran_start = int(i[4])
        tran_end = int(i[5])
        chrom_start = int(chrom_start)
        chr_on_gene = int(chr_on_gene)

        if tran_start <= chrom_start <= tran_end and current_chrom == chr_on_gene:
            return i[3], i[2], i[12], tran_start, tran_end
        if tran_start > chrom_start and current_chrom == chr_on_gene:
            return ' ', ' ', ' '

    # not sure if it ever gets here, didn't try to test it.
    # probably don't need this return
    return ' ', ' ', ' '


def build(ref_fileName):
    ref = open(ref_fileName, "r")
    reference_genome_list = createList(ref)

    Collection = []
    count = 0
    print "Crossing patient files with reference genome\nSaving files in directory 'referenced_patient_data'..."
    directory_size = len([name for name in os.listdir('./patient_data/')])
    Line_index_reference = chromosomeLineIndexReference(reference_genome_list)
    os.system("mkdir referenced_patient_data")
    percent_complete = []  # track print statements
    for name in glob.glob("./patient_data/*bedpe"):
        perComplete = int(((float(count) / directory_size) * 100))
        if perComplete % 5 == 0 and perComplete != 0:
            if perComplete not in percent_complete:
                print perComplete, "% complete"
                percent_complete.append(perComplete)

        filename = open(name, "r")
        filename.readline()
        with open("./referenced_patient_data/" + name.split('patient_data/')[1].split('.')[0] + '.tsv', "wb") as file:
            file.write(
                "geneA\tgeneA orientation\tgeneB\tgeneB orientation\tchr1\tposition\tchr1 "
                "orientation\tchr2\tposition\tchr2 orientation\n")
            for i, line in enumerate(filename):

                columns = line.split('\t')
                # orientation of chromosomes
                oriA = columns[8]
                oriB = columns[9]

                # convert chromosomes (can't use previous function because they're stored as digits in string form,
                # not "chr_")
                if columns[0] == "X" or columns[0] == "Y":
                    chr1 = 23
                else:
                    chr1 = int(columns[0])

                if columns[3] == "X" or columns[3] == "Y":
                    chr2 = 24
                else:
                    chr2 = int(columns[3])

                chr1_pos = int(columns[1])
                chr2_pos = int(columns[4])
                posA = Line_index_reference[chr1]
                posB = Line_index_reference[chr2]

                # return the required information for A and B genes
                Avalues = run_from_position(reference_genome_list, chr1_pos, chr1, posA[0], posA[1], file)

                Bvalues = run_from_position(reference_genome_list, chr2_pos, chr2, posB[0], posB[1], file)

                Aname = Avalues[2]
                Bname = Bvalues[2]

                # filter out valid gene fusions by

                if Aname != Bname and (Aname != ' ' and Bname != ' '):

                    if oriA != oriB and Avalues[0] == Bvalues[0]:
                        Collection.append((Aname, Avalues[0], Bname, Bvalues[0], Avalues[1], chr1_pos, oriA, Bvalues[1],
                                           chr2_pos, oriB))

                    if oriA == oriB and Avalues[0] != Bvalues[0]:
                        Collection.append((Aname, Avalues[0], Bname, Bvalues[0], Avalues[1], chr1_pos, oriA, Bvalues[1],
                                           chr2_pos, oriB))
            count += 1

            file.write(name)
            file.write('\n')

            for k in Collection:
                for i in k:
                    file.write(str(i) + '\t')
                file.write('\n')
            Collection = []


if __name__ == '__main__':
    start_time = time.time()

    print"Enter path to reference genome file\nenter file name if within working directory"
    path_to_ref = raw_input()
    build(path_to_ref)

    totalTime = time.time() - start_time
    secs = int(totalTime)
    mins = int(secs / 60)
    hours = int(mins / 60)
    print "--- ", hours, "hours", mins, "minutes", secs, "seconds ---"
