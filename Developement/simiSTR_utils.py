

import re

class Reader:
    def __init__(self, inputBedFile):
        self.inputBedFile = inputBedFile
        self.bedfile_d = dict()

    # This reader collects all bed file entrances in a dictionary for easier access
    def getBedFile(self, inputBedFile):
        bedfile_d = dict()
        with open(inputBedFile, 'r') as inBedFile:
            for bedline in inBedFile:
                """ 
                Here you can first filter out the header information that can contain
                '#', 'track' or 'browser' (http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1)
                """
                if not self.isBedHeader(bedline):
                    bedElements = bedline.split("\t")
                    # First three are mandatory, fourth is the sequence (name), we need at least 4
                    if len(bedElements) >= 4:
                        [chromosome, startPosition, endPosition, motifLength, motif] = bedElements[:4]
                        bedElements[4] = bedElements[4].strip() #this is needed as motif should be a string without "\n" or else
                        chromosomeNr = self.getChromosomeNumber(chromosome)
                        # sequence in CAPs
                        sequence = sequence.upper()
                        if chromosomeNr not in self.bedfile_d:
                            self.bedfile_d[chromosomeNr] = [bedElements]
                        else:
                            self.bedfile_d[chromosomeNr].append(bedElements)  # chr    from Pos    to Pos      lenMotif    motif
            return self.bedfile_d

    def isBedHeader(self, line):
        # from https://en.wikipedia.org/wiki/BED_(file_format)#Header
        # and http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1
        if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            return True

    def getChromosomeNumber(self, chromosome):
        chrNr = 0
        if type(chromosome) is not int:
            chr = re.search("(\d*)", chromosome)
            if chr is not None:
                chrNr = chr
            else:
                print("Check the first column in your assigned input bed file!")
        else:
            chrNr= chromosome
        return chrNr


class Writer:
    # write out the new coordinates of the adapted bed file.
    def printBedModifications(self, bedfile_l_copy, newBedFile):
        with open(newBedFile, 'w') as outBedfile:
            count = 0
            for chr in bedfile_l_copy:
                for line in chr:  # not recording bad mathces in new bedfile
                    if line[1] != 0 and line[2] != 0:  # thats the ones where pattern was not found
                        lines = str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" \
                                + str(line[4]) + "\t" + str(line[5] + "\n")
                        outBedfile.write(lines)
                    else:
                        count += 1
            print(count)