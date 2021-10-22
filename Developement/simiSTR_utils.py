

import re

class SimiSTR_Reader:
    def __init__(self,inputBedFile):
        self.inputBedFile = inputBedFile
        #self.inputFasta = inputFasta

    # This reader collects all bed file entrances in a dictionary for easier access
    def getBedFile(self):
        bedfile_d = dict()
        with open(self.inputBedFile, 'r') as inBedFile:
            for bedline in inBedFile:
                """ 
                Here you can first filter out the header information that can contain
                '#', 'track' or 'browser' (http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1)
                """
                if not self.__isBedHeader(bedline):
                    bedElements = bedline.split("\t")
                    # First three are mandatory, fourth is the sequence (name), we need at least 4
                    if len(bedElements) >= 5:
                        [chromosome, startPosition, endPosition, motifLength, motif] = bedElements[:5]
                        # sequence in CAPs
                        motif = motif.strip().upper() #this is needed as motif should be a string without "\n" or else
                        chromosomeNr = self.__getChromosomeNumber(chromosome)
                        if chromosomeNr not in bedfile_d:
                            bedfile_d[chromosomeNr] = [chromosome, startPosition, endPosition, motifLength, motif]
                        else:
                            bedfile_d[chromosomeNr].append([chromosome, startPosition, endPosition, motifLength, motif])
            return bedfile_d

    def __isBedHeader(self, line):
        # from https://en.wikipedia.org/wiki/BED_(file_format)#Header
        # and http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1
        if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            return True

    def __getChromosomeNumber(self, chromosome):
        chrNr = 0
        if type(chromosome) is not int:
            print(chromosome)
            chr = [int(s) for s in re.findall(r'\d+',chromosome)]
            if chr is not []:
                chrNr = chr[0]
            else:
                print("Check the first column in your assigned input bed file!")
        else:
            chrNr= chromosome
        return chrNr


class SimiSTR_Writer:
    def __init__(self, outputBedfile):
        self.outputBed = outputBedfile
        # self.outputFa = outputFasta

    # write out the new coordinates of the adapted bed file.
    def printBedModifications(self, bedfile_l_copy):
        with open(self.outputBed, 'w') as outBedfile:
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

