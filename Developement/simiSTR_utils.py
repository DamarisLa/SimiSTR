

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
                        expBaseNrchange = 0  # nr of bases changed through expansion change
                        noOfSubstitution = 0
                        nrOfIndels = [0,0]
                        if chromosomeNr not in bedfile_d:
                            bedfile_d[chromosomeNr] = []
                            bedfile_d[chromosomeNr].append([chromosomeNr, startPosition, endPosition, motifLength, motif, expBaseNrchange,noOfSubstitution, nrOfIndels])
                        else:
                            bedfile_d[chromosomeNr].append([chromosomeNr, startPosition, endPosition, motifLength, motif, expBaseNrchange,noOfSubstitution, nrOfIndels])
            return bedfile_d

    def __isBedHeader(self, line):
        # from https://en.wikipedia.org/wiki/BED_(file_format)#Header
        # and http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1
        if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
            return True

    def __getChromosomeNumber(self, chromosome):
        chrNr = 0
        if type(chromosome) is not int:
            #print(chromosome)
            chr = [int(s) for s in re.findall(r'\d+',chromosome)]

            # for X and Y chromosome
            if chr is None or len(chr) <= 0:
                chr = [str(s) for s in re.findall(r'\w+', chromosome)]
            if chr is not None and chr[0] is not None:
                #print(chr)
                if type(chr[0]) is str and not chr[0].isnumeric():
                    if "X" in chr[0] or "x" in chr[0]:
                        chrNr = "X"
                    elif "Y" in chr[0] or "y" in chr[0]:
                        chrNr = "Y"
                else:
                    chrNr = chr[0]
                print(chrNr)
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
            for entrance in bedfile_l_copy:
                for line in entrance:  # not recording bad mathces in new bedfile
                    if len(line)>= 8:
                        [chromosomeNr, startPosition, endPosition, motifLength, motif, expBaseNrchange,noOfSubstitution, nrOfIndels] = [str(s) for s in line[0:8]]
                        # a regions start and end was found, prepare line for new bedfile
                        if startPosition != 0 and endPosition != 0:
                            printString = chromosomeNr + "\t" + startPosition + "\t" + endPosition + "\t" \
                                          + motifLength + "\t" + motif + "\t" + expBaseNrchange + "\t" \
                                          + noOfSubstitution + "\t" + nrOfIndels + "\n"
                            outBedfile.write(printString)
                        else:
                            count += 1
            print("Not found regions (unchanged): ",count)

