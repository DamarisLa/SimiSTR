


# This reader collects all bed file entrances in a dictionary for easier access
def getBedFile(oldBedFile):
    bedfile_d = dict()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            splitline = line.split("\t")
            if len(splitline) > 3:
                splitline[4] = splitline[4].strip().upper()
                id = re.search("(\d*)",splitline[0]) #this line should enable to find the 1 in the "chr1" annotation
                if id[0] is not None:
                    if id not in bedfile_d.keys():
                        bedfile_d[str(id[0])] = list()
                        bedfile_d[str(id[0])].append(splitline)
                    else:
                        bedfile_d[str(id[0])].append(splitline)  # chr    from Pos    to Pos      lenMotif    motif
                else:
                    print("Bedfile reader does not recognize the chromosome id")
    return bedfile_d

def getBedFile(inputBedFile):
    bedfile_d = dict()
    with open(inputBedFile, 'r') as inBedFile:
        for bedline in inBedFile:
            """ 
            Here you can first filter out the header information that can contain
            '#', 'track' or 'browser' (http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1)
            """
            if not isBedHeader(bedline):
                bedElements = bedline.split("\t")
                # First three are mandatory, fourth is the sequence (name), we need at least 4
                if len(bedElements) >= 4:
                    [chromosome, startPosition, endPosition, motifLength, motif] = bedElements[:4]
                    bedElements[4] = bedElements[4].strip() #this is needed as motif should be a string without "\n" or else
                    chromosomeNr = getChromosomeNumber(chromosome)
                    # sequence in CAPs
                    sequence = sequence.upper()
                    if chromosomeNr not in bedfile_d:
                        bedfile_d[chromosomeNr] = [bedElements]
                    else:
                        bedfile_d[chromosomeNr].append(bedElements)  # chr    from Pos    to Pos      lenMotif    motif
        return bedfile_d

def isBedHeader(line):
    # from https://en.wikipedia.org/wiki/BED_(file_format)#Header
    # and http://genome.cse.ucsc.edu/FAQ/FAQformat.html#format1
    if line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return True

def getChromosomeNumber(chromosome):
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
