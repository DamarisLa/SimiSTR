


oldBedFile = "hs37_ver8.20short.bed"
newBedFile = "hs37_ver8.20short.adapt.bed"
newFastaFile = "InternshipFiles\FilteredViewed\\hs37d5.rand_adapt.fa"
oldFastaFile = "InternshipFiles\FilteredViewed\\hs37d5.fa"
#refFile =""

with open(newBedFile, 'w') as outBedfile:
    with open(oldFastaFile, 'w') as outBamfile:
        #with open(refFile, 'r') as inRefFile:
        with open(oldBedFile, 'r') as inBedFile:
            with open(newFastaFile, 'r') as inBamFile:
                for line in inBamFile:
                    print(line)




# read the regions from bedfile.
# region 1 :
    # navigate to
    # make an random


# GRCh37:1
# from Bio import SeqIO
#
# inFile = open('path to your fasta file','r')
#
# for record in SeqIO.parse(inFile,'fasta'):
#     if record.id == "sequence id you want to extract from":
#         print str(record.seq)[startCoordinates:endCoordinates + 1]