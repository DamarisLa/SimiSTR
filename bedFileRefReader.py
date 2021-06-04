import copy
from Bio import SeqIO
import random
import sys
import os

newBedFile = "..\FilteredViewed\\hs37_ver8.chr22.adapt.bed"
newFastaFile = "..\FilteredViewed\\hs37d5.chr22.rand_adapt.fa"
oldBedFile = "..\FilteredViewed\\hs37_ver8.chr22.bed "


def getBedFile(oldBedFile):
    bedfile_l = list()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            splitline = line.split("\t")
            if len(splitline)>3:
                bedfile_l.append(splitline) #chr    from Pos    to Pos      lenMotif    motif
    return bedfile_l


def main_manipulation(newFastaFile, newBedFile, oldBedfile):

    bedfile_l = getBedFile(newBedFile)
    bedfile_l_old =getBedFile(oldBedfile)
    # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates

    with open(newFastaFile, 'r') as inFastaFile:
        for record in SeqIO.parse(inFastaFile, "fasta"):
            # for i in range(0,len(record.seq)):
            sequence= record.seq
            recordLen = len(sequence)
            print(recordLen)
            # writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
            for i in range(0,len(bedfile_l)):
                shortTR = bedfile_l[i]
                shortTRold = bedfile_l_old[i]
                chrnr = shortTR[0]
                patternStart = int(shortTR[1])-1
                patternEnd = int(shortTR[2])
                patternLen = int(shortTR[3])
                pattern = shortTR[4].strip()
                patternStartold = int(shortTRold[1])-1
                patternEndold = int(shortTRold[2])
                patternLenold = int(shortTRold[3])
                patternold = shortTRold[4].strip()
                if record.id == chrnr:
                    seq_len = len(sequence)
                    partOfSeq = sequence[patternStart:patternEnd]

                    print(partOfSeq)
                    print("new: ",shortTR," ", (patternEnd-patternStart)/patternLen)
                    print("old: ",shortTRold, " ", (patternEndold-patternStartold)/patternLenold)


main_manipulation(newFastaFile,newBedFile,oldBedFile)