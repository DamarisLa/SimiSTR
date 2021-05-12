import copy
from Bio import SeqIO
import random

oldBedFile = "..\FilteredViewed\\hs37_ver8.chr22.bed"
newBedFile = "..\FilteredViewed\\hs37_ver8.chr22.adapt.bed"
newFastaFile = "..\FilteredViewed\\hs37d5.chr22.rand_adapt.fa"    #ref-file before STR simulations
oldFastaFile = "..\FilteredViewed\\hs37d5.chr22.fa"         #ref-file after STR simulations



#read the original bedfile. safe all coordinates maybe into a dict
# ?? whats the key ?? the coordinate? or should it be a list? => as iteration from first to last. List
# makes more sense than dict


def getBedFile(oldBedFile):
    bedfile_l = list()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            splitline = line.split("\t")
            if len(splitline)>3:
                bedfile_l.append(splitline) #chr    from Pos    to Pos      lenMotif    motif
    return bedfile_l


bedfile_l = getBedFile(oldBedFile)





#main function manipulate Fasta
def main_manipulation(bedfile_l, newFastaFile, oldFastaFile):
    # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
    bedfile_l_copy = copy.deepcopy(bedfile_l)
    # have an offset, that tells how much all following coordinates will be later or earlier
    offset = 0  # original in the beginning
    randomModification = 0
    with open(newFastaFile, 'w') as outFastaFile:
        with open(oldFastaFile, 'r') as inFastaFile:
            for record in SeqIO.parse(oldFastaFile, "fasta"):
                # writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
                for i in range(0,len(bedfile_l)):
                    shortTR = bedfile_l[i]
                    chrnr = shortTR[0]
                    patternStart = int(shortTR[1])-1
                    patternEnd = int(shortTR[2])-1
                    patternLen = int(shortTR[3])
                    pattern = shortTR[4].strip()
                    if record.id == chrnr:
                        #for i in range(0,len(record.seq)):
                        sequence = record.seq
                        seq_len = len(sequence)
                        partOfSeq = sequence[patternStart:patternEnd]
                        #print(pattern, "; ",partOfSeq)

                        #mabye check if startposition really matches...or adjust it to the left or right
                        partOfSeq_1 = sequence[patternStart:patternStart+patternLen]
                        startpoint = patternStart                           #if start is 1000
                        startPointCorrect = False
                        lower = False
                        number = -1                 #-1
                        while not startPointCorrect:
                            if partOfSeq_1 == pattern:
                                startPointCorrect = True #found startpoint
                            else: #will always go one further away from current start. trying both directions parallel
                                startpoint = startpoint + number            # start is 999, 1001, 998, 1002 usw.
                                partOfSeq_1 = sequence[startpoint:startpoint + patternLen]
                                if lower:   #true
                                    number += 1     #3
                                    number *= (-1)  #-3
                                    lower = False   #false
                                else:       #false
                                    number *= (-1)  #1 .... 3
                                    number += 1     #2 .... 4
                                    lower = True    #true


                        # find out if area in reference of this STR is acutally longer or shorter than according to
                        # the bed say? basically count if pattern goes left an right for longer.
                        # maybe later change to another comparison algorithm that takes costs into account (Blast)
                        ### Leftside check ###
                        check_left = patternStart - patternLen
                        goLeft = True
                        while check_left >= 0 and goLeft:
                            partOfSeq_2 = sequence[check_left:check_left+patternLen]
                            if partOfSeq_2 == pattern:                  #still matching
                                check_left = check_left - patternLen    #maybe your can even go further left
                            else:
                                goLeft = False                          #stop checking further right
                                check_left = check_left + patternLen    #because last length did still match
                        ### Rightside check ###
                        check_right = patternEnd + patternLen
                        goRight = True
                        while check_right+patternLen <= seq_len and goRight:
                            partOfSeq_3 = sequence[check_right:check_right+patternLen]
                            if partOfSeq_3 == pattern:                  #still matching
                                check_right = check_right + patternLen  #maybe your can even go further left
                            else:
                                goRight = False                         #stop checking further right

                        partOfSeq_4 = sequence[check_left:check_right]
                        numberOfRepeats = len(partOfSeq_4)/patternLen
                        #if you want to simulate a reduction you cannot reduce more than the available number of repeats.

                        #if you want to simulate an increase of repeats, do anything between

                        ### here starts the fun ###


main_manipulation(bedfile_l,newFastaFile,oldFastaFile)


#write out the new coordinates. The adapted bedfile. To find the regions and to
def printBedModifications(bedfile_l_copy):
    with open(newBedFile, 'w') as outBedfile:
        for line in bedfile_l_copy():
            outBedfile.write(line)




# fai format:
#chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
#22	    51304566	    2876892038	            60	    61
#when reading fasta

#bedfile format
#chr    from Pos    to Pos      lenMotif    motif
#22	    16052168	16052199	4	        AAAC    ==> 68:AAACAAACAAACAAACAAACAAACAAACAAAC:99