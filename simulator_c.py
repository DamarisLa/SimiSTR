import copy
from Bio import SeqIO
import random
import os




#read the original bedfile. safe all coordinates maybe into a dict
# ?? whats the key ?? the coordinate? or should it be a list? => as iteration from first to last. List
# makes more sense than dict
class simulator:
    def __init__(self,oldFastaFile,newFastaFile,oldBedFile,newBedFile,chanceOfChange):
        if os.path.isfile(oldFastaFile) and os.path.isfile(oldBedFile):
            print("running")
            self.__main_manipulation(newFastaFile, oldFastaFile, newBedFile, oldBedFile, chanceOfChange)
        else :
            if not os.path.isfile(oldBedFile):
                print(oldBedFile, " is not the wright filedir.")
            elif not os.path.isfile(oldFastaFile):
                print(oldFastaFile, " is not the wright filedir")


    def __getBedFile(self, oldBedFile):
        bedfile_l = list()
        with open(oldBedFile, 'r') as inBedFile:
            for line in inBedFile:
                splitline = line.split("\t")
                if len(splitline)>3:
                    bedfile_l.append(splitline) #chr    from Pos    to Pos      lenMotif    motif
        return bedfile_l


    #write out the new coordinates. The adapted bedfile. To find the regions and to
    def __printBedModifications(self, bedfile_l_copy, newBedFile):
        with open(newBedFile, 'w') as outBedfile:
            for line in bedfile_l_copy: #not recording bad mathces in new bedfile
                if line[1] != 0 and line[2] !=0:
                    lines = line[0]+"\t"+str(line[1])+"\t"+str(line[2])+"\t"+line[3]+"\t"+line[4]
                    outBedfile.write(lines)

    #bedfile_l = getBedFile(oldBedFile)
    #printBedModifications(bedfile_l,newBedFile)


    #main function manipulate Fasta
    def __main_manipulation(self, newFastaFile, oldFastaFile, newBedFile, oldBedFile, chanceOfChange):
        bedfile_l = self.__getBedFile(oldBedFile)
        # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
        bedfile_l_copy = copy.deepcopy(bedfile_l)
        # have an offset, that tells how much all following coordinates will be later or earlier
        offset = 0  # original in the beginning
        sequence2 = ""
        with open(newFastaFile, 'w') as outFastaFile:
            writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
            with open(oldFastaFile, 'r') as inFastaFile:
                for record in SeqIO.parse(inFastaFile, "fasta"):
                    record2 = copy.deepcopy(record)
                    # for i in range(0,len(record.seq)):
                    sequence= record.seq
                    sequence2 = copy.deepcopy(sequence)

                    # writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
                    for i in range(0,len(bedfile_l)):
                        shortTR = bedfile_l[i]
                        chrnr = shortTR[0]
                        patternStart = int(shortTR[1])-1 + offset
                        patternEnd = int(shortTR[2]) + offset
                        patternLen = int(shortTR[3])
                        pattern = shortTR[4].strip()

                        if record.id == chrnr:
                            seq_len = len(sequence2)
                            partOfSeq = sequence2[patternStart:patternEnd]
                            print(pattern, "; ",partOfSeq)

                            #mabye check if startposition really matches...or adjust it to the left or right
                            partOfSeq_1 = sequence2[patternStart:patternStart+patternLen]
                            startpoint = patternStart                           #if start is 1000
                            startPointCorrect = False
                            lower = False
                            number = -1                 #-1
                            noFit = False
                            while not startPointCorrect and not noFit:
                                if partOfSeq_1 == pattern:
                                    startPointCorrect = True #found startpoint
                                    patternStart = startpoint
                                    # find out if area in reference of this STR is acutally longer or shorter than according to
                                    # the bed say? basically count if pattern goes left an right for longer.
                                    # maybe later change to another comparison algorithm that takes costs into account (Blast)
                                    ### Leftside check ###
                                    check_left = patternStart - patternLen
                                    goLeft = True
                                    while check_left >= 0 and goLeft:
                                        partOfSeq_2 = sequence2[check_left:check_left + patternLen]
                                        if partOfSeq_2 == pattern:  # still matching
                                            check_left = check_left - patternLen  # maybe your can even go further left
                                        else:
                                            goLeft = False  # stop checking further right
                                            check_left = check_left + patternLen  # because last length did still match

                                    ### Rightside check ###
                                    check_right = check_left
                                    goRight = True
                                    while check_right + patternLen <= seq_len and goRight:
                                        partOfSeq_3 = sequence2[check_right:check_right + patternLen]
                                        if partOfSeq_3 == pattern:  # still matching
                                            check_right = check_right + patternLen  # maybe your can even go further left
                                        else:
                                            goRight = False  # stop checking further right

                                    patternStart = check_left
                                    patternEnd = check_right

                                    if patternEnd - patternStart == patternLen:
                                        startPointCorrect = False
                                        startpoint += number
                                        partOfSeq_1 = sequence2[startpoint:startpoint + patternLen]
                                        if lower:  # true
                                            number += 1  # 3
                                            number *= (-1)  # -3
                                            lower = False  # false
                                        else:  # false
                                            number *= (-1)  # 1 .... 3
                                            number += 1  # 2 .... 4
                                            lower = True  # true
                                        if abs(number) > 20:
                                            noFit = True


                                else: #will always go one further away from current start. trying both directions parallel
                                    startpoint = startpoint + number            # start is 999, 1001, 998, 1002 usw.
                                    partOfSeq_1 = sequence2[startpoint:startpoint + patternLen]
                                    if lower:   #true
                                        number += 1     #3
                                        number *= (-1)  #-3
                                        lower = False   #false
                                    else:       #false
                                        number *= (-1)  #1 .... 3
                                        number += 1     #2 .... 4
                                        lower = True    #true
                                    if abs(number) > 20:
                                        noFit = True

                            partOfSeq_4 = sequence2[patternStart:patternEnd]
                            numberOfRepeats = int(len(partOfSeq_4)/patternLen)

                            entrance = bedfile_l_copy[i]
                            entrance[1] = patternStart + 1
                            entrance[2] = patternEnd
                            bedfile_l_copy[i] = entrance
                            #do you want to increase or decrease?
                            if random.random()<=chanceOfChange:
                                #if you want to simulate a reduction you cannot reduce more than the available number of repeats.
                                #if you want to simulate an increase of repeats, do anything between
                                #manipulation = random.randint(0,10) if random.random()<0.5 else (-1)*(random.randint(0,numberOfRepeats))
                                manipulation = 10
                                # just to see in debugging that string replacement works
                                debugHelpPartOfSeq = sequence2[patternStart - 10:patternEnd + 10]
                                #print(debugHelpPartOfSeq)
                                numberOfRepeatsNew = numberOfRepeats+ manipulation      #total new number of repeats
                                patternEndNew = patternStart + numberOfRepeatsNew * patternLen  # current end
                                partOfSeqNew = pattern * numberOfRepeatsNew

                                offset += (patternEndNew - patternEnd) #remember for following
                                #replace
                                sequence2 = sequence2[:patternStart] + "" + sequence2[patternEnd:]
                                debugHelpPartOfSeq = sequence2[patternStart - 10:patternEnd + 10 ]
                                #print(debugHelpPartOfSeq)
                                sequence2 = sequence2[:patternStart] + partOfSeqNew + sequence2[patternStart:]
                                debugHelpPartOfSeq = sequence2[patternStart - 10:patternEnd + 10 + manipulation*patternLen]
                                #print(debugHelpPartOfSeq)
                                # just to see in debugging that string replacement works
                                debugHelpPartOfSeq2 = sequence2[patternStart - 10:patternEnd + 10 + -1*manipulation*patternLen]
                                #print(debugHelpPartOfSeq2)
                            if noFit:# no fit didnot match in 10 positions or is no STR anymore therefore is not a good coordinate for a STR
                                entrance = bedfile_l_copy[i]
                                entrance[1] = 0 #mark it as 0 later don't put it in new bedfile
                                entrance[2] = 0
                                bedfile_l_copy[i] = entrance


                    record2.seq = sequence2
                    writer.write_header()
                    writer.write_record(record2)
                self.__printBedModifications(bedfile_l_copy,newBedFile)










# fai format:
#chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
#22	    51304566	    2876892038	            60	    61
#when reading fasta

#bedfile format
#chr    from Pos    to Pos      lenMotif    motif
#22	    16052168	16052199	4	        AAAC    ==> 68:AAACAAACAAACAAACAAACAAACAAACAAAC:99