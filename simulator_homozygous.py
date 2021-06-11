# simulates STR changes compared to a reference genome
# hd stands for haploid and diploid, as in this version haploid or diploid genomes can be simulated

import copy
import random
import sys
from Bio import SeqIO
import os

# oldBedFile = "..\FilteredViewed\\hs37_ver8.chr22.bed "
# newBedFile = "..\FilteredViewed\\hs37_ver8.chr22.adapt.bed"
# newFastaFile = "..\FilteredViewed\\hs37d5.chr22.rand_adapt.fa"
# oldFastaFile = "..\FilteredViewed\\hs37d5.chr22.fa"


def getBedFile(oldBedFile):
    bedfile_l = list()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            splitline = line.split("\t")
            if len(splitline) > 3:
                bedfile_l.append(splitline) #chr    from Pos    to Pos      lenMotif    motif
    return bedfile_l


#write out the new coordinates. The adapted bedfile. To find the regions and to
def printBedModifications(bedfile_l_copy, newBedFile):
    with open(newBedFile, 'w') as outBedfile:
        for line in bedfile_l_copy: #not recording bad mathces in new bedfile
            if line[1] != 0 and line[2] !=0:
                lines = line[0]+"\t"+str(line[1])+"\t"+str(line[2])+"\t"+line[3]+"\t"+line[4]
                outBedfile.write(lines)

def mutate(sequence, chanceOfMutation, indelsLessMutation):
    sequence_2 = ""
    offset2= 0
    for i in sequence:
        chance2 = random.random()
        if chance2 < chanceOfMutation:
            chance = random.random()
            if i == "A":
                if chance < 0.85:
                    sequence_2 += "G"
                elif chance < 0.95:
                    sequence_2 += "T"
                else:
                    sequence_2 += "C"
            elif i == "T":
                if chance < 0.85:
                    sequence_2 += "C"
                elif chance < 0.95:
                    sequence_2 += "A"
                else:
                    sequence_2 += "G"
            elif i == "C":
                if chance < 0.85:
                    sequence_2 += "T"
                elif chance < 0.95:
                    sequence_2 += "G"
                else:
                    sequence_2 += "A"
            else: # i == "G":
                if chance < 0.85:
                    sequence_2 += "A"
                elif chance < 0.95:
                    sequence_2 += "C"
                else:
                    sequence_2 += "T"
        else:
            sequence_2+= i
    return sequence_2,offset2

#main function manipulate Fasta
def main_manipulation(newFastaFile, oldFastaFile, newBedFile, oldBedFile, chanceOfChange, nrOfChr, chanceOfMutationPerBase, indelsLessMutation, homozygousity):
    #cast inputs
    chanceOfChange = float(chanceOfChange)
    chanceOfMutationPerBase = float(chanceOfMutationPerBase)
    #safe bedfile as lists.
    bedfile_l = getBedFile(oldBedFile)
    bedfile_l_length = len(bedfile_l)
    # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
    bedfile_l_copy = copy.deepcopy(bedfile_l)
    sequence2 = ""
    with open(newFastaFile, 'w') as outFastaFile:
        writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
        with open(oldFastaFile, 'r') as inFastaFile: #read fastaFile
            for record in SeqIO.parse(inFastaFile, "fasta"): #every Record. 1 .... 2 .... 3..... .... 22 .... x...
                sequence= record.seq
                recordLen = len(sequence)               # old length
                print("old length", recordLen)

                for chr in range(0,nrOfChr): # per new to be created chromosome
                    record2 = copy.deepcopy(record)  # changes only on deep copies.
                    sequence2 = copy.deepcopy(sequence)  # changes only on deep copies.
                    # have an offset, that tells how much all following coordinates will be later or earlier
                    offset = 0  # original in the beginning
                    #naming of haploid and or diploid chromosome entrances.
                    nameOfChr = record2.name
                    idOfChr = record2.id
                    allele = 1
                    if chr == 1:
                        allele+=1
                    #else :
                    #    allele = allele
                    nameOfChr = nameOfChr + "_" + str(allele)
                    idOfChr = idOfChr + "_" + str(allele)
                    record2.name = nameOfChr
                    record2.id = idOfChr

                    for i in range(0,bedfile_l_length): #go through all coordinates in the bedfile.
                        shortTR = bedfile_l[i]  # entrance on i
                        chrnr = shortTR[0]      # which chromosome
                        chrnr_w = chrnr+"_"+str(allele)         # this number will be written down
                        patternStart = int(shortTR[1])-1 + offset
                        patternEnd = int(shortTR[2]) + offset
                        patternLen = int(shortTR[3])
                        pattern = shortTR[4].strip()
                        if record.id == chrnr: #recored id must be same as chrnr in the line of the bedfile.
                            seq_len = len(sequence2)
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

                            #general information
                            partOfSeq_4 = sequence2[patternStart:patternEnd]
                            numberOfRepeats = int(len(partOfSeq_4)/patternLen)

                            #preparation to change bedfile
                            entrance = bedfile_l[i]
                            entrance_c = copy.deepcopy(entrance)  # only change the copy
                            entrance_c[0] = chrnr_w
                            entrance_c[1] = patternStart + 1
                            entrance_c[2] = patternEnd

                            partOfSeqNew = pattern * numberOfRepeats
                            ##if STRs should be changed
                            if random.random()<=chanceOfChange:
                                #if you want to simulate a reduction you cannot reduce more than the available number of repeats.
                                #if you want to simulate an increase of repeats, do anything between
                                manipulation = random.randint(0,10) if random.random()<0.5 else (-1)*(random.randint(0,numberOfRepeats))
                                numberOfRepeatsNew = numberOfRepeats+ manipulation      #total new number of repeats
                                patternEndNew = patternStart + numberOfRepeatsNew * patternLen  # current end
                                entrance_c[2] = patternEndNew
                                offset += (patternEndNew - patternEnd)  # remember for following
                                partOfSeqNew = pattern * numberOfRepeatsNew # new middle sequence


                            #cut current sequence replace
                            sequence2 = sequence2[:patternStart] + "" + sequence2[patternEnd:] #cuts part inbetween

                            #mutate new sequence
                            partOfSeqNew, offset2 = mutate(partOfSeqNew, chanceOfMutationPerBase,indelsLessMutation)

                            #insert new sequence
                            sequence2 = sequence2[:patternStart] + partOfSeqNew + sequence2[patternStart:] #fills part inbetween

                            #changes in bedfile
                            if allele == 1: # only change the first entrance of the copy (first chromosome)
                                bedfile_l_copy[i] = entrance_c
                            else:           # append the entrances of the second chromosome
                                bedfile_l_copy.append(entrance_c)

                            #to many errors in pattern in general
                            if noFit:# no fit didnot match in 10 positions or is no STR anymore therefore is not a good coordinate for a STR
                                entrance = bedfile_l[i]
                                entrance_c = copy.deepcopy(entrance)
                                entrance[0] = chrnr_w
                                entrance[1] = 0  # mark it as 0 later don't put it in new bedfile
                                entrance[2] = 0
                                if allele == 1:
                                    bedfile_l_copy[i] = entrance_c
                                else:
                                    bedfile_l_copy[-1] = entrance_c

                    record2.seq = sequence2
                    record2.id = idOfChr
                    record2.name = nameOfChr
                    newrecordlen = len(sequence2)
                    print("new length sequence: ", newrecordlen)
                    writer.write_header()
                    writer.write_record(record2)

            printBedModifications(bedfile_l_copy,newBedFile)

#main_manipulation(newFastaFile,oldFastaFile,newBedFile,oldBedFile, 0.99, 2, 0.05)

if len(sys.argv) < 6:
    print("Please give fastafilenames+dir to the ref [1] and new Fasta that will be procues [2], the old bedfile [3], "
          "the  name and dir where the new bedfile should be [4], "
          "and a number between >0.0 && <1.0 indicating the changes of STR to the original ref file [5]")
    print("Additional: 'd' oder 'h' for diploid or haploid approach [6], chance of mutation >0.0 && <1.0 [7], how much rarer an indel is than a mutation (0-100, comma possible) [8], homozygousity >0.0 && <1.0 [9]")
else:
    for i in sys.argv:
        print(i)
    if len(sys.argv) > 6:
        print("more than 6 parameters")
        if sys.argv[6].isalpha():
            print("sys.argv[6].isalpha()")
            run = 1
            if sys.argv[6] == "d":
                print("sys.arg[6]: you assigned 'd' therefore run is diploid.")
                run = 2
            else:
                print("sys.arg[6]: run is haploid. 'h' or sth else was given.")
            if len(sys.argv) >7:
                print("more than 7 parameters")
                canBeCast = True
                try:
                    float(sys.argv[7])
                except ValueError:
                    print("Not a float")
                    canBeCast = False
                if canBeCast:
                    print("sys.argv[7].isfloat")
                    chanceOfMutation= float(sys.argv[7])
                    if len(sys.argv) > 8:
                        print("more than 8 parameters")
                        canBeCast = True
                        try:
                            float(sys.argv[8])
                        except ValueError:
                            print("Not a float")
                            canBeCast = False
                        if canBeCast: #can [8] be cast
                            print("sys.argv[8].isfloat")
                            indelLessTmutation = float(sys.argv[8])
                            if len(sys.argv) > 9: #len > 9
                                print("more than 9 parameters")
                                canBeCast = True
                                try:
                                    float(sys.argv[9])
                                except ValueError:
                                    print("Not a float")
                                    canBeCast = False
                                if canBeCast:
                                    print("sys.argv[8].isfloat")
                                    homozygous = float(sys.argv[8])
                                    if homozygous >= 0.0 and homozygous <= 1.0:
                                        print("in field 9 homozygous possibilty assigned was: ", homozygous)
                                        print("How much rarer is an Indel than Mutation, according to assignment: ",
                                              indelLessTmutation)
                                        if chanceOfMutation >= 0.0 and chanceOfMutation <= 1.0:
                                            print("your assigned chance of mutation rate ", chanceOfMutation)
                                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3],
                                                              sys.argv[5],
                                                              run,
                                                              chanceOfMutation, indelLessTmutation, homozygous)
                                        else:
                                            print(
                                                "Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3],
                                                              sys.argv[5],
                                                              run,
                                                              0.01, indelLessTmutation, homozygous)

                                    else:
                                        print(
                                            "wrong format. homozygous possibilty needs a value > 0.0 and smaller <=1.0. default is 0")
                                        print("How much rarer is an Indel than Mutation, according to assignment: ",
                                              indelLessTmutation)
                                        if chanceOfMutation >= 0.0 and chanceOfMutation <= 1.0:
                                            print("your assigned chance of mutation rate ", chanceOfMutation)
                                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3],
                                                              sys.argv[5],
                                                              run,
                                                              chanceOfMutation, indelLessTmutation, 0)
                                        else:
                                            print(
                                                "Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3],
                                                              sys.argv[5],
                                                              run,
                                                              0.01, indelLessTmutation, 0)
                                else: #9 was not cast.
                                    print("wrong format. homozygous possibilty needs a value > 0.0 and smaller <=1.0. default is 0")
                                    print("How much rarer is an Indel than Mutation, according to assignment: ",
                                          indelLessTmutation)
                                    if chanceOfMutation >= 0.0 and chanceOfMutation <= 1.0:
                                        print("your assigned chance of mutation rate ", chanceOfMutation)
                                        main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3],
                                                          sys.argv[5],
                                                          run,
                                                          chanceOfMutation, indelLessTmutation, 0)
                                    else:
                                        print(
                                            "Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                                        main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3],
                                                          sys.argv[5],
                                                          run,
                                                          0.01, indelLessTmutation, 0)
                            else: #len not 0 9
                                print(
                                    "in field 9 homozygous possibilty was missing. here the default is 0")
                                print("How much rarer is an Indel than Mutation, according to assignment: ",
                                      indelLessTmutation)
                                if chanceOfMutation >= 0.0 and chanceOfMutation <= 1.0:
                                    print("your assigned chance of mutation rate ", chanceOfMutation)
                                    main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5],
                                                      run,
                                                      chanceOfMutation, indelLessTmutation, 0)
                                else:
                                    print(
                                        "Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                                    main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5],
                                                      run,
                                                      0.01, indelLessTmutation, 0)
                        else:#[8] can not be cast. therefore default. And 7 fields ok and can be run. but check of range nessesary
                                print("wrong format. in field 8 a number is needed to express how much rarer an indel"
                                      "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                                if chanceOfMutation >= 0.0 and chanceOfMutation <= 1.0:
                                    print("your assigned chance of mutation rate ", chanceOfMutation)
                                    print("Default is 10 times rarer")
                                    print(
                                        "in field 9 homozygous possibilty was missing. if 'd' was chosen, here the default is 0")
                                    main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5],
                                                      run, chanceOfMutation, 10, 0)
                                else:
                                    print(
                                        "Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                                    print("Default is 10 times rarer")
                                    print(
                                        "in field 9 homozygous possibilty was missing. if 'd' was chosen, here the default is 0")
                                    main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5],
                                                      run, 0.01, 10, 0)
                    else: # 7 fields. but and chance of mutation was casted float. But is it in range?
                        if chanceOfMutation >= 0.0 and chanceOfMutation <= 1.0:
                            print("your assigned chance of mutation rate ", chanceOfMutation)
                            print("in field 8 a number is needed to express how much rarer an indel"
                                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                            print("Default is 10 times rarer")
                            print("in field 9 homozygous possibilty was missing. if 'd' was chosen, here the default is 0")
                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run, chanceOfMutation, 10, 0)
                        else :
                            print("Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                            print("in field 8 a number is needed to express how much rarer an indel"
                                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                            print("Default is 10 times rarer")
                            print(
                                "in field 9 homozygous possibilty was missing. if 'd' was chosen, here the default is 0")
                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run, 0.01, 10, 0)
                else: # 7 cannot be cast
                    print("sys.argv[7] is not float, therefore default mutation rate is at 0.01 == 1%")
                    print("in field [8] a number is needed to express how much rarer an indel"
                          "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                    print("Default is 10 times rarer")
                    print("in field [9] homozygous possibilty was missing. here the default is 0")
                    main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run, 0.01, 10, 0)
            else: # 7 fielts
                print("len < 7: 0.01 (1%) chance of mutation, as you did not give a parameter for chance of mutation[7]")
                print("in field [8] a number is needed to express how much rarer an indel"
                      "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                print("Default is 10 times rarer")
                print("in field [9] homozygous possibilty was missing. here the default is 0")
                main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run, 0.01,10, 0)
        else:
            print("sys.argv[6] is not the right format (string: 'd' or 'h'), therefore run is default haploid and with 0.01 == 1% mutationrate [7]")
            print("in field [8] a number is needed to express how much rarer an indel"
                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
            print("Default is 10 times rarer")
            print("in field [9] homozygous possibilty was missing. here the default is 0")
            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], 1, 0.01, 10, 0)
    else: #default
        if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[3]):
            print("given 5 parameters, therefore run is default haploid [6] and with 0.01 == 1% mutationrate [7]")
            print("in field [8] a number is needed to express how much rarer an indel"
                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
            print("Default is 10 times rarer")
            print("in field [9] homozygous possibilty was missing. here the default is 0")
            main_manipulation(sys.argv[2],sys.argv[1],sys.argv[4],sys.argv[3],sys.argv[5], 1, 0.01, 10, 0)
            #[0]./simulator2.py [1]../../reference/hs37d5.chr22.fa  [2]hs37d5.chr22.new1.fa [3]../bedfiles/hs37_ver8.chr22.bed [4]hs37_ver8.chr22.new1.bed [5]0.20
        else:
            if not os.path.isfile(sys.argv[1]):
                print(sys.argv[1], " is not a file")
            elif not os.path.isfile(sys.argv[3]):
                print(sys.argv[3], " is not a file")




# fai format:
#chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
#22	    51304566	    2876892038	            60	    61
#when reading fasta

#bedfile format
#chr    from Pos    to Pos      lenMotif    motif
#22	    16052168	16052199	4	        AAAC    ==> 68:AAACAAACAAACAAACAAACAAACAAACAAAC:99