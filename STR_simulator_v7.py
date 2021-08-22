# simulates STR changes compared to a reference genome
# hd stands for haploid and diploid, as in this version haploid or diploid genomes can be simulated

import copy
import os
import random
import sys
from Bio import SeqIO

def getBedFile(oldBedFile):
    #bedfile_l = list()
    bedfile_d = dict()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            splitline = line.split("\t")
            if len(splitline) > 3:
                splitline[4] = splitline[4].strip()
                if splitline[0] not in bedfile_d.keys():
                    bedfile_d[str(splitline[0])] = list()
                    bedfile_d[str(splitline[0])].append(splitline)
                else:
                    bedfile_d[str(splitline[0])].append(splitline)  # chr    from Pos    to Pos      lenMotif    motif
    return bedfile_d

# write out the new coordinates. The adapted bedfile. To find the regions and to
def printBedModifications(bedfile_l_copy, newBedFile):
    with open(newBedFile, 'w') as outBedfile:
        count = 0
        for chr in bedfile_l_copy:
            for line in chr:  # not recording bad mathces in new bedfile
                if line[1] != 0 and line[2] != 0:  # thats the ones where pattern was not found
                    lines = str(line[0]) + "\t" + str(line[1]) + "\t" + str(line[2]) + "\t" + str(line[3]) + "\t" + str(line[4] + "\n")
                    outBedfile.write(lines)
                else:
                    count += 1
        print(count)


def mutate(sequence, chanceOfMutation, indelsLessMutation):
    #sequence2 = sequence2[:patternStart] + partOfSeqNew + sequence2[patternEnd:]
    offset2 = 0
    chanceForIndels = 0
    if indelsLessMutation != 0:  # to avoid division through zero
        chanceForIndels = chanceOfMutation / indelsLessMutation
    for base in range(0,len(sequence)):
        chance2 = random.random()
        chance3 = random.random()
        if chance2 < chanceOfMutation:
            chance = random.random()
            if base+1 < len(sequence):
                if sequence[base] == "A":
                    if chance < 0.85:
                        #sequence_2 += "G"                                  #replace
                        sequence = sequence[:base] + 'G' + sequence[base+1:] #[inclusive:exclusive] => [:base]=>exclusive+ base + [base+1:]=>
                    elif chance < 0.95:
                        #sequence_2 += "T"
                        sequence = sequence[:base] + 'T' + sequence[base+1:]
                    else:
                        #sequence_2 += "C"
                        sequence = sequence[:base] + 'C' + sequence[base+1:]
                elif sequence[base]  == "T":
                    if chance < 0.85:
                        #sequence_2 += "C"
                        sequence = sequence[:base] + 'C' + sequence[base+1:]
                    elif chance < 0.95:
                        #sequence_2 += "A"
                        sequence = sequence[:base] + 'A' + sequence[base+1:]
                    else:
                        #sequence_2 += "G"
                        sequence = sequence[:base] + 'G' + sequence[base+1:]
                elif sequence[base]  == "C":
                    if chance < 0.85:
                        #sequence_2 += "T"
                        sequence = sequence[:base] + 'T' + sequence[base+1:]
                    elif chance < 0.95:
                        #sequence_2 += "G"
                        sequence = sequence[:base] + 'G' + sequence[base+1:]
                    else:
                        #sequence_2 += "A"
                        sequence = sequence[:base] + 'A' + sequence[base+1:]
                else:  # i == "G":
                    if chance < 0.85:
                        #sequence_2 += "A"
                        sequence = sequence[:base] + 'A' + sequence[base+1:]
                    elif chance < 0.95:
                        #sequence_2 += "C"
                        sequence = sequence[:base] + 'C' + sequence[base+1:]
                    else:
                        #sequence_2 += "T"
                        sequence = sequence[:base] + 'T' + sequence[base+1:]
        elif chance3 < chanceForIndels:  # zb 1% / 10  = 0.1%
            if base+1 < len(sequence):
                if random.random() < 0.5:  # insertion
                    #sequence_2 += base
                    if random.random() < 0.25:
                        #sequence_2 += "A"
                        sequence = sequence[:base] + 'G' + sequence[base:] #insertion
                    elif random.random() < 0.5:
                        #sequence_2 += "C"
                        sequence = sequence[:base] + 'C' + sequence[base:]
                    elif random.random() < 0.75:
                        #sequence_2 += "G"
                        sequence = sequence[:base] + 'G' + sequence[base:]
                    else:
                        #sequence_2 += "T"
                        sequence = sequence[:base] + 'T' + sequence[base:]
                    offset2 += 1
                else:  # deletion
                    #sequence_2 += ""
                    sequence = sequence[:base] + sequence[base+1:]
                    offset2 -= 1
        #else:
        #    sequence_2 += base
    return sequence, offset2


def findStartPoint(seq, start, pattern, patternLen):
    seq_len = len(seq)
    partOfSeq_1 = seq[start:start + patternLen]
    startpoint = start  # if start is 1000
    endpoint = start + patternLen
    startpoint_mem = startpoint
    endpoint_mem = endpoint
    startPointCorrect = False
    lower = False
    number = -1  # -1
    noFit = False
    times = False
    while not startPointCorrect and not noFit:

        if pattern != partOfSeq_1:  # pattern was not found on initial start point
            # will always go one further away from current start. trying both directions parallel
            startpoint = startpoint + number  # start is 999, 1001, 998, 1002 usw.
            partOfSeq_1 = seq[startpoint:startpoint + patternLen]
            if partOfSeq_1 == pattern:  # found
                number = -1  # reset and leave
            else:
                if lower:  # true
                    number += 1  # 3
                    number *= (-1)  # -3
                    lower = False  # false
                else:  # false
                    number *= (-1)  # 1 .... 3
                    number += 1  # 2 .... 4
                    lower = True  # true
                if abs(number) > 30:  # no possible patternstart or STR within +/- 10 positions distance of bedfile-patternstart
                    noFit = True
                    number = -1
                    startpoint = startpoint_mem
                    endpoint = endpoint_mem

        else:  # if a pattern is found
            startPointCorrect = True  # found startpoint
            patternStart = startpoint
            # find out if area in reference of this STR is acutally longer or shorter than according to
            # the bed say? basically count if pattern goes left an right for longer.
            # maybe later change to another comparison algorithm that takes costs into account (Blast)

            # __ Leftside check __ #
            check_left = patternStart - patternLen
            goLeft = True
            while check_left >= 0 and goLeft:
                partOfSeq_2 = seq[check_left:check_left + patternLen]
                if partOfSeq_2 == pattern:  # still matching
                    check_left = check_left - patternLen  # maybe your can even go further left
                else:
                    goLeft = False  # stop checking further left
                    check_left = check_left + patternLen  # because last position did not match, so go back to last
                    # one working

            # __ Rightside check __ #
            check_right = check_left
            goRight = True
            while check_right + patternLen <= seq_len and goRight:
                partOfSeq_3 = seq[check_right:check_right + patternLen]
                if partOfSeq_3 == pattern:  # still matching
                    check_right = check_right + patternLen  # maybe your can even go further left
                else:
                    goRight = False  # stop checking further right

            startpoint = check_left
            endpoint = check_right

            if endpoint - startpoint == patternLen:  # length only one time the pattern
                endpoint_mem = endpoint
                startpoint_mem = startpoint
                if not times:
                    if number < 0:
                        number += 9
                    else:
                        number -= 9
                    startPointCorrect = False
                    startpoint += number
                    partOfSeq_1 = seq[startpoint:startpoint + patternLen]
                    times = True
                else:
                    startPointCorrect = True

    return startpoint, endpoint, noFit


# main function manipulate Fasta
def main_manipulation(newFastaFile, oldFastaFile, newBedFile, oldBedFile, chanceOfChange, nrOfChr,
                      chanceOfMutationPerBase, indelsLessMutation, homozygousity_rate):
    # cast inputs
    chanceOfChange = float(chanceOfChange)
    chanceOfMutationPerBase = float(chanceOfMutationPerBase)
    # safe bedfile as lists.
    bedfile_d = getBedFile(oldBedFile)
    # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
    bedfile_total = list()

    with open(newFastaFile, 'w') as outFastaFile:
        writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
        with open(oldFastaFile, 'r') as inFastaFile:  # read fastaFile
            for record in SeqIO.parse(inFastaFile, "fasta"):  # every Record. 1 .... 2 .... 3..... .... 22 .... x...
                sequence = record.seq
                recordLen = len(sequence)  # old length
                print("old length", recordLen)
                homozygousity_d = dict()
                allele = 1
                for chr in range(0,
                                 nrOfChr):  # per chromosome a new to be created chromosome #eigther only 0 oder 0 and 1
                    record2 = copy.deepcopy(record)  # changes only on deep copies.
                    sequence2 = copy.deepcopy(sequence)  # changes only on deep copies.
                    # have an offset, that tells how much all following coordinates will be later or earlier
                    offset = 0  # original in the beginning
                    # naming of haploid and or diploid chromosome entrances.
                    nameOfChr = record2.name
                    idOfChr = record2.id

                    nameOfChr = nameOfChr + "_" + str(allele)
                    idOfChr = idOfChr + "_" + str(allele)
                    record2.name = nameOfChr
                    record2.id = idOfChr

                    id = str(record.id)
                    if id in bedfile_d.keys():
                        bedfile_l = bedfile_d[id]
                        bedfile_l_copy = copy.deepcopy(bedfile_l)
                        bedfile_l_length = len(bedfile_l)

                        for bedfEntrance in range(0,
                                                  bedfile_l_length):  # go through all coordinates in the bedfile.
                            shortTR = bedfile_l[bedfEntrance]
                            chrnr = shortTR[0]  # which chromosome
                            chrnr_w = chrnr + "_" + str(allele)  # this number will be written down
                            patternStart = offset + int(shortTR[1])  # - 1 #-1 for gangstr bedfiles
                            patternEnd = offset + int(shortTR[2])
                            patternLen = int(shortTR[3])
                            pattern = shortTR[4].strip()

                            if allele > 1 and random.random() < homozygousity_rate:  # should be the second allele and inside the homozygosity rate
                                # if below chance then copy allele 1
                                # homozy = True
                                chrnr2 = chrnr
                                ps = shortTR[1]
                                partOfSeq = homozygousity_d[chrnr2, ps]

                                # cut current sequence replace
                                sequence2 = sequence2[:patternStart] + "" + sequence2[
                                                                            patternEnd:]  # cuts part inbetween
                                # insert new sequence
                                sequence2 = sequence2[:patternStart] + partOfSeq + sequence2[
                                                                                   patternStart:]  # fills part inbetween
                                entrance_c = bedfile_l_copy[bedfEntrance]
                                entrance_cn = copy.deepcopy(entrance_c)
                                entrance_cn[0] = chrnr_w
                                entrance_cn[1] = patternStart  # + 1
                                entrance_cn[2] = patternStart + len(partOfSeq)
                                oldSeqLen = (patternEnd - patternStart)
                                change3 = len(partOfSeq) - oldSeqLen
                                # print("Change: ",change3)
                                offset += change3
                                # print("offset: ", offset)
                                bedfile_l_copy[bedfEntrance] = entrance_cn
                            else:  #
                                chrnr2 = chrnr
                                ps = shortTR[1]
                                homozygousity_d[chrnr2, ps] = ""
                                # find correct startpoit or if bedfile-entrance does not fit to sequence on that position
                                correctStart, correctEnd, noFit = findStartPoint(sequence2, patternStart, pattern,
                                                                                 patternLen)  # seq,start,pattern,patternLen
                                # if noFit: patternStart and patternEnd stay as they are, but will not be noted in new
                                #           bedfile. And offset most not be changed.

                                # debugHelpPartOfSeq = sequence2[patternStart - 30:patternEnd + 30]
                                # debugHelpPartOfSeq1 = sequence2[patternStart:patternEnd]
                                # debugHelpPartOfSeq2 = sequence2[correctStart - 30:correctEnd + 30]
                                # debugHelpPartOfSeq3 = sequence2[correctStart:correctEnd]
                                # print(helppattern,"\n",debugHelpPartOfSeq,"\n",debugHelpPartOfSeq1,"\n",debugHelpPartOfSeq2,"\n",debugHelpPartOfSeq3,"\n")

                                # preparation to change bedfile
                                # entrance = bedfile_l[bedfEntrance]
                                # entrance_c = copy.deepcopy(entrance)  # only change the copy

                                # partOfSeq_4 = sequence2[patternStart:patternEnd]  # just for debugging
                                if not noFit:  # it fits and has a start
                                    # assign correct ends
                                    patternStart = correctStart
                                    patternEnd = correctEnd
                                    patternEndNew = correctEnd
                                    # general information
                                    partOfSeq_4 = sequence2[patternStart:patternEnd]
                                    numberOfRepeats = int(len(partOfSeq_4) / patternLen)

                                    partOfSeqNew = pattern * numberOfRepeats
                                    # if STRs should be changed
                                    if random.random() <= chanceOfChange:
                                        # if you want to simulate a reduction you cannot reduce more than the available number of repeats.
                                        # if you want to simulate an increase of repeats, do anything between

                                        manipulation = random.randint(0, 5) if random.random() < 0.5 else (-1) * (
                                            random.randint(0, numberOfRepeats))
                                        numberOfRepeatsNew = numberOfRepeats + manipulation  # total new number of repeats

                                        # patternEndNew = patternStart + numberOfRepeatsNew * patternLen  # current end
                                        # offset += (patternEndNew - patternEnd)  # remember for following
                                        change = (manipulation * patternLen)
                                        # offset += change
                                        partOfSeqNew = pattern * numberOfRepeatsNew  # new middle sequence
                                        # patternEndNew = patternEnd + change
                                        # entrance_c[2] = patternEndNew

                                    # mutate new sequence
                                    partOfSeqNew, offset2 = mutate(partOfSeqNew, chanceOfMutationPerBase,
                                                                   indelsLessMutation)

                                    # debugHelpPartOfSeq9 = sequence2[patternStart - 3:patternEnd + 3]
                                    # cut current sequence replace
                                    patternEndNew += offset2
                                    # sequence2 = sequence2[:patternStart] + "" + sequence2[
                                    #                                            patternEnd:]  # cuts part inbetween
                                    # debugHelpPartOfSeq8 = sequence2[patternStart - 3:patternEnd + 3]
                                    # insert new sequence
                                    # sequence2 = sequence2[:patternStart] + partOfSeqNew + sequence2[
                                    #                                                      patternStart:]  # fills part inbetween
                                    sequence2 = sequence2[:patternStart] + partOfSeqNew + sequence2[patternEnd:]
                                    patternEndNew = patternEnd + (len(partOfSeqNew) - (patternEnd - patternStart))
                                    # debugHelpPartOfSeq6 = sequence2[patternStart - 3:patternEndNew + 3]

                                    # just to see in debugging that string replacement works
                                    # debugHelpPartOfSeq7 = sequence2[patternStart:patternEndNew]
                                    offsettemp = (len(partOfSeqNew) - (patternEnd - patternStart))
                                    offset += offsettemp
                                if noFit:  # no fit didnot match in 10 positions or is no STR anymore therefore is not a good coordinate for a STR
                                    entrance = bedfile_l[bedfEntrance]
                                    print("no fit, entrance: ", entrance)
                                    entrance_c = copy.deepcopy(entrance)
                                    entrance_c[0] = -1
                                    entrance_c[1] = 0  # mark it as 0 later don't put it in new bedfile
                                    entrance_c[2] = 0
                                else:
                                    if allele == 1:
                                        homozygousity_d[chrnr2, ps] = partOfSeqNew
                                    # preparation to change bedfile
                                    entrance = bedfile_l[bedfEntrance]
                                    entrance_c = copy.deepcopy(entrance)  # only change the copy
                                    entrance_c[0] = chrnr_w
                                    entrance_c[1] = patternStart  # + 1
                                    entrance_c[2] = patternEndNew
                                    # if patternEndNew < patternStart+1:
                                    #    entrance_c[2] = patternStart+1

                                # changes in bedfile
                                #if allele == 1:  # only change the first entrance of the copy (first chromosome)
                                bedfile_l_copy[bedfEntrance] = entrance_c
                                #else:  # append the entrances of the second chromosome
                                #    bedfile_l_copy.append(entrance_c)

                        if allele == 1:  # will only be visited once
                            allele += 1
                        record2.seq = sequence2
                        record2.id = idOfChr
                        record2.name = nameOfChr
                        newrecordlen = len(sequence2)
                        print("new length sequence: ", newrecordlen)
                        writer.write_header()
                        writer.write_record(record2)
                        bedfile_total.append(bedfile_l_copy)

        printBedModifications(bedfile_total, newBedFile)



if len(sys.argv) < 6:
    print("Please give fastafilename +dir to the ref [1] and new Fastafilename that will be created [2], the old bedfile [3], "
          "the  name and dir where the new bedfile should be [4], "
          "and a number between >0.0 && <1.0 indicating the changes of STR to the original ref file [5]")
    print(
        "Additional: 'd' oder 'h' for diploid or haploid approach [6], chance of mutation >0.0 && <1.0 [7], how much rarer an indel is than a mutation (0-100, comma possible) [8], homozygousity_rate >0.0 && <1.0 [9]")
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
            if len(sys.argv) > 7:
                print("more than 7 parameters")
                canBeCast = True
                try:
                    float(sys.argv[7])
                except ValueError:
                    print("Not a float")
                    canBeCast = False
                if canBeCast:
                    print("sys.argv[7].isfloat")
                    chanceOfMutation = float(sys.argv[7])
                    if len(sys.argv) > 8:
                        print("more than 8 parameters")
                        canBeCast = True
                        try:
                            float(sys.argv[8])
                        except ValueError:
                            print("Not a float")
                            canBeCast = False
                        if canBeCast:  # can [8] be cast
                            print("sys.argv[8].isfloat")
                            indelLessTmutation = float(sys.argv[8])
                            if len(sys.argv) > 9:  # len > 9
                                print("more than 9 parameters")
                                canBeCast = True
                                try:
                                    float(sys.argv[9])
                                except ValueError:
                                    print("Not a float")
                                    canBeCast = False
                                if canBeCast:
                                    print("sys.argv[9].isfloat")
                                    homozygous = float(sys.argv[9])
                                    if 0.0 <= homozygous <= 1.0:
                                        print("in field 9 homozygous possibilty assigned was: ", homozygous)
                                        print("How much rarer is an Indel than Mutation, according to assignment: ",
                                              indelLessTmutation)
                                        if 0.0 <= chanceOfMutation <= 1.0:
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
                                        if 0.0 <= chanceOfMutation <= 1.0:
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
                                else:  # 9 was not cast.
                                    print(
                                        "wrong format. homozygous possibilty needs a value > 0.0 and smaller <=1.0. default is 0")
                                    print("How much rarer is an Indel than Mutation, according to assignment: ",
                                          indelLessTmutation)
                                    if 0.0 <= chanceOfMutation <= 1.0:
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
                            else:  # len not 0 9
                                print(
                                    "in field 9 homozygous possibilty was missing. here the default is 0")
                                print("How much rarer is an Indel than Mutation, according to assignment: ",
                                      indelLessTmutation)
                                if 0.0 <= chanceOfMutation <= 1.0:
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
                        else:  # [8] can not be cast. therefore default. And 7 fields ok and can be run. but check of range nessesary
                            print("wrong format. in field 8 a number is needed to express how much rarer an indel"
                                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                            if 0.0 <= chanceOfMutation <= 1.0:
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
                    else:  # 7 fields. but and chance of mutation was casted float. But is it in range?
                        if 0.0 <= chanceOfMutation <= 1.0:
                            print("your assigned chance of mutation rate ", chanceOfMutation)
                            print("in field 8 a number is needed to express how much rarer an indel"
                                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                            print("Default is 10 times rarer")
                            print(
                                "in field 9 homozygous possibilty was missing. if 'd' was chosen, here the default is 0")
                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run,
                                              chanceOfMutation, 10, 0)
                        else:
                            print(
                                "Your chance of mutation was not a value between 0 and 1, therefore mutation rate is now at 0.01 == 1%")
                            print("in field 8 a number is needed to express how much rarer an indel"
                                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                            print("Default is 10 times rarer")
                            print(
                                "in field 9 homozygous possibilty was missing. if 'd' was chosen, here the default is 0")
                            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run,
                                              0.01, 10, 0)
                else:  # 7 cannot be cast
                    print("sys.argv[7] is not float, therefore default mutation rate is at 0.01 == 1%")
                    print("in field [8] a number is needed to express how much rarer an indel"
                          "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                    print("Default is 10 times rarer")
                    print("in field [9] homozygous possibilty was missing. here the default is 0")
                    main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run, 0.01, 10, 0)
            else:  # 7 fielts
                print(
                    "len < 7: 0.01 (1%) chance of mutation, as you did not give a parameter for chance of mutation[7]")
                print("in field [8] a number is needed to express how much rarer an indel"
                      "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
                print("Default is 10 times rarer")
                print("in field [9] homozygous possibilty was missing. here the default is 0")
                main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], run, 0.01, 10, 0)
        else:
            print(
                "sys.argv[6] is not the right format (string: 'd' or 'h'), therefore run is default haploid and with 0.01 == 1% mutationrate [7]")
            print("in field [8] a number is needed to express how much rarer an indel"
                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
            print("Default is 10 times rarer")
            print("in field [9] homozygous possibilty was missing. here the default is 0")
            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], 1, 0.01, 10, 0)
    else:  # default
        if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[3]):
            print("given 5 parameters, therefore run is default haploid [6] and with 0.01 == 1% mutationrate [7]")
            print("in field [8] a number is needed to express how much rarer an indel"
                  "should be than a mutation. like mutation = 1% but indel 10 or 9.9 times rarer.")
            print("Default is 10 times rarer")
            print("in field [9] homozygous possibilty was missing. here the default is 0")
            main_manipulation(sys.argv[2], sys.argv[1], sys.argv[4], sys.argv[3], sys.argv[5], 1, 0.01, 10, 0)
            # [0]./simulator2.py [1]../../reference/hs37d5.chr22.fa  [2]hs37d5.chr22.new1.fa [3]../bedfiles/hs37_ver8.chr22.bed [4]hs37_ver8.chr22.new1.bed [5]0.20
        else:
            if not os.path.isfile(sys.argv[1]):
                print(sys.argv[1], " is not a file")
            elif not os.path.isfile(sys.argv[3]):
                print(sys.argv[3], " is not a file")


