# simulates STR changes compared to a reference genome
# hd stands for haploid and diploid, as in this version haploid or diploid genomes can be simulated

import copy
import random
from Bio import SeqIO
import re
import argparse



from simiSTR_utils import SimiSTR_Writer
from simiSTR_utils import SimiSTR_Reader



class SimiSTR:

    def __init__(self, input_fasta, output_fasta, input_bedfile, output_bedfile,
                  expansion_possibility, diploidity, snv_chance, less_indels,
                  homozygousity, max_add=5, max_reduction=-1):
        self.input_fasta = input_fasta
        self.output_fasta = output_fasta
        self.input_bedfile = input_bedfile
        self.output_bedfile = output_bedfile
        self.expansion_possibility = expansion_possibility
        self.max_add = max_add
        self.max_reduction = max_reduction
        self.diploidity = diploidity
        self.snv_chance = snv_chance
        self.less_indels= less_indels
        self.homozygousity = homozygousity

    # mutation of sequence by chance and less likely mutate by insertion or deletion
    def snv_mutation(self, sequence, chanceOfSubstitution, indelsLessMutation):
        #this counts how much longer or shorter the sequence gets through the indels
        offsetIn = 0
        offsetOut = 0
        #this calculations avoid dividsion through zero
        chanceForIndels = self.calculateChanceForIndels(chanceOfSubstitution,indelsLessMutation)

        # here every base in the given region gets called +
        # and per random decision manipulated by substitution, insertion, deletion
        for base_idx in range(0,len(sequence)):
            randomSubstitution = random.random() # random number for substitituion
            randomIndel = random.random() # random number for insertions / deletions
            if randomSubstitution < chanceOfSubstitution:
                #substitute a base on this position
                sequence = self.substitutions(sequence, base_idx)
            elif randomIndel < chanceForIndels:  # zb 1% / 10  = 0.1%
                #insert or deletion following base in sequence and calculate the resulting offset
                sequence, offsetOut = self.indeletion(sequence, base_idx, offsetIn)
                offsetOut += offsetOut
        return sequence, offsetOut
    # needed be mutations function
    def calculateChanceForIndels(self, chanceOfSubstitution, indelsLessMutation):
        chanceForIndels = 0
        if indelsLessMutation != 0:  # to avoid division through zero
            chanceForIndels = chanceOfSubstitution / indelsLessMutation
        return chanceForIndels
    # certain amount of snv are substitutions
    def substitutions(self, sequence, base_idx):
        rand_num = random.random()
        chanceForTransition = 0.85
        sequence = self.transition(sequence, base_idx) if rand_num < chanceForTransition else self.transversion(sequence, base_idx)
        return sequence
    # certain amount of snv are transisitons
    def transition(self, sequence, base_idx):
        base = sequence[base_idx]
        replacementBase = ''
        if base == 'A':
            replacementBase = 'G'
        elif base == 'C':
            replacementBase = 'T'
        elif base == 'G':
            replacementBase = 'A'
        else: #base == 'T'
            replacementBase = 'C'
        sequence = sequence[:base] + replacementBase + sequence[base + 1:]
        return sequence
    # cretain amount of snv are transversions
    def transversion(self, sequence, base_idx):
        rand_num = random.random()
        base = sequence[base_idx]
        replacementBase = ''
        if base == 'A':
            replacementBase = 'T' if rand_num < 0.5 else 'C'
        elif base == 'C':
            replacementBase = 'G' if rand_num < 0.5 else 'A'
        elif base == 'G':
            replacementBase = 'C' if rand_num < 0.5 else 'T'
        else: #base == 'T'
            replacementBase = 'A' if rand_num < 0.5 else 'G'
        sequence = sequence[:base] + replacementBase + sequence[base + 1:]
        return sequence
    # certain mutations are insertions or deletions
    def indeletion(self, sequence, base_idx, offset):
        if base_idx + 1 < len(sequence):
            if random.random() < 0.5:  # insertion occurs in 50% of indel cases
                sequence = self.insertion(sequence,base_idx)
                offset += 1
            else:  # deletion occurs in the other 50% of indel cases
                sequence = self.deletion(sequence, base_idx)
                offset -= 1
        return sequence, offset
    # insert a certain base
    def insertion(self, sequence, base_idx):
        chanceForCertainBase = random.random()
        insertionBase = ''
        if chanceForCertainBase < 0.25:    # insertion of G
            insertionBase = 'G'
        elif chanceForCertainBase < 0.5:   # insertion of C
            insertionBase = 'C'
        elif chanceForCertainBase < 0.75:  # insertion of A
            insertionBase = 'A'
        else:                              # insertion of T
            insertionBase = 'T'
        sequence = sequence[:base_idx] + insertionBase + sequence[base_idx:]  # insertion
        return sequence
    # delete a certain base
    def deletion(self, sequence, base):
        sequence = sequence[:base] + sequence[base + 1:]
        return sequence


    # look up if bedfile coordinates match position & pattern on position in reference/ assigned genome
    # of find correct position near by (+/- 20)
    def findStartPoint(self, seq, start, pattern, patternLen):
        seq_len = len(seq)
        partOfSeq_1 = (seq[start:start + patternLen])
        startpoint = start  # if start is 1000
        endpoint = start + patternLen
        startpoint_memory = startpoint
        endpoint_memory = endpoint
        startPointCorrect = False
        lower = False
        number = -1  # -1
        noFit = False
        times = False
        while not startPointCorrect and not noFit:

            #  => replace this start-position search by sth linear or with a smith waterman.

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
                        startpoint = startpoint_memory
                        endpoint = endpoint_memory

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
                    endpoint_memory = endpoint
                    startpoint_memory = startpoint
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

    def STRexpansion(self, numberOfRepeats, patternLen, pattern):
        # if you want to simulate a reduction you cannot reduce more than the available number of repeats.
        # if you want to simulate an increase of repeats, do anything between
        expansion_factor_minus = 0
        if self.max_reduction == -1:
            expansion_factor_minus = numberOfRepeats
        else:
            if self.max_reduction < numberOfRepeats:
                expansion_factor_minus = self.max_reduction
            else:
                expansion_factor_minus = numberOfRepeats

        # here the mutations gets randomly calculated.
        manipulation = random.randint(0, self.max_add) if random.random() < 0.5 else (-1) * (
            random.randint(0, expansion_factor_minus))

        numberOfRepeatsNew = numberOfRepeats + manipulation  # total new number of repeats

        # patternEndNew = patternStart + numberOfRepeatsNew * patternLen  # current end
        # offset += (patternEndNew - patternEnd)  # remember for following
        change = (manipulation * patternLen)
        # offset += change
        partOfSeqNew = pattern * numberOfRepeatsNew  # new middle sequence
        # patternEndNew = patternEnd + change
        # entrance_c[2] = patternEndNew
        return partOfSeqNew

    # main function manipulate Fasta
    def main_manipulation(self):
        # Initialize Reader & Writer
        sReader = SimiSTR_Reader(self.input_bedfile)
        sWriter = SimiSTR_Writer(self.output_bedfile)

        # cast inputs
        #expansion_chance = float(self.expansion_possibility) # maybe not necessary anymore should be float due to argparse
        #snv_chance = float(self.snv_chance) # maybe not necessary anymore should be float due to argparse

        # safe bedfile as lists.
        bedfile_d = sReader.getBedFile()
        # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
        bedfile_total = list()

        with open(self.output_fasta, 'w') as outFastaFile:
            writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
            with open(self.input_fasta, 'r') as inFastaFile:  # read fastaFile
                for record in SeqIO.parse(inFastaFile, "fasta"):  # every Record. 1 .... 2 .... 3..... .... 22 .... x...
                    sequence = (record.seq).upper()
                    recordLen = len(sequence)  # old length
                    print("old length", recordLen)
                    homozygousity_d = dict()
                    allele = 1
                    for chr in range(0, self.diploidity):  # per chromosome a new to be created chromosome #eigther only "0"(haploid) oder "0 and 1" (diploid)
                        record2 = copy.deepcopy(record)  # changes only on deep copies.
                        sequence2 = copy.deepcopy(sequence)  # changes only on deep copies.
                        # have an offset, that tells how much all following coordinates will be later or earlier
                        offset = 0  # original in the beginning
                        # naming of haploid and or diploid chromosome entrances.
                        nameOfChr = record2.name
                        idOfChr = record2.id
                        #id = re.search("(\d*)",idOfChr)  # this line should enable to find the 1 in the "chr1" annotation

                        id = [int(s) for s in re.findall(r'\d+', idOfChr)]
                        if id is not []:
                            chrNr = id[0]

                            nameOfChr = nameOfChr + "_" + str(allele)
                            idOfChr = idOfChr + "_" + str(allele)
                            record2.name = nameOfChr
                            record2.id = id[0]

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

                                    if allele > 1 and random.random() < self.homozygousity:  # should be the second allele and inside the homozygosity rate
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
                                        entrance_cn[1] = patternStart  # + 1 #+1 when working with gangstr files
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
                                        correctStart, correctEnd, noFit = self.findStartPoint(sequence2, patternStart, pattern,
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
                                            # region-sequence if no STR Expansion will occur
                                            partOfSeqNew = pattern * numberOfRepeats

                                            # STRs EXPANSION
                                            chanceForExpansion = random.random()
                                            if chanceForExpansion <= self.expansion_possibility:
                                                # region-sequence get recalculated with new length
                                                partOfSeqNew = self.STRexpansion(numberOfRepeats,patternLen,pattern)

                                            # mutate new sequence
                                            partOfSeqNew, offset2 = self.snv_mutation(partOfSeqNew, self.snv_chance, self.less_indels)

                                            # debugHelpPartOfSeq9 = sequence2[patternStart - 3:patternEnd + 3]
                                            # cut current sequence replace
                                            #patternEndNew += offset2
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
                                            entrance_c[1] = patternStart  # + 1 when working
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
                        else:
                            print("Check the first column in your assigned input bed file!")

            sWriter.printBedModifications(bedfile_total)








parser = argparse.ArgumentParser(description="Run SimiSTR to change Expansionlength of STRs.")

parser.add_argument('-inf', '--input_fasta', type=str, metavar='', required=True, help="Path+Name to Fasta File that is template that needs STR expansion changes")
parser.add_argument('-outf', '--output_fasta', type=str, metavar='', required=True, help="Path+Name for newly generated Fasta File with expasion changes")
parser.add_argument('-ibf', '--input_bedfile', type=str, metavar='', required=True, help="Path+Name to Bedfile containing regions of known STRs in given Input Fasta")
parser.add_argument('-obf', '--output_bedfile', type=str, metavar='', required=True, help="Path+Name to Bedfile containing information about applied changes in given STR regions")
parser.add_argument('-expp', '--expansion_possibility',type=float, metavar='', required=True,  help="How many regions should be STR expansion length manipulated")
parser.add_argument('-dip', '--diploidity', type=int, metavar='', required=True, help="Diploid= 2 , Haploid= 1. Multiploid is not yet implemented" )
parser.add_argument('-snv', '--snv_chance', type=float, metavar='', required=True, help="[0.000-1.000] is the chance of a SNV.")
parser.add_argument('-lid', '--less_indels', type=int, metavar='', required=True, help="How much rarer should a insertion/deletion occur than a substitution.")
parser.add_argument('-ho', '--homozygousity', type=float, metavar='', required=True, help="How many regions should be homzygous. The rest will be heterozygous.")
parser.add_argument('-ma', '--max_add', type=int, metavar='', required=False, help="How many repeats can maximum be added [default: 5]")
parser.add_argument('-mr', '--max_reduction', type=int,metavar='', required=False, help="How many repeats can maximum be removed. [default: full length]")
args = parser.parse_args()


if __name__ == '__main__':
    sim = SimiSTR(args.input_fasta, args.output_fasta, args.input_bedfile, args.output_bedfile,
                  args.expansion_possibility, args.diploidity, args.snv_chance, args.less_indels,
                  args.homozygousity, args.max_add, args.max_reduction)
    sim.main_manipulation()
