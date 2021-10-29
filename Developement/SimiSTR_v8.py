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
                  homozygousity, max_add=5, max_reduction=-1, gangSTRflag=0):
        self.input_fasta = input_fasta
        self.output_fasta = output_fasta
        self.input_bedfile = input_bedfile
        self.output_bedfile = output_bedfile
        self.expansion_possibility = expansion_possibility
        self.max_add = max_add
        self.max_reduction = max_reduction
        self.diploidity = diploidity if diploidity is not None else 1
        self.snv_chance = snv_chance
        self.less_indels= less_indels
        self.homozygousity = homozygousity
        self.gangstr_flag = gangSTRflag

    # mutation of sequence by chance and less likely mutate by insertion or deletion
    def snv_mutation(self, sequence, chanceOfSubstitution, indelsLessMutation):
        #this counts how much longer or shorter the sequence gets through the indels
        nrOfIndels = [0,0] # [insertions, deletions]
        nrOfSubstitutions = 0
        #this calculations avoid dividsion through zero
        chanceForIndels = self.calculateChanceForIndels(chanceOfSubstitution,indelsLessMutation)

        # here every base in the given region gets called +
        # and per random decision manipulated by substitution, insertion, deletion
        for base_idx in range(0,len(sequence)):
            randomSubstitution = random.random() # random number for substitituion
            randomIndel = random.random() # random number for insertions / deletions
            if randomSubstitution < chanceOfSubstitution and base_idx+1 < len(sequence):
                #substitute a base on this position
                sequence, nrOfSubstitutions = self.substitutions(sequence, base_idx, nrOfSubstitutions)
            if randomIndel < chanceForIndels and base_idx+1 < len(sequence):  # zb 1% / 10  = 0.1%
                #insert or deletion following base in sequence and calculate the resulting offset
                sequence, nrOfIndels = self.indeletion(sequence, base_idx, nrOfIndels)

        return sequence, nrOfSubstitutions, nrOfIndels


    # needed be mutations function
    def calculateChanceForIndels(self, chanceOfSubstitution, indelsLessMutation):
        chanceForIndels = 0
        if indelsLessMutation != 0:  # to avoid division through zero
            chanceForIndels = chanceOfSubstitution / indelsLessMutation
        return chanceForIndels

    # certain amount of snv are substitutions
    def substitutions(self, sequence, base_idx,  nrOfSubstitutions):
        rand_num = random.random()
        chanceForTransition = 0.85
        sequence, nrOfSubstitutions = self.transition(sequence, base_idx,nrOfSubstitutions) if rand_num < chanceForTransition else self.transversion(sequence, base_idx,nrOfSubstitutions)
        return sequence, nrOfSubstitutions

    # certain amount of snv are transisitons
    def transition(self, sequence, base_idx, nrOfSubstitutions):
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
        sequence_new = sequence[:base_idx] + replacementBase + sequence[base_idx + 1:]
        return sequence_new, (nrOfSubstitutions+1)
    # cretain amount of snv are transversions

    def transversion(self, sequence, base_idx, nrOfSubstitutions):
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
        sequence = sequence[:base_idx] + replacementBase + sequence[base_idx + 1:]
        return sequence, (nrOfSubstitutions+1)

    # certain mutations are insertions or deletions
    def indeletion(self, sequence, base_idx, nrOfIndel):
        nrOfInsertions,nrOfDeletion = nrOfIndel[0],nrOfIndel[1]
        if base_idx + 1 < len(sequence):
            if random.random() < 0.5:  # insertion occurs in 50% of indel cases
                sequence = self.insertion(sequence,base_idx)
                nrOfInsertions += 1
            else:  # deletion occurs in the other 50% of indel cases
                sequence = self.deletion(sequence, base_idx)
                nrOfDeletion += 1
        return sequence, [nrOfInsertions,nrOfDeletion]

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
    def deletion(self, sequence, base_idx):
        sequence = sequence[:base_idx] + sequence[base_idx + 1:]
        return sequence


    # look up if bedfile coordinates match position & pattern on position in reference/ assigned genome
    # of find correct position near by (+/- 20)
    def findStartPoint(self, seq, start, pattern, patternLen):
        seq_len = len(seq)
        # initial sequence, startpoint till endpoint at end of patternlength
        partOfSeq_1 = (seq[start:start + patternLen])
        startpoint = start  # if start is 1000
        endpoint = start + patternLen
        # in case it was not found
        startpoint_memory = startpoint
        endpoint_memory = endpoint
        # flags that help find pattern at correct start with correct length
        startPointCorrect = False
        upstream = False
        offset = -1  # -1 # offset from initial start
        noFit = False
        stopTrying = False

        while not startPointCorrect and not noFit:
            if pattern != partOfSeq_1:  # pattern was not found on initial start point
                # will always go one further away from current start. trying both directions parallel
                startpoint = startpoint + offset  # start is 999, 1001, 998, 1002 usw.
                partOfSeq_1 = seq[startpoint:startpoint + patternLen]
                if partOfSeq_1 == pattern:  # found
                    offset = -1  # reset and leave
                else:
                    if upstream:
                        offset += 1  # 3
                        offset *= (-1)  # -3
                        upstream = False  # switch to downstream for next round
                    else:  #downstream
                        offset *= (-1)  # 1 .... 3
                        offset += 1  # 2 .... 4
                        upstream = True  # switch to upstream for next round
                    if abs(offset) > 30:  # no possible patternstart or STR within +/- 10 positions distance of bedfile-patternstart
                        noFit = True
                        offset = -1 #reset
                        startpoint = startpoint_memory
                        endpoint = endpoint_memory

            else:  # if a pattern is found
                startPointCorrect = True  # found startpoint
                patternStart = startpoint
                # find out if area in reference of this STR is acutally longer or shorter than according to
                # the bed say? basically count if pattern goes left an right for longer.
                # maybe later change to another comparison algorithm that takes costs into account (Blast)

                # __ Leftside check __ #
                check_upstream = patternStart - patternLen
                goUpstream = True
                while check_upstream >= 0 and goUpstream:
                    partOfSeq_2 = seq[check_upstream:check_upstream + patternLen]
                    if partOfSeq_2 == pattern:  # still matching
                        check_upstream = check_upstream - patternLen  # maybe your can even go further left
                    else:
                        goUpstream = False  # stop checking further left
                        check_upstream = check_upstream + patternLen  # because last position did not match, so go back to last
                        # one working

                # __ Rightside check __ #
                check_downstream = check_upstream
                goDownstream = True
                while check_downstream + patternLen <= seq_len and goDownstream:
                    partOfSeq_3 = seq[check_downstream:check_downstream + patternLen]
                    if partOfSeq_3 == pattern:  # still matching
                        check_downstream = check_downstream + patternLen  # maybe your can even go further right
                    else:
                        goDownstream = False  # stop checking further right

                startpoint = check_upstream
                endpoint = check_downstream

                if endpoint - startpoint == patternLen:  # length only one time the pattern
                    endpoint_memory = endpoint
                    startpoint_memory = startpoint
                    if not stopTrying:
                        #set the offset further off the regions searched already
                        if offset < 0:
                            offset += 9
                        else:
                            offset -= 9
                        startPointCorrect = False
                        startpoint += offset
                        partOfSeq_1 = seq[startpoint:startpoint + patternLen]
                        stopTrying = True
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
            if self.max_reduction <= numberOfRepeats:
                expansion_factor_minus = self.max_reduction
            else:
                expansion_factor_minus = numberOfRepeats

        # here the expansion length gets randomly calculated.
        manipulation = random.randint(0, self.max_add) if random.random() < 0.5 else (-1) * (random.randint(0, expansion_factor_minus))

        numberOfRepeatsNew = numberOfRepeats + manipulation  # total new number of repeats
        baseNrchange = (manipulation * patternLen)
        partOfSeqNew = pattern * numberOfRepeatsNew  # new middle sequence

        return partOfSeqNew, baseNrchange

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
                    for allele in range(1, (self.diploidity+1)):  # per allele create a chromosome #eigther only "1"(haploid) oder "1 and 2" (diploid)
                        record2 = copy.deepcopy(record)  # changes only on deep copies.
                        sequence2 = copy.deepcopy(sequence)  # changes only on deep copies.
                        # have an offset, that tells how much all following coordinates will be later or earlier
                        offset = 0  # original in the beginning
                        # naming of haploid and or diploid chromosome entrances.
                        nameOfChr = record2.name
                        idOfChr = record2.id
                        #id = re.search("(\d*)",idOfChr)  # this line should enable to find the 1 in the "chr1" annotation

                        id = [int(s) for s in re.findall(r'\d+', idOfChr)]
                        print(id)
                        if id is not []:
                            chrNr = id[0] # id is a list the size 1 from the regex

                            nameOfChr = nameOfChr + "_" + str(allele)
                            idOfChr = idOfChr + "_" + str(allele)
                            record2.name = nameOfChr
                            record2.id = chrNr

                            id = record2.id
                            if id in bedfile_d.keys():
                                bedfile_l = bedfile_d[id]
                                bedfile_l_copy = copy.deepcopy(bedfile_l)
                                bedfile_l_length = len(bedfile_l)

                                for bedfile_idx in range(0,
                                                          bedfile_l_length):  # go through all coordinates in the bedfile.
                                    shortTR = bedfile_l[bedfile_idx]
                                    chrnr = shortTR[0]  # which chromosome
                                    chrnr_w = str(chrnr) + "_" + str(allele)  # this number will be written down in the out-bedfile

                                    patternStart = offset + int(shortTR[1])
                                    patternEnd = offset + int(shortTR[2])
                                    patternLen = int(shortTR[3])
                                    pattern = shortTR[4].strip()

                                    """Cause Gangstr files have a different start point"""
                                    if self.gangstr_flag:
                                        patternStart = offset + int(shortTR[1]) - 1 #-1 for gangstr bedfiles

                                    if allele > 1 and random.random() < self.homozygousity:
                                        # should be the second allele and inside the homozygosity rate
                                        # if below chance then copy allele 1 == homozygous
                                        chrnr2 = chrnr  #original chromosome
                                        ps = shortTR[1] #patternstart
                                        # dictionary on position chromosome+ps gives back the sequence from the first allele
                                        partOfSeq = homozygousity_d[chrnr2, ps]

                                        # cut current sequence replace
                                        sequence2 = sequence2[:patternStart] + "" + sequence2[
                                                                                    patternEnd:]  # cuts part inbetween
                                        # insert new sequence
                                        sequence2 = sequence2[:patternStart] + partOfSeq + sequence2[
                                                                                           patternStart:]  # fills part inbetween

                                        entrance_allele1 = bedfile_l_copy[bedfile_idx]
                                        entrance_allele2 = copy.deepcopy(entrance_allele1)
                                        entrance_allele2[0] = chrnr_w   # chrNr
                                        entrance_allele2[1] = patternStart # allele2 start
                                        """allele2 pattern start for gangstr"""
                                        if self.gangstr_flag:
                                            entrance_allele2[1] = patternStart + 1 #+1 when working with gangstr files

                                        entrance_allele2[2] = patternStart + len(partOfSeq) # allele 2 pattern end new
                                        # entrance_allele2[3] == patternlength stays the same
                                        # entrance_allele2[4] == pattern stays the same

                                        #calculate new offset
                                        oldSeqLen = (patternEnd - patternStart)
                                        offset_allele2 = len(partOfSeq) - oldSeqLen
                                        offset += offset_allele2

                                        entrance_allele2[5] = offset_allele2
                                        #entrance_allele2[6] = snvs dont exist more than on allele 1 cause homozygous
                                        bedfile_l_copy[bedfile_idx] = entrance_allele2
                                    else:
                                        # this is the first or first and only allele getting mutated
                                        chrnr2 = chrnr
                                        ps = shortTR[1]  #ps patternstart
                                        # create entrance to safe the allele for second allele if necessary
                                        homozygousity_d[chrnr2, ps] = ""
                                        # find correct startpoit or if bedfile-entrance does not fit to sequence on that position
                                        correctStart, correctEnd, noFit = self.findStartPoint(sequence2, patternStart, pattern,
                                                                                         patternLen)  # seq,start,pattern,patternLen

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
                                                partOfSeqNew, expBaseNrchange = self.STRexpansion(numberOfRepeats,patternLen,pattern)

                                            # mutate new sequence
                                            partOfSeqNew, noOfSubstitution, nrOfIndels = self.snv_mutation(partOfSeqNew, self.snv_chance, self.less_indels)

                                            sequence2 = sequence2[:patternStart] + partOfSeqNew + sequence2[patternEnd:]
                                            patternEndNew = patternEnd + (len(partOfSeqNew) - (patternEnd - patternStart))
                                            offsettemp = (len(partOfSeqNew) - (patternEnd - patternStart))
                                            offset += offsettemp
                                        if noFit:  # no fit didnot match in 10 positions or is no STR anymore therefore is not a good coordinate for a STR
                                            entrance = bedfile_l[bedfile_idx]
                                            print("no fit, entrance: ", entrance)
                                            entrance_allele1 = copy.deepcopy(entrance)
                                            entrance_allele1[0] = -1
                                            entrance_allele1[1] = 0  # mark it as 0 later don't put it in new bedfile
                                            entrance_allele1[2] = 0
                                            # NOT FOUND SO NOT MANIPULATED
                                            # entrance_allele1[3] == patternlength does not exist
                                            # entrance_allele1[4] == pattern does not exist
                                            # entrance_allele1[5] = offset does not exist
                                            # entrance_allele1[6] = snvs dont exist
                                        else:
                                            if allele == 1:
                                                homozygousity_d[chrnr2, ps] = partOfSeqNew
                                            # preparation to change bedfile
                                            entrance = bedfile_l[bedfile_idx]
                                            entrance_allele1 = copy.deepcopy(entrance)  # only change the copy
                                            entrance_allele1[0] = chrnr_w
                                            if self.gangstr_flag:
                                                entrance_allele1[1] = patternStart + 1 # when working with one-based files
                                            else:
                                                entrance_allele1[1] = patternStart
                                            entrance_allele1[2] = patternEndNew

                                            # entrance_allele1[3] == patternLen stays the same
                                            # entrance_allele1[4] == pattern stays the same
                                            entrance_allele1[5] = expBaseNrchange #nr of bases changed through expansion change
                                            entrance_allele1[6] = noOfSubstitution
                                            entrance_allele1[7] = nrOfIndels #insertion, deletion


                                        # changes in bedfile
                                        #if allele == 1:  # only change the first entrance of the copy (first chromosome)
                                        bedfile_l_copy[bedfile_idx] = entrance_allele1
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

parser.add_argument('-inf', '--input_fasta', type=str,  required=True, help="Path+Name to Fasta File that is template that needs STR expansion changes")
parser.add_argument('-outf', '--output_fasta', type=str,  required=True, help="Path+Name for newly generated Fasta File with expansion changes")
parser.add_argument('-ibf', '--input_bedfile', type=str,  required=True, help="Path+Name to Bedfile containing regions of known STRs of given Input Fasta")
parser.add_argument('-obf', '--output_bedfile', type=str,  required=True, help="Path+Name to Bedfile containing information about applied changes in given STR regions")
parser.add_argument('-expp', '--expansion_possibility',type=float,  required=True,  help="[0.000-1.000] How many regions should be STR expansion length manipulated")
parser.add_argument('-dip', '--diploidity', type=int,choices=range(1, 3),  required=True, help="[1-2] Diploid= 2 , Haploid= 1. Multiploid is not yet implemented" )
parser.add_argument('-snv', '--snv_chance', type=float,  required=True, help="[0.000-1.000] is the chance of a SNV.")
parser.add_argument('-lid', '--less_indels', type=int,  required=True, help="[int] How much rarer should a insertion/deletion occur than a substitution.")
parser.add_argument('-ho', '--homozygousity', type=float,  required=True, help="[0.000-1.000] How many regions should be homzygous. The rest will be heterozygous.")
parser.add_argument('-ma', '--max_add', type=int,  required=False, help="[int] How many repeats per STR can maximum be added [default: 5]")
parser.add_argument('-mr', '--max_reduction', type=int, required=False, help="[int] How many repeats per STR can maximum be removed. [default: full length]")
parser.add_argument('-g', '--gangstr_flag', type=int,choices=range(0, 2),  required=False, help="[0-1] GangstrFile=1, else=0 [default: 0]")
args = parser.parse_args()


if __name__ == '__main__':
    sim = SimiSTR(args.input_fasta, args.output_fasta, args.input_bedfile, args.output_bedfile,
                  args.expansion_possibility, args.diploidity, args.snv_chance, args.less_indels,
                  args.homozygousity, args.max_add, args.max_reduction, args.gangstr_flag)
    sim.main_manipulation()
