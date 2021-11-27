# simulates STR changes compared to a reference genome
# hd stands for haploid and diploid, as in this version haploid or diploid genomes can be simulated

import copy
import random
from Bio import SeqIO
import re
import argparse

# for future quality of life improvement
# from pybedtools import BedTool # conda install --channel conda-forge --channel bioconda pybedtools
                               # pip3 install pybedtools

from simiSTR_utils import SimiSTR_Writer
from simiSTR_utils import SimiSTR_Reader
from simiSTR_utils import SimiSTR_bedline
from simiSTR_utils import logger



class SimiSTR:

    def __init__(self, user_args):
        # user_args is of class argparse
        self.input_fasta = user_args.input_fasta
        self.output_fasta = user_args.output_fasta
        self.input_bedfile = user_args.input_bedfile
        self.output_bedfile = user_args.output_bedfile
        self.expansion_possibility = user_args.expansion_possibility
        self.max_add = user_args.max_add
        self.max_reduction = user_args.max_reduction
        self.snv_chance = user_args.snv_chance
        # Luis: I changed it to a probability
        self.indels_chance = self.calculateChanceForIndels(user_args.snv_chance, user_args.less_indels)
        self.homozygousity = user_args.homozygousity
        self.gangstr_flag = user_args.gangstr_flag
        self.ploidy = int(user_args.ploidy) if user_args.ploidy is not None and (int(user_args.ploidy) == 1 or int(user_args.ploidy) == 2) else 1

    # mutation of sequence by chance and less likely mutate by insertion or deletion
    def snv_mutation(self, sequence):
        # uses self.snv_chance and self.indels_chance
        #this counts how much longer or shorter the sequence gets through the indels
        nrOfIndels = [0,0] # [insertions, deletions]
        nrOfSubstitutions = 0
        # here every base in the given region gets called +
        # and per random decision manipulated by substitution, insertion, deletion
        for base_idx in range(0,len(sequence)):
            randomSubstitution = random.random() # random number for substitituion
            randomIndel = random.random() # random number for insertions / deletions
            if randomSubstitution < self.snv_chance and base_idx+1 < len(sequence):
                #substitute a base on this position
                sequence, nrOfSubstitutions = self.substitutions(sequence, base_idx, nrOfSubstitutions)
            if randomIndel < self.indels_chance and base_idx+1 < len(sequence):  # zb 1% / 10  = 0.1%
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
        elif chanceForCertainBase >= 0.25 and chanceForCertainBase < 0.5:   # insertion of C
            insertionBase = 'C'
        elif chanceForCertainBase >= 0.5 and chanceForCertainBase < 0.75:  # insertion of A
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

        # if no parameter transmitted for max possible reduction of STR repeats.
        if self.max_reduction == -1:
            #then the maximum possible reduction is
            expansion_factor_minus = numberOfRepeats
        else:
            #if a parameter was assigned, assigned max_reduction has to be smaller then the the current repeatlength of the STR
            if self.max_reduction < numberOfRepeats:
                expansion_factor_minus = self.max_reduction
            else: #else reduce maximum the full number of repeats
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
        # bedfile_d = all the position in the bed file, id == chromomsome name
        bedfile_d = sReader.getBedFile()
        # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
        bedfile_total = list()

        outFastaFile = open(self.output_fasta, 'w')
        inFastaFile = open(self.input_fasta, 'r')  # read fastaFile
        writer = SeqIO.FastaIO.FastaWriter(outFastaFile)

        for record in SeqIO.parse(inFastaFile, "fasta"):  # every Record. 1 .... 2 .... 3..... .... 22 .... x...
        # @Damaris the original is in 'record' so 'sequence' can modified wirh no need of new
            recordLen = len(record.seq)  # old length

            logger(f'{record.id}', "info")
            logger(f'Original length {recordLen}\n', "info")

            homozygousity_d = dict()
            allele = 1 # defaul is haploid
            while allele <= self.ploidy:
                # FIXME: time consuming ??
                sequence = (record.seq).upper()
                record2 = copy.deepcopy(record)

                # have an offset, that tells how much all following coordinates will be later or earlier
                offset = 0  # original in the beginning
                # naming of haploid and or diploid chromosome entrances.

                id = record.id # same as in bedfile
                allele_id = record.id + "_" + str(allele)
                name_allele = record.name + "_" + str(allele)

                # update record2 id and name
                record2.id = allele_id
                record2.name = name_allele

                # use only "known" chromosomes (id)
                if id in bedfile_d.keys():
                    # bedfile_l is the list of all entrances in the bedfile for a ginve chromomsome (id)
                    # e.g chr1 start end length motif ... expansion_len nSEV nINDEL

                    bedfile_l = bedfile_d[id]
                    bedfile_l_output = copy.deepcopy(bedfile_l)

                    bedfile_idx = 0
                    for bedfile_entry in bedfile_l:
                        shortTR = SimiSTR_bedline(bedfile_entry)
                        chrnr = shortTR.chromosome
                        chrnr_w = str(chrnr) + "_" + str(allele)
                        ps = shortTR.start  #ps patternstart , it does NOT HAVE the gangstr mod. We look up for originals in list (for ploidy > 1)

                        """Cause Gangstr files have a different start point (gangstr docu)"""
                        patternStart = offset + shortTR.start - 1  if self.gangstr_flag  else offset + shortTR.start
                        patternEnd = offset + shortTR.end # see ganstr documentation
                        patternLen = shortTR.motiflength
                        pattern = shortTR.motif

                        is_homozygous = self.homozygousity >= random.random()
                        log_zygosity = 'NA' if allele == 1 else format(is_homozygous)
                        logger(f'allele -> {allele}; position -> {patternStart} | is homozygous: {log_zygosity}', "info")
                        # default is haploid, do this first and also if not homozygours
                        if (allele == 1) or (allele > 1 and not is_homozygous):
                            # this is the first or first and only allele getting mutated
                            chrnr2 = chrnr
                            # create entrance to safe the allele for second allele if necessary
                            homozygousity_d[chrnr2, ps] = ""
                            # find correct startpoit or if bedfile-entrance does not fit to sequence on that position
                            # seq,start,pattern,patternLen
                            # @Damaris why patternStart and not ps
                            # FIXMI: passing the whole sequence is an overkill ??
                            correctStart, correctEnd, noFit = self.findStartPoint(sequence, patternStart, pattern, patternLen)

                            # no fit didnot match in 10 positions or is no STR anymore therefore is not a good coordinate for a STR
                            if noFit:
                                shortTR_allele = SimiSTR_bedline(bedfile_entry)  # entrance
                                shortTR_allele1 = SimiSTR_bedline(bedfile_entry)  # entrance_allele1

                                logger("Entrance: {0} from your bedfile couldn't be located in the given genome.".format(shortTR_allele), "warning")
                                logger(f'no fit, entrance: {entrance}', "warning")

                                shortTR_allele1.chromosome = "NA" # -1 ???
                                shortTR_allele1.start = 0  # mark it as 0 later don't put it in new bedfile
                                shortTR_allele1.end = 0
                            else:
                                # it fits and has a start
                                # assign correct ends
                                patternStart = correctStart
                                patternEnd = correctEnd
                                patternEndNew = correctEnd
                                # general information
                                partOfSeq_4 = sequence[patternStart:patternEnd]
                                numberOfRepeats = int(len(partOfSeq_4) / patternLen)
                                # region-sequence if no STR Expansion will occur
                                partOfSeqNew = pattern * numberOfRepeats

                                # missing default value for false
                                expBaseNrchange = 0

                                # STRs EXPANSION
                                chanceForExpansion = random.random()
                                if chanceForExpansion <= self.expansion_possibility:
                                    # region-sequence get recalculated with new length
                                    partOfSeqNew, expBaseNrchange = self.STRexpansion(numberOfRepeats,patternLen,pattern)

                                # mutate new sequence
                                partOfSeqNew, noOfSubstitution, nrOfIndels = self.snv_mutation(partOfSeqNew)

                                sequence = sequence[:patternStart] + partOfSeqNew + sequence[patternEnd:]
                                patternEndNew = patternEnd + (len(partOfSeqNew) - (patternEnd - patternStart))
                                offsettemp = (len(partOfSeqNew) - (patternEnd - patternStart))
                                offset += offsettemp

                                if allele == 1:
                                    homozygousity_d[chrnr2, ps] = partOfSeqNew
                                # preparation to change bedfile
                                shortTR_allele = SimiSTR_bedline(bedfile_entry)  # entrance
                                shortTR_allele1 = SimiSTR_bedline(bedfile_entry)  # entrance_allele1
                                shortTR_allele1.chromosome = chrnr_w
                                if self.gangstr_flag:
                                    shortTR_allele1. start = patternStart + 1 # when working with one-based files
                                else:
                                    shortTR_allele1.start = patternStart
                                shortTR_allele1.end = patternEndNew
                                shortTR_allele1.expansion_length = expBaseNrchange #nr of bases changed through expansion change
                                shortTR_allele1.nSNV = noOfSubstitution
                                shortTR_allele1.nINDEL = nrOfIndels #insertion, deletion

                            # changes in bedfile
                            #if allele == 1:  # only change the first entrance of the copy (first chromosome)
                            bedfile_l_output[bedfile_idx] = shortTR_allele1.listBedline()
                            #else:  # append the entrances of the second chromosome
                            #    bedfile_l_output.append(entrance_c)

                        # second allele if needed
                        # should be the second allele and inside the homozygosity rate
                        # if below chance then copy allele 1 == homozygous
                        if allele > 1 and is_homozygous:
                            chrnr2 = chrnr  #original chromosome
                            # dictionary on position chromosome+ps gives back the sequence from the first allele
                            try:
                                partOfSeq = homozygousity_d[chrnr2, ps]
                            except:
                                logger(f'query: {format((chrnr2, ps))}\tdata: {format(homozygousity_d.keys())} ', "error")
                                partOfSeq = ""

                            # cut current sequence replace
                            sequence = sequence[:patternStart] + "" + sequence[patternEnd:]  # cuts part inbetween
                            # insert new sequence
                            sequence = sequence[:patternStart] + partOfSeq + sequence[patternStart:]  # fills part inbetween

                            # shortTR_allele1 = SimiSTR_bedline(bedfile_entry)
                            shortTR_allele2 = SimiSTR_bedline(bedfile_entry)
                            shortTR_allele2.chromosome = chrnr_w
                            shortTR_allele2.start = patternStart # allele2 start
                            shortTR_allele2.end = patternStart + len(partOfSeq) # allele 2 pattern end new

                            #calculate new offset
                            oldSeqLen = (patternEnd - patternStart)
                            offset_allele2 = len(partOfSeq) - oldSeqLen
                            offset += offset_allele2

                            shortTR_allele2.expansion_length = offset_allele2
                            #entrance_allele2[6] = snvs dont exist more than on allele 1 cause homozygous
                            bedfile_l_output[bedfile_idx] = shortTR_allele2.listBedline()

                        newrecordlen = len(sequence)
                        logger(f'Current sequence length: {newrecordlen}\n', "info")
                        bedfile_idx += 1
                else:
                    logger(f'Check the first column in your assigned input bed file!', "warning")

                record2.seq = sequence
                writer.write_header()
                writer.write_record(record2) 
                bedfile_total.append(bedfile_l_output)
                allele += 1 # next allele for diplod/polyploids

        sWriter.printBedModifications(bedfile_total)




def get_args():
    parser = argparse.ArgumentParser(description="Run SimiSTR to change Expansionlength of STRs.")
    parser.add_argument('-inf', '--input_fasta', type=str,  required=True, help="Path+Name to Fasta File that is template that needs STR expansion changes")
    parser.add_argument('-outf', '--output_fasta', type=str,  required=True, help="Path+Name for newly generated Fasta File with expansion changes")
    parser.add_argument('-ibf', '--input_bedfile', type=str,  required=True, help="Path+Name to Bedfile containing regions of known STRs of given Input Fasta")
    parser.add_argument('-obf', '--output_bedfile', type=str,  required=True, help="Path+Name to Bedfile containing information about applied changes in given STR regions")
    parser.add_argument('-expp', '--expansion_possibility',type=float,  required=True,  help="[0.000-1.000] How many regions should be STR expansion length manipulated")
    parser.add_argument('-dip', '--ploidy', type=int, choices=range(1, 3),  required=True, help="[1-2] Diploid= 2 , Haploid= 1. Multiploid is not yet implemented" )
    parser.add_argument('-snv', '--snv_chance', type=float,  required=True, help="[0.000-1.000] is the chance of a SNV.")
    parser.add_argument('-lid', '--less_indels', type=int,  required=True, help="[int] How much rarer should a insertion/deletion occur than a substitution.")
    parser.add_argument('-ho', '--homozygousity', type=float,  required=True, help="[0.000-1.000] How many regions should be homzygous. The rest will be heterozygous.")
    parser.add_argument('-ma', '--max_add', type=int,  required=False,default=5, help="[int] How many repeats per STR can maximum be added [default: 5]")
    parser.add_argument('-mr', '--max_reduction', type=int, required=False,default=-1, help="[int] How many repeats per STR can maximum be removed. [default: full length]")
    parser.add_argument('-g', '--gangstr_flag', type=int, choices=range(0, 2),default=0,  required=False, help="[0-1] GangstrFile=1, else=0 [default: 0]")
    return(parser.parse_args())


if __name__ == '__main__':
    args = get_args()
    sim = SimiSTR(args)
    sim.main_manipulation()
