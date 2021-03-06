from Bio import SeqIO
import sys
import os
import re

def getBedFile(oldBedFile):
    bedfile_l = list()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            #print(line)
            splitline = line.split("\t")
            #print(splitline)
            if len(splitline)>3:
                txt = re.search("chr(\d*)", splitline[1])
                chr = -1
                if txt is not None:
                    chr = txt[1]
                if chr != -1 and chr != '':
                    #if int(chr) > 10:
                    #    print(chr)
                    important_fields = [chr,splitline[2],splitline[3],splitline[7],splitline[-1].strip()]
                    #print(important_fields)
                    bedfile_l.append(important_fields) #chr    from Pos    to Pos      lenMotif    motif
    return bedfile_l


def main_reader(newFastaFile, newBedFile, oldBedfile):
    bedfile_l = getBedFile(newBedFile)
    print("Successful first newBedFile read")
    #bedfile_l_old =getBedFile(oldBedfile)
    #print("Successful first oldBedFile read")
    # copy that list or dict and always safe the changes of the offset, to memorize the new coordinates

    with open(newFastaFile, 'r') as inFastaFile:
        for record in SeqIO.parse(inFastaFile, "fasta"):
            # for i in range(0,len(record.seq)):
            sequence= record.seq
            recordLen = len(sequence)
            print(recordLen)
            # writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
            count = 0
            for i in range(0,len(bedfile_l)):
                shortTR = bedfile_l[i]
                #shortTRold = bedfile_l_old[i%(len(bedfile_l_old))]
                chrnr = shortTR[0]
                patternStart = int(shortTR[1])
                patternEnd = int(shortTR[2])
                patternLen = int(shortTR[3])
                pattern = shortTR[4].strip()
                # patternStartold = int(shortTRold[1])
                # patternEndold = int(shortTRold[2])
                # patternLenold = int(shortTRold[3])
                # patternold = shortTRold[4].strip()


                if record.id == chrnr:
                    seq_len = len(sequence)
                    partOfSeq = sequence[patternStart:patternEnd]
                    partOfSeq2 = sequence[patternStart-6:patternEnd+6]

                    #print("old: ", shortTRold, " ", (patternEndold - patternStartold) / patternLenold)

                    print(partOfSeq2)
                    print("     ",partOfSeq)
                    print("new: ",shortTR," ", (patternEnd-patternStart)/patternLen,"\n")
                    patterntotal = int((patternEnd-patternStart)/patternLen)*pattern
                    if patterntotal != partOfSeq and not patternStart == patternEnd:
                        count +=1


            print("Unequal pairs: ", count)


main_reader("..\\FilteredViewed\\hs37d5.fa","..\\FilteredViewed\\simplerepeats_37_min3bp_max10bp.bed","..\\FilteredViewed\\simplerepeats_37_min3bp_max10bp.bed")



# if len(sys.argv) < 4:
#     print("Please give a fastafile, the name and dir where the new dir has to be, the old bedfile, "
#           "the  name and dir where the new bedfile should be, "
#           "and a number between 0 and 1 indicating the change with which the original ref file")
# else:
#     for i in sys.argv:
#         print(i)
#     if os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[3]):
#         print("Enter main reader.")
#         main_reader(sys.argv[1],sys.argv[3],sys.argv[2])
#         #[0]./bedFileRefReader.py [1]../../hs37d5.chr22.new1.fa  [2]../bedfiles/hs37_ver8.chr22.bed [3]hs37_ver8.chr22.new1.bed
#     else:
#         if not os.path.isfile(sys.argv[1]):
#             print(sys.argv[1], " is not a file")
#         elif not os.path.isfile(sys.argv[2]):
#             print(sys.argv[2], " is not a file")
#         elif not os.path.isfile(sys.argv[3]):
#             print(sys.argv[3], " is not a file")
