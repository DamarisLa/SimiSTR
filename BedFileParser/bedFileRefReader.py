from Bio import SeqIO
import sys
import os

def getBedFile(oldBedFile):
    bedfile_l = list()
    with open(oldBedFile, 'r') as inBedFile:
        for line in inBedFile:
            splitline = line.split("\t")
            if len(splitline)>3:
                bedfile_l.append(splitline) #chr    from Pos    to Pos      lenMotif    motif
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
                patternStart = int(shortTR[1])-1
                patternEnd = int(shortTR[2])
                patternLen = int(shortTR[3])
                pattern = shortTR[4].strip()
                #patternStartold = int(shortTRold[1])-1
                #patternEndold = int(shortTRold[2])
                #patternLenold = int(shortTRold[3])
                #patternold = shortTRold[4].strip()


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


newBedFile = "..\\FilteredViewed\\Grch38\\GangstrBedfiles\\randomSubset.hg38_ver13.sorted_noXY.bed"
oldBedFile = "..\\FilteredViewed\\Grch38\\GangstrBedfiles\\randomSubset.hg38_ver13.sorted_noXY.bed"
newFastaFile = "..\\FilteredViewed\\Grch38\\grch38_minchrs_rnamed.fa"
oldFastaFile = "..\\FilteredViewed\\Grch38\\grch38_minchrs_rnamed.fa"






main_reader(newFastaFile,newBedFile,oldBedFile)



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
