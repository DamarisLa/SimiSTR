from Bio import SeqIO


oldFastaFile = "..\\FilteredViewed\\hs37d5.fa"
newFastaFile = "..\\FilteredViewed\\hs37d5.chr22.fa"


def getOneChromosomeFasta(chromosmeNr_str):
    with open(newFastaFile, 'w') as outFastaFile:
        with open(oldFastaFile, 'r') as inFastaFile:
            for record in SeqIO.parse(oldFastaFile, "fasta"):
                if record.id == chromosmeNr_str:
                    writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
                    writer.write_header()
                    writer.write_record(record)

#chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
#22	    51304566	    2876892038	            60	    61

#mmap und pigz => reader vielleicht schneller

# new fasta => new fai


#samtools faidx ref.fasta
#test faidx
# def checkOffsetFromFai(number):
#     f = open(newFastaFile, "r")
#     count = 0
#     for y in f.readlines():
#         print(y)
#     for x in f.readline():
#         #for y in x:
#         if count == (number-1) or (count == number) or count == (number+1):
#             print(x)
#         count +=1
#         print(x, count)
#
#
# checkOffsetFromFai(54)