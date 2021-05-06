import copy
from Bio import SeqIO

oldBedFile = "..\FilteredViewed\\hs37_ver8.20short.bed"
newBedFile = "..\FilteredViewed\\hs37_ver8.20short.adapt.bed"
newFastaFile = "..\FilteredViewed\\hs37d5.rand_adapt.fa"    #ref-file before STR simulations
oldFastaFile = "..\FilteredViewed\\hs37d5.chr22.fa"         #ref-file after STR simulations
fastaFai = "..\FilteredView\\"


#read the original bedfile. safe all coordinates maybe into a dict
# ?? whats the key ?? the coordinate? or should it be a list? => as iteration from first to last. List
# makes more sense than dict
bedfile_l = list()
with open(oldBedFile, 'r') as inBedFile:
    for line in inBedFile:
        print("hello")


#have an offset, that tells how much all following coordinates will be later or earlier
offset = 0 #original in the beginning

#copy that list or dict and always safe the changes of the offset, to memorize the new coordinates
bedfile_l_copy = copy.deepcopy(bedfile_l)


# get all coordnates and offsets of the ref-file to know to which position need to be jumped
# as according to pos from bedfile the exact ref-position needs to be calculated
# here a dict could be more helpfull => further this needs to get a function that gets called
# by the main function
fai_info_d = dict()
with open(fastaFai, 'r') as inFaiFile:
    for line in inFaiFile:
        print("hello")



#main function manipulate Fasta
with open(newFastaFile, 'w') as outFastaFile:
    with open(oldFastaFile, 'r') as inFastaFile:
        # maybe needed later
        #for record in SeqIO.parse(oldFastaFile, "fasta"):
        #    writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
        #    writer.write_header()
        #    writer.write_record(record)
        # where is the first STR
        for str in bedfile_l:
            chrnr = str[0]
            chrpos = str[1]
            strPattern = str[2]



#write out the new coordinates. The adapted bedfile. To find the regions and to
with open(newBedFile, 'w') as outBedfile:
    print("Hello")








# fai format:
#chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
#22	    51304566	    2876892038	            60	    61
#when reading fasta

#bedfile format
#chr    from Pos    to Pos      lenMotif    motif
#22	    16052168	16052199	4	        AAAC    ==> 68:AAACAAACAAACAAACAAACAAACAAACAAAC:99