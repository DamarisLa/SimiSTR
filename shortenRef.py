from Bio import SeqIO

oldFastaFile = "..\\FilteredViewed\\hs37d5.fa"
newFastaFile = "..\\FilteredViewed\\hs37d5.chr22.fa"

with open(newFastaFile, 'w') as outFastaFile:
    with open(oldFastaFile, 'r') as inFastaFile:
        for record in SeqIO.parse(oldFastaFile, "fasta"):
            if record.id == '22':
                writer = SeqIO.FastaIO.FastaWriter(outFastaFile)
                writer.write_header()
                writer.write_record(record)

#chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
#22	    51304566	    2876892038	            60	    61
