from Bio import SeqIO
import re

def splitHaploid(fastaFile):

    with open(fastaFile, 'r') as inFastaFile:
        fastaName = fastaFile[0:-3]
        fastaFile_h1 = fastaName + "_hap01.fa"
        fastaFile_h2 = fastaName + "_hap02.fa"
        with open(fastaFile_h1, 'w') as fh1:
            writer_fh1 = SeqIO.FastaIO.FastaWriter(fastaFile_h1)
            with open(fastaFile_h2, 'w') as fh2:
                writer_fh2 = SeqIO.FastaIO.FastaWriter(fastaFile_h2)


                for record in SeqIO.parse(inFastaFile, "fasta"):
                    id = record.id
                    print(id)
                    reg1 = re.search('\d*\_(\d)', id)
                    if reg1 != None:
                        hap = reg1[1]
                        if hap=='1':
                            writer_fh1.write_record(record)
                        else :
                            writer_fh2.write_record(record)




newFastaFile = "..\\FilteredViewed\\Grch38\\grch38.rand_adapt.fa"
splitHaploid(newFastaFile)

