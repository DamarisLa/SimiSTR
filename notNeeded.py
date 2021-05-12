# get all coordnates and offsets of the ref-file to know to which position need to be jumped
# as according to pos from bedfile the exact ref-position needs to be calculated
# here a dict could be more helpfull => further this needs to get a function that gets called
# by the main function



def getFai(fastaFai, fai_info_d):
    with open(fastaFai, 'r') as inFaiFile:
        for line in inFaiFile:
            splittedline = line.strip().split("\t")
            fai_info_d[(splittedline[0], splittedline[1])] = splittedline
            #chr    length of Chr   offset(CoordinateStart) bases each line (60without 61with \n)
    return fai_info_d

faispec_info_d = dict()
fastaFai = "..\FilteredViewed\\hs37d5.chr22.fa.fai"
faispec_info_d = getFai(fastaFai,faispec_info_d)
# faigeneral_info_d =dict()
# fastaFaiAll = "..\FilteredViewed\\hs37d5.fa.fai"
# faigeneral_info_d = getFai(fastaFaiAll,faigeneral_info_d)
#
# for id in faigeneral_info_d.keys():
#     entrance = faigeneral_info_d[id]
#     start = entrance[2]
#
#     print("Start: ", )
#     print("hello")