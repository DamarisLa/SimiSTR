from bedFileRefReader import getBedFile
import json

bedfile = "..\\FilteredViewed\\Grch38\\GangstrBedfiles\\randomSubset.hg38_ver13.sorted_noXY.bed"
jsonfile = "..\\FilteredViewed\\Grch38\\GangstrBedfiles\\randomSubset.hg38_ver13.sorted_noXY.json"


def catalogJsonWriter (bedfile, jsonfile):
    bf = getBedFile(bedfile)
    bf_l = list()
    bf_d = dict()
    for line in bf:
        print(line)
        referenceRegion = str(line[0])+":"+str(line[1])+"-"+str(line[2])
        str_string = '{"LocusID": "'+ str(line[-1]).strip() + '", "LocusStructure": "('+ line[-1].strip() + ')*", '\
                     + '"ReferenceRegion": "' + referenceRegion + '", '\
                     +'"VariantType": "Repeat"}'
        bf_l.append(json.loads(str_string))
    with open(jsonfile, 'w') as json_file:
        #for i in bf_l:
        #    json.dump(i, json_file)
        json.dump(bf_l,json_file)




catalogJsonWriter(bedfile,jsonfile)