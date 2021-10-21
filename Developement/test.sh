import SimiSTR_v8


inputBedFile="..\\..\\FilteredViewed\\Grch38\\GangstrBedfiles\\randomSubset.hg38_ver13.sorted_noXY.bed"
inputFastaFile="..\\..\\FilteredViewed\\Grch38\\grch38_minchrs_rnamed.fa"
outputBedFile="..\\..\\FilteredViewed\\Grch38\\grch38.adapt.bed"
outputFasta="..\\..\FilteredViewed\\Grch38\\grch38.rand_adapt.fa"

SimiSTR_v8 -of ${outputFasta} -if ${inputFastaFile} -obf ${outputBedFile} -ibf ${inputBedFile}  -expp 1.0 -dip 2 -snv 0.5 -lid 10 -ho 0.5