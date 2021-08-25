# STRsimulator
Manipulation of a Referencefile in order to simulate STR!!

The simulator takes a haploid file as reference(.fasta) and a region file (.bed) containing information about known STR-regions as input.
All of the supplied regions can be modified in 
  - expansion (% of regions that will randomly be positive or negative expanded [0.00-1.00]), 
  - mutation (% chance for a base to be substituted [0.00-1.00]), 
  - number of indels (X times less likely than chance for mutation to insert or delete a base [0.00-1.00]). 
Further can the simulation file (.fasta) be created as 
  - haploid [h] or 
  - diploid  [d]. If diploid is chosen, 
    - the percentage  of regions that should get homozygous can be set [0.00-1.00].

## Installation
- Implementation and testing occurred in python=3.8.
- Python modules used: os, sys, random, copy and Bio (SeqIO)

## Usage
-      python  ./STRsimulator_v7.py \
>    - “ref.fasta” \ ________________________# **obligation** (Reference Fasta)
>    - “name_for_out.fasta” \ ____________# **obligation** (Result Fata)
>    - “in_regionfile.bed” \ _______________# **obligation** (! Important: need to be sorted !) (Bedfile)
>    - “name_for_out_regionfile.bed” \ ___# **obligation** (Result Bedfile – contains info about changes)
>    - [0.00-1.00] \ 	________________________# **obligation** (as in 0-100% chance a STR expansion gets modified)
>    - "[d|h]" \ ____________________________# h = haploid | d = diploid
>    - [0.00-1.00] \ _______________________# chance for each base in region to be mutated (as in 0-100%)
>    - [0-100 or more] \ __________________# Indels times more unlikely than substitution (float, everything possible)
>    - [0.00-1.00] \ _______________________# percentage of homozygous regions (as in 0-100%)

#### Important: Input Bedfile, needs to be sorted in region-start and region-end ( f.e. bedtools sort -i myfile.bed > myfile.sorted.bed )
- Current Version expects fasta and bedfiles with Chromosome-names without "chr". 
There is a folder attached with other versions, that might contain readers dealing with "chr" naming. [OldVersions](https://github.com/DamarisLa/STRsimulator/tree/main/OldVersions) 
- Bedfile should not have a header
- Further Input fasta should have the columns:
- ChromosomeNr; RegionStart; RegionEnd; MotifLength; Motif
- f.e. :
    -       22	20348371	20348390	4	TTTA
    -       22	20353556	20353575	4	AAAC
    -       22	20354654	20354669	4	ATTT
    -       22	20374713	20374727	3	TTG
    
#### Important: [GangSTR Bedfile](https://github.com/gymreklab/GangSTR#tr-regions---regions) 
Currently the main function has a calculation startposition-1 when reading and start position +1 when working with GangSTR bedfiles. 
line nr 248 #-1 , when not working with gangstr-bedfiles. 
line nr 248 and 360 -1  , when working with gangstr-bedfiles.
line nr 269 and 360 # +1, when not working with gangstr-bedfiles. 
!It is important to be aware of the meaning of the start-position in the bedfile one uses, and adapt the code if nessesary!

### Future Improvements
- Implementation is linear. Fasta reader is a bottle neck. Future improvement will be threading before fasta-reader, that chromosome will be run parallel, than after each other. 
- Improvements in structure and function delegation
- Improvements in parameter handling
- Parameter to submit if a bedfile-start-position is +1 oder not, and automatically deals with it
- Parameter to submit min/max possible change of expansion of an STR
- Parameters in Config-File
- Implementation will be changed. Position search like until know, however, if region cannot be found in certain window, or only 1 repeat can be found in that windoe, try to do a smith-waterman for a slightly increased region. (?) (it is a vage idea to solve the issue with overlapping regions)

#### [BedFileParser](https://github.com/DamarisLa/STRsimulator/tree/main/BedFileParser) 
Contains Parser that reformat Bedfiles into Tool-specific Files needed by several STR-Tools compared in 
STR Benchmarkpaper [here will be a link to our paper](https://)







 
 



