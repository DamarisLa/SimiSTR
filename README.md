# STRsimulator
Manipulation of a Referencefile in order to simulate STR!!

The simulator takes a haploid file as reference(.fasta) and a region file (.bed) containing information about known STR-regions as input.
All of supplied regions can be modified in 
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

#### Important: Input Bedfile, needs to be sorted in regionstart and regionend ( f.e. bedtools sort -i myfile.bed > myfile.sorted.bed )
- Current Version expects fasta and bedfiles with Chromosome-names without "chr". 
There is a folder attached with other versions, that might contain readers dealing with "chr" naming. [OldVersions](https://github.com/DamarisLa/STRsimulator/tree/main/OldVersions) 

- Further Input fasta should have the columns:
- ChromosomeNr; RegionStart; RegionEnd; MotifLength; Motif
- f.e. :
    -       22	20348371	20348390	4	TTTA
    -       22	20353556	20353575	4	AAAC
    -       22	20354654	20354669	4	ATTT
    -       22	20374713	20374727	3	TTG
    

### Future Improvements
- Implementation is linear. Fasta reader is a bottle neck. Future improvement will be threading before fasta-reader, that chromosome will be run parallel, than after each other. 
- Improvements in structure and function delegation
- Improvements in parameter handling


#### [BedFileParser](https://github.com/DamarisLa/STRsimulator/tree/main/BedFileParser) 
Contains Parser that reformat Bedfiles into Tool-specific Files needed by several STR-Tools compared in 
STR Benchmarkpaper [here will be a link to our paper](https://)







 
 



