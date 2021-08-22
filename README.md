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
>- python  ./STRsimulator_v7.py \
>    - “ref.fasta” \ __________________# *obligation* (Reference Fasta)
>    - “name_for_out.fasta” \ _________# *obligation* (Result Fata)
>    - “in_regionfile.bed” \ ___________# *obligation* (! Important: need to be sorted !) (Bedfile)
>    - “name_for_out_regionfile.bed” \ _# *obligation* (Result Bedfile – contains info about changes)
>    - [0.00-1.00] \ 	_________________# *obligation* (as in 0-100% chance a STR expansion gets modified)
>    - "[d|h]" \ _____________________# h = haploid | d = diploid
>    - [0.00-1.00] \ _________________# chance for each base in region to be mutated (as in 0-100%)
>    - [0-100 or more] \ _____________# Indels times more unlikely than substitution (float, everything possible)
>    - [0.00-1.00] \ _________________# percentage of homozygous regions (as in 0-100%)

#### Important: Input Bedfile, needs to be sorted ( bedtools sort -i myfile.bed > myfile.sorted.bed )

### Future Improvements:
Implementation is linear. Fasta reader is a bottle neck. Future improvement will be threading before fasta-reader, that chromosome will be run parallel, than after each other. 





 
 



