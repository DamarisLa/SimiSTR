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

Installation
- Implementation and testing occurred in python=3.8.
- Modules used: os, sys, random, copy and Bio (SeqIO) 



 
 



