
# ![grafik](https://user-images.githubusercontent.com/62024017/137172832-8453c9de-573a-4ff1-8ad2-8cf056a8cf01.png)

Manipulation of a Referencefile in order to simulate STR!!

The simulator takes a haploid file as reference(.fasta) and a region file (.bed) containing information about known STR-regions as input.
All of the supplied regions can be modified in 
  - expansion (% of regions that will randomly be positive or negative expanded [0.00-1.00]), 
  - mutation (SNV) (% chance for a base to be substituted [0.00-1.00]), 
  - number of indels (single base insertion/ deletions) (X times less likely than chance for mutation to insert or delete a base [0.00-1.00]). 
Further can the simulation file (.fasta) be created as 
  - haploid [h] or 
  - diploid  [d]. If diploid is chosen, 
    - the percentage of regions that should be homozygous can be set [0.00-1.00].

Expansion is randomly adding or removing repeats of a pattern in a STR- regions. 
- It removes  by default [0 - "all available"]  repeats of the pattern 
- It adds     by default [0 - 5]                repeats of the pattern
- The parameters -mr (--max_reduction) and -ma (--max_add) allow to adapt the maximum change. However its never possible to remove more than the whole STRregion. 

The simulator works on assembled genomes, as well as on only one or more assembled chromosomes, if the bed-file contains such entrances likewise (anything else could run errorfree, but will not manipulate anything, as manipulations only occur in the known regions). 


## Installation
- Implementation and testing occurred in python=3.8.
- Python modules used: os, sys, random, copy and Bio (SeqIO)

## Usage
```
usage: SimiSTR_v8.py [-h] -inf INPUT_FASTA -outf OUTPUT_FASTA -ibf INPUT_BEDFILE -obf OUTPUT_BEDFILE  
                     -expp EXPANSION_POSSIBILITY -dip {1,2} DIPLOIDITY -snv SNV_CHANCE -lid LESS_INDELS 
                     -ho HOMOZYGOUSITY [-ma MAX_ADD] [-mr MAX_REDUCTION] [-g {0,1}]

Run SimiSTR to change Expansionlength of STRs.

optional arguments:
  -h, --help            show this help message and exit
  -inf INPUT_FASTA, --input_fasta INPUT_FASTA
                        Path+Name to Fasta File that is template that needs
                        STR expansion changes
  -outf OUTPUT_FASTA, --output_fasta OUTPUT_FASTA
                        Path+Name for newly generated Fasta File with
                        expansion changes
  -ibf INPUT_BEDFILE, --input_bedfile INPUT_BEDFILE
                        Path+Name to Bedfile containing regions of known STRs
                        of given Input Fasta
  -obf OUTPUT_BEDFILE, --output_bedfile OUTPUT_BEDFILE
                        Path+Name to Bedfile containing information about
                        applied changes in given STR regions
  -expp EXPANSION_POSSIBILITY, --expansion_possibility EXPANSION_POSSIBILITY
                        [0.000-1.000] How many regions should be STR expansion
                        length manipulated
  -dip {1,2}, --diploidity {1,2}
                        [1-2] Diploid= 2 , Haploid= 1. Multiploid is not yet
                        implemented
  -snv SNV_CHANCE, --snv_chance SNV_CHANCE
                        [0.000-1.000] is the chance of a SNV.
  -lid LESS_INDELS, --less_indels LESS_INDELS
                        [int] How much rarer should a insertion/deletion occur
                        than a substitution.
  -ho HOMOZYGOUSITY, --homozygousity HOMOZYGOUSITY
                        [0.000-1.000] How many regions should be homzygous.
                        The rest will be heterozygous.
  -ma MAX_ADD, --max_add MAX_ADD
                        [int] How many repeats per STR can maximum be added
                        [default: 5]
  -mr MAX_REDUCTION, --max_reduction MAX_REDUCTION
                        [int] How many repeats per STR can maximum be removed.
                        [default: full length]
  -g {0,1}, --gangstr_flag {0,1}
                        [0-1] GangstrFile=1, else=0 [default: 0]
```                        
                        

#### Important: Input Bedfile, needs to be sorted in region-start and region-end ( f.e. bedtools sort -i myfile.bed > myfile.sorted.bed )

- Further Input fasta should have the columns:
- ChromosomeNr; RegionStart; RegionEnd; MotifLength; Motif
- f.e. :
    -       22	20348371	20348390	4	TTTA
    -       22	20353556	20353575	4	AAAC
    -       22	20354654	20354669	4	ATTT
    -       22	20374713	20374727	3	TTG
    
#### Important: [GangSTR Bedfile](https://github.com/gymreklab/GangSTR#tr-regions---regions) 
The main function need a calculation startposition-1 when reading and start position +1 when working with GangSTR bedfiles. 
Please set the GangSTR flag as 1 in case a gangstr-bedfile is used, else it is assumed that the startpositions in the bedfile are 0-based.
(The correct start usually will be found even if the 0-flag is not set. However it increases chances for errors)
!It is important to be aware of the meaning of the start-position in the bedfile one uses, and adapt the code if nessesary!

### Future Improvements
- Implementation is linear. Fasta reader is a bottle neck. Future improvement will be threading before fasta-reader, that chromosome will be run parallel, than after each other. 




[BedFileParser](https://github.com/DamarisLa/STRsimulator/tree/main/BedFileParser) 
Contains Parser that reformat Bedfiles into Tool-specific Files needed by several STR-Tools compared in 










 
 



