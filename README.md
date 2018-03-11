# Computational Biology RNA-Protein Interactions

## Introduction
The objective of this repository is to develop a computational program to study RNA-Protein Interactions
The entire code base is designed to be a standalone container and is written in Python 2 with UNIX commands executed on MacOSx and executables (HBplus and DSSR).

RNA and Protein interactions are of interest because RNA is part of the central dogma of biology with respect to DNA,RNA, protein and transcription/translation.

RNA and Proteins depend upon hydrogen bonding and secondary forces.

## Code Architecture design

This code base utilizes PDB files and hydogen and second force calculations are determined through HBplus to calculate physical/chemical bonding. 
The secondary structure of RNA is constructed utilizing DSSR.

### Computational Biology and Code Design
This code base is designed to be a self contained application where each individual step outlined can be ran as single a operation:

- Create a file called structures of interest (proteins interested i.e. 1MNB)
- Download the pdb files from the Protein Database Bank (PDB).
- Calculate hydrogen bonding of the RNA-Protein complex utilizing HBplus based on each pdb file (3.9A).
- Determine RNA secondary structures from DSSR (i.e. stems, hairpines)
- Bin the nucleotide base involed in bonding to it's corresponding secondary structure.

### Source of information
The rna-protein structure pdb files were retrieved from the Protein Database Bank (PDB)
http://www.rcsb.org/

### Dependencies
*DSSR http://x3dna.org/

*HBplus  http://www.ebi.ac.uk/thornton-srv/software/HBPLUS/

### To run this program

```
python master.py <command>
python master.py help
```

### Commands available
```
$ python master.py help

Welcome to the Computatinal Biology RNA-Protein Interaction help section


The following commands are available:


dssrcli
fileretriever
hbplushbcli
hbplusvdwcli
dssrprocessedreader
hbplusprocessedreader
```

### To run HBplus

```
hbplus.exe [options] [cleaned filename] [uncleaned filename]
```

### To run DSSR
```
x3dna -dssr --input=1msy.pdb --output=1msy.out
```
