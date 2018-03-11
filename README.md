# Computational Biology RNA-Protein Interactions

## Introduction
The objective of this repository is to develop a program to study RNA-Protein Interactions
The entire code base is written in Python 2 with UNIX commands executed on MacOSx.

RNA and Protein interactions are of interest because RNA is part of the central dogma of biology with respect to DNA,RNA, protein and transcription/translation.

RNA and Proteins depend upon hydrogen bonding and secondary forces.
This code base utilizes PDB files and runs it through HBplus to calculate physical/chemical bonding. 
The secondary structure of RNA is constructed utilizing DSSR.


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
### To run HBplus

```
hbplus.exe [options] [cleaned filename] [uncleaned filename]
```

### To run DSSR
```
x3dna -dssr --input=1msy.pdb --output=1msy.out
```
