# R-looper 2

## Synopsis
R-looper is an application designed to enable theoretical analysis of R-loops and R-loop forming regions. The suite contains different applications that share a common code base. Each application in the suite is described below R-looper is compatible only with UNIX based architectures (Linux, OSX, etc), and must be compiled on the user's machine.
 
## Installation
Prerequisites:
- Unix based system (Ubuntu, OSX, etc)
- Install git if required: https://git-scm.com/book/en/v2/Getting-Started-Installing-Git. Use git to pull the repository. 
- OR unzip the downloaded repository

To install follow these steps:
1. Open a terminal, and use the `cd` command to navigate into the directory where the software has been downloaded and/or unzipped. 
2. Build the software using `make all`
4. If the build was successful (no error messages), a new directory should have been created within your rlooper directory called `bin`. The executable `rlooper` should be within that directory, and can now be moved anywhere on your computer and executed from the terminal with the arguments described below. 

## Arguments
### Required
The **first argument** is the path of the input file, which should be a fasta formatted DNA sequence. The header should contain the name of the sequence, coordinates, and which is the template strand. The sequence will be automatically handled correctly based on the provided orientation of the template strand. TODO: Make the header parser more flexible. 

As an example: `>HG19_AIRN_PFC53_REVERSE_dna range=PFC53FIXED:1-3908 5'pad=0 3'pad=0 strand=- repeatMasking=none`
This header specifies the gene name (HG19_AIRN_PFC53_REVERSE), the coordinates (chromosome PFC53FIXED, bases 1-3908), and the strand ('-'). 

When downloading a sequence from the genome browser, make sure to download the sequence with no repeat masking.  

The **second argument** is the name that will be used as the base for the output .bed files and genome browser tracks. Specify the name with no file extension, appropriate file extensions will be appended automatically.

### Optional
#### Biophysical Parameters
`--N` followed by an integer specifies the size of the superhelical domain in base-pairs. This is the number of bases over which DNA is not in its relaxed state. Alternatively, `--N auto` can be used to set the size of the superhelical domain to the size of the provided sequence. The default if unspecified is 1500bp, (Kouzine and Levins, Myc C FUSE work)

`--a` followed by an floating point value specifies the energetic disfavorability in kcal/mol associated with the two junctions at each end of the R-loop structure. This quantity is assumed to be length independent. The default value if unspecified is 10 kcal/mol (5 per junction), consistent with the experimentally measured value for other types of superhelically driven DNA structural transition. (Levins, Benham)

`--minlength` followed by an integer specifies the size of R-loops in base-pairs below which should be excluded from the biophysical ensemble. The default if unspecified is 2bp, which is the full ensemble. 

`--reverse` reverses the correct orientation of the provided input sequence.

`--complement` complements the correct orientation of the provided input sequence. 

`--invert` transforms the provided input sequence to its reverse complement. Useful when seeking to examine the non-template  strand of the provided sequence.

`--homopolymer` followed by a floating point value overrides the base-pairing energetics with that value in kcal/mol. In this case, the content of the input sequence does not matter, as all dinucleotide pairs will be treated as having the same provided energy.

`--bedfile` indicates to the software that peaks corresponding with regions of high r-loop favorability should be called and outputted to a .bed formatted file. 

`--sensitivity` is used only when `--bedfile` is speficied and specifies the sensitivity of the threshold-based peak caller. 10 seems to be a good default value when calling peaks. 

`--residuals` computes and outputs residual quantities to stdout. These quantities describe the expected amount of remaining superhlicity in a molecule consisting of the provided sequence, devided between twist and writhe. Can be easily used to determine how much relaxation would be expected on this molecule if single R-loops were allowed to form. 

`--circular` treats the provided of sequence as circular by adding an additional set of r-loop structures to the ensemble that would span the boundry between the beginning and end of the sequence if the sequence is circular. Useful in a case where a small circular molecule is being transcribed and a region of energetic favorability lies at the end of the provided sequence.

`--top` followed by some integer n, indicates that information about the top n most favorable single R-loop structures should be outputted to stdout. Useful to get a sense of how the most favorabl structures over the sequence are distributed in position and energy.

`--dump` dumps the full statistical ensemble (every possible single R-loop structure) to a file for analysis. It provides the start location, stop location, energy, and probability of each structure, and the structures are sorted by decreasing probability. 

`--localaverageenergy` adds an aditional signal to the output as a .wig file. This signal is a measure of local average energy, and can supplement the base pair involvement probability. Disabled by default because it is a computationally expensive signal to compute.

## Example Usage
Suppose that you wanted to run R-looper on a FASTA formatted sequence file, "pfc53.fa", and output the results to a set of files using `pfc53_output` as the base name. Suppose also that you want to assume negative 7% superhelical density over the 1500nt (default superhelical domain size). You would use the following:
`rlooper pFC53_full.fa pfc53_output --sigma -0.07`

Now suppose you wanted to alter the cost of the DNA:RNA DNA:DNA junctions at each end of the R-loop from the default value (10 kcal/mol) to a higher value (20kcal/mol) and use a superhelical domain that is scaled to the length of the sequence to hold the superhelical density constant at -7%. Suppose you also wanted to compute the more costly "local average energy" signal over the sequence. 
`rlooper pFC53_full.fa pfc53_output --sigma -0.07 --a 20 --N auto --localaverageenergy`

## Contributors


This software repository has been created and maintained by Robert Stolz (rstolzATucdavis.edu). 

Development of this software was funded in part by NIH grant XXX and NSF grants XXX and XXX. 

## License
All rights reserved 2018
