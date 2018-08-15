# rlooper

## Synopsis
R-looper is an application designed to enable theoretical analysis of R-loops and R-loop forming regions. The suite contains different applications that share a common code base. Each application in the suite is described below R-looper is compatible only with UNIX based architectures (Linux, OSX, etc), and must be compiled on the user's machine.
 
## Installation
TBD

## Arguments
### Required
The first argument is the path of the input file, which should be a fasta formatted DNA sequence. The header should contain the name of the sequence, coordinates, and which is the template strand. TODO: Make the header parser more flexible. 

As an example: `>HG19_AIRN_PFC53_REVERSE_dna range=PFC53FIXED:1-3908 5'pad=0 3'pad=0 strand=- repeatMasking=none`
This header specifies the gene name (HG19_AIRN_PFC53_REVERSE), the coordinates (chromosome PFC53FIXED, bases 1-3908), and the strand ('-'). 

The second argument is the name that will be used as the base for the output .bed files and genome browser tracks. Specify the name with no file extension, appropriate file extensions will be appended automatically.

### Optional
#### Biophysical Parameters
`--N` followed by an integer specifies the size of the superhelical domain, 

## Usage

## Contributors

## License
