# Identify the N-terminal Domain Type of a VSG Sequence

The N-terminal domain of a VSG sequence can be classified as one of the following types/subtypes:

Type A:
* A1
* A2
* A3

Type B:
* B1
* B2

## Usage

```
python find_VSG_Nterm.py file path
```

For more details about arguments see [Input](#input).

## Required Software
* [Biopython](https://anaconda.org/anaconda/biopython) (written using version 1.68)
* [pandas](https://anaconda.org/anaconda/pandas) (written using version 0.24.2)
* [HMMER](http://hmmer.org) (written using version 3.1b2)
	- HMMER must be installed in PATH.

## Input
* (1) FASTA file containing one or more VSG protein sequences. - (`file`)
* (1) Path to directory containing HMM profiles to run sequences against. - (`path`)

### Required files for HMMER hmmscan
* VSG N-terminal TypeA hmm profile ([Wickstead et al.](https://www.sciencedirect.com/science/article/pii/S0166685114000772)) - (`VSG-N-TypeA.hmm`)
* VSG N-terminal TypeB hmm profile([Wickstead et al.](https://www.sciencedirect.com/science/article/pii/S0166685114000772)) - (`VSG-N-TypeB.hmm`)
* VSG N-terminal TypeSubtype hmm profile - (`VSG-N-TypeSubtype.hmm`)

HMM profile files are located in [HMMprofiles](HMMprofiles).  
The files must be in a sinlge directory (and the path to this directory provided as `path`).  
HMM profile files must first be compressed and indexed using HMMER hmmpress.

For example:
```
hmmpress profile.hmm
```

## Output
* (1) Text file containing progress of the executed script.
* (1) FASTA file containing the most probable N-terminal domain(s) of the input VSG(s) and the corresponding typesubtype.
* (1) CSV file containing a summary of the identified N-terminal domain(s).
* (1) Directory containing:
	- (1) FASTA file containing trimmed variants of the input VSG(s).
	- (2) HMMER hmmscan table output files containing data on each trimmed variant of the input VSG(s) (TypeA and TypeB).
	- (2) HMMER hmmscan standard output files containing annotation and statistical information about each trimmed variant of the input VSG(s) (TypeA and TypeB).
	- (1) FASTA file containing the N-terminal domain sequence(s) and the corresponding type (w/o subtype).
	- (1) HMMER hmmscan table output file containing data on the N-terminal domain sequence(s) of the input VSG(s) (TypeA1-3/TypeB1-2).
	- (1) HMMER hmmscan standard output file containing annotation and statistical information about the N-terminal domain(s) of the input VSG(s) (TypeA1-3/TypeB1-2).
