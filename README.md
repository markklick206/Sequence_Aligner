Sequence_Aligner
================

Globally aligns sequences

NW_ALIGN v0.1.1

Forrest Ireland
Nov. 29, 2014

Uses the Needleman-Wunsch algorithm to globally align to sequences.

REQUIRED FILES
*******************************
input.txt - REQUIRED - Contains the two imput sequences to be aligned by the program. The first sequence is on line 1, then there needs to be a blank line, followed by a second sequence. Currently the program accepts a full alphabet of characters.


OUTPUT FILES
*******************************
matrices.txt - PRODUCED - Used for debug purposes, displays the alignment scoring matrix and the traceback matrix used to produce the alignment.

AlignedSequences.txt - PRODUCED - An output file that displays the aligned sequences with gaps introduced.


CHANGE LOG
*******************************
v0.1
	- First release
	- Aligns two sequences properly
	- Read an input file containing two sequences
	- Writes aligned sequences to file
	- Writes traceback and score matrices to file

v0.1.1
	- Fixed issue with release/debug builds
	- Statically links some c++ dll crap so it should work on all windows machines. Had issues with computers not having the c++ redistributables installed
