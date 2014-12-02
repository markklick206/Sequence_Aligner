/*
*	Implements the Needleman-Wunsch Algorithm for aligning two sequences
*	Input: input.txt - contains the two sequences to be aligned. Each one on a separate line of the file.
*	Output: output.txt - Stores the sequence alignments
*	Also outputs a score for the alignment.
*/
// Forrest Ireland
//
#define HERE cout << "At line " << __LINE__ << endl;

#include "NW_Align.h"

using namespace std;

int main() {
	NWAlignment* NWA = new NWAlignment();

	if (NWA->ReadInputSequenceFromFASTA("sequence1.fasta", 1) && NWA->ReadInputSequenceFromFASTA("sequence.fasta", 2))
		NWA->AlignSequences();
    NWA->WriteAlignedSequencesToFile("AlignedSequences.txt");
	//NWA->PrintSequenceToConsole(1, false);
	//NWA->PrintSequenceToConsole(2, false);
	//NWA->PrintSequenceToConsole(1, true);
	//NWA->PrintSequenceToConsole(2, true);
    NWA->CloseNWAlign();
    
	return 0;
}