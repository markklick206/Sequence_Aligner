/*
*	Implements the Needleman-Wunsch Algorithm for aligning two sequences
*	Input: input.txt - contains the two sequences to be aligned. Each one on a separate line of the file.
*	Output: output.txt - Stores the sequence alignments
*	Also outputs a score for the alignment.
*/
// Forrest Ireland
//
#define HERE cout << "At line " << __LINE__ << endl;

#include "NW_MultiAlign.h"
#include "Multi_Sequence.h"

using namespace std;

int main() {
	MultiSequence* MS1 = new MultiSequence();
	MultiSequence* MS2 = new MultiSequence();
	MultiSequence* MS3 = new MultiSequence();
	MultiSequence* MS4 = new MultiSequence();

	Sequence S1, S2, S3, S4;
	S1.setSequenceFromTextFile("seq1.txt");
	S2.setSequenceFromTextFile("seq2.txt");
	S3.setSequenceFromTextFile("seq3.txt");
	S4.setSequenceFromTextFile("seq4.txt");

	MS1->setFirstSequence(S1);
	MS2->setFirstSequence(S2);
	MS3->setFirstSequence(S3);
	MS4->setFirstSequence(S4);
	
    NWMultiAlign A1, A2, A3;
	
    A1.SetMultiSequence(MS1, 1);
    A1.SetMultiSequence(MS2, 2);
	
	A2.SetMultiSequence(MS3, 1);
	A2.SetMultiSequence(MS4, 2);
	
	A1.AlignMultiSequences();
	
	A2.AlignMultiSequences();
    
	A3.SetMultiSequence(A1.GetAlignedMultiSequence(), 1);
	A3.SetMultiSequence(A2.GetAlignedMultiSequence(), 2);
	
	A3.AlignMultiSequences();
    
	A1.WriteAlignedMultiSequenceToFile("A1.txt");
	A2.WriteAlignedMultiSequenceToFile("A2.txt");
	A3.WriteAlignedMultiSequenceToFile("A3.txt");
	
	return 0;
}