/*
*	Implements the Needleman-Wunsch Algorithm for aligning two sequences
*	Input: input.txt - contains the two sequences to be aligned. Each one on a separate line of the file.
*	Output: output.txt - Stores the sequence alignments
*	Also outputs a score for the alignment.
*/
// Forrest Ireland Mark Klick John Talbot
//
#define HERE cout << "At line " << __LINE__ << endl;

#include "NW_MultiAlign.h"
#include "Multi_Sequence.h"

using namespace std;

int main() {
	MultiSequence* MS1 = new MultiSequence();

	char* c = new char[3];
	MS1->setNumSequences(3);

	c[0] = 'A';
	c[1] = 'A';
	c[2] = 'A';
	MS1->push_back(c);

	c[0] = 'T';
	c[1] = 'T';
	c[2] = 'T';
	MS1->push_back(c);

	c[0] = 'C';
	c[1] = '-';
	c[2] = 'G';
	MS1->push_back(c);

	c[0] = 'G';
	c[1] = 'G';
	c[2] = 'G';
	MS1->push_back(c);

    
    
	MultiSequence* MS2 = new MultiSequence();
	
	MS2->setNumSequences(3);

	c[0] = 'A';
	c[1] = 'A';
	c[2] = 'A';
	MS2->push_back(c);

	c[0] = 'T';
	c[1] = 'T';
	c[2] = 'T';
	MS2->push_back(c);

	c[0] = 'C';
	c[1] = '-';
	c[2] = 'G';
	MS2->push_back(c);

	c[0] = 'G';
	c[1] = 'G';
	c[2] = 'G';
	MS2->push_back(c);
    
    NWMultiAlign NWMA;
    
    NWMA.SetMultiSequence(MS1, 1);
    NWMA.SetMultiSequence(MS2, 2);
    
    NWMA.AlignMultiSequences();
    
    NWMA.WriteAlignedMultiSequenceToFile("MS.txt");
    NWMA.OutputTraceAndScoringMatrices("Mat.txt");
    
	return 0;
}
