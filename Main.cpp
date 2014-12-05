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
    string str = "AATAGA";
    Sequence S;
    S.setSequence(str);
    MS1->setFirstSequence(S);
    
    
	MultiSequence* MS2 = new MultiSequence();
    string str1 = "AATTGA";
    Sequence S1;
    S1.setSequence(str1);
    MS2->setFirstSequence(S1);
    
    NWMultiAlign NWMA;
    
    NWMA.SetMultiSequence(MS1, 1);
    NWMA.SetMultiSequence(MS2, 2);
    
    NWMA.AlignMultiSequences();
    
    NWMA.WriteAlignedMultiSequenceToFile("Sequence.txt");
    NWMA.OutputTraceAndScoringMatrices("Mat.txt");
    
	return 0;
}