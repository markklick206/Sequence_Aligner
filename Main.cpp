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
#include "NJ.h"

using namespace std;

int main() {
	vector<MultiSequence*>* SeqSet;
	SeqSet = new vector<MultiSequence*>;

	ifstream input;
	input.open("InputSequences.txt");
	int i = 0;
	while (!input.eof()) {
		string str;
		getline(input, str);
		Sequence seq;
		seq.setSequence(str);
		seq.setAccessionNum(i);
		SeqSet->push_back(new MultiSequence(seq));
		i++;
	}
	input.close();

	NeighborJoin(SeqSet);

	(*SeqSet)[0]->WriteMultiSequenceToFile("ALIGNEDSEQUENCESWEGNIURNU.txt");

	return 0;
}
