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
#include "NJ.h"

using namespace std;

int main() {
	vector<MultiSequence*>* SeqSet;
	SeqSet = new vector<MultiSequence*>;

	ifstream input;
	input.open("InputSequences.txt");
	int i = 0;
	while (!input.eof()) {
		string filename;
		getline(input, filename);
		Sequence seq;
		seq.setSequenceFromFASTAFile(filename);
		seq.setAccessionNum(i);
		SeqSet->push_back(new MultiSequence(seq));
		i++;
	}
	input.close();

	// Start Timer
	std::clock_t startTime = clock();

	NeighborJoin(SeqSet);

	//End Timer
	std::clock_t endTime = clock();

	cout << "Run Time: " << double(endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

	(*SeqSet)[0]->WriteMultiSequenceToFile("ALIGNMENT.txt");

	return 0;
}
