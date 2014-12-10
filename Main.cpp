/*
*	Implements the Needleman-Wunsch Algorithm for aligning two sequences
*	Input: input.txt - contains the two sequences to be aligned. Each one on a separate line of the file.
*	Output: output.txt - Stores the sequence alignments
*	Also outputs a score for the alignment.
*/
// Forrest Ireland
//
#define HERE std::cout << "At line " << __LINE__ << std::endl;

#include "NW_MultiAlign.h"
#include "Multi_Sequence.h"
#include "NJ.h"

using namespace std;

int main() {
	vector<MultiSequence> SeqSet;
	ifstream input;
    input.open("InputFilenames.txt");
	int i = 0;
	while (!input.eof()) {
		Sequence seq;
		string filename;
		getline(input, filename);
		if (!seq.setSequenceFromFASTAFile(filename))
            break;
		seq.setAccessionNum(i);
		SeqSet.push_back(MultiSequence(seq));
		i++;
	}
	input.close();

	// Start Timer
	std::clock_t startTime = clock();

	NeighborJoin(SeqSet);

	SeqSet[0].WriteMultiSequenceToFile("ALIGNMENT.txt");

	//End Timer
	std::clock_t endTime = clock();

	cout << "Run Time: " << double(endTime - startTime) / (double)CLOCKS_PER_SEC << " seconds." << endl;

	string str;
	getline(cin, str);

	return 0;
}
