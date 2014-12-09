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

//#define READ_FROM_FASTA

using namespace std;

int main() {
	vector<MultiSequence> SeqSet;

	ifstream input;
#ifdef READ_FROM_FASTA
	input.open("InputFASTAFiles.txt");
#else
    input.open("InputSequences.txt");
#endif
	int i = 0;
	while (!input.eof()) {
		Sequence seq;
#ifdef READ_FROM_FASTA
		string filename;
		getline(input, filename);
		seq.setSequenceFromFASTAFile(filename);
#else
        string str;
        getline(input, str);
        seq.setSequence(str);
#endif
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

	return 0;
}
