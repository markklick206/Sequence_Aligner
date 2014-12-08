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

//#define READ_FROM_FASTA
#define READ_FROM_INPUT_SEQ_TXT

using namespace std;

int main() {
	vector<MultiSequence*>* SeqSet;
	SeqSet = new vector<MultiSequence*>;

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
    
    for (int i = 0; i < SeqSet->size(); i++)
        delete (*SeqSet)[i];
    
    delete SeqSet;

	return 0;
}
