/*	Forrest Ireland
*
*	Needleman-Wunsch Sequence Alignment
*	
*	Class that allows for use of the Needleman-Wunsch Algorithm for sequence alignment
*	
*	GOALS FOR THE CLASS IMPLEMENTATION
	_	Add ability to read sequences from more than text files
	_	Write sequences to text files
	_	Implement the NW algorithm so it works with 2 sequences and simple scoring
		-	Simple scoring: MATCH = +1, MISMATCH = -1, GAP = -1
	_	Add a more enhanced scoring matrix for alignments
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <stack>

#include "Sequence.h"

class NWAlignment {
public:	
// Variables

// Functions
	NWAlignment();
	~NWAlignment();

	/************************************************/
	/* CONTROL FUNCTIONS							*/
	/************************************************/

	// Starts the algorithm for aligning two sequences. Needs to have two sequences already loaded into the class members
	bool AlignSequences();

	// Deletes all memory associated with class. Call before destructor
	void CloseNWAlign();

	// Deletes or clears the memory associated with sequences
	void ClearAlignmentSequences();

	// Creates the alignment score matrix. Requires Seq1 and Seq2 to contain sequences
	bool CreateAlignScoreMatrix();

	// Creates the traceback matrix. Requires Seq1 and Seq2 to contain sequences
	bool CreateTracebackMatrix();
    
    // Initialize alignment score matrix with gap penalty multiples of the row/column number for the first row and column
    void InitializeAlignScoreMatrix();
    
    // Initialize the traceback matrix to show left movement in the top row, and up movement in the first column
    void InitializeTracebackMatrix();

	// Deletes the memory associated with the Alignment score matrix
	void DeleteAlignScoreMatrix();

	// Deletes memory associated with the traceback matrix
	void DeleteTracebackMatrix();

	/************************************************/
	/* INPUT FUNCTIONS								*/
	/************************************************/

	// ## Reads a sequence from a FASTA file, sets the sequence seqNum to this input sequence. IDK how this works?
	bool ReadInputSequenceFromFASTA(std::string filename, int seqNum);

	// Reads a pair of sequences from a text file
	bool ReadInputSequencesFromFile(std::string filename);

	// Allows direct setting of the pre-aligned sequences
	void SetInputSequence(Sequence &seq, int seqNum);

	/************************************************/
	/* OUTPUT FUNCTIONS								*/
	/************************************************/

	// Copies aligned sequence seqNum (1 or 2) to the string seq
	void GetAlignedSequence(Sequence &seq, int seqNum);

	// Writes the new aligned sequences to a text file
	bool WriteAlignedSequencesToFile(std::string filename);
    
    // Prints a sequence to the console
    void PrintSequenceToConsole(int seqNum, bool aligned);

	// Writes the traceback and scoring matrices to a file
	void OutputTraceAndScoringMatrices(std::string filename);

	/************************************************/
	/* CONFIGURATION FUNCTIONS						*/
	/************************************************/

	// Loads a mismatch scoring matrix from a text file
	bool SetMismatchScoringMatrix(std::string filename);

private:
// Variables
	Sequence Seq1, Seq2, Seq1_Aligned, Seq2_Aligned;
	int Seq1Length, Seq2Length;
	int alignmentScore;
	int MATCH = 2;
	int MISMATCH = -1;
	int GAP = -2;
	bool SeqAligned;
	std::vector<std::vector<int>> mismatchScoringMatrix; // Future mismatch scoring matrix
	int** alignmentScoreMatrix;
	int** traceBackMatrix;

// Functions

	// Returns largest of 3 ints
	int max3(int A, int B, int C);

	// Converts a string into a wstring
	std::wstring s2ws(const std::string& s);

};