// Class: MultiSequence Header
// Forrest Ireland
// Dec. 3, 2014

/*
*   Implements a MultiSequence class
*
*   Acts as a container for multiple sequences that have been aligned. Allows alignment do be performed on a set of aligned sequences
*/

#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "Sequence.h"

typedef std::vector<char> VC;

class MultiSequence {
public:
	MultiSequence() { }
	MultiSequence(Sequence s);
	~MultiSequence();

	/************************************************/
	/* OPERATOR FUNCTIONS                           */
	/************************************************/

	// Returns an array of characters at a specific index of the multiSequence
	VC operator[](int i);
	void operator()(int i, VC& c);

	// Pushes an array of char onto the end of the multiSequence sequences. Size of char array needs to be equal to the number of sequences in the multiSequence
	void push_back(VC& c);

	// Set number of sequences in multiSequence. Must be done before adding any characters
	void setNumSequences(int num);

	// Adds a single sequence to produce a multiSequence of only 1 sequence. Can only add a sequence to an empty multiSequence
	void setFirstSequence(Sequence& seq);

	// Sets the IDs of all sequences in multiSequence
	void setSequenceIDs(std::vector<int>& id);

	/************************************************/
	/* GET FUNCTIONS                                */
	/************************************************/

	// Writes all sequences contained into a file
	void WriteMultiSequenceToFile(std::string filename);

	// Returns the number of Sequences contained in multiSequence
	int numSequences();

	// Returns length of sequences in multiSequence vector
	int Length();

	// Returns the sequence ID of the input index
	int getSequenceIDs(int index);

private:
	void push_back(char c);

	std::vector<Sequence> multiSequence;
	int length;
};
