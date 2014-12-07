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
#include <windows.h>
#include <sstream>

#include "Sequence.h"

class MultiSequence {
public:
	MultiSequence() { }
	MultiSequence(Sequence s);
	~MultiSequence();

	/************************************************/
	/* OPERATOR FUNCTIONS                           */
	/************************************************/

	// Returns an array of characters at a specific index of the multiSequence
	char* operator[](int i);

	// Pushes an array of char onto the end of the multiSequence sequences. Size of char array needs to be equal to the number of sequences in the multiSequence
	void push_back(char* c);

	// Set number of sequences in multiSequence. Must be done before adding any characters
	void setNumSequences(int num);

	// Adds a single sequence to produce a multiSequence of only 1 sequence. Can only add a sequence to an empty multiSequence
	void setFirstSequence(Sequence &seq);

	/************************************************/
	/* GET FUNCTIONS                                */
	/************************************************/

	// Writes all sequences contained into a file
	void WriteMultiSequenceToFile(std::string filename);

	// Returns the number of Sequences contained in multiSequence
	int numSequences();

	// Returns length of sequences in multiSequence vector
	int Length();

private:
	std::vector<Sequence> multiSequence;
	int length;
};
