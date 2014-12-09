#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <istream>
#include <vector>
#include <stack>
#include <algorithm>

#include "Multi_Sequence.h"

typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<double> VD;
typedef std::vector<VD> VVD;
typedef std::vector<char> VC;

class NWMultiAlign {
public:
	NWMultiAlign();
	~NWMultiAlign();

	bool AlignMultiSequences();
	bool CreateAlignScoreMatrix();
	bool CreateTracebackMatrix();
	void InitializeAlignScoreMatrix();
	void InitializeTracebackMatrix();
	void DeleteAlignScoreMatrix();
	void DeleteTracebackMatrix();
  
	void SetMultiSequence(MultiSequence& MSIn, int seqNum);
	MultiSequence& GetAlignedMultiSequence();

	void ClearMSOut();
  
	bool WriteAlignedMultiSequenceToFile(std::string filename);
	void OutputTraceAndScoringMatrices(std::string filename);

	double alignmentScore;
	double SequenceDistance();
	int LevenshteinDistance();

private:
// Variables
	MultiSequence MS1, MS2, MSOut;
  
	int MATCH = 2;
	int MISMATCH = -1;
	int GAP = -2;
	VVI alignmentScoreMatrix;
	VVI traceBackMatrix;
// Functions

	// Returns largest of 3 ints
	int max3(int A, int B, int C);
};
