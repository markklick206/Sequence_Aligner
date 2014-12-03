#define HERE std::cout << "At line " << __LINE__ << std::endl;

#include "NW_Align.h"
#include <Windows.h>

NWAlignment::NWAlignment() {
	alignmentScoreMatrix = 0;
	traceBackMatrix = 0;
}

NWAlignment::~NWAlignment() { }

/************************************************/
/* CONTROL FUNCTIONS							*/
/************************************************/

bool NWAlignment::AlignSequences() {
	
	std::stack<char> S1A, S2A;
	if (!CreateAlignScoreMatrix())
		return false;
	if (!CreateTracebackMatrix())
		return false;

	alignmentScore = 0;

	InitializeAlignScoreMatrix();
	InitializeTracebackMatrix();

	int rows = Seq1.Length();
	int cols = Seq2.Length();
	int score;
	int choice[3];
	
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			if (Seq1[i - 1] == Seq2[j - 1])
				score = MATCH;
			else
				score = MISMATCH;

			choice[0] = alignmentScoreMatrix[i - 1][j - 1] + score;
			choice[1] = alignmentScoreMatrix[i - 1][j] + GAP;
			choice[2] = alignmentScoreMatrix[i][j - 1] + GAP;
			alignmentScoreMatrix[i][j] = max3(choice[0], choice[1], choice[2]);

			if (alignmentScoreMatrix[i][j] == choice[0])
				traceBackMatrix[i][j] = 0;
			else if (alignmentScoreMatrix[i][j] == choice[1])
				traceBackMatrix[i][j] = 1;
			else
				traceBackMatrix[i][j] = -1;
		}
	}
	
	int i = Seq1.Length();
	int j = Seq2.Length();
	int x = 0;
	int y = 0;
	while (i > 0 || j > 0) {
		if (traceBackMatrix[i][j] == 0) {
			alignmentScore += alignmentScoreMatrix[i][j];
			S1A.push(Seq1[i - 1]);
			x++;
			S2A.push(Seq2[j - 1]);
			y++;
			i--;
			j--;
		}
		else if (traceBackMatrix[i][j] == 1) {
			alignmentScore += alignmentScoreMatrix[i][j];
			S1A.push(Seq1[i - 1]);
			x++;
			S2A.push('-');
			y++;
			i--;
		}
		else if (traceBackMatrix[i][j] == -1) {
			alignmentScore += alignmentScoreMatrix[i][j];
			S1A.push('-');
			x++;
			S2A.push(Seq2[j - 1]);
			y++;
			j--;
		}
	}
	
	while (!S1A.empty()) {
		Seq1_Aligned.push_back(S1A.top());
		S1A.pop();
	}
	while (!S2A.empty()) {
		Seq2_Aligned.push_back(S2A.top());
		S2A.pop();
	}

	return true;
}

void NWAlignment::CloseNWAlign() {
	DeleteAlignScoreMatrix();
	DeleteTracebackMatrix();
    ClearAlignmentSequences(); 
}

void NWAlignment::ClearAlignmentSequences() {

}

bool NWAlignment::CreateAlignScoreMatrix() {
	if (Seq1.Length() == 0 || Seq2.Length() == 0)
		return false;


	alignmentScoreMatrix = new int*[Seq1.Length() + 1];
	for (int i = 0; i < Seq1.Length() + 1; i++)
		alignmentScoreMatrix[i] = new int[Seq2.Length() + 1];


	for (int i = 0; i <= Seq1.Length(); i++) {
		for (int j = 0; j <= Seq2.Length(); j++) {
			alignmentScoreMatrix[i][j] = 0;
		}
	}

	return true;
}

bool NWAlignment::CreateTracebackMatrix() {
	if (Seq1.Length() == 0 || Seq2.Length() == 0)
		return false;


	traceBackMatrix = new int*[Seq1.Length() + 1];
	for (int i = 0; i < Seq1.Length() + 1; i++)
		traceBackMatrix[i] = new int[Seq2.Length() + 1];


	for (int i = 0; i <= Seq1.Length(); i++) {
		for (int j = 0; j <= Seq2.Length(); j++) {
			traceBackMatrix[i][j] = 0;
		}
	}

	return true;
}

void NWAlignment::InitializeAlignScoreMatrix() {
	for (int i = 0; i <= Seq1.Length(); i++) {
        alignmentScoreMatrix[i][0] = MISMATCH * i;
    }


	for (int i = 0; i <= Seq2.Length(); i++) {
        alignmentScoreMatrix[0][i] = MISMATCH * i;
    }
}

void NWAlignment::InitializeTracebackMatrix() {
	for (int i = 0; i <= Seq1.Length(); i++) {
        traceBackMatrix[i][0] = 1;
    }


	for (int i = 0; i <= Seq2.Length(); i++) {
        traceBackMatrix[0][i] = -1;
    }
}

void NWAlignment::DeleteAlignScoreMatrix() {
	if (alignmentScoreMatrix) {
		for (int i = 0; i < Seq1.Length() + 1; i++)
			delete[] alignmentScoreMatrix[i];

		delete[] alignmentScoreMatrix;
	}
}

void NWAlignment::DeleteTracebackMatrix() {
	if (traceBackMatrix) {
		for (int i = 0; i < Seq1.Length() + 1; i++)
			delete[] traceBackMatrix[i];

		delete[] traceBackMatrix;
	}
}

/************************************************/
/* INPUT FUNCTIONS								*/
/************************************************/

void NWAlignment::SetInputSequence(Sequence seq, int seqNum) {
	if (seqNum == 1) {
		Seq1 = seq;
	}

	if (seqNum == 2) {
		Seq2 = seq;
	}
}

/************************************************/
/* OUTPUT FUNCTIONS								*/
/************************************************/

void NWAlignment::GetAlignedSequence(Sequence &seq, int seqNum) {
    if (seqNum == 1)
        seq = Seq1_Aligned;
    if (seqNum == 2)
        seq = Seq2_Aligned;
}

bool NWAlignment::WriteAlignedSequencesToFile(std::string filename) {
	std::ofstream output;
    output.open(filename);
    if (output.is_open()) {
		HERE
		output << "Alignment Score: " << alignmentScore / static_cast<double>(Seq1_Aligned.Length()) << std::endl << std::endl;
        for (int i = 0; i < Seq1_Aligned.Length(); i++)
            output << Seq1_Aligned[i];
		output << std::endl;
		
        for (int i = 0; i < Seq2_Aligned.Length(); i++)
            output << Seq2_Aligned[i];
		output << std::endl;
        output.close();
		
        return true;
    }
    else {
		std::string error = "Error: Unable to open - ";
		error.append(filename);
		std::wstring temp = s2ws(error);
		LPCWSTR err = temp.c_str();
		MessageBox(NULL, err, L"Error", MB_OK);
    }
    return false;
}
 
void NWAlignment::PrintSequenceToConsole(int seqNum, bool aligned) {
    if (!aligned) {
		if (seqNum == 1) {
			std::cout << "Seq1: ";
			for (int i = 0; i < Seq1.Length(); i++)
				std::cout << Seq1[i];
			std::cout << std::endl;
		}
		if (seqNum == 2) {
			std::cout << "Seq2: ";
			for (int i = 0; i < Seq2.Length(); i++)
				std::cout << Seq2[i];
			std::cout << std::endl;
		}
    }
    else {
        if (seqNum == 1) {
			std::cout << "Aligned SEQ1: ";
			for (int i = 0; i < Seq1_Aligned.Length(); i++)
				std::cout << Seq1_Aligned[i];
			std::cout << std::endl;
        }
        if (seqNum == 2) {
			std::cout << "Aligned SEQ2: ";
			for (int i = 0; i < Seq2_Aligned.Length(); i++)
				std::cout << Seq2_Aligned[i];
			std::cout << std::endl;
        }
    }
}

void NWAlignment::OutputTraceAndScoringMatrices(std::string filename) {
	std::ofstream file;
	file.open(filename);
	file << "Alignment Score Matrix" << std::endl;

	file << "      ";
	for (int i = 0; i < Seq2.Length(); i++)
		file << std::setw(3) << Seq2[i];
	file << std::endl;
	for (int i = 0; i <= Seq1.Length(); i++) {
		if (i != 0)
			file << Seq1[i - 1] << "  ";
		else
			file << "   ";
		for (int j = 0; j <= Seq2.Length(); j++) {
			file << std::setw(3) << alignmentScoreMatrix[i][j];
		}
		file << std::endl;
	}
	file << std::endl << std::endl;


	file << "Traceback Matrix" << std::endl;
	file << "      ";
	for ( int i = 0; i < Seq2.Length(); i++)
		file << std::setw(3) << Seq2[i];
	file << std::endl;
	for (int i = 0; i <= Seq1.Length(); i++) {
		if (i != 0)
			file << Seq1[i - 1] << "  ";
		else
			file << "   ";
		for (int j = 0; j <= Seq2.Length(); j++) {
			file << std::setw(3) << traceBackMatrix[i][j];
		}
		file << std::endl;
	}
}

/************************************************/
/* CONFIGURATION FUNCTIONS						*/
/************************************************/

// NEED TO IMPLEMENT
bool NWAlignment::SetMismatchScoringMatrix(std::string filename) {
    return true;
}

/************************************************/
/* PRIVATE: UTILITY FUNCTIONS					*/
/************************************************/

// Returns largest of 3 ints
int NWAlignment::max3(int A, int B, int C) {
	return max(A, max(B, C));
}

// Converts a string into a wstring
std::wstring NWAlignment::s2ws(const std::string& s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring r(buf);
	delete[] buf;
	return r;
}