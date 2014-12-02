#define HERE cout << "At line " << __LINE__ << endl;

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
	if (!CreateAlignScoreMatrix())
		return false;
	if (!CreateTracebackMatrix())
		return false;

	alignmentScore = 0;

	InitializeAlignScoreMatrix();
	InitializeTracebackMatrix();

	int rows = Seq1Length;
	int cols = Seq2Length;
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
	OutputTraceAndScoringMatrices("matrices.txt");
	int i = Seq1Length;
	int j = Seq2Length;
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
    Seq1.clear();
    Seq2.clear();
    Seq1_Aligned.clear();
    Seq2_Aligned.clear();
}

bool NWAlignment::CreateAlignScoreMatrix() {
	if (Seq1Length == 0 || Seq2Length == 0)
		return false;


	alignmentScoreMatrix = new int*[Seq1Length + 1];
	for (int i = 0; i < Seq1Length + 1; i++)
		alignmentScoreMatrix[i] = new int[Seq2Length + 1];


	for (int i = 0; i <= Seq1Length; i++) {
		for (int j = 0; j <= Seq2Length; j++) {
			alignmentScoreMatrix[i][j] = 0;
		}
	}

	return true;
}

bool NWAlignment::CreateTracebackMatrix() {
	if (Seq1Length == 0 || Seq2Length == 0)
		return false;


	traceBackMatrix = new int*[Seq1Length + 1];
	for (int i = 0; i < Seq1Length + 1; i++)
		traceBackMatrix[i] = new int[Seq2Length + 1];


	for (int i = 0; i <= Seq1Length; i++) {
		for (int j = 0; j <= Seq2Length; j++) {
			traceBackMatrix[i][j] = 0;
		}
	}

	return true;
}

void NWAlignment::InitializeAlignScoreMatrix() {
	for (int i = 0; i <= Seq1Length; i++) {
        alignmentScoreMatrix[i][0] = MISMATCH * i;
    }


	for (int i = 0; i <= Seq2Length; i++) {
        alignmentScoreMatrix[0][i] = MISMATCH * i;
    }
}

void NWAlignment::InitializeTracebackMatrix() {
	for (int i = 0; i <= Seq1Length; i++) {
        traceBackMatrix[i][0] = 1;
    }


	for (int i = 0; i <= Seq2Length; i++) {
        traceBackMatrix[0][i] = -1;
    }
}

void NWAlignment::DeleteAlignScoreMatrix() {
	if (alignmentScoreMatrix) {
		for (int i = 0; i < Seq1Length + 1; i++)
			delete[] alignmentScoreMatrix[i];

		delete[] alignmentScoreMatrix;
	}
}

void NWAlignment::DeleteTracebackMatrix() {
	if (traceBackMatrix) {
		for (int i = 0; i < Seq1Length + 1; i++)
			delete[] traceBackMatrix[i];

		delete[] traceBackMatrix;
	}
}

/************************************************/
/* INPUT FUNCTIONS								*/
/************************************************/

bool NWAlignment::ReadInputSequenceFromFASTA(string filename, int seqNum) {
	ifstream fasta;
	string line, s;
	vector<char> seq;
	fasta.open(filename);
	if (fasta.is_open()) {
		getline(fasta, line);
		while (!fasta.eof()) {
			getline(fasta, line);
			s.append(line);
		}
		for (unsigned int i = 0; i < s.size(); i++)
			seq.push_back(s[i]);
		SetInputSequence(seq, seqNum);
		return true;
	}
	else {
		string error = "Error: Unable to open - ";
		error.append(filename);
		wstring temp = s2ws(error);
		LPCWSTR err = temp.c_str();
		MessageBox(NULL, err, L"Error", MB_OK);
		return false;
	}
    return false;
}

bool NWAlignment::ReadInputSequencesFromFile(string filename) {
    ifstream input;
	string str;
    input.open(filename);
    if (input.is_open()) {
        getline(input, (str));
		for (unsigned int i = 0; i < str.size(); i++)
			Seq1.push_back(str[i]);
		Seq1Length = static_cast<int>(Seq1.size());
		cout << "Seq1Length: " << Seq1Length << endl;
        input.ignore();
        getline(input, (str));
		for (unsigned int i = 0; i < str.size(); i++)
			Seq2.push_back(str[i]);
		Seq2Length = static_cast<int>(Seq2.size());
		cout << "Seq2Length: " << Seq2Length << endl;
        input.close();
        if (!Seq1.empty() && !Seq2.empty())
            return true;
    }
	else {
		string error = "Error: Unable to open - ";
		error.append(filename);
		wstring temp = s2ws(error);
		LPCWSTR err = temp.c_str();
		MessageBox(NULL, err, L"Error", MB_OK);
		return false;
    }
	return false;
}

void NWAlignment::SetInputSequence(vector<char> &seq, int seqNum) {
	if (seqNum == 1) {
		Seq1 = seq;
		Seq1Length = static_cast<int>(Seq1.size());
		cout << "Seq1Length: " << Seq1Length << endl;
	}

	if (seqNum == 2) {
		Seq2 = seq;
		Seq2Length = static_cast<int>(Seq2.size());
		cout << "Seq2Length: " << Seq2Length << endl;
	}
}

/************************************************/
/* OUTPUT FUNCTIONS								*/
/************************************************/

void NWAlignment::GetAlignedSequence(vector<char> &seq, int seqNum) {
    if (seqNum == 1)
        seq = Seq1_Aligned;
    if (seqNum == 2)
        seq = Seq2_Aligned;
}

bool NWAlignment::WriteAlignedSequencesToFile(string filename) {
    ofstream output;
    output.open(filename);
    if (output.is_open()) {
		output << "Alignment Score: " << alignmentScore << endl << endl;
        for (unsigned int i = 0; i < Seq1_Aligned.size(); i++)
            output << Seq1_Aligned[i];
        output << endl;
        for (unsigned int i = 0; i < Seq2_Aligned.size(); i++)
            output << Seq2_Aligned[i];
		output << endl;
        output.close();
        return true;
    }
    else {
		string error = "Error: Unable to open - ";
		error.append(filename);
		wstring temp = s2ws(error);
		LPCWSTR err = temp.c_str();
		MessageBox(NULL, err, L"Error", MB_OK);
    }
    return false;
}
 
void NWAlignment::PrintSequenceToConsole(int seqNum, bool aligned) {
    if (!aligned) {
		if (seqNum == 1) {
			cout << "Seq1: ";
			for (unsigned int i = 0; i < Seq1.size(); i++)
				cout << Seq1[i];
			cout << endl;
		}
		if (seqNum == 2) {
			cout << "Seq2: ";
			for (unsigned int i = 0; i < Seq2.size(); i++)
				cout << Seq2[i];
			cout << endl;
		}
    }
    else {
        if (seqNum == 1) {
            cout << "Aligned SEQ1: ";
            for (unsigned int i = 0; i < Seq1_Aligned.size(); i++)
                cout << Seq1_Aligned[i];
            cout << endl;
        }
        if (seqNum == 2) {
            cout << "Aligned SEQ2: ";
            for (unsigned int i = 0; i < Seq2_Aligned.size(); i++)
                cout << Seq2_Aligned[i];
            cout << endl;
        }
    }
}

void NWAlignment::OutputTraceAndScoringMatrices(string filename) {
	ofstream file;
	file.open(filename);
	file << "Alignment Score Matrix" << endl;

	file << "      ";
	for (unsigned int i = 0; i < Seq2.size(); i++)
		file << setw(3) << Seq2[i];
	file << endl;
	for (unsigned int i = 0; i <= Seq1.size(); i++) {
		if (i != 0)
			file << Seq1[i - 1] << "  ";
		else
			file << "   ";
		for (unsigned int j = 0; j <= Seq2.size(); j++) {
			file << setw(3) << alignmentScoreMatrix[i][j];
		}
		file << endl;
	}
	file << endl << endl;


	file << "Traceback Matrix" << endl;
	file << "      ";
	for (unsigned int i = 0; i < Seq2.size(); i++)
		file << setw(3) << Seq2[i];
	file << endl;
	for (unsigned int i = 0; i <= Seq1.size(); i++) {
		if (i != 0)
			file << Seq1[i - 1] << "  ";
		else
			file << "   ";
		for (unsigned int j = 0; j <= Seq2.size(); j++) {
			file << setw(3) << traceBackMatrix[i][j];
		}
		file << endl;
	}
}

/************************************************/
/* CONFIGURATION FUNCTIONS						*/
/************************************************/

// NEED TO IMPLEMENT
bool NWAlignment::SetMismatchScoringMatrix(string filename) {
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
wstring NWAlignment::s2ws(const string& s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	wstring r(buf);
	delete[] buf;
	return r;
}