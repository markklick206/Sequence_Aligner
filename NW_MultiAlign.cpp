#include "NW_MultiAlign.h"

#define HERE std::cout << "At line " << __LINE__ << std::endl;

NWMultiAlign::NWMultiAlign() {
	alignmentScoreMatrix = 0;
	traceBackMatrix = 0;
	MS1 = 0;
	MS2 = 0;
	MSOut = 0;
}

NWMultiAlign::~NWMultiAlign() {
	DeleteAlignScoreMatrix();
	DeleteTracebackMatrix();
}

bool NWMultiAlign::AlignMultiSequences() {
	std::stack<char*> MSStack;
	if (!CreateAlignScoreMatrix())
		return false;
	if (!CreateTracebackMatrix())
		return false;
    
	alignmentScore = 0;

	InitializeAlignScoreMatrix();
	InitializeTracebackMatrix();
	
	int rows = MS1->Length();
	int cols = MS2->Length();
	int x = MS1->numSequences();
	int y = MS2->numSequences();
	int score;
	int choice[3];
	
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			char* ci = (*MS1)[i - 1];
			char* cj = (*MS2)[j - 1];
			int pairWiseScore = 0;
			
			for (int k = 0; k < x; k++) {
				for (int l = 0; l < y; l++) {
					if (ci[k] == cj[l])
						pairWiseScore += MATCH;
					else
						pairWiseScore += MISMATCH;
				}
			}
			
			score = pairWiseScore;
			//double score = pairWiseScore \ (double) (x * y);

			choice[0] = alignmentScoreMatrix[i - 1][j - 1] + score;
			choice[1] = alignmentScoreMatrix[i - 1][j] + (GAP * y);
			choice[2] = alignmentScoreMatrix[i][j - 1] + (GAP * x);
			alignmentScoreMatrix[i][j] = max3(choice[0], choice[1], choice[2]);
			
			if (alignmentScoreMatrix[i][j] == choice[0])
				traceBackMatrix[i][j] = 0;
			else if (alignmentScoreMatrix[i][j] == choice[1])
				traceBackMatrix[i][j] = 1;
			else
				traceBackMatrix[i][j] = -1;
			
		}
	}
	
	int i1 = MS1->Length();
	int j1 = MS2->Length();
	char* c;
	char* ms;
	alignmentScore = alignmentScoreMatrix[i1][j1];
	while (i1 > 0 || j1 > 0) {
		if (traceBackMatrix[i1][j1] == 0) {
			c = new char[x + y];
			ms = (*MS1)[i1 - 1];
			for (int k = 0; k < x; k++)
				c[k] = ms[k];
			delete [] ms;
			ms = (*MS2)[j1 - 1];
			for (int k = x; k < (x + y); k++)
				c[k] = ms[k - x];
			delete [] ms;
			MSStack.push(c);
			i1--;
			j1--;
		}
		else if (traceBackMatrix[i1][j1] == 1) {
			c = new char[x + y];
			ms = (*MS1)[i1 - 1];
			for (int k = 0; k < x; k++)
				c[k] = ms[k];
			delete [] ms;
			for (int k = x; k < (x + y); k++)
				c[k] = '-';
			MSStack.push(c);
			i1--;
		}
		else if (traceBackMatrix[i1][j1] == -1) {
			c = new char[x + y];
			for (int k = 0; k < x; k++)
				c[k] = '-';
			ms = (*MS2)[j1 - 1];
			for (int k = x; k < (x + y); k++)
				c[k] = ms[k - x];
			delete [] ms;
			MSStack.push(c);
			j1--;
		}
	}
    
	std::vector<int> IDS;
	for (unsigned int i = 0; i < MS1->numSequences(); i++)
		IDS.push_back(MS1->getSequenceIDs(i));
	for (unsigned int i = 0; i < MS2->numSequences(); i++)
		IDS.push_back(MS2->getSequenceIDs(i));

	MSOut = new MultiSequence();
	MSOut->setNumSequences(x + y);
	while (!MSStack.empty()) {
		MSOut->push_back(MSStack.top());
		MSStack.pop();
	}

	MSOut->setSequenceIDs(IDS);

	alignmentScore = alignmentScore / (double)MSOut->Length();
	return true;
}

bool NWMultiAlign::CreateAlignScoreMatrix() {
	if (MS1->Length() == 0 || MS2->Length() == 0)
		return false;


	alignmentScoreMatrix = new int*[MS1->Length() + 1];
	for (int i = 0; i < MS1->Length() + 1; i++)
		alignmentScoreMatrix[i] = new int[MS2->Length() + 1];


	for (int i = 0; i <= MS1->Length(); i++) {
		for (int j = 0; j <= MS2->Length(); j++) {
			alignmentScoreMatrix[i][j] = 0;
		}
	}

	return true;
}

bool NWMultiAlign::CreateTracebackMatrix() {
	if (MS1->Length() == 0 || MS2->Length() == 0)
		return false;


	traceBackMatrix = new int*[MS1->Length() + 1];
	for (int i = 0; i < MS1->Length() + 1; i++)
		traceBackMatrix[i] = new int[MS2->Length() + 1];


	for (int i = 0; i <= MS1->Length(); i++) {
		for (int j = 0; j <= MS2->Length(); j++) {
			traceBackMatrix[i][j] = 0;
		}
	}

	return true;
}

void NWMultiAlign::InitializeAlignScoreMatrix() {
	for (int i = 0; i <= MS1->Length(); i++) {
        alignmentScoreMatrix[i][0] = MISMATCH * i;
    }


	for (int i = 0; i <= MS2->Length(); i++) {
        alignmentScoreMatrix[0][i] = MISMATCH * i;
    }
}

void NWMultiAlign::InitializeTracebackMatrix() {
	for (int i = 0; i <= MS1->Length(); i++) {
        traceBackMatrix[i][0] = 1;
    }


	for (int i = 0; i <= MS2->Length(); i++) {
        traceBackMatrix[0][i] = -1;
    }
}

void NWMultiAlign::DeleteAlignScoreMatrix() {
	if (alignmentScoreMatrix) {
		for (int i = 0; i < MS1->Length() + 1; i++)
			delete[] alignmentScoreMatrix[i];

		delete[] alignmentScoreMatrix;
	}
	alignmentScoreMatrix = 0;
}

void NWMultiAlign::DeleteTracebackMatrix() {
	if (traceBackMatrix) {
		for (int i = 0; i < MS1->Length() + 1; i++)
			delete[] traceBackMatrix[i];

		delete[] traceBackMatrix;
	}
	traceBackMatrix = 0;
}

void NWMultiAlign::SetMultiSequence(MultiSequence* MSIn, int seqNum){
	if (seqNum == 1){
		MS1 = MSIn;
	}
	if (seqNum == 2){
		MS2 = MSIn;
	}
}

MultiSequence* NWMultiAlign::GetAlignedMultiSequence() {
	return MSOut;
}

bool NWMultiAlign::WriteAlignedMultiSequenceToFile(std::string filename) {
		MSOut->WriteMultiSequenceToFile(filename);
		return true;
}

void NWMultiAlign::OutputTraceAndScoringMatrices(std::string filename) {
	std::ofstream file;
	file.open(filename);
	file << "Alignment Score Matrix" << std::endl;

	file << "      ";
	file << std::endl;
	for (int i = 0; i <= MS1->Length(); i++) {
		for (int j = 0; j <= MS2->Length(); j++) {
			file << std::setw(3) << alignmentScoreMatrix[i][j];
		}
		file << std::endl;
	}
	file << std::endl << std::endl;


	file << "Traceback Matrix" << std::endl;
	file << "      ";
	file << std::endl;
	for (int i = 0; i <= MS1->Length(); i++) {
		for (int j = 0; j <= MS2->Length(); j++) {
			file << std::setw(3) << traceBackMatrix[i][j];
		}
		file << std::endl;
	}
}

int NWMultiAlign::max3(int A, int B, int C) {
    return max(A, max(B, C));
}

double NWMultiAlign::SequenceDistance() {
	if (MSOut->numSequences() != 2) {
		std::cout << "NO" << std::endl;
		return -3.14159;
	}

	int x = MSOut->Length();
	int matches = 0, l = 0;

	for (int i = 0; i < x; i++){
		if ((*MSOut)[0][i] != '-' && (*MSOut)[1][i] != '-') {
			if (((*MSOut)[i][0] == (*MSOut)[i][1]))
				matches++;
			l++;
		}
	}

	return 1 - (matches / (double)l);
}

int NWMultiAlign::LevenshteinDistance() {
	if (MSOut->numSequences() != 2)
		return -314159;

	int x = MSOut->Length();

	int** Mat = new int*[x];
	for (int i = 0; i < x; i++)
		Mat[i] = new int[x];

	for (int i = 0; i < x; i++) {
		Mat[i][0] = i;//JOHN
		Mat[0][i] = i;
	}

	for (int i = 1; i < x; i++) {
		for (int j = 1; j < x; j++) {
			if ((*MSOut)[i - 1][0] == (*MSOut)[j - 1][1])
				Mat[i][j] = Mat[i - 1][j - 1];
			else {
				Mat[i][j] = min(Mat[i - 1][j] + 1, min(Mat[i][j - 1] + 1, Mat[i - 1][j - 1] + 1));
			}
		}
	}

	int a = Mat[x - 1][x - 1];

	for (int i = 0; i < x; i++)
		delete [] Mat[i];
	delete [] Mat;
	return a;
}
