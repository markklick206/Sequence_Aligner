#include "NW_MultiAlign.h"

bool NWMultiAlign::AlignMultiSequences() {
  std::stack<char*> MSStack;
  if (!CreateAlignScoreMatrix())
  	return false;
  if (!CreateTraceBackMatrix())
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
			// Pair wise sum scoring here
			char* ci = MS1[i];
			char* cj = MS2[j];
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
	
	int i = MS1->Length();
	int j = MS2->Length();
	alignmentScore = alignmentScoreMatrix[i][j];
	while (i > 0 || j > 0) {
		if (traceBackMatrix[i][j] == 0) {
			char* c = new char[x + y];
			char* ms = MS1[i];
			for (int k = 0; k < x; k++)
				c[k] = ms[k];
			ms = MS2[j];
			for (int k = x; k < (x + y); k++)
				c[k] = ms[k - x];
			MSStack.push(c)
			i--;
			j--;
		}
		else if (traceBackMatrix[i][j] == 1) {
			char* c = new char[x + y];
			char* ms = MS1[i];
			for (int k = 0; k < x; k++)
				c[k] = ms[k];
			ms = MS2[j];
			for (int k = x; k < (x + y); k++)
				c[k] = '-';
			MSStack.push(c)
			i--;
		}
		else if (traceBackMatrix[i][j] == -1) {
			char* c = new char[x + y];
			char* ms = MS1[i];
			for (int k = 0; k < x; k++)
				c[k] = '-';
			ms = MS2[j];
			for (int k = x; k < (x + y); k++)
				c[k] = ms[k - x];
			MSStack.push(c)
			j--;
		}
	}
	
	MSOut = new MultiSequence();
	MSOut->setNumSequences(x + y);
	while (!MSStack.empty()) {
		MSOut->push_back(MSStack.top());
		MSStack.pop();
	}
	
	return true;
}

bool NWMultiAlign::CreateAlignScoreMatrix() {
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

bool NWMultiAlign::CreateTracebackMatrix() {
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

void NWMultiAlign::InitializeAlignScoreMatrix() {
	for (int i = 0; i <= Seq1.Length(); i++) {
        alignmentScoreMatrix[i][0] = MISMATCH * i;
    }


	for (int i = 0; i <= Seq2.Length(); i++) {
        alignmentScoreMatrix[0][i] = MISMATCH * i;
    }
}

void NWMultiAlign::InitializeTracebackMatrix() {
	for (int i = 0; i <= Seq1.Length(); i++) {
        traceBackMatrix[i][0] = 1;
    }


	for (int i = 0; i <= Seq2.Length(); i++) {
        traceBackMatrix[0][i] = -1;
    }
}

void NWMultiAlign::DeleteAlignScoreMatrix() {
	if (alignmentScoreMatrix) {
		for (int i = 0; i < Seq1.Length() + 1; i++)
			delete[] alignmentScoreMatrix[i];

		delete[] alignmentScoreMatrix;
	}
}

void NWMultiAlign::DeleteTracebackMatrix() {
	if (traceBackMatrix) {
		for (int i = 0; i < Seq1.Length() + 1; i++)
			delete[] traceBackMatrix[i];

		delete[] traceBackMatrix;
	}
}

void NWMultiAlign::SetMultiSequence(MultiSequence* MSIn, int seqNum){
	if (seqNum == 1){
		MS1 = MSIn;
	}
	if (seqNum == 2){
		MS2 = MSIn;
	}
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
