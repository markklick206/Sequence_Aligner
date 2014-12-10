// Class: MultiSequence Code
// Forrest Ireland
// Dec. 3, 2014

#include "Multi_Sequence.h"

MultiSequence::MultiSequence(Sequence s) {
	setFirstSequence(s);
}

MultiSequence::~MultiSequence() {

}

/************************************************/
/* OPERATOR FUNCTIONS                           */
/************************************************/

VC MultiSequence::operator[](int a) {
	VC c;
	for (unsigned int i = 0; i < multiSequence.size(); i++) {
		c.push_back(multiSequence[i][a]);
	}
	return c;
}


void MultiSequence::operator()(int a, VC& c) {
	c.clear();
	for (unsigned int i = 0; i < multiSequence.size(); i++) {
		c.push_back(multiSequence[i][a]);
	}
}

void MultiSequence::push_back(VC& c) {
	for (unsigned int i = 0; i < multiSequence.size(); i++) {
		multiSequence[i].push_back(c[i]);
	}
}

void MultiSequence::push_back(char c) {
	for (unsigned int i = 0; i < multiSequence.size(); i++) {
		multiSequence[i].push_back(c);
	}
}

void MultiSequence::setNumSequences(int num) {
	for (int i = 0; i < num; i++) {
		multiSequence.push_back(Sequence());
	}
}

void MultiSequence::setFirstSequence(Sequence &seq) {
	if (!multiSequence.empty()) {
		std::cout << "Sequence has been set. Cannot set another one." << std::endl;
		return;
	}
	setNumSequences(1);
	char c;
	for (int i = 0; i < seq.Length(); i++) {
		c = seq[i];
		push_back(c);
	}
	std::vector<int> ID;
	ID.push_back(seq.getAccessionNum());
	setSequenceIDs(ID);
}

void MultiSequence::setSequenceIDs(std::vector<int>& id) {
	for (unsigned int i = 0; i < multiSequence.size(); i++) {
		multiSequence[i].setAccessionNum(id[i]);
	}
}

/************************************************/
/* GET FUNCTIONS                                */
/************************************************/

void MultiSequence::WriteMultiSequenceToFile(std::string filename) {
    std::vector<std::string> vOut(multiSequence.size());
	std::ofstream output;
    int L = Length();
    int charPerLine = 100;
	output.open(filename);
	if (output.is_open()) {
        
        
        for (int k = 0; k < L / charPerLine; k++) {
            output << "\t" << std::setw(5) << std::left << (k * charPerLine) + 1 << std::setw(charPerLine - 5) << std::right << ((k + 1) * charPerLine) << std::endl;
            for (unsigned int i = 0; i < multiSequence.size(); i++) {
                output << multiSequence[i].getAccessionNum() << "\t";
                for (int j = (k * charPerLine); j < (k + 1) * charPerLine; j++) {
                    output << multiSequence[i][j];
                }
                output << std::endl;
            }
            output << std::endl;
        }
        
        
        int k = L / charPerLine;
        output << "\t" << std::setw(5) << std::left << (k * charPerLine) + 1 << std::setw(L % charPerLine - 5) << std::right << L << std::endl;
        for (unsigned int i = 0; i < multiSequence.size(); i++) {
            output << multiSequence[i].getAccessionNum() << "\t";
            for (int j = ((L / charPerLine) * charPerLine); j < L; j++) {
                output << multiSequence[i][j];
            }
            output << std::endl;
        }
        output << std::endl;
        
		output.close();
	}
	else {
		std::cout << "Error: Unable to open - " << filename << std::endl;
	}
}

int MultiSequence::numSequences() {
	return static_cast<int>(multiSequence.size());
}

int MultiSequence::Length() {
	return multiSequence[0].Length();
}

int MultiSequence::getSequenceIDs(int index) {
	return multiSequence[index].getAccessionNum();
}
