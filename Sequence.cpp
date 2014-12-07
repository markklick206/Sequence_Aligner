// Class: Sequence Code
// Forrest Ireland
// Dec. 1, 2014

#include "Sequence.h"

Sequence::Sequence() {
    
}

Sequence::Sequence(const Sequence &s) {
    this->sequence = s.sequence;
    this->organismName = s.organismName;
    this->accessionNum = s.accessionNum;
    this->file = s.file;
    this->length = s.length;
}

Sequence::~Sequence() {
    
}

/************************************************/
/* OPERATOR FUNCTIONS                           */
/************************************************/

char Sequence::operator[](int i) {
	return sequence[i];
}

void Sequence::push_back(char c) {
	sequence.push_back(c);
	this->length++;
}

/************************************************/
/* INPUT FUNCTIONS                              */
/************************************************/

bool Sequence::setSequenceFromTextFile(std::string filename) {
    std::ifstream input;
	std::string str;
    input.open(filename);
    if (input.is_open()) {
        getline(input, (str));
		for (unsigned int i = 0; i < str.size(); i++) {
			sequence.push_back(str[i]);
		}
		length = sequence.size();
        return true;
    }
	else {
        std::cout << "Error: Unable to open - " << filename << std::endl;
		return false;
    }
	return false;
}

bool Sequence::setSequenceFromFASTAFile(std::string filename) {
	std::ifstream fasta;
	std::string line, s;
	std::vector<char> seq;
	fasta.open(filename);
	if (fasta.is_open()) {
		getline(fasta, line);
		while (!fasta.eof()) {
			getline(fasta, line);
			s.append(line);
		}
		for (unsigned int i = 0; i < s.size(); i++) {
			sequence.push_back(s[i]);
		}
		return true;
	}
	else {
        std::cout << "Error: Unable to open - " << filename << std::endl;
		return false;
	}
    return false;
}

void Sequence::setSequence(std::vector<char> &s) {
    sequence = s;
}

void Sequence::setSequence(std::string &s) {
    for (unsigned int i = 0; i < s.size(); i++)
        sequence.push_back(s[i]);
}


void Sequence::setAccessionNum(int aNum) {
	accessionNum = aNum;
}

/************************************************/
/* GET FUNCTIONS                                */
/************************************************/

void Sequence::getSequence(std::vector<char> &R) {
    R = sequence;
}

int Sequence::getAccessionNum() {
    return accessionNum;
}

std::string Sequence::getOrganismName() {
    return organismName;
}

std::string Sequence::getFilename() {
    return file;
}

int Sequence::Length() {
    return sequence.size();
}
