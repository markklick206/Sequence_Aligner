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
		std::string error = "Error: Unable to open - ";
		error.append(filename);
		std::wstring temp = s2ws(error);
		LPCWSTR err = temp.c_str();
		MessageBox(NULL, err, L"Error", MB_OK);
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
		std::string error = "Error: Unable to open - ";
		error.append(filename);
		std::wstring temp = s2ws(error);
		LPCWSTR err = temp.c_str();
		MessageBox(NULL, err, L"Error", MB_OK);
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

/************************************************/
/* GET FUNCTIONS                                */
/************************************************/

void Sequence::getSequence(std::vector<char> &R) {
    R = sequence;
}

std::string Sequence::getAccessionNum() {
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

std::wstring Sequence::s2ws(const std::string &s)
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