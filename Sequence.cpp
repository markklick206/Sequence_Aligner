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
/* INPUT FUNCTIONS                                */
/************************************************/

bool Sequence::setSequenceFromTextFile(std::string filename) {
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

bool Sequence::setSequenceFromFASTAFile(std::string filename) {
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

void Sequence::setSequence(std::vector<char> &s) {
    sequence = s;
}

void Sequence::setSequence(std::string &s) {
    for (unsigned int = 0; i < s.size(); i++)
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
    return length;
}


