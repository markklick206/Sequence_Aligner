// Class: Sequence Header
// Forrest Ireland
// Dec. 1, 2014

/*
*   Implements a Sequence class
*
*   Holds information about a sequence, such as its sequence, organism, accession numbers, etc.
*/

#pragma once

#include <vector>
#include <iostream>
#include <string>
#include <fstream>

class Sequence {
public:
    Sequence();
    Sequence(const Sequence &s);
    ~Sequence();

	/************************************************/
	/* OPERATOR FUNCTIONS                           */
	/************************************************/

	char operator[](int i);
	void push_back(char c);
    
	/************************************************/
	/* INPUT FUNCTIONS                              */
	/************************************************/
    
    bool setSequenceFromTextFile(std::string filename);
    bool setSequenceFromFASTAFile(std::string filename);
    void setSequence(std::vector<char> &s);
    void setSequence(std::string &s);

	void setOrganismName(std::string oName);
	void setAccessionNum(std::string aNum);
	void setFile(std::string file);
    
	/************************************************/
	/* GET FUNCTIONS                                */
	/************************************************/
    
    void getSequence(std::vector<char> &R);
    std::string getAccessionNum();
    std::string getOrganismName();
    std::string getFilename();
    int Length();
    
private:
    std::vector<char> sequence;
    std::string organismName;
    std::string accessionNum;
    std::string file;
    int length;
};