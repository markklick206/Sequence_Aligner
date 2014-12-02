// Class: Sequence Header
// Forrest Ireland
// Dec. 1, 2014

/*
*   Implements a Sequence class
*
*   Holds information about a sequence, such as its sequence, organism, accession numbers, etc.
*/

#include <vector>
#include <string>

class Sequence {
public:
    Sequence();
    Sequence(const Sequence &s);
    ~Sequence();
    
	/************************************************/
	/* INPUT FUNCTIONS                                */
	/************************************************/
    
    bool setSequenceFromTextFile(std::string filename);
    bool setSequenceFromFASTAFile(std::string filename);
    void setSequence(std::vector<char> &s);
    void setSequence(std::string &s);
    
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