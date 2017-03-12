#ifndef PRSUITE_INC_RE_H_
#define PRSUITE_INC_RE_H_

#include "ProbeDataStructs.h"
#include "OutStream.h"

class RESitesClass{
public :
	int span;
	std::vector <std::string> chr_names;
	std::unordered_map< std::string, int > chroffsets_indexfile;
	std::unordered_map<std::string, int > chr_starts; // first RE site
	std::unordered_map<std::string, int > chr_ends; // last RE site
	std::vector < PrDes::REindexes > indexes;
	std::vector < int > posvector;
	
	RESitesClass(OutStream& rlog) : rLog (rlog) {}
	
	void InitialiseVars(std::string);
	bool GettheREPositions(std::string, int, int*);
	void CleanClass();
	
	
	
private:

	OutStream& rLog;
};


#endif // PRSUITE_INC_RE_H_
