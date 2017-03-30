/*** 
   HiCapTools.
   Copyright (c) 2017 Pelin Sahlén <pelin.akan@scilifelab.se>

	Permission is hereby granted, free of charge, to any person obtaining a 
	copy of this software and associated documentation files (the "Software"), 
	to deal in the Software with some restriction, including without limitation 
	the rights to use, copy, modify, merge, publish, distribute the Software, 
	and to permit persons to whom the Software is furnished to do so, subject to
	the following conditions:

	The above copyright notice and this permission notice shall be included in all 
	copies or substantial portions of the Software. The Software shall not be used 
	for commercial purposes.

	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
	INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
	PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
	HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
	CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
	OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
***/

//
//  NegativeDesignProbes.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROBDES_INC_NEGPRDES_H_
#define HCT_PROBDES_INC_NEGPRDES_H_

#include "DesignProbes.h"

class NegativeProbeDesign:protected DesignClass{ //Negative probe design class
    friend class ProbeFeatureClass;
public:
   
	std::map< std::string, IntervalTree <std::string> > geneIntTree;
    std::map< std::string, IntervalTree <std::string> > exonIntTree;
    std::map< std::string, IntervalTree <std::string> > regulRegIntTree;
    
    std::vector<PrDes::FeatureStruct> intergenicPool;
    std::vector<PrDes::FeatureStruct> exonPool;
    std::vector<PrDes::FeatureStruct> intronPool;

    void InitialiseDesign(ProbeFeatureClass&, std::string, std::string, bool, int, std::string, std::string, int, PrDes::RENFileInfo&, int, std::string );
    int ConstructNegativeControlProbes(int, std::string,  Repeats, int, int);
    void WritetoFile(bioioMod&);
    NegativeProbeDesign(OutStream& dlog) : DesignClass (dlog) {}
    
private:
	struct negProbe{
		int start;
		int end;
		std::string name;
		negProbe(int s, int e, std::string n): start(s), end(e), name(n) {}
	};
	std::map <std::string, std::vector<int>> chrToIndex;
	
	std::vector<negProbe> toWriteSorted;
	
	int minREfragLen, distforbidIntergenic, distforbidReg;
	bool ifExistRegRegionFile;
    std::string fileName, summaryFileName, designName, genAssem, fasFileName;
    bool ifRep, ifMap;
    
    void ConstructPools(std::string, ProbeFeatureClass&);
    std::string FindOverlaps( std::map< std::string, IntervalTree <std::string> >& , std::string, unsigned long int, unsigned long int);
    void chooseRandomProbesFromPool(int, std::vector<PrDes::FeatureStruct>&, Repeats&, std::ofstream &, std::string, std::ofstream &);
    bool CheckRepeatOverlaps(std::string, int&, bool, Repeats&);
    
};

#endif //HCT_PROBDES_INC_NEGPRDES_H_
