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
//  DesignProbes.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROBDES_INC_DESPR_H_
#define HCT_PROBDES_INC_DESPR_H_

#include "ProbeFeatureClass.h"
#include "RepeatOverlaps.h"
#include "bioioMod.h"

class DesignClass{ //Design probe class
    friend class ProbeFeatureClass;
public:
    std::unordered_map< std::string, int > chrIndex; // key = chrname, value = index of REfragments
    std::vector < PrDes::DesignLayers > OneDesign; // TSS is on plus strand
    
    void InitialiseDesign(ProbeFeatureClass&, std::vector< PrDes::REPosMap >&);
    void DesignProbes(ProbeFeatureClass&, RESitesClass&, Repeats&, bioioMod&, std::string, std::string, std::string, int, int, PrDes::RENFileInfo&, int, int, int);
    void MergeAllChrOutputs(ProbeFeatureClass&, PrDes::RENFileInfo&);
    bool ConstructSeq(PrDes::RENFileInfo&, bioioMod&, std::string);
    DesignClass(OutStream& dlog) : dLog (dlog) {}
    
protected:
	int ProbeLen;
	std::string bigwigsummarybinary, mappabilityfile;
	int reLeftCut, reRightCut;
	OutStream& dLog;
	double mapThreshold;
	int repOverlapExtent;
	int BUFSIZE;
	int MaxDistancetoTSS;
	int minDistToTSS;
	struct BaseProbe{
		std::string chr;
		int start;
		int end;
		std::string strand;
		std::string feature;
	};
	std::vector<BaseProbe> probeList;
	
	bool CheckFragmentSize(RESitesClass &, std::string, int, int);
    double BigWigSummary(std::string, int, int);   
    bool CheckRESite(RESitesClass&, bioioMod&,std::string, int, int&, bool, Repeats&, bool, bool, bool);
    bool createNewEntry(std::unordered_map<int, PrDes::REposStruct >&, std::unordered_map<int, PrDes::REposStruct >&, int, std::string, int,bool);
    int CheckREsiteAroundProbe(RESitesClass&, std::string, int, int);
    bool WritetoFile(std::ofstream&, std::string, int, int, std::vector< std::string >, bool, std::string, ProbeFeatureClass&);
    int CheckDistanceofProbetoTSS(RESitesClass&, std::string, int, int, int);
    bool CheckRepeats(Repeats&, std::string, int, int, bool);
    bool CheckMappability(std::string, int, int, bool);
	bool CheckCountofNs(bioioMod&, std::string, int, int,bool);
    
    
};

#endif //HCT_PROBDES_INC_DESPR_H_
