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
//  Probes.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROXDET_INC_PROBES_H_
#define HCT_PROXDET_INC_PROBES_H_

#include <map>
#include "Features.h"
#include "Data_Structures.h"

class ProbeSet{
    friend class FeatureClass;
public:
   
	void ReadProbeCoordinates(std::string, std::map <std::string, std::string>&, int, bool, PrDes::RENFileInfo&);
    
	int FindOverlaps(std::string, unsigned long int, unsigned long int, std::string);
    int FindOverlaps_NegCtrls(std::string, unsigned long int, unsigned long int, std::string);
    
    ProbeSet(OutStream& prlog, int fCount, int fRCount) : prLog (prlog), fileCount(fCount), filesReadCount (fRCount) {}
    
private:
	OutStream& prLog;
	std::map< std::string, std::string> probeTypeMap;
	int fileCount;
	int filesReadCount;
	
	void GetProbeFeats(std::stringstream&, CaptureProbes&, std::string&);
    int FindClosestFeature(int, std::vector<int>&, std::string);
    int AddTotheIntervalTree(std::vector<Interval< int > >&, std::vector<Interval< int > >&, std::string, std::string, bool); 
    void ProcessProbeLine(std::map< std::string, std::vector < std::string > >::iterator, std::vector < int>&, CaptureProbes, int&, std::map< std::string, std::vector<Interval<int>>>&, std::map< std::string, std::vector<Interval<int>>>&, int);
   
};


#endif //HCT_PROXDET_INC_PROBES_H_
