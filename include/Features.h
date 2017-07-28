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
//  Features.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROXDET_INC_PROM_H_
#define HCT_PROXDET_INC_PROM_H_

#include <sstream>
#include <map>
#include "RESitesCount.h"
#include "IntervalTree.h"

struct temppars{
    std::string chr;
    int start;
    int end;
    
    std::string strand;
    std::string name;
    std::string tr_id;
    std::string probe_id;
    std::string probetarget;
    int FeatureType;
};




class FeatureClass{ //Probe Clusters Associated with a Promoter
public:
    
    temppars *tp;
    void InitialiseData(int, int, int);
    
    std::map< std::string, IntervalTree <std::string> > promIntTree; //Promoter Interval Tree
    
    void ReadFeatureAnnotation(RESitesClass&, std::string, std::string);
    std::string FindOverlaps(std::string, unsigned long int, unsigned long int);
    
    FeatureClass(OutStream& plog) : pLog (plog) {}
       
private:
	OutStream& pLog;
	int ClusterPromoters; 
	int fileCount ;
	int filesReadCount;
	int fOverlapPad;
	std::map< std::string, std::vector < Interval < std::string > > > chrIntervals;
	
    void GetTrFeats(std::stringstream&, temppars&, std::string);
    void ClusterIsoformPromoters(std::vector<int>&, std::vector<std::string>&, std::vector<int>&, std::vector<std::string>&, std::string, int);
    int FindLeftMostCoord(std::vector<int>);
    int FindRightMostCoord(std::vector<int>);
};

#endif //HCT_PROXDET_INC_PROM_H_




