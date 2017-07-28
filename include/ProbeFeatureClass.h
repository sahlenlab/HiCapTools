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
//  ProbeFeatureClass.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROBDES_INC_PROMOTER_H_
#define HCT_PROBDES_INC_PROMOTER_H_

#include "RESitesCount.h"
#include "IntervalTree.h"
#include <map>

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

struct cCounter{
	int count;
	bool probesDesigned;
};

class ProbeFeatureClass{ //Probe Clusters Associated with a Promoter
public:
	int promPadding;//padding for Promoter Interval Tree
	int NofPromoters;
	
	std::map<std::string, PrDes::PromoterStruct> promFeatures; //featid map to 

	std::vector <std::string> ChrNames_proms;
	std::map< std::string, IntervalTree <std::string> > promIntTree; //Promoter Interval Tree

	void InitialiseData(int, int, int);
	void ReadFeatureAnnotation(RESitesClass&, std::string, std::string);
	ProbeFeatureClass(OutStream& plog) : pLog (plog) {}

private:
	int ClusterPromoters;
	std::map< std::string, std::vector < Interval < std::string > > > chrIntervals;
	int fileCount, fileReadCount;
	OutStream& pLog;
	std::map<std::string, cCounter> transcriptCounter;
	
	void GetTrFeats(std::stringstream&,temppars&, std::string);
    void ClusterIsoformPromoters(std::vector<int>, std::vector<std::string>, std::vector<std::string>, std::vector<int>&, std::vector<std::string>&, std::vector<std::string>&, std::string);
    int FindLeftMostCoord(std::vector<int>);
    int FindRightMostCoord(std::vector<int>);
	void DealwithSharedPromoters();
};
#endif // HCT_PROBDES_INC_PROMOTER_H_
