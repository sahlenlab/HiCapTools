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
//  ProbeDataStructs.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROBDES_INC_DS_H_
#define HCT_PROBDES_INC_DS_H_

#include <vector>
#include <string>
#include <unordered_map>

namespace PrDes{
	
	struct REposStruct{
		std::vector< std::string > prom_indexes; //feat_id
		bool processed;
		bool whichside; //direction, [0 = upstream - left (L)] // [1 = downstream - right (R)]
	};

	struct REPosMap{ // This struct keeps the numbers of restriction fragments where there is at least one pair coming from a feature, one struct per chr
		std::unordered_map<int, REposStruct > repmap; //key:REsite, [1..n]:  // [1] : promindex
	};

	struct DesignLayers{
		std::vector < REPosMap  > Layer; // Probes that are close to each other will move to another design.
		//each element of this vector will be a new layer of design. Until there is no probes close to each other, new design layers will be created.
	};

	struct FeatureStruct{
		std::vector <std::string> ProbeID;
		int FeatureType; //Promoter, NegCtrl, etc..
		std::string chr;
		int nofRESites; // how many restriction sites the core promoter contains, will be used to normalise the signal  
		int closestREsitenums[2];  
		int start;
		int end;
	};

	struct PromoterStruct: FeatureStruct{
		std::vector < std::string > transcripts;
		std::vector < std::string > genes;
		std::string strand;
		int TSS;
		bool sharedpromoter; 
		std::vector < std::string > genes_sharingproms;
		bool probesSkip;
	//restriction fragments will be numbered with respect to the closest site upstream [0] or downstream [1].
	};


//For RE Sites
	struct REindexes{
		std::vector <int> binstart;
		std::vector <int> binend;
		std::vector <int> offset;
		std::vector <int> count;
	}; // each index struct will keep a chromosome

	struct Mappindexes{
		std::vector <int> binstart;
		std::vector <int> binend;
		std::vector <int> offset;
		std::vector <int> count;
	}; // each index struct will keep a chromosome

	struct RENFileInfo{
		std::string REName;
		std::string desName;
		int leftOfCut;
		int rightOfCut;
		char currTime[100];
		double mappabilityThreshold;
		int repeatOverlapExtent;
		bool ifRepeatAvail;
		bool ifMapAvail;
		std::string genomeAssembly;
	};
}

#endif // HCT_PROBDES_INC_DS_H_
