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
//  Proximities.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROXDET_INC_PROX_H_
#define HCT_PROXDET_INC_PROX_H_

#include "Probes.h"

typedef struct{
    std::string chr1, chr2;
    int probeid1, probeid2;
    int *resites1, *resites2;
}Alignment;

typedef struct{
    std::string chr1, chr2;
    std::string enid1, enid2;
    int *resites1, *resites2;
}enAlignment;


class ProximityClass{
	friend class ProbeSet;
    friend class ProcessBAM;
    friend class FeatureClass;

public:

	void RecordProximities(Alignment, std::string, std::string, int);
	void AnnotateDistalInteractor(std::string, std::string, std::string, int*, int);
	void AnnotateFeatFeatInteraction(std::string, std::string, int);
	void PopulateInteractions(std::map<int, Junction >&, int*, int);
	void CountProximities(ProbeSet, int);
	ProximityClass(int nOfExp) : NOFEXPERIMENTS (nOfExp) {}

private:

	bool chrfound, foundbefore;
	int n;
	int NOFEXPERIMENTS;
	void createtable(std::vector< int >);
	

};

#endif //HCT_PROXDET_INC_PROX_H_
