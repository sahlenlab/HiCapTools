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
//  Find_Interactions.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PROXDET_INC_FI_H_
#define HCT_PROXDET_INC_FI_H_

#include "BackgroundInteractionFrequency.h"

class DetectInteractions{
public:
	
	void CalculatePvalAndPrintInteractionsProbeDistal(ProbeSet&, std::vector<DetermineBackgroundLevels>, std::string, int, std::vector < std::string >&, std::string, int, PrDes::RENFileInfo&);
	void CalculatePvalAndPrintInteractionsProbeProbe(ProbeSet&, std::vector<DetermineBackgroundLevels>, std::string, int,std::vector < std::string >&, std::string, int, PrDes::RENFileInfo&);
	
	void CalculatePvalAndPrintInteractionsProbeDistal_NegCtrls(ProbeSet&, std::vector<DetermineBackgroundLevels>, std::string, int, std::vector < std::string >&, std::string, int, PrDes::RENFileInfo&);
	void CalculatePvalAndPrintInteractionsProbeProbe_NegCtrls(ProbeSet&, std::vector<DetermineBackgroundLevels>, std::string, int, std::vector < std::string >&, std::string, int, PrDes::RENFileInfo&);

	DetectInteractions(OutStream& flog, int minNSuppPair, bool p_val, int minJDist) : fLog (flog), MinNumberofSupportingPairs (minNSuppPair), CALCULATE_P_VALUES(p_val), MinimumJunctionDistance (minJDist) {}

private:
	int MinNumberofSupportingPairs;
	bool CALCULATE_P_VALUES;
	int MinimumJunctionDistance;
	OutStream& fLog;
	
	bool CheckSupportingPairs(int*, int);
	double CalculatepVal(std::map< int, double >,std::map< int, double >, int,int);
    int FindClosestTranscriptTSS(int, std::vector<int>);
};

typedef struct{
    std::string ExperimentName;
    ProbeSet probes;
    DetermineBackgroundLevels background;
    DetectInteractions interactions;
    int ExperimentNo;
    
}intparams;


#endif // HCT_PROXDET_INC_FI_H_
