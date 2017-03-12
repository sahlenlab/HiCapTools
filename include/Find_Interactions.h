#ifndef INTSUITE_INC_FI_H_
#define INTSUITE_INC_FI_H_

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


#endif // INTSUITE_INC_FI_H_
