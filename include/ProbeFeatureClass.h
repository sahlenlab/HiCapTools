#ifndef PRSUITE_INC_PROMOTER_H_
#define PRSUITE_INC_PROMOTER_H_

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
	
	void GetTrFeats(std::stringstream&,temppars&, std::string);
    void ClusterIsoformPromoters(std::vector<int>, std::vector<std::string>, std::vector<std::string>, std::vector<int>&, std::vector<std::string>&, std::vector<std::string>&, std::string);
    int FindLeftMostCoord(std::vector<int>);
    int FindRightMostCoord(std::vector<int>);
	void DealwithSharedPromoters();
};
#endif // PRSUITE_INC_PROMOTER_H_
