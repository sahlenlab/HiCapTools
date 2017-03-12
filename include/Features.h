#ifndef INTSUITE_INC_PROM_H_
#define INTSUITE_INC_PROM_H_

#include <sstream>
#include "RESitesCount.h"

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
    void InitialiseData(int, int);
    
    void ReadFeatureAnnotation(RESitesClass&, std::string, std::string);
    FeatureClass(OutStream& plog) : pLog (plog) {}
       
private:
	OutStream& pLog;
	int ClusterPromoters; 
	int fileCount ;
	int filesReadCount;
	
    void GetTrFeats(std::stringstream&, temppars&, std::string);
    void ClusterIsoformPromoters(std::vector<int>&, std::vector<std::string>&, std::vector<int>&, std::vector<std::string>&, std::string, int);
    int FindLeftMostCoord(std::vector<int>);
    int FindRightMostCoord(std::vector<int>);
};

#endif //INTSUITE_INC_PROM_H_




