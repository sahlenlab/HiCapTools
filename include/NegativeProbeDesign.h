//
//  DesignProbes.h
//  HiCap_ProbeDesigner
//
//  Created by Pelin Sahlen on 12/03/2015.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
//
#ifndef PRSUITE_INC_NEGPRDES_H_
#define PRSUITE_INC_NEGPRDES_H_

#include "DesignProbes.h"

class NegativeProbeDesign:protected DesignClass{ //Negative probe design class
    friend class ProbeFeatureClass;
public:
   
	std::map< std::string, IntervalTree <std::string> > geneIntTree;
    std::map< std::string, IntervalTree <std::string> > exonIntTree;
    std::map< std::string, IntervalTree <std::string> > regulRegIntTree;
    
    std::vector<PrDes::FeatureStruct> intergenicPool;
    std::vector<PrDes::FeatureStruct> exonPool;
    std::vector<PrDes::FeatureStruct> intronPool;

    void InitialiseDesign(ProbeFeatureClass&, std::string, std::string, bool, int, std::string, std::string, int, PrDes::RENFileInfo&, int, std::string );
    int ConstructNegativeControlProbes(int, std::string,  Repeats, int, int);
    void WritetoFile();
    NegativeProbeDesign(OutStream& dlog) : DesignClass (dlog) {}
    
private:
	struct negProbe{
		int start;
		int end;
		std::string name;
		negProbe(int s, int e, std::string n): start(s), end(e), name(n) {}
	};
	std::map <std::string, std::vector<int>> chrToIndex;
	
	std::vector<negProbe> toWriteSorted;
	
	int minREfragLen, distforbidIntergenic, distforbidReg;
	bool ifExistRegRegionFile;
    std::string fileName, summaryFileName, designName, genAssem;
    bool ifRep, ifMap;
    
    void ConstructPools(std::string, ProbeFeatureClass&);
    std::string FindOverlaps( std::map< std::string, IntervalTree <std::string> >& , std::string, unsigned long int, unsigned long int);
    void chooseRandomProbesFromPool(int, std::vector<PrDes::FeatureStruct>&, Repeats&, std::ofstream &, std::string, std::ofstream &);
    bool CheckRepeatOverlaps(std::string, int&, bool, Repeats&);
   // void WritetoFile(std::ofstream &, std::string, int, int, std::string, int);
    
};

#endif //PRSUITE_INC_NEGPRDES_H_
