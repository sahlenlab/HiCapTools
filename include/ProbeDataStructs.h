#ifndef PRSUITE_INC_DS_H_
#define PRSUITE_INC_DS_H_

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
		//int *closestREsitenums;  
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
		bool ifRepeatAvail;
		bool ifMapAvail;
		std::string genomeAssembly;
	};
}

#endif // PRSUITE_INC_DS_H_
