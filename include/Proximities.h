#ifndef INTSUITE_INC_PROX_H_
#define INTSUITE_INC_PROX_H_

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

#endif //INTSUITE_INC_PROX_H_
