#ifndef INTSUITE_INC_DS_H_
#define INTSUITE_INC_DS_H_

#include <map>
#include <vector>
#include "IntervalTree.h"



struct Junction{
    int* paircount; // [0..n]: number of times observed for experiments 0..n,
    int* strandcombination; //count of forward-forward = [0], forward-reverse = [1], reverse-forward = [2], reverse-reverse = [3] //calculate index as ExperimentNo*4 + strcomb
    std::vector< double > p_val;
    int refragend;
    bool reportit;
    int distance;
};

struct SignalStruct{
    std::map<int, Junction > junctions; //key:RE frag start, [0]:REfrag end [1..n]: number of times observed for experiments 1..n, before: frag starts before
};

struct SignalStruct_CTX{ // This struct keeps the numbers restriction fragments where there is at least one pair coming from a feature, one struct per chr
	std::map<int, Junction > junctions_ctx;
	std::string maptochrname;
};

struct FeattoFeatSignalStruct{
    std::string interacting_feature_id; // chr_start
    std::vector < int > signal; // The number of reads within the core promoter of the interactor promoter
    int *strandcombination; //count of forward-forward = [0], forward-reverse = [1], reverse-forward = [2], reverse-reverse = [3]
};

struct FeatureStruct{
    std::string feature_id; // Index of the feature it belongs in FeatureStruct "chr_start"
	std::string probe_name;
    std::string probe_target; // the probe  targeting this feature will be annotated with this value. Read from transcript file and compared with value from probe file
	int FeatureType; // 0 = not annotated, 1 = annotated with a promoter, 2 = annotated with a cad_gwas, 3 = annotated with a negative control
	
	std::string chr;
	std::string Name;
	std::string TranscriptName; //NM..
	int start; //Transcription Start Site in case of genes
    int end;
	int *closestREsitenums; // offset of the closest RE sites on each side of the negctrl this will be used to index interactions
	std::string strand;

 	SignalStruct proximities; // Intra-chromosomal interactions
    std::vector < SignalStruct_CTX > proximities_ctx; // Each element of this vector will represent a chromosome
    std::vector < FeattoFeatSignalStruct > Inter_feature_ints; // Each element represent a promoter

};


//DECLARE FEATURES
extern std::map< std::string, FeatureStruct > Features; //key = featureid, (chr_start); value = feature struct
extern std::map< std::string, std::vector < std::string > > MetaFeatures; // key = name (A1BG, rs1101 etc); value = list of feature_id (e.g. chr1_199202)

struct CaptureProbes{
	std::string name_of_design;
    std::string chr;
    std::vector < std::string > probe_feats;
    std::string target_id;
    std::string side; //which side of the feature it is placed, left or right -> left is upstream design, right is downstream design
    int repeat_overlap;
    double mappability;
    int start;
    int end;
    int closestTSS_index;
    int annotated;// 0 = not annotated, 1 = annotated with a promoter, 2 = annotated with a cad_gwas, 3 = annotated with a negative control
    bool conflicting_annotations; // if a probe is associated with more than one feature
    std::string feature_id; // Index of the feature it belongs in FeatureStruct "chr_start"
};

struct Probe_Design{
    std::vector < CaptureProbes > Probes;
    std::map< std::string, IntervalTree < int > > Probe_Tree;  //chr, probes
};

//DECLARE PROBES
extern std::map < std::string, Probe_Design > Design;
extern std::map < std::string, Probe_Design > Design_NegCtrl;


//For RE Sites
struct REindexes{
	std::vector <int> binstart;
	std::vector <int> binend;
	std::vector <int> offset;
	std::vector <int> count;
}; // each index struct will keep a chromosome

struct IndexMap{
    int ChrRowStartIndex;
    int ChrRowEndIndex;
    std::string chrname;
};

struct BG_signals{

	std::map< int, double > mean;
    std::map< int, double > smoothed;
	std::map< int, double > stdev;
    std::map< int, double > smoothed_stdev;
    std::map< int, int > samplesize;
	
	double a; 
	double b; 
};

#endif  // INTSUITE_INC_DS_H_
