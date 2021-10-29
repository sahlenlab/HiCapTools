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
//  ProcessPairs.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "ProcessPairs.h"
#include "Global.h"

#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
#include "api/BamIndex.h"
using namespace BamTools;

void ProcessBAM::Initialize(std::string bamfilename, int nOfExp, int padd, int readLen){
    
    NOFEXPERIMENTS = nOfExp;
    padding=padd;
    ReadLen=readLen;
    BamReader reader;
    
    if ( !reader.Open(bamfilename.c_str()) )
        std::cerr << "Could not open input BAM file." << std::endl;
    bLog << bamfilename << "  opened" << std::endl;
    
    // retrieve 'metadata' from BAM files, these are required by BamWriter
    const RefVector references = reader.GetReferenceData();
    
    // Make a map of chr names to RefIDs
    RefVector::const_iterator chrit;
    for (chrit = references.begin(); chrit != references.end(); ++chrit){
        int key = reader.GetReferenceID(chrit->RefName);
        RefIDtoChrNames[key] = chrit->RefName;
        RefIDtoChrLength[key]= chrit->RefLength;
    }
    
    // Make a map of RefIDs to chrnames
    for (chrit = references.begin(); chrit != references.end(); ++chrit){
        int key = reader.GetReferenceID(chrit->RefName);
        ChrNamestoRefID[chrit->RefName] = key;
    }
    
    reader.Close();
}


void ProcessBAM::ProcessSortedBamFile_NegCtrls(ProbeSet& ProbeClass, RESitesClass& dpnII, ProximityClass& proximities, std::string BAMFILENAME, int ExperimentNo, std::string DesignName, std::string StatsOption, FeatureClass& proms){
    
    BamAlignment al, almate;
    BamReader reader;
    BamRegion probeRegion;
    int probeRefID;
    int invalidRECoord=0;
      
    if ( !reader.Open(BAMFILENAME.c_str()) )
        std::cerr << "Could not open input BAM file." << std::endl;
    else{
        bLog << ExperimentNo+1 << " out of " << NOFEXPERIMENTS << " will be read to create negative controls" << std::endl;

	if(!(reader.LocateIndex())){
                reader.CreateIndex();
                bLog<<"Index created, is located:"<<reader.LocateIndex()<<std::endl;
        }
     }
    
    Alignment pairinfo;
    int strandcomb = 0, sc_index=0;
    bool read1_strand, read2_strand;
    std::string feature_id1, feature_id2;
    bool re1f, re1r;
    pairinfo.resites1 = new int [2];
    pairinfo.resites2 = new int [2];
    
    while (reader.GetNextAlignment(al)) {

       totalNumberofPairs=totalNumberofPairs+1;
    }
    
	
	for(int i = 0; i < Design_NegCtrl[DesignName].Probes.size(); i++){
		probeRefID = ChrNamestoRefID[Design_NegCtrl[DesignName].Probes[i].chr];
		probeRegion.LeftRefID = probeRefID;
		probeRegion.RightRefID = probeRefID;
	
		//left and right
		if(Design_NegCtrl[DesignName].Probes[i].side=="L"){
			probeRegion.LeftPosition =  Design_NegCtrl[DesignName].Probes[i].start;
			probeRegion.RightPosition =  Design_NegCtrl[DesignName].Probes[i].end + padding;
		}
		else if(Design_NegCtrl[DesignName].Probes[i].side=="R"){
			probeRegion.LeftPosition =  Design_NegCtrl[DesignName].Probes[i].start- padding;;
			probeRegion.RightPosition =  Design_NegCtrl[DesignName].Probes[i].end;
		}
		
		reader.SetRegion(probeRegion);

		while(reader.GetNextAlignmentCore(al)){
			++NumberofPairs;
			pairinfo.chr1 = RefIDtoChrNames[al.RefID];
			pairinfo.chr2 = RefIDtoChrNames[al.MateRefID];
			pairinfo.probeid1 = i;
			pairinfo.probeid2 = ProbeClass.FindOverlaps_NegCtrls(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen), DesignName);
			feature_id1 = Design_NegCtrl[DesignName].Probes[pairinfo.probeid1].feature_id;
                 //Find the closest RE
            re1r = dpnII.GettheREPositions(pairinfo.chr2, (al.MatePosition), pairinfo.resites2, invalidRECoord);
            if ( !re1r){
                pairinfo.resites2[0] = al.MatePosition;
                pairinfo.resites2[1] = al.MatePosition;
            }
     
                
           if((pairinfo.probeid2 != -1)){
			   feature_id2 = Design_NegCtrl[DesignName].Probes[pairinfo.probeid2].feature_id;
			   ++NofPairs_Both_on_Probe;
               }
               else{
				   std::string checkFeat = proms.FindOverlaps(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen));
				   if(checkFeat!="null"){
					   feature_id2 = checkFeat;
						++NofPairs_Both_on_Probe;
					}
					else{
						feature_id2 = "null";
						++NofPairs_One_on_Probe;
					}
			   }
			   
			   if(StatsOption!="ComputeStatsOnly"){
				   read1_strand = al.IsReverseStrand();
				   read2_strand = al.IsMateReverseStrand();
				   if(read1_strand == 0){
					   if(read2_strand == 0)
						   strandcomb = 0; // forward-forward
						else
							strandcomb = 1; //forward-reverse
					}
					else{
						if(read2_strand == 0)
							strandcomb = 2; // reverse-forward
						if(read2_strand == 1)
							strandcomb = 3; // reverse-reverse
					}
					sc_index = (ExperimentNo*4) + strandcomb;
					proximities.RecordProximities(pairinfo, feature_id1, feature_id2, sc_index, ExperimentNo);
				}
        if(NumberofPairs % 20000000 == 0)
           bLog << NumberofPairs/1000000 << " million pairs are processed  " << BAMFILENAME << std::endl;  
    }
}
	NofPairsNoAnn = totalNumberofPairs/2 - (NumberofPairs);
	
	if(invalidRECoord>0)
		bLog<<"!!Warning!! Some invalid coordinates were encountered with coordinates out of chromosome boundaries and were therefore skipped. Number of Skipped coordinates:"<<invalidRECoord<<std::endl;
    bLog << BAMFILENAME <<  "    BAM file finished" << std::endl;
}



void ProcessBAM::ProcessSortedBAMFile(ProbeSet& ProbeClass, RESitesClass& dpnII, ProximityClass& proximities, std::string BAMFILENAME, int ExperimentNo, std::string whichchr, std::string DesignName, std::string StatsOption, FeatureClass& proms){
   
	BamAlignment al, almate;
	BamReader reader;
	BamRegion probeRegion;
	int probeRefID;
	int invalidRECoord=0;

	std::string indexFilename;
    
	indexFilename = BAMFILENAME;
    
	if ( !reader.Open(BAMFILENAME.c_str()) )
		std::cerr << "Could not open input BAM file." << std::endl;
	else{
		bLog << ExperimentNo+1 << " out of " << NOFEXPERIMENTS << " will be read" << std::endl;
		if(!(reader.LocateIndex())){
                reader.CreateIndex();
                bLog<<"Index created, is located:"<<reader.LocateIndex()<<std::endl;
        }
	}

	Alignment pairinfo;
	int strandcomb = 0, sc_index=0;
	std::string feature_id1, feature_id2;
	bool read1_strand, read2_strand;
	bool re1f, re1r;
	pairinfo.resites1 = new int [2];
	pairinfo.resites2 = new int [2];

	std::map<int, std::string>::const_iterator it1;
	int chrindex = 0;
	if(whichchr.find("All")==std::string::npos){
		for(it1 = RefIDtoChrNames.begin(); it1 != RefIDtoChrNames.end(); ++it1){
			if(it1->second == whichchr){
				chrindex = it1->first;
				break;
			}
		}
	}

	while (reader.GetNextAlignment(al)) {

        totalNumberofPairs=totalNumberofPairs+1;
    }
    
	
	for(int i = 0; i < Design[DesignName].Probes.size(); i++){
		probeRefID = ChrNamestoRefID[Design[DesignName].Probes[i].chr];
		probeRegion.LeftRefID = probeRefID;
		probeRegion.RightRefID = probeRefID;
		if(Design[DesignName].Probes[i].side=="L"){
			probeRegion.LeftPosition =  Design[DesignName].Probes[i].start;
			probeRegion.RightPosition =  Design[DesignName].Probes[i].end + padding;
		}
		else if(Design[DesignName].Probes[i].side=="R"){
			probeRegion.LeftPosition =  Design[DesignName].Probes[i].start- padding;;
			probeRegion.RightPosition =  Design[DesignName].Probes[i].end;
		}
		
		reader.SetRegion(probeRegion);
		while(reader.GetNextAlignmentCore(al)){
		
			++NumberofPairs;
		
			if(whichchr.find("All")==std::string::npos){
				if(al.RefID == chrindex){
					pairinfo.chr1 = RefIDtoChrNames[al.RefID];
					pairinfo.chr2 = RefIDtoChrNames[al.MateRefID];
					pairinfo.probeid1 = i;
					pairinfo.probeid2 = ProbeClass.FindOverlaps(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen), DesignName);
					feature_id1 = Design[DesignName].Probes[pairinfo.probeid1].feature_id;
        
					re1r = dpnII.GettheREPositions(pairinfo.chr2, (al.MatePosition), pairinfo.resites2, invalidRECoord);
					if ( !re1r){
						pairinfo.resites2[0] = al.MatePosition;
						pairinfo.resites2[1] = al.MatePosition;
					}
                
					if(( pairinfo.probeid2 != -1)){
                
						feature_id2 = Design[DesignName].Probes[pairinfo.probeid2].feature_id;
						++NofPairs_Both_on_Probe;
					}
					else{
						std::string checkFeat = proms.FindOverlaps(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen));
						if(checkFeat!="null"){
							feature_id2 = checkFeat;
							++NofPairs_Both_on_Probe;
						}
						else{
							feature_id2 = "null";
							++NofPairs_One_on_Probe;
						}
					}
            
					if(StatsOption!="ComputeStatsOnly"){
						read1_strand = al.IsReverseStrand();
						read2_strand = al.IsMateReverseStrand();
						if(read1_strand == 0){
							if(read2_strand == 0)
								strandcomb = 0; // forward-forward
							else
								strandcomb = 1; //forward-reverse
						}
						else{
							if(read2_strand == 0)
								strandcomb = 2; // reverse-forward
							if(read2_strand == 1)
								strandcomb = 3; // reverse-reverse
						}
						sc_index = (ExperimentNo*4) + strandcomb;
						proximities.RecordProximities(pairinfo, feature_id1, feature_id2, sc_index, ExperimentNo);
					}
       
					if(NumberofPairs % 20000000 == 0)
						bLog << NumberofPairs/1000000 << " million pairs are processed  " << BAMFILENAME << std::endl;
				}
			}
			else{
			
				pairinfo.chr1 = RefIDtoChrNames[al.RefID];
				pairinfo.chr2 = RefIDtoChrNames[al.MateRefID];
				pairinfo.probeid1 = i;
				pairinfo.probeid2 = ProbeClass.FindOverlaps(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen), DesignName);
				feature_id1 = Design[DesignName].Probes[pairinfo.probeid1].feature_id;
        
            	re1r = dpnII.GettheREPositions(pairinfo.chr2, (al.MatePosition), pairinfo.resites2, invalidRECoord);
				if ( !re1r){
					pairinfo.resites2[0] = al.MatePosition;
					pairinfo.resites2[1] = al.MatePosition;
				}
        		if(( pairinfo.probeid2 != -1)){
                	feature_id2 = Design[DesignName].Probes[pairinfo.probeid2].feature_id;
					++NofPairs_Both_on_Probe;
				}
				else{
					std::string checkFeat = proms.FindOverlaps(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen));
					if(checkFeat!="null"){
							feature_id2 = checkFeat;
							++NofPairs_Both_on_Probe;
					}
					else{
						feature_id2 = "null";
						++NofPairs_One_on_Probe;
					}
				}
            
				if(StatsOption!="ComputeStatsOnly"){
					read1_strand = al.IsReverseStrand();
					read2_strand = al.IsMateReverseStrand();
					if(read1_strand == 0){
						if(read2_strand == 0)
							strandcomb = 0; // forward-forward
						else
							strandcomb = 1; //forward-reverse
					}
					else{
						if(read2_strand == 0)
							strandcomb = 2; // reverse-forward
						if(read2_strand == 1)
							strandcomb = 3; // reverse-reverse
					}
					sc_index = (ExperimentNo*4) + strandcomb;
					proximities.RecordProximities(pairinfo, feature_id1, feature_id2, sc_index, ExperimentNo);
				}
       
				if(NumberofPairs % 20000000 == 0)
					bLog << NumberofPairs/1000000 << " million pairs are processed  " << BAMFILENAME << std::endl;
			
			}
		}
	}
	//------------------
	if(StatsOption=="ComputeStatsOnly"){
		
		for(int i = 0; i < Design_NegCtrl[DesignName].Probes.size(); i++){
		probeRefID = ChrNamestoRefID[Design_NegCtrl[DesignName].Probes[i].chr];
		probeRegion.LeftRefID = probeRefID;
		probeRegion.RightRefID = probeRefID;
	
		//left and right
		if(Design_NegCtrl[DesignName].Probes[i].side=="L"){
			probeRegion.LeftPosition =  Design_NegCtrl[DesignName].Probes[i].start;
			probeRegion.RightPosition =  Design_NegCtrl[DesignName].Probes[i].end + padding;
		}
		else if(Design_NegCtrl[DesignName].Probes[i].side=="R"){
			probeRegion.LeftPosition =  Design_NegCtrl[DesignName].Probes[i].start- padding;;
			probeRegion.RightPosition =  Design_NegCtrl[DesignName].Probes[i].end;
		}
		
		reader.SetRegion(probeRegion);

		while(reader.GetNextAlignmentCore(al)){
			++NumberofPairs;
			pairinfo.chr1 = RefIDtoChrNames[al.RefID];
			pairinfo.chr2 = RefIDtoChrNames[al.MateRefID];
			pairinfo.probeid1 = i;
			pairinfo.probeid2 = ProbeClass.FindOverlaps_NegCtrls(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen), DesignName);
			feature_id1 = Design_NegCtrl[DesignName].Probes[pairinfo.probeid1].feature_id;
                 //Find the closest RE
            re1r = dpnII.GettheREPositions(pairinfo.chr2, (al.MatePosition), pairinfo.resites2, invalidRECoord);
            if ( !re1r){
                pairinfo.resites2[0] = al.MatePosition;
                pairinfo.resites2[1] = al.MatePosition;
            }
     
                
           if((pairinfo.probeid2 != -1)){
			   feature_id2 = Design_NegCtrl[DesignName].Probes[pairinfo.probeid2].feature_id;
			   ++NofPairs_Both_on_Probe;
               }
               else{
				   std::string checkFeat = proms.FindOverlaps(pairinfo.chr2, al.MatePosition, (al.MatePosition + ReadLen));
				   if(checkFeat!="null"){
					   feature_id2 = checkFeat;
						++NofPairs_Both_on_Probe;
					}
					else{
						feature_id2 = "null";
						++NofPairs_One_on_Probe;
					}
			   }
        if(NumberofPairs % 20000000 == 0)
           bLog << NumberofPairs/1000000 << " million pairs are processed  " << BAMFILENAME << std::endl;  
		}
	}
	//=-------------------------------------------	
	}
	if(invalidRECoord>0)
		bLog<<"!!Warning!! Some invalid coordinates were encountered with coordinates out of chromosome boundaries and were therefore skipped. Number of Skipped coordinates:"<<invalidRECoord<<std::endl;
	bLog << BAMFILENAME <<  "    BAM file finished" << std::endl;
	NofPairsNoAnn = totalNumberofPairs/2 - (NumberofPairs);
}
