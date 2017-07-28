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
//  ProbeFeatureClass.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "ProbeFeatureClass.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

temppars *tp;


void ProbeFeatureClass::InitialiseData(int clustProm, int fCount, int pPadding){

	promPadding = pPadding, //padding for Promoter Interval Tree
	tp=new temppars [2];
	
	NofPromoters=0;
	
	fileCount=fCount;
	fileReadCount=0;
	
	ClusterPromoters = clustProm;

    pLog << "Promoter Class Initialised" << std::endl;

}


void ProbeFeatureClass::GetTrFeats(std::stringstream &trx, temppars &tpars, std::string option){
	
	std::string field, start, end;
    
    //Transcript Line Format
    //name2	 name	chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds
    //#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
    
    //BED detail Line Format for SNP
    //chrom	Start	End	name	ID	Description
	
	if(option=="transcript"){
		getline(trx,field,'\t'); 
		getline(trx,tpars.tr_id,'\t'); 
    }
    
    getline(trx,tpars.chr,'\t'); 
    
    if(option=="transcript"){
		getline(trx,tpars.strand,'\t');
    }
    
    getline(trx,start,'\t');
    getline(trx,end,'\t');
    
    if(option=="transcript"){
		getline(trx,field,'\t'); //cdsStart
		getline(trx,field,'\t'); //cdsEnd
		getline(trx,field,'\t'); //exonCount
		getline(trx,field,'\t'); //exonStarts
		getline(trx,field,'\t'); //exonEnds
		getline(trx,field,'\t'); //score
		
		getline(trx,tpars.name,'\t'); 
		
		getline(trx,field,'\t'); //cdsStartStat
		getline(trx,field,'\t'); //cdsEndStat
		getline(trx,field,'\t'); //exonFrames
    
		if(tpars.strand=="+"){
			tpars.start=std::stoi(start);
			tpars.end=std::stoi(end);
			tpars.probe_id=tpars.tr_id+"."+start;
		}
		else{
			tpars.start=std::stoi(end);
			tpars.end=std::stoi(start);
			tpars.probe_id=tpars.tr_id+"."+end;
		}
		 tpars.FeatureType = 1;
    }
    
   if(option=="SNV"){
	   
	   tpars.start=std::stoi(start);
	   tpars.end=std::stoi(end);
	   tpars.strand="+";
	   
	   getline(trx,tpars.name,'\t');
	   tpars.tr_id=tpars.name;
	   tpars.FeatureType = 2;
	   tpars.probe_id=tpars.tr_id+"."+start;
	   getline(trx,field,'\t');
	   getline(trx,tpars.tr_id,'\t');
   }  
}

int ProbeFeatureClass::FindLeftMostCoord(std::vector<int> coords){
    int leftmosttss;
    
    leftmosttss = coords[0];
    for (int i = 0; i < coords.size();++i){
        if (coords[i] <= leftmosttss )
            leftmosttss = coords[i];
    }
    return leftmosttss;
}

int ProbeFeatureClass::FindRightMostCoord(std::vector<int> coords){
    int rightmosttss;

    rightmosttss = coords[0];
    for (int i = 0; i < coords.size();++i){
        if (coords[i] >= rightmosttss )
            rightmosttss = coords[i];
    }
    return rightmosttss;

}

																																																										/////////////
void ProbeFeatureClass::ClusterIsoformPromoters(std::vector <int> isoformprs, std::vector<std::string> tr_ids, std::vector<std::string> probe_ids, std::vector <int>& clusteredtrcoords, std::vector<std::string>& clustered_ids, std::vector<std::string>& clustered_probe_ids, std::string strand){
    
    unsigned int j,k,l,cluster_idx=0;
    
    //cluster_idx
    
    std::vector<int> clustercoords;
    std::vector<std::string> clustered_tr_ids;
    std::vector<std::string> clustered_pr_ids;
    
    std::vector<bool> clustered;
    
    for(j = 0;j < isoformprs.size(); ++j)
        clustered.push_back(0);
    
    for(j = 0; j < isoformprs.size(); ++j){
        l=0;
        if(!clustered[j]){
            clustercoords.push_back(isoformprs[j]);
            clustered_tr_ids.push_back(tr_ids[j]);
            clustered_pr_ids.push_back(probe_ids[j]);
            if(j+1 < isoformprs.size()){
                for(k = j+1; k < isoformprs.size();++k){
                    if(!clustered[j]){
                        if(abs(isoformprs[j] - isoformprs[k]) < ClusterPromoters){
                            ++l;
                            clustered[k] = 1;
                            clustercoords.push_back(isoformprs[k]);
                            clustered_tr_ids.push_back(tr_ids[k]);
                            clustered_pr_ids.push_back(probe_ids[k]);
                        }
                    }
                }
                if (strand == "+") {
                    int leftmosttss = FindLeftMostCoord(clustercoords);
                    clusteredtrcoords.push_back(leftmosttss);
                    clustered_ids.push_back(clustered_tr_ids[0]);
                    clustered_probe_ids.push_back(clustered_pr_ids[0]);
                    ++cluster_idx;
                    clustercoords.clear();
                    clustered_tr_ids.clear();
                    clustered_pr_ids.clear();
                }
                else{
                    int rightmosttss = FindRightMostCoord(clustercoords);
                    clusteredtrcoords.push_back(rightmosttss);
                    clustered_ids.push_back(clustered_tr_ids[0]);
                    clustered_probe_ids.push_back(clustered_pr_ids[0]);
                    ++cluster_idx;
                    clustercoords.clear();
                    clustered_tr_ids.clear();
                    clustered_pr_ids.clear();
                }
            }
            else{
                clusteredtrcoords.push_back(isoformprs[j]);
                clustered_ids.push_back(tr_ids[j]);
                clustered_probe_ids.push_back(probe_ids[j]);
            }
        }
    }
}

void ProbeFeatureClass::ReadFeatureAnnotation(RESitesClass& dpnIIsites, std::string transcriptfile, std::string option)
{
	std::string temp,tr1,tr2; 
	//std::string feature_id;
	int promindex = 0; 
	int invalidRECoordinate=0;
	
	bool flag = true; 
	
	

	std::string RefSeqfilename;

	RefSeqfilename.append(transcriptfile);
	    
    pLog << option<< "File is " << RefSeqfilename << std::endl;
    
	std::ifstream RefSeq_file(RefSeqfilename.c_str());

    std::vector < int > isoformprs; // to keep isoform promoters
    std::vector < std::string > tr_ids;
    std::vector < std::string > probe_ids;
    
    std::vector < int > clusteredcoords;
	std::vector < std::string > ids_of_clustered;
    std::vector < std::string > probe_ids_of_clustered;
    
	getline(RefSeq_file, temp);

    getline(RefSeq_file, tr1);

	std::stringstream trx1 ( tr1 ); 
	GetTrFeats(trx1, tp[0], option);

	isoformprs.push_back(tp[0].start);
	tr_ids.push_back(tp[0].tr_id);
    probe_ids.push_back(tp[0].probe_id);
     
	while(getline(RefSeq_file,tr2)){
		std::stringstream trx2 ( tr2);
		GetTrFeats(trx2, tp[1], option);
		while(tp[0].name == tp[1].name){
            isoformprs.push_back(tp[1].start);
			tr_ids.push_back(tp[1].tr_id);
            probe_ids.push_back(tp[1].probe_id);
			if(getline(RefSeq_file,tr2)){
				std::stringstream trx2 ( tr2);
				GetTrFeats(trx2, tp[1], option);
			}
			else
				break;
		};
        
        if(isoformprs.size() > 1){
            ClusterIsoformPromoters(isoformprs, tr_ids, probe_ids, clusteredcoords, ids_of_clustered, probe_ids_of_clustered, tp[0].strand); //Promoters that are within "ClusterPromoters" of each other are clustered
            // Fill promoter struct each unclustered isoform will have its promoter
            for(int y = 0; y < clusteredcoords.size();++y){
				std::string feature_id;
				feature_id.append(tp[0].name);
                feature_id.append("_");
                std::string tr_st = std::to_string(clusteredcoords[y]);
                feature_id.append(tr_st);
                
                promFeatures[feature_id].genes.push_back(tp[0].name);
                promFeatures[feature_id].transcripts.push_back(ids_of_clustered[y]);
                promFeatures[feature_id].strand = tp[0].strand;
                promFeatures[feature_id].chr = tp[0].chr;
                
                promFeatures[feature_id].TSS = clusteredcoords[y];
                
                //////////////Add intervals to vector clusteredcoords[y]
				if(chrIntervals.find(tp[0].chr)!=chrIntervals.end()){
					chrIntervals[tp[0].chr].push_back(Interval <std::string>((clusteredcoords[y] - promPadding),(clusteredcoords[y] + promPadding), feature_id));
					ChrNames_proms.push_back(tp[0].chr);
				}
				else{
					std::vector < Interval < std::string > > tempvector ;
					tempvector.push_back(Interval <std::string>((clusteredcoords[y] - promPadding),(clusteredcoords[y] + promPadding), feature_id));
					chrIntervals.emplace(tp[0].chr, tempvector);
				}
				
				//deal with shared promoters
				if(transcriptCounter.find(tp[0].name+"_"+tp[0].chr)!=transcriptCounter.end()){
					transcriptCounter[tp[0].name+"_"+tp[0].chr].count = transcriptCounter[tp[0].name].count + 1;
				}
				else{
					cCounter tmp;
					tmp.count=1;
					tmp.probesDesigned=false;
					transcriptCounter.emplace(tp[0].name+"_"+tp[0].chr, tmp);
				}
				
				promFeatures[feature_id].sharedpromoter = false;
				promFeatures[feature_id].probesSkip = false;
                
                promFeatures[feature_id].ProbeID.push_back(probe_ids_of_clustered[y]);
                promFeatures[feature_id].FeatureType = tp[0].FeatureType;
                promFeatures[feature_id].end = tp[0].end;
                                
                dpnIIsites.GettheREPositions(promFeatures[feature_id].chr, promFeatures[feature_id].TSS, promFeatures[feature_id].closestREsitenums, invalidRECoordinate);
            }
        }
        else{
			std::string feature_id;
            feature_id.append(tp[0].name);
            feature_id.append("_");
            
            std::string tr_st = std::to_string(isoformprs[0]);
            
            feature_id.append(tr_st);
            
            promFeatures[feature_id].genes.push_back(tp[0].name);
            promFeatures[feature_id].transcripts.push_back(tr_ids[0]);
            promFeatures[feature_id].strand = tp[0].strand;
            promFeatures[feature_id].chr.append(tp[0].chr);
            promFeatures[feature_id].TSS = isoformprs[0];
            promFeatures[feature_id].ProbeID.push_back(probe_ids[0]);
            promFeatures[feature_id].FeatureType = tp[0].FeatureType;
            promFeatures[feature_id].end = tp[0].end;
            
            
			if(chrIntervals.find(tp[0].chr)!=chrIntervals.end()){
				chrIntervals[tp[0].chr].push_back(Interval <std::string>((isoformprs[0] - promPadding),(isoformprs[0] + promPadding), feature_id));
				ChrNames_proms.push_back(tp[0].chr);
			}
			else{
				std::vector < Interval < std::string > > tempvector ;
				tempvector.push_back(Interval <std::string>((isoformprs[0] - promPadding),(isoformprs[0] + promPadding), feature_id));
				chrIntervals.emplace(tp[0].chr, tempvector);
			}
			
			//deal with shared promoters
			if(transcriptCounter.find(tp[0].name+"_"+tp[0].chr)!=transcriptCounter.end()){
				transcriptCounter[tp[0].name+"_"+tp[0].chr].count = transcriptCounter[tp[0].name+"_"+tp[0].chr].count + 1;
			}
			else{
				cCounter tmp;
				tmp.count=1;
				tmp.probesDesigned=false;
				transcriptCounter.emplace(tp[0].name+"_"+tp[0].chr, tmp);
			}
				
			promFeatures[feature_id].sharedpromoter = false;
			promFeatures[feature_id].probesSkip = false;
			
		        
            dpnIIsites.GettheREPositions(promFeatures[feature_id].chr, promFeatures[feature_id].TSS, promFeatures[feature_id].closestREsitenums, invalidRECoordinate);
        }
        
        clusteredcoords.clear();
        ids_of_clustered.clear();
        probe_ids_of_clustered.clear();
		isoformprs.clear();
        probe_ids.clear();
		isoformprs.push_back(tp[1].start);
		tr_ids.clear();
		tr_ids.push_back(tp[1].tr_id);
        probe_ids.push_back(tp[1].probe_id);
        
	//SWAP TP[0] and TP[1]
		tp[0].name.clear();
		tp[0].name.append(tp[1].name);
		tp[0].tr_id.clear();
		tp[0].tr_id.append(tp[1].tr_id);
        tp[0].probe_id.clear();
        tp[0].probe_id.append(tp[1].probe_id);
		tp[0].chr.clear();
		tp[0].chr.append(tp[1].chr);
		tp[0].end = tp[1].end;
		tp[0].start = tp[1].start;
		tp[0].strand = tp[1].strand;
		++promindex;
		
		if(promindex%10000 == 0){
		  pLog << "..."<< promindex<<" Features Annotated from "<< option <<" file"<< std::endl;
	  }
    }
	
	pLog << "Number of Features Annotated from "<< option <<" file: "<<promindex<< std::endl;

	NofPromoters += promindex; 
	++fileReadCount; 
	
	if(fileReadCount==fileCount){ 
		pLog<<"Total Number of Features annotated: "<< NofPromoters<<std::endl;
		DealwithSharedPromoters();
		pLog << "Shared Promoters Determined" << std::endl;
		for(auto it = chrIntervals.begin(); it != chrIntervals.end(); ++it){
				std::vector< Interval < std::string > > temp;
				promIntTree[it->first] = IntervalTree< std::string >(it->second);
				it->second.swap(temp);
		}
	}
}


void ProbeFeatureClass::DealwithSharedPromoters(){ // If promoters are too close to each other
	for (auto it = promFeatures.begin(); it != promFeatures.end(); ++it) {		
		for (auto itj = std::next(it , 1); itj != promFeatures.end(); ++itj) {
			if (it->second.chr == itj->second.chr){
				if(!(it->second.probesSkip)){
					if(!(itj->second.probesSkip)){
						std::string tmpTrans1 = (it->second).genes[0]+"_"+it->second.chr;
						std::string tmpTrans2 = (itj->second).genes[0]+"_"+itj->second.chr;;
						if(abs(it->second.TSS - itj->second.TSS) < ClusterPromoters){		
							if(transcriptCounter[tmpTrans1].count < transcriptCounter[tmpTrans2].count){
								itj->second.probesSkip=true;
								transcriptCounter[tmpTrans2].count = transcriptCounter[tmpTrans2].count - 1;
								transcriptCounter[tmpTrans1].probesDesigned=true;
							}
							else if(transcriptCounter[tmpTrans1].count > transcriptCounter[tmpTrans2].count){
								it->second.probesSkip=true;
								transcriptCounter[tmpTrans1].count = transcriptCounter[tmpTrans1].count - 1;
								transcriptCounter[tmpTrans2].probesDesigned=true;
									
							}
							else if(transcriptCounter[tmpTrans1].count == transcriptCounter[tmpTrans2].count){
								if(transcriptCounter[tmpTrans1].probesDesigned && !(transcriptCounter[tmpTrans2].probesDesigned)){
									it->second.probesSkip=true;
									transcriptCounter[tmpTrans1].count = transcriptCounter[tmpTrans1].count - 1;	
									transcriptCounter[tmpTrans2].probesDesigned=true;
								}
								else if(!(transcriptCounter[tmpTrans1].probesDesigned) && (transcriptCounter[tmpTrans2].probesDesigned)){
									itj->second.probesSkip=true;
									transcriptCounter[tmpTrans2].count = transcriptCounter[tmpTrans2].count - 1;
									transcriptCounter[tmpTrans1].probesDesigned=true;
									
								}
								else{
									if(abs(it->second.TSS-it->second.end) > abs(itj->second.TSS-itj->second.end)){
										itj->second.probesSkip=true;
										transcriptCounter[tmpTrans2].count = transcriptCounter[tmpTrans2].count - 1;
										transcriptCounter[tmpTrans1].probesDesigned=true;
									}
									else if(abs(it->second.TSS-it->second.end) < abs(itj->second.TSS-itj->second.end)){
										it->second.probesSkip=true;
										transcriptCounter[tmpTrans1].count = transcriptCounter[tmpTrans1].count - 1;
										transcriptCounter[tmpTrans2].probesDesigned=true;
									}
								}
								
							}
						}
						else{
							transcriptCounter[tmpTrans1].probesDesigned=true;
						}
					}
				}
			}
		}
	}
}
