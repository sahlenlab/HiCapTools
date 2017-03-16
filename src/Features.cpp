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
//  Features.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "Features.h"
#include "Data_Structures.h"
#include <map>
#include <fstream>
#include <algorithm>

std::map< std::string, FeatureStruct > Features;
std::map< std::string, std::vector < std::string > > MetaFeatures;

void FeatureClass::InitialiseData(int clustProm, int fCount){
	
	ClusterPromoters = clustProm;
	fileCount = fCount;
	filesReadCount=0;
	tp = new temppars [2];
    pLog << "Promoter Class Initialised" << std::endl;

}

void FeatureClass::GetTrFeats(std::stringstream &trx, temppars &tpars, std::string option){
	
	std::string field, start, end;
    
    //Transcript Line Format
    //name2	 name	chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds
    
    //BED Line Format for SNP and Negative Control Region
    //chrom	Start	End	name
	
	if(option=="transcript"){
		getline(trx,tpars.name,'\t'); 
		getline(trx,tpars.tr_id,'\t'); 
    }
    
    getline(trx,tpars.chr,'\t'); 
    
    if(option=="transcript"){
		getline(trx,tpars.strand,'\t');
    }
    
    getline(trx,start,'\t');
    getline(trx,end,'\t');
    
    if(option=="transcript"){
		getline(trx,field,'\t');
		getline(trx,field,'\t');
		getline(trx,field,'\t');
    
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
    
   if(option=="SNV" || option=="neg_ctrl"){
	   
	   tpars.start=std::stoi(start);
	   tpars.end=std::stoi(end);
	   tpars.strand="+";
	   
	   getline(trx,tpars.name,'\t');
	   for (auto & c: tpars.name) c = toupper(c);
	   tpars.tr_id=tpars.name;
	   if(option=="SNV")
		tpars.FeatureType = 2;
	   if(option=="neg_ctrl")
		tpars.FeatureType = 3;
	   tpars.probe_id=tpars.tr_id+"."+start;
   }
}



int FeatureClass::FindLeftMostCoord(std::vector<int> coords){
    
    int leftmosttss;
    
    leftmosttss = coords[0];
    for (int i = 0; i < coords.size();++i){
        if (coords[i] <= leftmosttss )
            leftmosttss = coords[i];
    }
    return leftmosttss;
}


int FeatureClass::FindRightMostCoord(std::vector<int> coords){
	
    int rightmosttss;
    
    rightmosttss = coords[0];
    for (int i = 0; i < coords.size();++i){
        if (coords[i] >= rightmosttss )
            rightmosttss = coords[i];
    }
    return rightmosttss;
}


void FeatureClass::ClusterIsoformPromoters(std::vector <int>& isoformprs, std::vector<std::string>& tr_ids, std::vector <int>& clusteredtrcoords, std::vector<std::string>& clustered_ids, std::string strand, int ClusterPromoters){
    
    int j,k,l,z;
    
    std::vector<int> clustercoords;
    std::vector<std::string> clustered_tr_ids;
    std::vector<bool> clustered;
    
    for(j = 0;j < isoformprs.size(); ++j)
        clustered.push_back(0);
    
    for(j = 0; j < isoformprs.size(); ++j){
        l = 0;
        if(!clustered[j]){
            clustercoords.push_back(isoformprs[j]);
            clustered_tr_ids.push_back(tr_ids[j]);
        
            if(j+1<isoformprs.size()){
                for(k = j+1;k < isoformprs.size();++k){
                    if(!clustered[j]){
                        if(abs(isoformprs[j] - isoformprs[k]) < ClusterPromoters){
                            ++l;
                            clustered[k] = 1;
                            clustercoords.push_back(isoformprs[k]);
                            clustered_tr_ids.push_back(tr_ids[k]);
                        }
                    }
                }
                if (strand == "+") {
                    int leftmosttss = FindLeftMostCoord(clustercoords);
                    clusteredtrcoords.push_back(leftmosttss);
                    clustered_ids.push_back(clustered_tr_ids[0]);
   
                    clustercoords.clear();
                    clustered_tr_ids.clear();
                    
                }
                else{
                    int rightmosttss = FindRightMostCoord(clustercoords);
                    clusteredtrcoords.push_back(rightmosttss);
                    clustered_ids.push_back(clustered_tr_ids[0]);
                    
                    clustercoords.clear();
                    clustered_tr_ids.clear();
                    
                }
            }
            else{
                clusteredtrcoords.push_back(isoformprs[j]);
                clustered_ids.push_back(tr_ids[j]);
            }
        }
    }	
}




void FeatureClass::ReadFeatureAnnotation(RESitesClass& dpnIIsites, std::string TranscriptListFileName, std::string option){
    std::string temp,tr1, tr2, feature_id;
    int nofgenes = 0, z = 0;
    std::locale l;
    
    std::string RefSeqfilename;
    
    RefSeqfilename.append(TranscriptListFileName);
    
    pLog<< option << " File is " << RefSeqfilename << std::endl;
    
    std::ifstream RefSeq_file(RefSeqfilename.c_str());
    
    std::vector< int > isoformprs;
    std::vector< std::string > tr_ids;
    
    
    std::vector < int > clusteredcoords;
    std::vector < std::string > ids_of_clustered;
    
    
    //Discard headers
    
	getline(RefSeq_file, temp);
	
	getline(RefSeq_file, tr1);
		
    
    std::stringstream trx1 ( tr1 );
    GetTrFeats(trx1,tp[0], option);
    
    isoformprs.push_back(tp[0].start);
	tr_ids.push_back(tp[0].tr_id);
    
    while(getline(RefSeq_file,tr2)){
		
        std::stringstream trx2 ( tr2);
        GetTrFeats(trx2,tp[1], option);
        while(tp[0].name == tp[1].name){
        
            isoformprs.push_back(tp[1].start);
            tr_ids.push_back(tp[1].tr_id);
            
            if(getline(RefSeq_file,tr2)){
				std::stringstream trx2 ( tr2);
				GetTrFeats(trx2,tp[1], option);
			}
			else
				break;
        };
        
       
       if(isoformprs.size() > 1){
            
            ClusterIsoformPromoters(isoformprs, tr_ids, clusteredcoords, ids_of_clustered, tp[0].strand, ClusterPromoters); //Promoters that are within "ClusterPromoters" of each other are clustered
            for(int y = 0; y < clusteredcoords.size();++y){
                feature_id.append(tp[0].name);
                feature_id.append("_");
                
                std::string tr_st = std::to_string(clusteredcoords[y]);
                
                feature_id.append(tr_st);
                
                Features[feature_id].Name = tp[0].name;
                Features[feature_id].TranscriptName = ids_of_clustered[y];
                Features[feature_id].probe_name = tp[0].probe_id;
                Features[feature_id].strand = tp[0].strand;
                Features[feature_id].chr = tp[0].chr;
                Features[feature_id].start = clusteredcoords[y];
                Features[feature_id].end = tp[0].end;
                Features[feature_id].feature_id = feature_id;
                Features[feature_id].probe_target = tp[0].probetarget;
                Features[feature_id].FeatureType = tp[0].FeatureType;
                
             
                MetaFeatures[Features[feature_id].Name].push_back(feature_id);
                feature_id = "";
            }
        }
        else{
            feature_id.append(tp[0].name);
            feature_id.append("_");
            
            std::string tr_st = std::to_string(isoformprs[0]);
            
            feature_id.append(tr_st);
            
            Features[feature_id].Name = tp[0].name;
            Features[feature_id].TranscriptName = tr_ids[0];
            Features[feature_id].probe_name = tp[0].probe_id;
            Features[feature_id].strand = tp[0].strand;
            Features[feature_id].chr = tp[0].chr;
            Features[feature_id].start = isoformprs[0];
            Features[feature_id].end = tp[0].end;
            Features[feature_id].feature_id = feature_id;
            Features[feature_id].probe_target = tp[0].probetarget;
            Features[feature_id].FeatureType = tp[0].FeatureType;
            
            MetaFeatures[Features[feature_id].Name].push_back(feature_id);
        }
        feature_id ="";
        isoformprs.clear();
        clusteredcoords.clear();
        ids_of_clustered.clear();
        tr_ids.clear();
        
        
		isoformprs.push_back(tp[1].start);
		tr_ids.push_back(tp[1].tr_id);
        
        //SWAP TP[0] and TP[1]
		tp[0].name.clear();
		tp[0].name.append(tp[1].name);
		tp[0].tr_id.clear();
		tp[0].tr_id.append(tp[1].tr_id);
		tp[0].chr.clear();
		tp[0].chr.append(tp[1].chr);
        tp[0].probe_id.clear();
        tp[0].probe_id.append(tp[1].probe_id);
        tp[0].probetarget.clear();
        tp[0].probetarget.append(tp[1].probetarget);
		tp[0].end = tp[1].end;
		tp[0].start = tp[1].start;
		tp[0].strand = tp[1].strand;
        
		++nofgenes;
		
		if(nofgenes%10000 == 0){
		  pLog << nofgenes << "  Features Annotated from " <<option<< "file"<< std::endl;
        } 
        
        
    }
    ++filesReadCount;
    pLog << nofgenes << "  lines read corresponding to " << Features.size() << " features from " <<option<<" file."<< std::endl;
    
    if(filesReadCount==fileCount){ 
		pLog<<"Total Number of Features annotated: "<< Features.size()<<std::endl;
	}
}



