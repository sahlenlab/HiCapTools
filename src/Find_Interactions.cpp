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
//  Find_Interactions.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "Find_Interactions.h"
#include "Global.h"
#include <fstream>
//ALGLIB PACKAGE HEADERS

#include "alglibmisc.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "dataanalysis.h"
#include "specialfunctions.h"
#include "solvers.h"
#include "optimization.h"
#include "diffequations.h"
#include "fasttransforms.h"
#include "integration.h"
#include "interpolation.h"
using namespace alglib_impl;

bool DetectInteractions::CheckSupportingPairs(int* supppairs, int nOfExp){
	bool recordit = 0;
	for (int t = 0; t < nOfExp; ++t){
		if (supppairs[t] >= MinNumberofSupportingPairs){
			recordit = 1;
			break;
		}
	}
	return recordit;
}

double DetectInteractions::CalculatepVal(std::map< int, double > mean, std::map< int, double > std, int bin,int sig){
	
    double q = 0, p_value = 1;
	if(mean.find(bin) == mean.end()){
		mean[bin] = 0;
		std[bin] = 1;
		p_value = 1;
	}
	else{
        if (std[bin] == 0.0)
            q = 0;
        else
            q = alglib::normaldistribution(( sig - mean[bin]) / std[bin]); //z-score
		p_value = 1- q;
	}
	return p_value;
}

void DetectInteractions::CalculatePvalAndPrintInteractionsProbeDistal(ProbeSet& prs, std::vector<DetermineBackgroundLevels> background, std::string BaseFileName, int NumberofExperiments, std::vector < std::string >& ExperimentNames, std::string whichchr, int BinSize, PrDes::RENFileInfo& reInfo){

    fLog << "will print proximities " << std::endl;
//Probe to Distal Interactions
	std::string FileName, FileName2;
    
    int enoughpairs, count = 0;
    double p_value;
    
    
	FileName.append(BaseFileName);
    FileName.append(".");
	FileName.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileName.append(".");
    FileName.append(whichchr);
    FileName.append(".Proximities.Probe_Distal.");
	FileName.append(reInfo.currTime);
	FileName.append(".txt");
	std::ofstream outf1(FileName.c_str());
	
	/**
	FileName2.append(BaseFileName); //washU ouput
	FileName2.append(".");
	FileName2.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileName2.append(".");
    FileName2.append(whichchr);
    FileName2.append(".Proximities.Probe_Distal.WashU.");
	FileName2.append(reInfo.currTime);
	FileName2.append(".txt");
	std::ofstream outf2(FileName2.c_str());
    **/
	
    outf1 << "RefSeqName" << '\t' << "TranscriptName" << '\t' << "Feature_ID" << '\t' 
		<< "Feature_Chr" << '\t' << "Feature_Start" << '\t' << "Feature_End" << '\t' << "Annotation" << '\t' <<  "Strand" << '\t';
    
    outf1 << "Interactor_Chr" << '\t' << "Interactor_Start" << '\t' << "Interactor_End" << '\t' << "distance" ;
    
	for (int e = 0; e < NumberofExperiments; ++e)
		outf1 << '\t'<< ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_p_value"  << '\t' << ExperimentNames[e] << "_StrandCombination" ;
    outf1 << std::endl;

    std::map<std::string, FeatureStruct>::iterator featiter;
    std::map<std::string, FeatureStruct>::iterator featiter2;
    
    std::map<int, Junction >::iterator it; //key:REpos
    std::vector < SignalStruct_CTX >::iterator itx;
    std::map<int, Junction >::iterator itt;
    //initialise memory
    for(int i=0; i< Features.size(); ++i){
		
		auto featiter = Features.begin() ;
		std::advance(featiter, i);
		
        for (it = featiter->second.proximities.junctions.begin(); it != featiter->second.proximities.junctions.end(); ++it){
            
            if(featiter->second.FeatureType != 3 && ((abs(it->first - featiter->second.start)) > MinimumJunctionDistance)){
                
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
            
            if(it->second.reportit){
                                             
               it->second.distance = it->first - featiter->second.start;
                 it->second.p_val.reserve(NumberofExperiments);
                 
            }
		  }
        }        
    }
    
    //calculate pvalue in parallel
   // #pragma omp parallel for
    for(int i=0; i< Features.size(); ++i){
		
		auto featiter = Features.begin() ;
		std::advance(featiter, i);
		
        for (it = featiter->second.proximities.junctions.begin(); it != featiter->second.proximities.junctions.end(); ++it){
            
            if(featiter->second.FeatureType != 3 && ((abs(it->first - featiter->second.start)) > MinimumJunctionDistance)){
                
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
            
            if(it->second.reportit){
                                             
               it->second.distance = it->first - featiter->second.start;
                
                int bin;
                bin = abs(it->second.distance) / BinSize;
                for (int e = 0; e < NumberofExperiments; ++e){
                    if (CALCULATE_P_VALUES && it->second.paircount[e] >= MinNumberofSupportingPairs){
						if(background[e].bglevels.smoothed.find(bin) == background[e].bglevels.smoothed.end()){
							int i=0;
							
							while(background[e].bglevels.smoothed.find(bin-i)==background[e].bglevels.smoothed.end()){
								if(bin==0)
									--i;
								else
									++i;
							}
							if(bin-i>=0 && bin-i<=background[e].bglevels.smoothed.rbegin()->first){
								background[e].bglevels.smoothed[bin] = background[e].bglevels.smoothed[bin-i];
								background[e].bglevels.smoothed_stdev[bin] = background[e].bglevels.smoothed_stdev[bin-i];
									
							}
							else{
								background[e].bglevels.smoothed[bin] = 0;
								background[e].bglevels.smoothed_stdev[bin] = 1;
							}
						}
     
						if (background[e].bglevels.smoothed[bin] == 0.0)
							it->second.p_val.push_back(1.0);
								
						else{
							double t= 1.0 - alglib::normaldistribution(( it->second.paircount[e] - background[e].bglevels.smoothed[bin]) / background[e].bglevels.smoothed_stdev[bin]); 
							it->second.p_val.push_back(t);
						}
				   }
                    else
                        it->second.p_val.push_back(1.0);
                }
                
            }
		  }
        }        
    }
	
    //output to file
    for(featiter = Features.begin() ; featiter != Features.end(); ++featiter){
		
		if (count % 2500 == 0)
            fLog << count << " probes processed in Probe-Distal" << std::endl;
        ++count;
            
        for (it = featiter->second.proximities.junctions.begin(); it != featiter->second.proximities.junctions.end(); ++it){
            
            if(featiter->second.FeatureType != 3 && ((abs(it->first - featiter->second.start)) > MinimumJunctionDistance)){
                
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
            
            if(it->second.reportit){
               
                outf1 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' 
                      << featiter->second.chr  << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                      << featiter->second.FeatureType << '\t' << featiter->second.strand << '\t';
                
                outf1 << featiter->second.chr << '\t' << it->first << '\t' << it->second.refragend << '\t';
               /**
               //WashU
				int wStart, wEnd;
               if(featiter->second.start > featiter->second.end){
					wStart=featiter->second.end - 100;
					wEnd= featiter->second.end + 100;	
			   }
               else if(featiter->second.start == featiter->second.end){
					wStart=featiter->second.start-100;
					wEnd= featiter->second.end + 100;	
			   }
               else if(featiter->second.start < featiter->second.end){
					wStart=featiter->second.start - 100;
					wEnd= featiter->second.start + 100;	
			   }
               
				outf2 << featiter->second.chr << ':' << wStart << '-' << wEnd<< '\t'<< featiter->second.chr<<':' << it->first << '-' << it->second.refragend << '\t';
               
               **/
                outf1 << it->second.distance;
                
                //WashU
               // double avgscore=0;
                
                for (int e = 0; e < NumberofExperiments; ++e){
                    
                    outf1 << '\t' << it->second.paircount[e] << '\t' << it->second.p_val[e] << '\t';
                    for(int s = 0; s < 3; ++s)
                        outf1 << it->second.strandcombination[(e*4)+s] << "_";
					outf1 << it->second.strandcombination[(e*4)+3];
                    
                    //calculate washU score
                    //avgscore = avgscore + (it->second.paircount[e]/double(NumberofExperiments));
                    
                    
                }
                outf1 << std::endl;
                
                //washU
               // outf2 << avgscore <<std::endl;
            }
		  }
        }
        
        for (itx = featiter->second.proximities_ctx.begin(); itx != featiter->second.proximities_ctx.end(); ++itx){
            if(featiter->second.FeatureType != 3 ){
				for(itt = itx->junctions_ctx.begin(); itt != itx->junctions_ctx.end(); ++itt){
					
					enoughpairs = CheckSupportingPairs(itt->second.paircount, NumberofExperiments);
					if(enoughpairs){
                    
						outf1 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t'
							<< featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
							<< featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                          
						outf1 << itx->maptochrname << '\t' << itt->first << '\t' << itt->second.refragend << '\t' << -1 ;
						
						/**
						//WashU
						int wStart, wEnd;
						if(featiter->second.start > featiter->second.end){
							wStart=featiter->second.end - 100;
							wEnd= featiter->second.end + 100;	
						}
						else if(featiter->second.start == featiter->second.end){
							wStart=featiter->second.start - 100;
							wEnd= featiter->second.end + 100;	
						}
						else if(featiter->second.start < featiter->second.end){
							wStart=featiter->second.start -100;
							wEnd= featiter->second.start + 100;	
						}
						
						int dStart, dEnd;
						if(itt->first == itt->second.refragend){
							dStart=itt->first;
							dEnd= itt->second.refragend + 1;	
						}
						else{
							dStart=itt->first;
							dEnd= itt->second.refragend;	
						}
						
						outf2 << featiter->second.chr << ':' << wStart << '-' << wEnd<< '\t'<< itx->maptochrname<<':' << dStart << '-' << dEnd ;
						
						//WashU
						double avgscore=0; **/
						
						for (int e = 0; e < NumberofExperiments; ++e){
							outf1 << '\t'<< itt->second.paircount[e] << '\t' << -1 << '\t';
							for(int s = 0; s < 3; ++s)
								outf1 << itt->second.strandcombination[(e*4) + s] << "_";
							outf1 << itt->second.strandcombination[(e*4) + 3];
							
							//calculate washU score
							//avgscore = avgscore + (itt->second.paircount[e]/double(NumberofExperiments));
							
						}
						outf1 << std::endl;
						
						//washU
						//outf2 << avgscore <<std::endl;
					}
				}
			}
        }
    }
    outf1.close();
    count = 0;
    
    
}  
    
void DetectInteractions::CalculatePvalAndPrintInteractionsProbeProbe(ProbeSet& prs, std::vector<DetermineBackgroundLevels> background, std::string BaseFileName, int NumberofExperiments, std::vector < std::string >& ExperimentNames, std::string whichchr, int BinSizeProbeProbe, PrDes::RENFileInfo& reInfo){    
//Promoter-promoter Interactions

	int bin;
	int flag;
	int count=0;
	std::string FileName3, FileNameWashU;
	
	//WashU
	//std::vector<std::string> seenPP;
	
	FileName3.append(BaseFileName);
    FileName3.append(".");
	FileName3.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileName3.append(".");
    FileName3.append(whichchr);
	FileName3.append(".Proximities.Probe_Probe.");
	FileName3.append(reInfo.currTime);
	FileName3.append(".txt");
	
	std::ofstream outf3(FileName3.c_str());
	/**
	//WashU
	FileNameWashU.append(BaseFileName);
    FileNameWashU.append(".");
	FileNameWashU.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileNameWashU.append(".");
    FileNameWashU.append(whichchr);
	FileNameWashU.append(".Proximities.Probe_Probe.WashU.");
	FileNameWashU.append(reInfo.currTime);
	FileNameWashU.append(".txt");
		
	std::ofstream outfwu(FileNameWashU.c_str());
	**/
    
	outf3 << "RefSeqName_1" << '\t' << "TranscriptName_1" << '\t' << "Feature_ID_1" << '\t' 
	<< "FeatureChr_1" << '\t' << "FeatureStart_1" << '\t' << "FeatureEnd_1"
    << '\t' << "Annotation_1" << '\t' <<  "Strand_1" << '\t';
    
    outf3 << "RefSeqName_2" << '\t' << "TranscriptName_2" << '\t' << "Feature_ID_2" << '\t' 
    << "FeatureChr_2" << '\t' << "FeatureStart_2" << '\t' << "FeatureEnd_2"
    << '\t' << "Annotation_2" << '\t' <<  "Strand_2" << '\t';
   
    outf3 << "abs(Distance)";

    for (int e = 0; e < NumberofExperiments; ++e)
        outf3 << '\t' << ExperimentNames[e] << "_SuppPairs"<< '\t' << ExperimentNames[e] << "_p_value" << '\t' << ExperimentNames[e] << "_StrandCombination";
    outf3 << std::endl;
	std::vector< FeattoFeatSignalStruct >::const_iterator itff; //first: REpos, second: signal
	std::string f;
	
	for(auto featiter = Features.begin() ; featiter != Features.end(); ++featiter){
		
        if (count % 2500 == 0)
            fLog << count << " probes processed in Probe-Probe" << std::endl;
        ++count;
        
        
        for (itff = featiter->second.Inter_feature_ints.begin(); itff != featiter->second.Inter_feature_ints.end(); ++itff){
			f = itff->interacting_feature_id;
            auto featiter2 = Features.find(f);
            if((featiter->second.FeatureType != 3 && featiter2->second.FeatureType != 3) && featiter->second.TranscriptName != featiter2->second.TranscriptName && (abs(featiter->second.start - featiter2->second.start) >= MinimumJunctionDistance)){
				
				
								
                outf3 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' 
                      << featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                      << featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                
                outf3 << featiter2->second.Name << '\t' << featiter2->second.TranscriptName << '\t' << featiter2->first << '\t' 
                      << featiter2->second.chr << '\t' << featiter2->second.start << '\t' << featiter2->second.end << '\t'
                      << featiter2->second.FeatureType << '\t' << featiter2->second.strand << '\t';
                /***
                //WashU
                bool seen=false;      
                std::string constructPP = featiter2->first+":"+featiter->first;      
                if(std::find(seenPP.begin(), seenPP.end(), constructPP)==seenPP.end()){ //If the interaction does not already have the duplicate written
					seenPP.push_back(constructPP);
					int start1, end1, start2, end2;
					if(featiter->second.start>featiter->second.end){
						start1=featiter->second.end - 100;
						end1=featiter->second.end + 100;
					}
					else if(featiter->second.start==featiter->second.end){
						start1=featiter->second.start -100;
						end1=featiter->second.end + 100;
					}
					else if(featiter->second.start<featiter->second.end){
						start1=featiter->second.start -100;
						end1=featiter->second.start + 100;
					}
					
					if(featiter2->second.start > featiter2->second.end){
						start2=featiter2->second.end -100;
						end2=featiter2->second.end + 100;
					}
					else if(featiter2->second.start==featiter2->second.end){
						start2=featiter2->second.start - 100;
						end2=featiter2->second.end + 100;
					}
					else if(featiter2->second.start < featiter2->second.end){
						start2=featiter2->second.start -100;
						end2=featiter2->second.start + 100;
					}
					
					outfwu<< featiter->second.chr << ':' << start1 << '-' << end1<< '\t'<< featiter2->second.chr<<':' << start2 << '-' << end2 << '\t';
				} 
				else 
					seen = true;
                ***/
                bin=0;
                flag=0;
                if(featiter->second.chr == featiter2->second.chr){
                    outf3 << abs(featiter->second.start - featiter2->second.start) ;
                    bin=abs(featiter->second.start - featiter2->second.start) / BinSizeProbeProbe;
                    flag=1;   
				}
                else
                    outf3 << -1 ;
                    
                //WashU
				//double avgscore=0;
				
                for (int e = 0; e < NumberofExperiments; ++e){
                    outf3 << '\t' << itff->signal[e] << '\t';
                    
                    //WashU
                  //  if(!seen)
					//	avgscore = avgscore + (itff->signal[e]/double(NumberofExperiments));
                    
                    ///////Pval calc begins
                                      					
					if (CALCULATE_P_VALUES && flag){
						if(background[e].bglevelsProbeProbe.smoothed.find(bin) == background[e].bglevelsProbeProbe.smoothed.end()){
							int i=0;
							while(background[e].bglevelsProbeProbe.smoothed.find(bin-i)==background[e].bglevelsProbeProbe.smoothed.end()){
								if(bin==0)
									--i;
								else
									++i;
							}
							if(bin-i>=0 && bin-i<=background[e].bglevelsProbeProbe.smoothed.rbegin()->first){
								background[e].bglevelsProbeProbe.smoothed[bin] = background[e].bglevelsProbeProbe.smoothed[bin-i];
								background[e].bglevelsProbeProbe.smoothed_stdev[bin] = background[e].bglevelsProbeProbe.smoothed_stdev[bin-i];
							}
							else{
								background[e].bglevelsProbeProbe.smoothed[bin] = 0;
								background[e].bglevelsProbeProbe.smoothed_stdev[bin] = 1;
							}
						}
						if (background[e].bglevelsProbeProbe.smoothed[bin] == 0.0)
							outf3 << "1.0"; 
						
						else{
							double t= 1.0 - alglib::normaldistribution(( itff->signal[e] - background[e].bglevelsProbeProbe.smoothed[bin]) / background[e].bglevelsProbeProbe.smoothed_stdev[bin]);
							outf3 << t ;
						}
					}
					else
                       outf3 << "1.0";
                       
                   /////Pval calc ends
                   outf3<<"\t";
                   for(int s=0; s<3; ++s){
					   outf3 << itff->strandcombination[(e*4) + s] << "_";
				   }
				   outf3 << itff->strandcombination[(e*4) + 3];
                }
                outf3 << std::endl;
                
                //WashU
               // if(!seen)
				//	outfwu<<  avgscore<<std::endl;
            }
  
        }
    }
	outf3.close();

}

void DetectInteractions::CalculatePvalAndPrintInteractionsProbeDistal_NegCtrls(ProbeSet& prs, std::vector<DetermineBackgroundLevels> background, std::string BaseFileName, int NumberofExperiments, std::vector<std::string>& ExperimentNames, std::string whichchr, int BinSize, PrDes::RENFileInfo& reInfo){

    fLog << "will print proximities of NegCtrls" << std::endl; 
    //Probe to Distal Interactions
    std::string FileName2, washUFileName;
    int enoughpairs, count = 0;
    double p_value;
    
    FileName2.append(BaseFileName);
    FileName2.append(".");
    FileName2.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	//FileName2.append(".");
    //FileName2.append(whichchr);
    FileName2.append(".Proximities.Probe_Distal.NegCtrls.");
    FileName2.append(reInfo.currTime);
    FileName2.append(".txt");
    std::ofstream outf2(FileName2.c_str());
    

    outf2 << "RefSeqName" << '\t' << "TranscriptName" << '\t' << "Feature_ID" << '\t'
    << "Feature_Chr" << '\t' << "Feature_Start" << '\t' << "Feature_End" << '\t' << "Annotation" << '\t' <<  "Strand" << '\t';
    
    outf2 << "Interactor_Chr" << '\t' << "Interactor_Start" << '\t' << "Interactor_End" << '\t' << "distance" ;
    for (int e = 0; e < NumberofExperiments; ++e)
        outf2 << '\t'<< ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_p_value" << '\t' << ExperimentNames[e] << "_StrandCombination";
    outf2 << std::endl;

    std::map<std::string, FeatureStruct>::iterator featiter;
    std::map<std::string, FeatureStruct>::iterator featiter2;
    std::map<int, Junction >::iterator it; //key:REpos
    std::vector < SignalStruct_CTX >::iterator itx;
    std::map<int, Junction >::iterator itt;
  //initialise memory
      for(int i=0 ; i< Features.size(); ++i){
		auto featiter = Features.begin();
		std::advance(featiter, i);
		for (it = featiter->second.proximities.junctions.begin(); it != featiter->second.proximities.junctions.end(); ++it){
            if(featiter->second.FeatureType == 3){
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
                if(it->second.reportit){
                    it->second.p_val.reserve(NumberofExperiments);
                }
            }
        }
    }
  
  
    //Parallelise p-val calculation
    //#pragma omp parallel for 
    for(int i=0 ; i< Features.size(); ++i){
		auto featiter = Features.begin() ;
		std::advance(featiter, i);
		
		for (it = featiter->second.proximities.junctions.begin(); it != featiter->second.proximities.junctions.end(); ++it){
            if(featiter->second.FeatureType == 3){
                
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
                
                if(it->second.reportit){
                    
                    it->second.distance = it->first - featiter->second.start;
                    
                    int bin;
                    bin = abs(it->second.distance) / BinSize;
                    for (int e = 0; e < NumberofExperiments; ++e){
                        if (CALCULATE_P_VALUES && it->second.paircount[e] >= MinNumberofSupportingPairs){
                            if(background[e].bglevels.smoothed.find(bin) == background[e].bglevels.smoothed.end()){
								int i=0;
								while(background[e].bglevels.smoothed.find(bin-i)==background[e].bglevels.smoothed.end()){
									if(bin==0)
										--i;
									else
										++i;
								}
								if(bin-i>=0 && bin-i<=background[e].bglevels.smoothed.rbegin()->first){
									background[e].bglevels.smoothed[bin] = background[e].bglevels.smoothed[bin-i];
									background[e].bglevels.smoothed_stdev[bin] = background[e].bglevels.smoothed_stdev[bin-i];
								}
								else{
									background[e].bglevels.smoothed[bin] = 0;
									background[e].bglevels.smoothed_stdev[bin] = 1;
								}
							}
                            
							if (background[e].bglevels.smoothed[bin] == 0.0)
								it->second.p_val.push_back(1.0);
						    else{
						         double t = 1.0 - alglib::normaldistribution(( it->second.paircount[e] - background[e].bglevels.smoothed[bin]) / background[e].bglevels.smoothed_stdev[bin]);
						         it->second.p_val.push_back(t); //z-score      
							}
							
						}
                        else
                            it->second.p_val.push_back(1.0);
                        
                    }
                    
                }
            }
        }
    }
    //write to file
    for(featiter = Features.begin() ; featiter != Features.end(); ++featiter){
       if (count % 2500 == 0)
            fLog << count << " probes processed in Probe-Distal Neg Ctrls" << std::endl;
        ++count;
        for (it = featiter->second.proximities.junctions.begin(); it != featiter->second.proximities.junctions.end(); ++it){
            if(featiter->second.FeatureType == 3){
                
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
                
                if(it->second.reportit){
                    outf2 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' 
                    << featiter->second.chr  << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                    << featiter->second.FeatureType << '\t' << featiter->second.strand << '\t';
                    
                    outf2 << featiter->second.chr << '\t' << it->first << '\t' << it->second.refragend << '\t';
                    
                    outf2 << it->second.distance ;
  
                    for (int e = 0; e < NumberofExperiments; ++e){
                                     
                        outf2 << '\t'<< it->second.paircount[e] << '\t' << it->second.p_val[e] << '\t';
                        for(int s = 0; s < 3; ++s)
							 outf2 << it->second.strandcombination[(e*4)+s] << "_";
						outf2 << it->second.strandcombination[(e*4)+3];  
                    }
                    outf2 << std::endl;
                    
                }
            }
        }
           
        for (itx = featiter->second.proximities_ctx.begin(); itx != featiter->second.proximities_ctx.end(); ++itx){
            if(featiter->second.FeatureType == 3) {
                for(itt = itx->junctions_ctx.begin(); itt != itx->junctions_ctx.end(); ++itt){
                    enoughpairs = CheckSupportingPairs(itt->second.paircount, NumberofExperiments);
                    if(enoughpairs){
                        outf2 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' 
                        << featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                        << featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                        
                        outf2 << itx->maptochrname << '\t' << itt->first << '\t' << itt->second.refragend << '\t' << -1 ;
                        
                        
                        
                        for (int e = 0; e < NumberofExperiments; ++e){
                            outf2 << '\t' << itt->second.paircount[e] << '\t' << -1 << '\t';
                            for(int s = 0; s < 3; ++s)
								outf2 << itt->second.strandcombination[(e*4) + s] << "_";
							outf2 << itt->second.strandcombination[(e*4) + 3];     
                        }
                        outf2 << std::endl;
                    }
                }
            }
        }
    }
    
    outf2.close();
    count = 0;
}



void DetectInteractions::CalculatePvalAndPrintInteractionsProbeProbe_NegCtrls(ProbeSet& prs, std::vector<DetermineBackgroundLevels> background, std::string BaseFileName, int NumberofExperiments, std::vector<std::string>& ExperimentNames, std::string whichchr, int BinSizeProbeProbe, PrDes::RENFileInfo& reInfo){
   
	int bin;
	int flag;
	int count = 0;
    std::string FileName4;
   
    FileName4.append(BaseFileName);
    FileName4.append(".");
    FileName4.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	//FileName4.append(".");
    //FileName4.append(whichchr);
    FileName4.append(".Proximities.Probe_Probe.NegCtrls.");
    FileName4.append(reInfo.currTime);
    FileName4.append(".txt");
    std::ofstream outf4(FileName4.c_str());
    
    outf4<< "RefSeqName_1" << '\t' << "TranscriptName_1" << '\t' << "Feature_ID_1" << '\t' << "FeatureChr_1" << '\t' << "FeatureStart_1" << '\t' << "FeatureEnd_1" << '\t' << "Annotation_1" << '\t' <<  "Strand_1" << '\t';
    
    outf4 << "RefSeqName_2" << '\t' << "TranscriptName_2" << '\t' << "Feature_ID_2" << '\t' << "FeatureChr_2" << '\t' << "FeatureStart_2" << '\t' << "FeatureEnd_2" << '\t' << "Annotation_2" << '\t' <<  "Strand_2" << '\t';
   
    outf4 << "abs(Distance)";
    for (int e = 0; e < NumberofExperiments; ++e)
        outf4 << '\t' << ExperimentNames[e] << "_SuppPairs"<< '\t' << ExperimentNames[e]<<"_pval" << '\t' << ExperimentNames[e]<<"_StrandCombination";
    outf4 << std::endl;
    
    std::vector< FeattoFeatSignalStruct >::const_iterator itff; //first: REpos, second: signal
    std::string f;
    for(auto featiter = Features.begin() ; featiter != Features.end(); ++featiter){
        if (count % 2500 == 0)
            fLog << count << " probes processed in Probe-Probe Neg Ctrls" << std::endl;
        ++count;
        for (auto itff = featiter->second.Inter_feature_ints.begin(); itff != featiter->second.Inter_feature_ints.end(); ++itff){
            f = itff->interacting_feature_id;
            auto featiter2 = Features.find(f);
            if((featiter->second.FeatureType == 3 || featiter2->second.FeatureType == 3) && featiter->second.TranscriptName != featiter2->second.TranscriptName && (abs(featiter->second.start - featiter2->second.start) >= MinimumJunctionDistance)){
                outf4 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t'
                << featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                << featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                
                outf4 << featiter2->second.Name << '\t' << featiter2->second.TranscriptName << '\t' << featiter2->first << '\t' 
                << featiter2->second.chr << '\t' << featiter2->second.start << '\t' << featiter2->second.end << '\t'
                << featiter2->second.FeatureType << '\t' << featiter2->second.strand << '\t';
                
                
                bin=0;
                flag=0;
                
                if(featiter->second.chr == featiter2->second.chr){
                    outf4 << abs(featiter->second.start - featiter2->second.start) << '\t';
                    bin = abs(featiter->second.start - featiter2->second.start) / BinSizeProbeProbe;   
                    flag=1;
				}
                else
                    outf4 << -1 ;
                    
                			    
                for (int e = 0; e < NumberofExperiments; ++e){
                    outf4 << '\t' << itff->signal[e] << '\t';
                    
                    ///////Pval calc begins
                                      					
					if (CALCULATE_P_VALUES && flag){
						if(background[e].bglevelsProbeProbe.smoothed.find(bin) == background[e].bglevelsProbeProbe.smoothed.end()){
							int i=0;
							while(background[e].bglevelsProbeProbe.smoothed.find(bin-i)==background[e].bglevelsProbeProbe.smoothed.end()){
								if(bin==0)
									--i;
								else
									++i;
							}
							if(bin-i>=0 && bin-i<=background[e].bglevelsProbeProbe.smoothed.rbegin()->first){
								background[e].bglevelsProbeProbe.smoothed[bin] = background[e].bglevelsProbeProbe.smoothed[bin-i];
								background[e].bglevelsProbeProbe.smoothed_stdev[bin] = background[e].bglevelsProbeProbe.smoothed_stdev[bin-i];
							}
							else{
								background[e].bglevelsProbeProbe.smoothed[bin] = 0;
								background[e].bglevelsProbeProbe.smoothed_stdev[bin] = 1;
							}
						}
						
						if (background[e].bglevelsProbeProbe.smoothed[bin] == 0.0)
							outf4 << "1.0";
							
						else{
							double t= 1.0 - alglib::normaldistribution(( itff->signal[e] - background[e].bglevelsProbeProbe.smoothed[bin]) / background[e].bglevelsProbeProbe.smoothed_stdev[bin]);
							outf4 << t ;
						}
					}
					else        
                       outf4 << "1.0" ;
                       
					
                   /////Pval calc ends
                   outf4 << "\t" ;
                   for(int s = 0; s < 3; ++s)
					   outf4 << itff->strandcombination[(e*4) + s] << "_";
				   outf4 << itff->strandcombination[(e*4) + 3];

                }
                outf4 << std::endl;
            }
        }
    }
 	outf4.close();
}
