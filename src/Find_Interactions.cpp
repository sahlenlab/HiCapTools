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

    fLog << "will print interactions " << std::endl;
//Probe to Distal Interactions
	std::string FileName, FileName2;
    
    int enoughpairs, count = 0;
    double p_value;
    
    
	FileName.append(BaseFileName);
    FileName.append(".");
	FileName.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileName.append(".");
    FileName.append(whichchr);
    FileName.append(".Interactions.Probe_Distal.");
	FileName.append(reInfo.currTime);
	FileName.append(".txt");
	std::ofstream outf1(FileName.c_str());
    
	
    outf1 << "RefSeqName" << '\t' << "TranscriptName" << '\t' << "Feature_ID" << '\t' << "Probe_ID" << '\t' << "Feature Chr" << '\t' << "Feature Start" << '\t' << "Feature End" << '\t' << "Annotation" << '\t' <<  "Strand" << '\t';
    
    outf1 << "Interactor Chr" << '\t' << "Interactor Start" << '\t' << "Interactor End" << '\t' << "distance" << '\t';
    
	for (int e = 0; e < NumberofExperiments; ++e)
		outf1 << ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_p_value" ;
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
						   background[e].bglevels.smoothed[bin] = 0;
						   background[e].bglevels.smoothed_stdev[bin] = 1;
						   it->second.p_val.push_back(1.0);
						}
						else{
							if (background[e].bglevels.smoothed_stdev[bin] == 0.0)
								it->second.p_val.push_back(1.0);
								
							else{
								double t= 1.0 - alglib::normaldistribution(( it->second.paircount[e] - background[e].bglevels.smoothed[bin]) / background[e].bglevels.smoothed_stdev[bin]); 
								it->second.p_val.push_back(t);
							}
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
               
                outf1 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' << featiter->second.probe_name << '\t'
                      << featiter->second.chr  << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                      << featiter->second.FeatureType << '\t' << featiter->second.strand << '\t';
                
                outf1 << featiter->second.chr << '\t' << it->first << '\t' << it->second.refragend << '\t';
               
               
                outf1 << it->second.distance << '\t';
                
                for (int e = 0; e < NumberofExperiments; ++e){
                    outf1 << it->second.paircount[e] << '\t' << it->second.p_val[e] ;
                    
                }
                outf1 << std::endl;
            }
		  }
        }
        
        for (itx = featiter->second.proximities_ctx.begin(); itx != featiter->second.proximities_ctx.end(); ++itx){
            if(featiter->second.FeatureType != 3 ){
				for(itt = itx->junctions_ctx.begin(); itt != itx->junctions_ctx.end(); ++itt){
					
					enoughpairs = CheckSupportingPairs(itt->second.paircount, NumberofExperiments);
					if(enoughpairs){
                    
						outf1 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' << featiter->second.probe_name << '\t'
							<< featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
							<< featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                          
						outf1 << itx->maptochrname << '\t' << itt->first << '\t' << itt->second.refragend << '\t' << -1 << '\t';
						for (int e = 0; e < NumberofExperiments; ++e){
							outf1 << itt->second.paircount[e] << '\t' << -1 ;
							
						}
						outf1 << std::endl;
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
	std::string FileName3;
	
	FileName3.append(BaseFileName);
    FileName3.append(".");
	FileName3.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileName3.append(".");
    FileName3.append(whichchr);
	FileName3.append(".Interactions.Probe_Probe.");
	FileName3.append(reInfo.currTime);
	FileName3.append(".txt");
	
		
	std::ofstream outf3(FileName3.c_str());

    
	outf3 << "RefSeqName_1" << '\t' << "TranscriptName_1" << '\t' << "Feature_ID_1" << '\t' << "Probe_ID_1" << '\t' << "FeatureChr_1" << '\t' << "FeatureStart_1" << '\t' << "FeatureEnd_1"
    << '\t' << "Annotation_1" << '\t' <<  "Strand_1" << '\t';
    
    outf3 << "RefSeqName_2" << '\t' << "TranscriptName_2" << '\t' << "Feature_ID_2" << '\t' << "Probe_ID_2" << '\t' << "FeatureChr_2" << '\t' << "FeatureStart_2" << '\t' << "FeatureEnd_2"
    << '\t' << "Annotation_2" << '\t' <<  "Strand_2" << '\t';
   
    outf3 << "abs(Distance)";

    for (int e = 0; e < NumberofExperiments; ++e)
        outf3 << '\t' << ExperimentNames[e] << "_SuppPairs" << ExperimentNames[e] << "_p_value";
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
				
				
								
                outf3 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' << featiter->second.probe_name << '\t'
                      << featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                      << featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                
                outf3 << featiter2->second.Name << '\t' << featiter2->second.TranscriptName << '\t' << featiter2->first << '\t' << featiter2->second.probe_name << '\t'
                      << featiter2->second.chr << '\t' << featiter2->second.start << '\t' << featiter2->second.end << '\t'
                      << featiter2->second.FeatureType << '\t' << featiter2->second.strand << '\t';
                
                bin=0;
                flag=0;
                if(featiter->second.chr == featiter2->second.chr){
                    outf3 << abs(featiter->second.start - featiter2->second.start) << '\t';
                    bin=abs(featiter->second.start - featiter2->second.start) / BinSizeProbeProbe;
                    flag=1;   
				}
                else
                    outf3 << -1 << '\t';
                    
                
                for (int e = 0; e < NumberofExperiments; ++e){
                    outf3 << itff->signal[e] << '\t';
                    ///////Pval calc begins
                                      					
					if (CALCULATE_P_VALUES && flag){
						if(background[e].bglevelsProbeProbe.smoothed.find(bin) == background[e].bglevelsProbeProbe.smoothed.end()){
							background[e].bglevelsProbeProbe.smoothed[bin] = 0;
							background[e].bglevelsProbeProbe.smoothed_stdev[bin] = 1;
							outf3 << "1.0" ;
						}
						else{
							if (background[e].bglevelsProbeProbe.smoothed_stdev[bin] == 0.0)
								outf3 << "1.0"; 
							
							else{
								double t= 1.0 - alglib::normaldistribution(( itff->signal[e] - background[e].bglevelsProbeProbe.smoothed[bin]) / background[e].bglevelsProbeProbe.smoothed_stdev[bin]);
								outf3 << t ;
							}
						}
					}
					else
                       outf3 << "1.0";
                       
					
                   /////Pval calc ends
                  
                  
                }
                outf3 << std::endl;
            }
  
        }
    }
	outf3.close();

}

void DetectInteractions::CalculatePvalAndPrintInteractionsProbeDistal_NegCtrls(ProbeSet& prs, std::vector<DetermineBackgroundLevels> background, std::string BaseFileName, int NumberofExperiments, std::vector<std::string>& ExperimentNames, std::string whichchr, int BinSize, PrDes::RENFileInfo& reInfo){

    fLog << "will print interactions of NegCtrls" << std::endl; ////////////////////
    //Probe to Distal Interactions
    std::string FileName, FileName2;
    int enoughpairs, count = 0;
    double p_value;
    
    FileName2.append(BaseFileName);
    FileName2.append(".");
    FileName2.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileName2.append(".");
    FileName2.append(whichchr);
    FileName2.append(".Interactions.Probe_Distal.NegCtrls.");///////////////////////////////////
    FileName2.append(reInfo.currTime);
    FileName2.append(".txt");
    std::ofstream outf2(FileName2.c_str());
    
    outf2 << "RefSeqName" << '\t' << "TranscriptName" << '\t' << "Feature_ID" << '\t' << "Probe_ID" << '\t' << "Feature Chr" << '\t' << "Feature Start" << '\t' << "Feature End" << '\t' << "Annotation" << '\t' <<  "Strand" << '\t';
    
    
    outf2 << "Interactor Chr" << '\t' << "Interactor Start" << '\t' << "Interactor End" << '\t' << "distance" << '\t';
    for (int e = 0; e < NumberofExperiments; ++e)
        outf2 << ExperimentNames[e] << "_SuppPairs" << '\t' << ExperimentNames[e] << "_p_value";
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
            if(featiter->second.FeatureType == 3){///////////////////
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
            if(featiter->second.FeatureType == 3){////////////////////////////////
                
                it->second.reportit=0;
                for (int t = 0; t < NumberofExperiments; ++t){
					if (it->second.paircount[t] >= MinNumberofSupportingPairs){
						it->second.reportit = 1;
						break;
					}
				}
                
                if(it->second.reportit){
                    
                    it->second.distance = it->first - featiter->second.start;/////////////////////////
                    
                    int bin;
                    bin = abs(it->second.distance) / BinSize;
                    for (int e = 0; e < NumberofExperiments; ++e){
                        if (CALCULATE_P_VALUES && it->second.paircount[e] >= MinNumberofSupportingPairs){
                            
                            if(background[e].bglevels.smoothed.find(bin) == background[e].bglevels.smoothed.end()){
								background[e].bglevels.smoothed[bin] = 0;
								background[e].bglevels.smoothed_stdev[bin] = 1;
								it->second.p_val.push_back(1.0);
							}
							else{
								if (background[e].bglevels.smoothed_stdev[bin] == 0.0)
									it->second.p_val.push_back(1.0);
						        else{
						            double t = 1.0 - alglib::normaldistribution(( it->second.paircount[e] - background[e].bglevels.smoothed[bin]) / background[e].bglevels.smoothed_stdev[bin]);
						            it->second.p_val.push_back(t); //z-score      
								}
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
                    outf2 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' << featiter->second.probe_name << '\t'
                    << featiter->second.chr  << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                    << featiter->second.FeatureType << '\t' << featiter->second.strand << '\t';
                    
                    outf2 << featiter->second.chr << '\t' << it->first << '\t' << it->second.refragend << '\t';
                    
                    outf2 << it->second.distance << '\t';
                    
                    for (int e = 0; e < NumberofExperiments; ++e){
                                     
                        outf2 << it->second.paircount[e] << '\t' << it->second.p_val[e] ;
                        
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
                        outf2 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' << featiter->second.probe_name << '\t'
                        << featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                        << featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                        
                        outf2 << itx->maptochrname << '\t' << itt->first << '\t' << itt->second.refragend << '\t' << -1 << '\t';
                        for (int e = 0; e < NumberofExperiments; ++e){
                            outf2 << itt->second.paircount[e] << '\t' << -1 ;
                            
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
	FileName4.append(".");
    FileName4.append(whichchr);
    FileName4.append(".Interactions.Probe_Probe.NegCtrls.");
    FileName4.append(reInfo.currTime);
    FileName4.append(".txt");
    std::ofstream outf4(FileName4.c_str());
    
    outf4<< "RefSeqName_1" << '\t' << "TranscriptName_1" << '\t' << "Feature_ID_1" << '\t' << "Probe_ID_1" << '\t' << "FeatureChr_1" << '\t' << "FeatureStart_1" << '\t' << "FeatureEnd_1" << '\t' << "Annotation_1" << '\t' <<  "Strand_1" << '\t';
    
    outf4 << "RefSeqName_2" << '\t' << "TranscriptName_2" << '\t' << "Feature_ID_2" << '\t' << "Probe_ID_2" << '\t' << "FeatureChr_2" << '\t' << "FeatureStart_2" << '\t' << "FeatureEnd_2" << '\t' << "Annotation_2" << '\t' <<  "Strand_2" << '\t';
   
    outf4 << "abs(Distance)";
    for (int e = 0; e < NumberofExperiments; ++e)
        outf4 << '\t' << ExperimentNames[e] << "_SuppPairs"<< '\t' << ExperimentNames[e]<<"_pval";
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
                outf4 << featiter->second.Name << '\t' << featiter->second.TranscriptName << '\t' << featiter->first << '\t' << featiter->second.probe_name << '\t'
                << featiter->second.chr << '\t' << featiter->second.start << '\t' << featiter->second.end << '\t'
                << featiter->second.FeatureType << '\t' <<  featiter->second.strand << '\t';
                
                outf4 << featiter2->second.Name << '\t' << featiter2->second.TranscriptName << '\t' << featiter2->first << '\t' << featiter2->second.probe_name << '\t'
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
                    outf4 << -1 << '\t';
                    
                			    
                for (int e = 0; e < NumberofExperiments; ++e){
                    outf4 << itff->signal[e] << '\t';
                    
                    ///////Pval calc begins
                                      					
					if (CALCULATE_P_VALUES && flag){
						if(background[e].bglevelsProbeProbe.smoothed.find(bin) == background[e].bglevelsProbeProbe.smoothed.end()){
							background[e].bglevelsProbeProbe.smoothed[bin] = 0;
							background[e].bglevelsProbeProbe.smoothed_stdev[bin] = 1;
							outf4 << "1.0"; 
						}
						else{
							if (background[e].bglevelsProbeProbe.smoothed_stdev[bin] == 0.0)
								outf4 << "1.0";
							
							else{
								double t= 1.0 - alglib::normaldistribution(( itff->signal[e] - background[e].bglevelsProbeProbe.smoothed[bin]) / background[e].bglevelsProbeProbe.smoothed_stdev[bin]);
								outf4 << t ;
							}
						}
					}
					else        
                       outf4 << "1.0" ;
                       
					
                   /////Pval calc ends

                }
                outf4 << std::endl;
            }
        }
    }
 	outf4.close();
}
