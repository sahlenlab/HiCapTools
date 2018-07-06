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
//  BackgroundInteractionFrequency.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/rolling_mean.hpp>
using namespace boost::accumulators;

#include "BackgroundInteractionFrequency.h"
#include "Global.h"
#include "linear.h"
#include <fstream>
#include <cmath>


///For Probe-Distal and Probe-Probe
void DetermineBackgroundLevels::CalculateMeanandStdRegress(std::string eName, int ExperimentNo, std::string DesignName, BG_signals& bglevelsloc, int binsize, std::string whichProx, int MinimumJunctionDistance, OutStream& log, int WindowSizeloc){


	expName= eName; //set experiment name
	
	std::map< int, int > nofentries_perBin; // required for calculating mean
	std::map< int, int > signal_square; // required for calculating stdev
    double sum_square, variance;
	int distance;
    std::string feature_id;

	
    std::map< int, Junction >::const_iterator iter;
	for(int i = 0; i < Design_NegCtrl[DesignName].Probes.size(); i++){
        feature_id = Design_NegCtrl[DesignName].Probes[i].feature_id;
        //////////Probe-distal
        if(whichProx=="ProbeDistal"){
			for (iter = Features[feature_id].proximities.junctions.begin(); iter != Features[feature_id].proximities.junctions.end(); ++iter){
				
				distance = iter->first - Design_NegCtrl[DesignName].Probes[i].end;
            
				int bin = abs(distance) / binsize;

				if(iter->second.paircount[ExperimentNo] > 0){
           
					if(bglevelsloc.mean.find(bin) == bglevelsloc.mean.end())
						bglevelsloc.mean[bin] = iter->second.paircount[ExperimentNo];
					else
						bglevelsloc.mean[bin] = bglevelsloc.mean[bin] + iter->second.paircount[ExperimentNo];
                
					if(nofentries_perBin.find(bin) == nofentries_perBin.end()){
						nofentries_perBin[bin] = 1;
						signal_square[bin] = (iter->second.paircount[ExperimentNo])*(iter->second.paircount[ExperimentNo]);
					}
					else{
						nofentries_perBin[bin] = nofentries_perBin[bin] + 1;
						signal_square[bin] = (signal_square[bin] + ((iter->second.paircount[ExperimentNo])*(iter->second.paircount[ExperimentNo])));
					}
				}
			}
		}
		///////////Probe-Probe
        else if(whichProx=="ProbeProbe"){
			for (auto iter = Features[feature_id].Inter_feature_ints.begin(); iter != Features[feature_id].Inter_feature_ints.end(); ++iter){
								
				distance = Features[feature_id].start - Features[(*iter).interacting_feature_id].start;
				
				if((Features[feature_id].FeatureType == 3 && Features[(*iter).interacting_feature_id].FeatureType == 3) && Features[feature_id].TranscriptName != Features[(*iter).interacting_feature_id].TranscriptName && (abs(distance) >= MinimumJunctionDistance)){
            
					int bin = abs(distance) / binsize; 
					if((*iter).signal[ExperimentNo] > 0){
						if(bglevelsloc.mean.find(bin) == bglevelsloc.mean.end())
							bglevelsloc.mean[bin] = (*iter).signal[ExperimentNo];
						else
							bglevelsloc.mean[bin] = bglevelsloc.mean[bin] + (*iter).signal[ExperimentNo];
                
						if(nofentries_perBin.find(bin) == nofentries_perBin.end()){
							nofentries_perBin[bin] = 1;
							signal_square[bin] = ((*iter).signal[ExperimentNo])*((*iter).signal[ExperimentNo]);
						}
						else{
							nofentries_perBin[bin] = nofentries_perBin[bin] + 1;
							signal_square[bin] = (signal_square[bin] + (((*iter).signal[ExperimentNo])*((*iter).signal[ExperimentNo])));
						}
					}
				}
			}
		}        
    }
    
// Calculate Mean and stdev
	std::map< int, double>::iterator it; // iterator for bin signals
	for (it = bglevelsloc.mean.begin(); it != bglevelsloc.mean.end(); ++it){
        bglevelsloc.samplesize[it->first] = nofentries_perBin[it->first];
        if(bglevelsloc.samplesize[it->first] == 0){
            it->second = 1; //Mean
        }
        else{
            it->second /= nofentries_perBin[it->first]; //Mean
        }
        if(bglevelsloc.samplesize[it->first] <= 1){
            bglevelsloc.stdev[it->first] = 0;
        }
        else{
            sum_square = ((it->second)*(it->second));
            sum_square /= (double(nofentries_perBin[it->first])); // Average Sum Squared
            variance = signal_square[it->first] - sum_square;
            variance /= (nofentries_perBin[it->first] -1);
            bglevelsloc.stdev[it->first] = sqrt(variance); //stdev
        }
	}
    
    int first25kb = ceil(25000/(double)binsize); // to exclude bins of first 25 kb from smoothing 
    
    boost::accumulators::accumulator_set<double, stats<tag::rolling_mean> > acc(tag::rolling_window::window_size = (WindowSizeloc*2));
    boost::accumulators::accumulator_set<double, stats<tag::rolling_mean> > acc2(tag::rolling_window::window_size = (WindowSizeloc*2));
    
    int w = 0, z = 0, s = 0, a=0;
    std::deque< double > dm, ds, db, dwm, dws;

    
    // Put values into a queue
    for (it = bglevelsloc.mean.begin(); it != bglevelsloc.mean.end(); ++it){
	
        dm.push_back(it->second);
        ds.push_back(bglevelsloc.stdev[it->first]);
        db.push_back(it->first);
	
    }
    //Push the first WindowSize-1 elements into a queue
	if(dm.empty())
		log<<"No Element in dm Outloop"<<std::endl;
		
	bool firstpass=true;
    for(w = 0; w < (dm.size()); ++w){
		if(w < first25kb){
			bglevelsloc.smoothed[db[w]] = dm[w];
			bglevelsloc.smoothed_stdev[db[w]] = ds[w];
		}
		else{
			boost::accumulators::accumulator_set<double, stats<tag::rolling_mean> > acc(tag::rolling_window::window_size = (WindowSizeloc));
			boost::accumulators::accumulator_set<double, stats<tag::rolling_mean> > acc2(tag::rolling_window::window_size = (WindowSizeloc));
			
			if(w >= first25kb && w < first25kb + WindowSizeloc - 1){
				
	
				/**
				if(dwm.size()==WindowSizeloc){
					for(z = 0; z < WindowSizeloc;++z){
						std::cout<<dwm[z]<<std::endl;
					}
					dwm.pop_front();
				}
			**/
				//change to rolling mean if needed
				bglevelsloc.smoothed[db[w]] = dm[w];
				bglevelsloc.smoothed_stdev[db[w]] = ds[w];
				

			}
			else {
			    if(firstpass){
			        for(auto it=w-(WindowSizeloc/2); it<w+WindowSizeloc/2;++it){ //for the first bin that is to be smoothed
			            dwm.push_back(dm[it]);
			            dws.push_back(ds[it]);
			        }
			        firstpass=false;
			    }
				bool ifrolling=false;
				
				//To check that the window stops at the last bin
				if(w+WindowSizeloc/2 < dm.size()){
				    dwm.push_back(dm.at(w+WindowSizeloc/2));
				    dws.push_back(ds.at(w+WindowSizeloc/2));
				}
				
			
				if(dwm.size()==WindowSizeloc){
					ifrolling=true;
					for(z = 0; z < WindowSizeloc;++z){
						acc(dwm[z]);
						acc2(dws[z]);
					}
				}
				if(ifrolling){
					bglevelsloc.smoothed[db[w]] = rolling_mean(acc);
					bglevelsloc.smoothed_stdev[db[w]] = rolling_mean(acc2); 
					
				}
				else{
					bglevelsloc.smoothed[db[w]] = dm[w];
					bglevelsloc.smoothed_stdev[db[w]] = ds[w]; 
					
				}
		
				dwm.pop_front();
				dws.pop_front();
				
			}
		}
	}
}

void DetermineBackgroundLevels::PrintBackgroundFrequency(int bSize, int bSizePP, PrDes::RENFileInfo& reInfo){
	
	std::string FileNamePE, FileNamePP;
	
	FileNamePE.append(expName);
	FileNamePE.append(".");
	FileNamePE.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	FileNamePE.append(".Probe_Distal.BackgroundLevels.");
	FileNamePE.append(reInfo.currTime);
	FileNamePE.append(".txt");
	
	FileNamePP.append(expName);
	FileNamePP.append(".");
	FileNamePP.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));	
	FileNamePP.append(".Probe_Probe.BackgroundLevels.");
	FileNamePP.append(reInfo.currTime);
	FileNamePP.append(".txt");
	
	std::ofstream outfPE(FileNamePE.c_str());
	std::ofstream outfPP(FileNamePP.c_str());
    outfPE << "Distance Bin" << '\t' << "Interaction Distance" << '\t' << "Mean " << '\t' << "Stdev" << '\t' << "Sample Size" << '\t' << "Mean (Smoothed)" << '\t' << "StDev (Smoothed)" <<  std::endl;
    
	for (auto it = bglevels.smoothed.begin(); it != bglevels.smoothed.end(); ++it){
        outfPE << (it->first) << '\t'<< (it->first)*bSize << '\t' << bglevels.mean[it->first] << '\t' << bglevels.stdev[it->first] << '\t' << bglevels.samplesize[it->first] << '\t' << it->second << '\t' << bglevels.smoothed_stdev[it->first] << std::endl;	
	}
	
	outfPP << "Distance Bin" << '\t' << "Interaction Distance" << '\t' << "Mean " << '\t' << "Stdev" << '\t' << "Sample Size" << '\t' << "Mean (Smoothed)" << '\t' << "StDev (Smoothed)" <<  std::endl;
	
	for (auto it = bglevelsProbeProbe.smoothed.begin(); it != bglevelsProbeProbe.smoothed.end(); ++it){
        outfPP << (it->first) << '\t' << (it->first)*bSizePP << '\t' << bglevelsProbeProbe.mean[it->first] << '\t' << bglevelsProbeProbe.stdev[it->first] << '\t' << bglevelsProbeProbe.samplesize[it->first] << '\t' << it->second << '\t' << bglevelsProbeProbe.smoothed_stdev[it->first] << std::endl;	
	}
	
}

