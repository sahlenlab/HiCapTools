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
//  main.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include <iostream>
#include <string>
#include <algorithm>

#include "PrintUsage.h"
#include "HiCapTools.h"


int main(int argc, const char * argv[]) {
	
	std::string whichMod;
	std::string whichChr;
	std::string proxType;
	std::string statsOption;
	std::string printOption;
	std::string extraConfig; 
	std::string probeOption;
	
	HiCapTools prde;
	
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    whichMod = argv[1];
    if (argc == 4 ||argc == 6 || argc == 8 || argc==10) {
		for(int i=2; i<argc; i+=2){
			if(std::string(argv[i])=="-o" || std::string(argv[i])=="--option"){
				probeOption=argv[i+1];
			}
			else if(std::string(argv[i])=="-config"){
				extraConfig=argv[i+1];
			}
			else if((std::string(argv[i])=="-c" || std::string(argv[i])=="--chr")){
				whichChr=argv[i+1];
			}
			else if(std::string(argv[i])=="-m" || std::string(argv[i])=="--outputmode"){
				statsOption=argv[i+1];
			}
			else if(std::string(argv[i])=="-p" || std::string(argv[i])=="--proximitytype"){
				printOption=argv[i+1];
			}
			else{
				print_usage();
				return 1;
			}
		}
	}
	else{
		print_usage();
		return 1;
	}
    
    if(whichMod=="ProbeDesigner"){
		if (!(argc == 4 || argc == 6 || argc == 8)){
			print_usage();
			return 1;
		}
		else{
			if(probeOption=="FeatureProbes" && argc <= 4){
				print_usage();
				return 1;
			}
			if(probeOption!="FeatureProbes" && probeOption!="NegativeControls" ){						
				std::cerr<<"!!Error!! Invalid argument for -o"<<std::endl;
				print_usage();
				return 1;
			}
		
			if(probeOption=="FeatureProbes" && (argc == 6 || argc == 8) ){										
				if(whichChr.empty()){
					std::cerr<<"!!Error!! option 'FeatureProbes' requires input chromosome!"<<std::endl;
					print_usage();
					return 1;
				}
			}
			else if(probeOption=="NegativeControls"){
				if(!whichChr.empty()){
					std::cerr<<"!!Warning!! option 'NegativeControls' does not require input chromosome! Negative Controls will be generated for all chromosomes"<<std::endl;
				}			
			}
			prde.ProbeDesignMain(whichChr, extraConfig, probeOption);
		}
		
	}
	
	else if(whichMod=="ProximityDetector"){
		if (!(argc == 4 ||argc == 6 || argc == 8 || argc==10)) {
			print_usage();
			return 1;
		}
		else{
				if(statsOption=="ComputeStatsOnly"){
					if(!whichChr.empty() || !printOption.empty()){
						std::cerr<<"!!Warning!! mode 'ComputeStatsOnly' does not require input chromosome or proximitytype! Stats will be generated for all chromosomes and proximity types"<<std::endl;
					}
				}
				else if(statsOption=="PrintProximities" && (argc == 6 || argc == 8 || argc == 10)){
					if(printOption=="Neg"){
						if(!whichChr.empty()){
							std::cerr<<"!!Warning!! proximitytype 'Neg' does not require input chromosome! Negative Control proximities will be generated for all chromosomes"<<std::endl;
						}
					}
					else if(printOption=="NonNeg" || printOption=="Both"){
						if(whichChr.empty()){
							std::cerr<<"!!Error!! proximitytype 'NonNeg' or 'Both' require input chromosome!"<<std::endl;
							print_usage();
							return 1;
						}
					}
				}
				else{
					print_usage();
					return 1;
				}
				
		}
		prde.ProxDetectMain(whichChr, statsOption, printOption, extraConfig);
			
	}
	else{
		print_usage();
		return 1;
	}
	
	return 0;		
}
