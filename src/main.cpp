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
	
	HiCapTools prde;
	
    if (argc < 2) {
        print_usage();
        return 1;
    }
    
    whichMod = argv[1];
    
    if(whichMod=="ProbeDesigner"){
		if (!(argc == 4 || argc == 6)){
			print_usage();
			return 1;
		}
		else{
			for(int i=2; i<argc; i+=2){
				if((std::string(argv[i])=="-c" || std::string(argv[i])=="--chr")){
					whichChr=argv[i+1];					
				}
				else if(argc==6 && (std::string(argv[i])=="-config")){
					extraConfig=argv[i+1];
				}
				else{
				print_usage();
				return 1;
				}
			}
			prde.ProbeDesignMain(whichChr, extraConfig);
		}
		
	}
	
	else if(whichMod=="ProximityDetector"){
		if (!(argc == 8 || argc==10)) {
			print_usage();
			return 1;
		}
		else{
			for(int i=2; i<argc; i+=2){
				if(std::string(argv[i])=="-c" || std::string(argv[i])=="--chr"){
					whichChr=argv[i+1];
				}
				else if(std::string(argv[i])=="-m" || std::string(argv[i])=="--outputmode"){
					if(std::string(argv[i+1])=="ComputeStatsOnly" || std::string(argv[i+1]) == "PrintProximities")
						statsOption=argv[i+1];
					else{
						std::cout<<"!!Error!! Invalid argument for -m"<<std::endl;
						print_usage();
						return 1;
					}
				}
				else if(std::string(argv[i])=="-p" || std::string(argv[i])=="--proximitytype"){
					if(std::string(argv[i+1])=="Neg" || std::string(argv[i+1]) == "NonNeg"||std::string(argv[i+1])=="Both")
						printOption=argv[i+1];
					else{
						std::cout<<"!!Error!! Invalid argument for -p"<<std::endl;
						print_usage();
						return 1;
					}
				}
				else if(argc==10 && (std::string(argv[i])=="-config")){
						extraConfig=argv[i+1];
				}
				else{
					print_usage();
					return 1;
				}
					
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
