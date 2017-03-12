//
//  main.cpp
//  PrDe
//
//  Created by Pelin Sahlen and Anandashankar Anil on 28/02/2017.
//  Copyright (c) 2015 Pelin Sahlen. All rights reserved.
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
	
	HiCapTools prde;
	
    if (argc < 2) {
        print_usage();
        return -1;
    }
    
    whichMod = argv[1];
    
    if(whichMod=="ProbeDesigner"){
		if (argc != 4) {
			print_usage();
			return -1;
		}
		else{
			if(std::string(argv[2])=="-c" || std::string(argv[2])=="--chr"){
				whichChr=argv[3];
				prde.ProbeDesignMain(whichChr);
			}
			else{
				print_usage();
				return -1;
			}
		}
		
	}
	
	else if(whichMod=="ProximityDetector"){
		if (argc != 8) {
			print_usage();
			return -1;
		}
		else{
			for(int i=2; i<argc; ++i){
				if(std::string(argv[i])=="-c" || std::string(argv[i])=="--chr"){
					whichChr=argv[i+1];
				}
				else if(std::string(argv[i])=="-m" || std::string(argv[i])=="--outputmode"){
					if(std::string(argv[i+1])=="ComputeStatsOnly" || std::string(argv[i+1]) == "PrintProximities")
						statsOption=argv[i+1];
					else{
						std::cout<<"Invalid argument for -m"<<std::endl;
						print_usage();
						return -1;
					}
				}
				else if(std::string(argv[i])=="-p" || std::string(argv[i])=="--proximitytype"){
					if(std::string(argv[i+1])=="Neg" || std::string(argv[i+1]) == "NonNeg"||std::string(argv[i+1])=="Both")
						printOption=argv[i+1];
						else{
						std::cout<<"Invalid argument for -p"<<std::endl;
						print_usage();
						return -1;
					}
				}
				else{
					print_usage();
					return -1;
				}
			}
			prde.ProxDetectMain(whichChr, statsOption, printOption);
			
		}
		
	}
	else{
		print_usage();
		return -1;
	}
}
