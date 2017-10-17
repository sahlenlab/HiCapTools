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
//  PrintUsage.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_PRINTUSE_INC_H_
#define HCT_PRINTUSE_INC_H_

#include <iostream>
#include <iomanip>


static int print_usage()
{
    std::cerr<<std::endl;
    std::cerr<<"Program: HiCapTools (A software suite for Probe Design and Proximity Detection for targeted chromosome conformation capture applications)"<<std::endl;
    std::cerr<<"Contact: Pelin Sahlén <pelin.akan@scilifelab.se>"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<<"Usage:  HiCapTools <option> [arguments]"<<std::endl;
    std::cerr<< "Options:  "<<std::endl;
    std::cerr<<std::endl;
    std::cerr<< "    ProbeDesigner"<<std::endl;
    std::cerr<< '\t'<< "Required Arguments:"<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<"-o or --option"<< "'FeatureProbes' to create probes targeting Features for chromosomes or 'NegativeControls' to create negative control probes for all chromosomes"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<< '\t'<< "Optional Arguments:"<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<"-c or --chr"<< "the Chromosome to process in the format chrN, where N can be the name/number of the chromosome or All if processing all available chromosomes."<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<" "<< "'-c' is required input for option 'FeatureProbes' and is not required for 'NegativeControls'"<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<"-config"<< "the path to the Probe config file, if changed from default in bin/config directory"<<std::endl;
    std::cerr<<std::endl;
    
    std::cerr<< "    ProximityDetector"<<std::endl;
    std::cerr<< '\t'<<"Required Arguments:"<<std::endl;
    std::cerr<< '\t'<<std::left <<std::setw(25)<<"-m or --outputmode"<<"'ComputeStatsOnly' to compute only the stats or 'PrintProximities' to also find and print proximities"<<std::endl;
    std::cerr<<std::endl;
    std::cerr<< '\t'<< "Optional Arguments:"<<std::endl;
    std::cerr<< '\t'<<std::left<< std::setw(25)<<"-c or --chr" <<"the Chromosome to process in the format chrN, where N can be the name/number of the chromosome or All if processing all available chromosomes"<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<" "<< "'-c' is required input for option 'FeatureProbes' and is not required for 'NegativeControls'"<<std::endl;
    std::cerr<< '\t'<<std::left<< std::setw(25)<<"-p or --proximitytype"<<"'Neg' to print only negative control Probe proximities or 'NonNeg' to print only Feature Probe proximities or 'Both' to print both "<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<" "<< "'-p' is required input for mode 'PrintProximities' and is not required for 'ComputeStatsOnly'"<<std::endl;
    std::cerr<< '\t'<<std::left<<std::setw(25)<<"-config"<< "the path to the config file, if changed from default in bin/config directory"<<'\n'<<'\t'<<std::setw(25)<<" "<<std::endl;
    return 0;
}

#endif // HCT_PRINTUSE_INC_H_
