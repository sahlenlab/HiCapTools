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
//  CallHiCUP.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//
#include <vector>
#include <dirent.h>
#include "CallHiCUP.h"

int CallHiCUP::GenerateRestrictionFile(std::string hicupdigesterPath, std::string enzymeMotif, std::string genomeName, std::string fastaFilepath, std::string& digestFile){
	
	std::string cmd;
	DIR *dr;
	struct dirent *dirStruct;
	dr=opendir(".");
	
	std::string nameComp="Digest_"+genomeName+"_"+enzymeMotif+"_None";
	
	 while ((dirStruct = readdir(dr)) != NULL) {
		std::string fname = dirStruct->d_name;
		if(fname.find(nameComp) != std::string::npos){
			cLog<<"A digest file already exists :"<<fname<<std::endl;
			digestFile=fname;
			return 1;
		}
	}
	
	cmd=hicupdigesterPath+" --re1 "+enzymeMotif+" --genome "+genomeName+" "+fastaFilepath;
	
	char *cmdchar = &cmd[0];
    
    char buf[512];
    double a = -1;
    FILE *fp;
    
    if ((fp = popen(cmdchar, "r")) == NULL) {
        cLog<<"Error opening pipe!"<<std::endl;
        return -1;
    }
    
    while (fgets(buf, 512, fp) != NULL) {        
        a = atof(buf);
    }
    
    if(pclose(fp))  {
        cLog<<"HiCUP not found in bin/hicup OR HiCUP returned error"<<std::endl;
        return -1;
    }
    
    
    std::vector<std::string> retNames;
    
    
    while ((dirStruct = readdir(dr)) != NULL) {
		std::string fname = dirStruct->d_name;
		int count=0;
		if(fname.find(nameComp) != std::string::npos){
			retNames.push_back(fname);
		}
	}
	if(!retNames.empty()){
		digestFile = retNames[0];
		cLog<<"Digest File generated :"<<digestFile<<std::endl;
		return 1;
	}
	else 
		return 0;
	
}
