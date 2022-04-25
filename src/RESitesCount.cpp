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
//  RESitesCount.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "RESitesCount.h"
#include <fstream>
#include <iostream>

void RESitesClass::InitialiseVars(std::string DigestedGenomeFileName){

    std::string s;
    s.append(DigestedGenomeFileName);
	std::ifstream RESitesf(s.c_str());
    
    rLog << "Digest File is " << s << std::endl;
    
	s.clear();

	std::string chrname, chrp, temp;
	int pos, chrstart; 
	int startpos;
	int span = 1000000; // Window size
	
	//For indexing
	rLog << "Initialising RE site Class..." << std::endl;

	bool f=0;

	std::getline(RESitesf,temp); //get the header row1
	std::getline(RESitesf,temp); //get the header row
	
	
	RESitesf >> chrp >> temp >> pos >> temp >> temp >> temp >> temp; // Read the first line
	
	chrname = chrp; // get the chr name outside the loop
	chr_names.push_back(chrp);
	chrstart = (pos ); // chromosome start
	indexes.push_back(PrDes::REindexes());
	chroffsets_indexfile[chrp] = (indexes.size()-1);
	chr_starts[chrp] = chrstart;
	
	while(!(RESitesf.eof())){
		int binstart = (pos);
		int binend = binstart + span;
		startpos = posvector.size();
		while( chrname == chrp && (pos >= binstart && pos <= binend)){
			posvector.push_back(pos );
			RESitesf >> chrp >> temp >> pos >> temp >> temp >> temp >> temp;
			if(RESitesf.eof()){
				f = 1; //End of file
				break;
			}
		}
		
		if(indexes.back().binend.empty())
			indexes.back().binstart.push_back(chrstart);
		else
			indexes.back().binstart.push_back((indexes.back().binend.back())+1);
		
		indexes.back().binend.push_back((posvector.back() + 1));
		indexes.back().offset.push_back(startpos);
		indexes.back().count.push_back((posvector.size()-startpos));
		
		if((chrname != chrp) || f ){
			chr_ends[chrname] = indexes.back().binend.back();
			indexes.push_back(PrDes::REindexes());
			chroffsets_indexfile[chrp] = (indexes.size()-1);
			chrname = chrp;
			chrstart = (pos);
			chr_names.push_back(chrp);
			chr_starts[chrp] = chrstart;
		}
	}
	rLog << "RE site class initialised " << std::endl;

}

bool RESitesClass::GettheREPositions(std::string chr, int pos, int* renums, int& invalidCounter){ // Returns closest RE sites to a position
	
    int HalfClusterDist = 5000; //in case the pos is at the end of a bin
    int rightchr;
    int starttosearch = 0;
    int bitcount = 0;
    int REposprev, REposat, REposnext;

    std::unordered_map< std::string, int >::iterator it = chroffsets_indexfile.find(chr);
    
    if(it == chroffsets_indexfile.end()){ // Chromosome is not in the list
        return 0;
    }
    
    rightchr = it->second; // Get the right index vector
    std::unordered_map< std::string, int >::iterator its = chr_starts.find(chr);
    std::unordered_map< std::string, int >::iterator ite = chr_ends.find(chr);
    
    if ((pos ) <= its->second || (pos ) >= ite->second){
		//rLog<<"!!Error!! : Encountered invalid coordinates. A coordinate is out of chromosome boundaries and is therefore skipped: chr "<< chr<<" Position "<< pos <<". Check if the correct genome assembly is being used"<< std::endl;
		invalidCounter=invalidCounter+1;
        return 0;
    }
    
    for(int i = 0; i < indexes[rightchr].binstart.size();++i){ // Iterate over the bins
        if (indexes[rightchr].binstart[i] <= (pos) && indexes[rightchr].binend[i] >= (pos)){ // If start-HalfClusterDist is within a bin
            if(i > 0){
                starttosearch = indexes[rightchr].offset[i-1]; // mark the index in the posvector to start to search
                bitcount = indexes[rightchr].count[i-1]; // this many elements of the posvector is contained within that bin
                bitcount += indexes[rightchr].count[i];
            }
            else{
                starttosearch = indexes[rightchr].offset[i]; // mark the index in the posvector to start to search
                bitcount = indexes[rightchr].count[i]; // this many elements of the posvector is contained within that bin
            }
            if (((pos + HalfClusterDist) > indexes[rightchr].binend[i]) && (i+1 < indexes[rightchr].binstart.size())) // if end is included in the next bin
                bitcount += indexes[rightchr].count[i+1]; // mark the number of elements in the next bin
            break;
        }
    }
    for (int i = (starttosearch + 1); i < starttosearch + bitcount; i++){
        REposprev = posvector[i - 1];
        REposat   = posvector[i];
        REposnext = posvector[i + 1];
        while (REposat <  pos){
            REposprev = REposat;
            ++i;
            REposat = posvector[i];
            REposnext = posvector[i + 1];
        }
        if (REposat == pos) {
            renums[0] = REposprev;
            renums[1] = REposnext;
            break;
        }
        else{
            renums[0] = REposprev;
            renums[1] = REposat;
            break;
        }
    }
    return 1;
}


void RESitesClass::CleanClass(){

	posvector.clear();
	chr_names.clear();
	indexes.clear();
	chroffsets_indexfile.clear();
}
