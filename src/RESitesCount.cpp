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
#include <algorithm>
#include <fstream>
#include <iostream>

void RESitesClass::InitialiseVars(std::string DigestedGenomeFileName){

    posvector.clear();
    chr_names.clear();
    indexes.clear();
    chroffsets_indexfile.clear();
    chr_starts.clear();
    chr_ends.clear();
    chr_ranges.clear();

    std::string s;
    s.append(DigestedGenomeFileName);
    std::ifstream RESitesf(s.c_str());

    rLog << "Digest File is " << s << std::endl;

    s.clear();

    std::string chrp, temp;
    int pos;
    span = 1000000; // Window size

    //For indexing
    rLog << "Initialising RE site Class..." << std::endl;

    std::getline(RESitesf,temp); //get the header row1
    std::getline(RESitesf,temp); //get the header row

    std::string currentChr;
    std::vector<int> positions;

    auto finalizeChromosome = [&](const std::string& chr, const std::vector<int>& chrPositions){
        if(chr.empty() || chrPositions.empty()){
            return;
        }

        chr_names.push_back(chr);
        chroffsets_indexfile[chr] = static_cast<int>(indexes.size());
        chr_starts[chr] = chrPositions.front();
        chr_ends[chr] = chrPositions.back();

        const size_t startOffset = posvector.size();
        posvector.insert(posvector.end(), chrPositions.begin(), chrPositions.end());
        chr_ranges[chr] = std::make_pair(startOffset, chrPositions.size());

        indexes.push_back(PrDes::REindexes());
        auto& idx = indexes.back();

        size_t binStartIndex = 0;
        while(binStartIndex < chrPositions.size()){
            const int binStartCoord = chrPositions[binStartIndex];
            const int binLimit = binStartCoord + span;

            size_t binEndIndex = binStartIndex;
            while(binEndIndex < chrPositions.size() && chrPositions[binEndIndex] <= binLimit){
                ++binEndIndex;
            }

            idx.binstart.push_back(binStartCoord);
            idx.binend.push_back(chrPositions[binEndIndex - 1] + 1);
            idx.offset.push_back(static_cast<int>(startOffset + binStartIndex));
            idx.count.push_back(static_cast<int>(binEndIndex - binStartIndex));

            binStartIndex = binEndIndex;
        }
    };

    while(RESitesf >> chrp >> temp >> pos >> temp >> temp >> temp >> temp){
        if(currentChr.empty()){
            currentChr = chrp;
        }

        if(chrp != currentChr){
            finalizeChromosome(currentChr, positions);
            positions.clear();
            currentChr = chrp;
        }

        positions.push_back(pos);
    }

    finalizeChromosome(currentChr, positions);

    rLog << "RE site class initialised " << std::endl;

}

bool RESitesClass::GettheREPositions(std::string chr, int pos, int* renums, int& invalidCounter){ // Returns closest RE sites to a position

    std::unordered_map<std::string, std::pair<size_t, size_t>>::iterator rangeIt = chr_ranges.find(chr);

    if(rangeIt == chr_ranges.end()){ // Chromosome is not in the list
        return 0;
    }

    std::unordered_map< std::string, int >::iterator its = chr_starts.find(chr);
    std::unordered_map< std::string, int >::iterator ite = chr_ends.find(chr);

    if (pos < its->second || pos > ite->second){
                //rLog<<"!!Error!! : Encountered invalid coordinates. A coordinate is out of chromosome boundaries and is therefore skipped: chr "<< chr<<" Position "<< pos <<". Check if the correct genome assembly is being used"<< std::endl;
                invalidCounter=invalidCounter+1;
        return 0;
    }

    const size_t startOffset = rangeIt->second.first;
    const size_t siteCount = rangeIt->second.second;

    const auto startIter = posvector.begin() + startOffset;
    const auto endIter = startIter + siteCount;

    auto lower = std::lower_bound(startIter, endIter, pos);

    int upstreamSite;
    int downstreamSite;

    if(lower == endIter){
        upstreamSite = downstreamSite = *(endIter - 1);
    }
    else if(*lower == pos){
        upstreamSite = (lower == startIter) ? *lower : *(lower - 1);
        downstreamSite = (std::next(lower) == endIter) ? *lower : *std::next(lower);
    }
    else{ // *lower > pos
        downstreamSite = *lower;
        upstreamSite = (lower == startIter) ? *lower : *(lower - 1);
    }

    renums[0] = upstreamSite;
    renums[1] = downstreamSite;

    return 1;
}


void RESitesClass::CleanClass(){

        posvector.clear();
        chr_names.clear();
        indexes.clear();
        chroffsets_indexfile.clear();
        chr_starts.clear();
        chr_ends.clear();
        chr_ranges.clear();
}
