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
//  RepeatOverlaps.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//
#include "RepeatOverlaps.h"
#include <fstream>

int Repeats::ReadRepeatIntervals(std::string repeatfile, OutStream& log){

    repeat_features feats;
    
    std::string chr_of_vector, chr, temp1, temp2, temp3, tmp;
    int start, end, index = 0;
    std::string fname;
    std::vector< Interval< repeat_features >  > repeat_intervals;
    
    fname.append(repeatfile);
    std::ifstream infile(fname.c_str());
    
    log << "Repeat File is " << fname << std::endl;
    
    
    while(infile >> tmp>> temp2 >> tmp>> tmp>>tmp>> chr>>start >> end >>tmp>> temp3>> temp1>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp){

        chr_of_vector = chr;
        ++index;
        bool endoffile = false;
        
        while(chr == chr_of_vector && !(endoffile)){
            feats.a = (temp1); feats.b = (temp2); feats.c = (temp3);
            repeat_intervals.push_back(Interval< repeat_features >(start, end, feats));
            
            if(infile >> tmp>> temp2 >> tmp>> tmp>>tmp>> chr>>start >> end >>tmp>> temp3>> temp1>>tmp>>tmp>>tmp>>tmp>>tmp>>tmp)
                endoffile = false;
            else
                endoffile = true;
            ++index;
        }
        
        if((chr_of_vector.length() < 6)){  //filter out haplotype assemblies/random chromosomes
            AddTotheIntervalTree(repeat_intervals, chr_of_vector);
        }
        chr_of_vector = chr;
    }

    return 0;
}


int Repeats::AddTotheIntervalTree(std::vector<Interval<repeat_features> >& intervals, std::string chr_of_vector){
    
    IntervalTree< repeat_features > tree;

    std::vector<Interval< repeat_features >  > temp;
    
    tree = IntervalTree< repeat_features >(intervals);
    repeat_map[chr_of_vector] = tree;
    intervals.swap(temp);
    
    return 1;
}



int Repeats::FindOverlaps(std::string chr, unsigned long int probestart, unsigned long int probeend){
    
    std::vector<Interval<repeat_features> > results;
    
    repeat_map[chr].findOverlapping(probestart, probeend, results);
    
    return results.size();
}
