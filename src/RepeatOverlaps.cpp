
#include "RepeatOverlaps.h"
#include <fstream>

int Repeats::ReadRepeatIntervals(std::string repeatfile, OutStream& log, int overlapExtent){

    repeat_features feats;
    
    std::string chr_of_vector, chr, temp1, temp2, temp3;
    int start, end, index = 0;
    std::string fname;
    std::vector< Interval< repeat_features >  > repeat_intervals;
    
    fname.append(repeatfile);
    std::ifstream infile(fname.c_str());
    
    log << "Repeat File is " << fname << std::endl;
    
    
    while(infile >> chr){
        infile >> start >> end >> temp1 >> temp2 >> temp3;
        chr_of_vector = chr;
        ++index;
        bool endoffile = false;
        
        while(chr == chr_of_vector && !(endoffile)){
            feats.a = (temp1); feats.b = (temp2); feats.c = (temp3);
            repeat_intervals.push_back(Interval< repeat_features >(start, end, feats));
            
            if(infile >> chr)
                infile >> start >> end >> temp1 >> temp2 >> temp3;
            else
                endoffile = true;
            ++index;
        }
        
        if((chr_of_vector.length() < overlapExtent)){ 
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
