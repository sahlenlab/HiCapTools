//
//  RepeatOverlaps.h
//  HiCap_ProbeDesigner
//
//  Created by Pelin Sahlen on 12/02/2016.
//  Copyright © 2016 Pelin Sahlen. All rights reserved.
//

#ifndef PRSUITE_INC_REPOVER_H_
#define PRSUITE_INC_REPOVER_H_

#include "IntervalTree.h"
#include "ProbeDataStructs.h"
#include "OutStream.h"


struct repeat_features{
    std::string a;
    std::string b;
    std::string c;
};

class Repeats{
    friend class Probes;
public:
    std::unordered_map< std::string, IntervalTree < repeat_features > > repeat_map;
    int ReadRepeatIntervals(std::string, OutStream&, int);
    int FindOverlaps(std::string, unsigned long int, unsigned long int);
    
private:
    int AddTotheIntervalTree(std::vector <Interval < repeat_features >  >&, std::string);
};

#endif // PRSUITE_INC_REPOVER_H_
