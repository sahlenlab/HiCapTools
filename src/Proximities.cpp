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
//  Proximities.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//
#include "Proximities.h"

void ProximityClass::AnnotateDistalInteractor(std::string feature_id, std::string anchored_chr, std::string interactor_chr, int *interactor_resite, int sc_index, int ExperimentNo){
    
    n = interactor_resite[0];
    // Intra chromosomal interaction
    if(anchored_chr.compare(interactor_chr) == 0){
        PopulateInteractions(Features[feature_id].proximities.junctions, interactor_resite, sc_index, ExperimentNo);
    }
    else{  // inter chromosomal interaction
        chrfound = 0;
        for(auto itx = Features[feature_id].proximities_ctx.begin() ; itx < Features[feature_id].proximities_ctx.end();++itx){
            if (interactor_chr.compare(itx->maptochrname) == 0){
                if(itx->junctions_ctx.find(n) == itx->junctions_ctx.end()){
                    itx->junctions_ctx[n].paircount = new int[NOFEXPERIMENTS]; // add a new entry
                    itx->junctions_ctx[n].strandcombination = new int[((NOFEXPERIMENTS)*4)];
                    for(int z = 0; z < (NOFEXPERIMENTS); ++z)
                        itx->junctions_ctx[n].paircount[z] = 0;
                    for (int z = 0; z < ((NOFEXPERIMENTS)*4); ++z) {
                        itx->junctions_ctx[n].strandcombination[z] = 0;
                    }
                    itx->junctions_ctx[n].refragend = interactor_resite[1];
                    itx->junctions_ctx[n].paircount[ExperimentNo] = 1;
                    itx->junctions_ctx[n].strandcombination[sc_index] = 1;
                    
                 }
                else{ // if inserted before
                    itx->junctions_ctx[n].paircount[ExperimentNo] += 1;
                    itx->junctions_ctx[n].strandcombination[sc_index] += 1;
                }
                
                chrfound = 1;
                break;
            }
        }
        if(!chrfound){
            Features[feature_id].proximities_ctx.push_back(SignalStruct_CTX());
            Features[feature_id].proximities_ctx.back().maptochrname.append(interactor_chr);
            Features[feature_id].proximities_ctx.back().junctions_ctx[n].paircount = new int[NOFEXPERIMENTS];
            Features[feature_id].proximities_ctx.back().junctions_ctx[n].strandcombination = new int[NOFEXPERIMENTS*4];
            for(int z = 0; z < (NOFEXPERIMENTS); ++z)
                Features[feature_id].proximities_ctx.back().junctions_ctx[n].paircount[z] = 0;
            Features[feature_id].proximities_ctx.back().junctions_ctx[n].strandcombination = new int[((NOFEXPERIMENTS)*4)];
            for (int z = 0; z < ((NOFEXPERIMENTS)*4); ++z) {
                Features[feature_id].proximities_ctx.back().junctions_ctx[n].strandcombination[z] = 0;
            }
            Features[feature_id].proximities_ctx.back().junctions_ctx[n].refragend = interactor_resite[1];
            Features[feature_id].proximities_ctx.back().junctions_ctx[n].paircount[ExperimentNo] = 1;
            Features[feature_id].proximities_ctx.back().junctions_ctx[n].strandcombination[sc_index] += 1;
           
        }
    }
}


void ProximityClass::AnnotateFeatFeatInteraction(std::string feature_id1, std::string feature_id2, int sc_index, int ExperimentNo){
    foundbefore = 0;
    for (auto it = Features[feature_id1].Inter_feature_ints.begin(); it < Features[feature_id1].Inter_feature_ints.end(); ++it){ // Check if the interaction with that promoter is seen before
        if (it->interacting_feature_id == feature_id2){
            it->signal[ExperimentNo] += 1.0;
            it->strandcombination[sc_index] += 1;
            foundbefore = 1;
            break;
        }
    }
    if(!foundbefore){ // Create a new entry for that promoter
        Features[feature_id1].Inter_feature_ints.push_back(FeattoFeatSignalStruct());
        Features[feature_id1].Inter_feature_ints.back().interacting_feature_id = feature_id2;
        Features[feature_id1].Inter_feature_ints.back().signal.resize(NOFEXPERIMENTS);
        Features[feature_id1].Inter_feature_ints.back().signal[ExperimentNo] = 1;
        Features[feature_id1].Inter_feature_ints.back().strandcombination = new int[((NOFEXPERIMENTS)*4)];
        for (int z = 0; z < ((NOFEXPERIMENTS)*4); ++z)
            Features[feature_id1].Inter_feature_ints.back().strandcombination[z] = 0;
        
        Features[feature_id1].Inter_feature_ints.back().strandcombination[sc_index] = 1;
    }
}


void ProximityClass::PopulateInteractions(std::map<int, Junction >& signals, int *interactor_resite, int sc_index, int ExperimentNo){
int n = interactor_resite[0];
    
    if(signals.find(n) == signals.end()){
        signals[n].paircount = new int[NOFEXPERIMENTS];
        for(int z = 0; z < (NOFEXPERIMENTS); ++z)
            signals[n].paircount[z] = 0; // Initialise
        
        signals[n].strandcombination = new int[4*(NOFEXPERIMENTS)];
        for (int z = 0; z < 4*(NOFEXPERIMENTS); ++z)
            signals[n].strandcombination[z] = 0;
        
        signals[n].paircount[ExperimentNo] = 1;
        signals[n].strandcombination[sc_index] = 1;
        signals[n].refragend = interactor_resite[1];
        
        
    }
    else{
        signals[n].paircount[ExperimentNo] += 1;
        signals[n].strandcombination[sc_index] += 1;
    }
    
}

void ProximityClass::RecordProximities(Alignment pair, std::string feature_id1, std::string feature_id2, int sc_index, int ExperimentNo) {
 
    if ((feature_id1.length() != 4 && feature_id2.length() != 4)) { //"null"
        AnnotateFeatFeatInteraction(feature_id1, feature_id2, sc_index, ExperimentNo);
        AnnotateFeatFeatInteraction(feature_id2, feature_id1, sc_index, ExperimentNo);
        
    }
    else{
        
        if (feature_id1.length() != 4 ) // If first read is annotated with a probe, second read is the interactor
            AnnotateDistalInteractor(feature_id1, pair.chr1, pair.chr2, pair.resites2, sc_index, ExperimentNo);
        else
            AnnotateDistalInteractor(feature_id2, pair.chr2, pair.chr1, pair.resites1, sc_index, ExperimentNo);
    }
}

