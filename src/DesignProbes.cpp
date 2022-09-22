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
//  DesignProbes.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "DesignProbes.h"
#include <fstream>
#include <cstdio>
#include <algorithm>

void DesignClass::InitialiseDesign(ProbeFeatureClass& proms, std::vector<PrDes::REPosMap>& fragmap ){
     
    for(int k = 0;k < proms.ChrNames_proms.size(); ++k){
        chrIndex[proms.ChrNames_proms[k]] = k;
        fragmap.push_back(PrDes::REPosMap());
    }
}


double DesignClass::BigWigSummary(std::string chr, int start, int end){

    std::string cmd;
    std::string s, e;
    
    cmd.append(bigwigsummarybinary);
    cmd.append(" ");
    cmd.append(mappabilityfile);
    cmd.append(" ");
    cmd.append(chr);
    cmd.append(" ");
    
    s = std::to_string(start);
    e = std::to_string(end);
    cmd.append(s);
    cmd.append(" ");
    cmd.append(e);
    cmd.append(" 1");
    
    char *cmdchar = &cmd[0];
    
    char buf[BUFSIZE];
    double a = -1;
    FILE *fp;
    
    if ((fp = popen(cmdchar, "r")) == NULL) {
        dLog<<"BigWigSummary: Error opening pipe!"<<std::endl;
        return -1;
    }
    
    while (fgets(buf, BUFSIZE, fp) != NULL) {        
        a = atof(buf);
    }
    
    if(pclose(fp))  {
        dLog<<"BigWigSummary: Command not found or exited with error status "<<std::endl;
        return -1;
    }
    return a;
    
}

bool DesignClass::CheckFragmentSize(RESitesClass &dpnII, std::string chr, int closest_re, int whichside){
	
	bool refound;
	int *resites;
    resites = new int[2];
    int invalidRECoordinate=0;
	
	refound = dpnII.GettheREPositions(chr,closest_re, resites, invalidRECoordinate);
	
	if(abs(closest_re - resites[1-whichside])> ProbeLen){ //
		return true;
	}
	else 
		return false;
	
}

bool DesignClass::CheckRepeats(Repeats& repeat_trees, std::string chr, int probe_start, int probe_end, bool ifRep){
	
	if(!ifRep)
		return true;
	
	int overlaps = repeat_trees.FindOverlaps(chr, probe_start ,probe_end);
	
	if(overlaps < repOverlapExtent)
		return true;
	else
		return false;
}

bool DesignClass::CheckMappability(std::string chr, int probe_start, int probe_end, bool ifMap){
	
	if(!ifMap)
		return true;
	
	double mappability = BigWigSummary(chr, probe_start, probe_end);
	
	if(mappability > mapThreshold)
		return true;
	else
		return false;
}
bool DesignClass::CheckCountofNs(bioioMod& getSeq, std::string chr, int probe_start, int probe_end, bool ifCountNs){
	if(!ifCountNs)
		return true;

	std::string toFas = getSeq.GetFasta(chr+":"+std::to_string(probe_start)+"-"+std::to_string(probe_end));
	int cnt = std::count(toFas.begin(), toFas.end(), 'N');
	if(cnt>0){
		return false;
	}
	else{
		return true;
	}
}
bool DesignClass::CheckRESite(RESitesClass & dpnII, bioioMod& getSeq, std::string chr, int tss, int& closest_re, bool whichside, Repeats& repeat_trees, bool ifRep, bool ifMap, bool ifCountNs){
    
    if(!ifRep && !ifMap && !ifCountNs){
		return true;
	}
	int probeStart, probeEnd;
	int *resites;
    resites = new int[2];
    bool refound = false, passed = false;
    int invalidRECoordinate=0;
	
    
    if (whichside) {
		probeStart = (closest_re - ProbeLen + reRightCut -1 );
		probeEnd = ( closest_re + reRightCut - 1 );
	}
	else{
		probeStart = (closest_re - reLeftCut - 1 );
		probeEnd = (closest_re + ProbeLen - reLeftCut - 1 );
	}
    
    bool checkRep = CheckRepeats(repeat_trees, chr, probeStart, probeEnd, ifRep);
    bool checkMap = CheckMappability(chr, probeStart, probeEnd, ifMap);
	bool checkNs = CheckCountofNs(getSeq, chr,probeStart,probeEnd,ifCountNs); // For genomes no mappability or repeat files available
    bool checkFrag = CheckFragmentSize(dpnII, chr, closest_re, whichside);
       
    if(checkRep && checkMap && checkNs && checkFrag && abs(closest_re - tss) <= MaxDistancetoTSS){
		passed = true;
	}
	else{
		while(abs(closest_re - tss) <= MaxDistancetoTSS){
			
			refound = dpnII.GettheREPositions(chr,closest_re, resites, invalidRECoordinate);
			if(refound){
				closest_re = resites[whichside];
			}
			else
				break;
			if(abs(closest_re - tss) <= MaxDistancetoTSS){	
				if(CheckFragmentSize(dpnII, chr, closest_re, whichside)){ //Frag size passes
					if (whichside == 1) {
						probeStart = closest_re - ProbeLen + reRightCut -1 ;
						probeEnd = closest_re + reRightCut - 1 ;
					}
					else{
						probeStart = closest_re - reLeftCut - 1 ;
						probeEnd = closest_re + ProbeLen - reLeftCut - 1 ; 
					}				
					checkRep = CheckRepeats(repeat_trees, chr, probeStart, probeEnd, ifRep);
					checkMap = CheckMappability(chr, probeStart, probeEnd, ifMap);
					checkNs = CheckCountofNs(getSeq, chr,probeStart,probeEnd,ifCountNs);
					if(checkRep && checkMap && checkNs){
						passed = true;
						break;
					}
				}
			}			
		}
	}
   return passed;
}


int DesignClass::CheckDistanceofProbetoTSS(RESitesClass& dpnII, std::string chr, int tss, int closest_re, int whichside){
   //return closest re site to TSS 
    int *resites;
    resites = new int[2];
    bool refound = 0;
    int invalidRECoordinate=0;

    while (abs(tss - closest_re) < minDistToTSS) {
		refound = dpnII.GettheREPositions(chr, closest_re, resites, invalidRECoordinate);
		if (refound){
			closest_re = resites[whichside];
		}
		else{	
			return closest_re;
		}
	}
    return closest_re;
}


bool DesignClass::createNewEntry(std::unordered_map<int, PrDes::REposStruct >& thismap, std::unordered_map<int, PrDes::REposStruct >& maptosearch, int chrind, std::string promind, int REpos, bool whichside){

    if ( REpos == 0 )
        return 0;
    
    if (thismap.find(REpos) == maptosearch.end()) {
        thismap[REpos].prom_indexes.push_back(promind);
        thismap[REpos].processed = 0;
        thismap[REpos].whichside = whichside; //This is to handle cases where a downstream probe is moved to upstream design.
        thismap[REpos].turnedoff = false; 
        
    }
    else{
		if(whichside && !thismap[REpos].whichside){
			auto it = thismap[REpos].prom_indexes.begin();
			thismap[REpos].prom_indexes.insert(it, promind);	
			thismap[REpos].whichside = whichside;
		}
		else
			thismap[REpos].prom_indexes.push_back(promind);
    }
    return 1;
}


void DesignClass::DesignProbes(ProbeFeatureClass & Feats, RESitesClass & dpnII, Repeats& repeat_trees, bioioMod& getSeq, std::string bgwgsumbin, std::string mappabilityfilepath, std::string whichchr, int mDisttoTSS, int prlen, PrDes::RENFileInfo& reInfo, int bufSize, int dFromTSS, int minDistBetProbes){
    
    // There is only one design layer which contains both upstream and downstream designs. Each feature will have two probes if possible.
    ProbeLen = prlen;
    MaxDistancetoTSS=mDisttoTSS;
    bigwigsummarybinary=bgwgsumbin;
    mappabilityfile=mappabilityfilepath;
	reLeftCut = reInfo.leftOfCut;
	reRightCut = reInfo.rightOfCut;
    mapThreshold = reInfo.mappabilityThreshold;
    repOverlapExtent = reInfo.repeatOverlapExtent;
    BUFSIZE=bufSize;
    minDistToTSS = dFromTSS;
    
    OneDesign.push_back(PrDes::DesignLayers());
    InitialiseDesign(Feats, OneDesign.back().Layer);
    
    dLog << "Design initialised in DesignProbes" << std::endl;
    
    std::unordered_map<std::string, int>::iterator it;
    int left_res, right_res; 
    std::string fn;
    
    fn.append(reInfo.desName);
    fn.append(".");
    fn.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
    fn.append(".ProbeDesignSummary_");
    fn.append(whichchr);
    fn.append(".");
    fn.append(reInfo.REName);
    fn.append(".");
    fn.append(reInfo.currTime);
    fn.append(".txt");
    std::ofstream summaryfile(fn.c_str());
    
    summaryfile<<"Feature_Name"<<'\t'<<"Feature_chr"<<'\t'<<"Feature_coordinate"<<'\t'<<"Dist between Left ProbeStart and Closest RESite to fragment"<<'\t'<<"Dist between Right ProbeEnd and Closest RESite to fragment"<<'\t'<<"If Upstream(left)/Downstream(right) probes created"<<std::endl;
    
    std::string fname = reInfo.desName;
    fname.append(".");
    fname.append(reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
    fname.append("_");
    fname.append(whichchr);
    fname.append(".");
    fname.append(reInfo.REName);
    fname.append(".");
    fname.append(reInfo.currTime);
    fname.append(".gff3");
    std::ofstream outfile;
    outfile.open(fname, std::fstream::out);

	outfile<<"##gff-version 3.2.1"<<std::endl;
	outfile<<"##genome-build "<<reInfo.genomeAssembly.substr(reInfo.genomeAssembly.find_first_of(',')+1)<<" "<<reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')) <<std::endl;
	////////add genome and sequence region if possible
		
    for (auto pIt = Feats.promFeatures.begin(); pIt != Feats.promFeatures.end(); ++pIt) {
        if(pIt->second.chr == whichchr && !(pIt->second.probesSkip)){ //if chosen chromosome and no skip of the transcript
            it = chrIndex.find(pIt->second.chr);
            bool passed_upstream = false, passed_downstream = false;
            
            left_res = CheckDistanceofProbetoTSS(dpnII, pIt->second.chr, pIt->second.TSS, pIt->second.closestREsitenums[0], 0);
            right_res = CheckDistanceofProbetoTSS(dpnII, pIt->second.chr, pIt->second.TSS, pIt->second.closestREsitenums[1], 1);
            //last parameter is design direction not the promoter (RE on the left or right side of TSS)

            passed_upstream = CheckRESite(dpnII, getSeq, pIt->second.chr, pIt->second.TSS, left_res, 0, repeat_trees, reInfo.ifRepeatAvail, reInfo.ifMapAvail,reInfo.ifCountNs);
            passed_downstream = CheckRESite(dpnII, getSeq, pIt->second.chr, pIt->second.TSS, right_res, 1, repeat_trees, reInfo.ifRepeatAvail, reInfo.ifMapAvail,reInfo.ifCountNs);
            
            if (passed_upstream && passed_downstream) { //check if fragment is same and longer than ProbeLen
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, left_res,0); 
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, right_res,1); 
            }
            if((passed_upstream) && (!passed_downstream)){
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, (left_res),0); 
            }
            if((!passed_upstream) && passed_downstream){
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, (right_res),1); 
            }
        }
    }
//filter out probes that are too close to each other   
    int rightlimit= minDistBetProbes+ProbeLen-reRightCut+1;
    for (it = chrIndex.begin(); it != chrIndex.end(); ++it) {
        for (auto itr = OneDesign.back().Layer[it->second].repmap.begin(); itr != OneDesign.back().Layer[it->second].repmap.end(); ++itr){
			if(!itr->second.turnedoff){ //if RE site has active probe
				for(auto its = std::next(itr , 1); its != OneDesign.back().Layer[it->second].repmap.end(); ++its){ 
					if(itr->second.prom_indexes[0]!=its->second.prom_indexes[0]){ //if the probes are not from same feature
						if(abs(itr->first-its->first) < rightlimit){ //if the probes are closer than the threshold
							if(itr->second.whichside && its->second.whichside){ //if both probes are downstream, deactivate the one which is upstream 
								if(itr->first<its->first){
									itr->second.turnedoff=true;
								}
								else{
									its->second.turnedoff=true;
								}
							}
							else{ //in other cases deactivate the probe downstream of the other
								if(itr->first<its->first){
									its->second.turnedoff=true;
								}
								else{
									itr->second.turnedoff=true;
								}
							}
						}
					}
				}
			}
		}
	}
//Write to Probe File
    for (it = chrIndex.begin(); it != chrIndex.end(); ++it) {
        for (auto itr = OneDesign.back().Layer[it->second].repmap.begin(); itr != OneDesign.back().Layer[it->second].repmap.end(); ++itr){
			if(!itr->second.turnedoff){
				WritetoFile(outfile, it->first, it->second, itr->first, itr->second.prom_indexes, itr->second.whichside, reInfo.desName, Feats);
				Feats.promFeatures[itr->second.prom_indexes[0]].probeREs[itr->second.whichside] = abs(Feats.promFeatures[itr->second.prom_indexes[0]].closestREsitenums[itr->second.whichside] - itr->first);
			}
		}
    }
 
//Write to summary file    
    for (auto pIt = Feats.promFeatures.begin(); pIt != Feats.promFeatures.end(); ++pIt) {
		if(pIt->second.chr == whichchr){
			std::string probesStatus;
		
			if(pIt->second.probeREs[0] > -1 && pIt->second.probeREs[1] > -1)
				probesStatus="both";
			else if(pIt->second.probeREs[0] > -1 && pIt->second.probeREs[1] == -1)
				probesStatus="only upstream";
			else if(pIt->second.probeREs[0] == -1 && pIt->second.probeREs[1] > -1)
				probesStatus="only downstream";
			else 
				probesStatus="none";	
			
			summaryfile << pIt->second.genes[0] << '\t' << pIt->second.chr  << '\t' << pIt->second.TSS << '\t' << pIt->second.probeREs[0]  << '\t' << pIt->second.probeREs[1] << '\t' << probesStatus << std::endl;
		}
	}
    
    
    
    dLog << "Design for "<< whichchr<<" written to Output files " << std::endl;

    
    outfile.close();
}



int DesignClass::CheckREsiteAroundProbe(RESitesClass& dpnII, std::string chr, int probe_re, int direction){
    // It checks if for neighbouring restriction enzyme positions any probes designed. if yes, it will move it to a new design layer.
    int *resites;
    resites = new int[2];
    bool refound = 0;
    int invalidRECoordinate=0;
    
    refound = dpnII.GettheREPositions(chr,probe_re,resites, invalidRECoordinate);
    
    if (refound) {
        return resites[direction];
    }
    else
        return 0;

}


bool DesignClass::WritetoFile(std::ofstream &outfile, std::string chr, int chrind, int repos, std::vector< std::string > values, bool direction, std::string design, ProbeFeatureClass& feats){
	
	std::string target, side;
	int disttotss, probestart, probeend;
	switch(feats.promFeatures[values[0]].FeatureType){
		case 1:
			target="promoter";
			break;
		case 2:
			target="SNV";
			break;
		case 3:
			target="neg_ctrl";
			break;
		default:
			target ="other";
	}
    
    if ((repos == 0 || repos == 120 || repos == -120 ))
        return 0;
    
    if (direction) {// DOWNSTREAM - RIGHT
		probestart=( repos - ProbeLen + reRightCut); //add
		probeend=( repos + reRightCut); //1 based
		
		side="R";
		disttotss=abs(probeend - feats.promFeatures[values[0]].TSS); 
	}
	
	else if(!direction){ // UPSTREAM - LEFT
		probestart = ( (repos+1) );
		probeend =( (repos+1) + ProbeLen ); //+1 to select the right fragment
		side = "L";
		disttotss = abs(probestart - feats.promFeatures[values[0]].TSS);  
	}   

    outfile  << chr  << '\t' << "." << '\t' <<"probe"<<'\t';
    
    outfile <<probestart + 1 << '\t' << probeend << '\t'<< "." << '\t'<< "." << '\t'<< "." << '\t'; // to adjust for 1-based coords
        
    outfile << "Name="<<feats.promFeatures[values[0]].genes[0]<<"; " <<"transcriptid="<< feats.promFeatures[values[0]].transcripts[0]<<"; "<< "side="<<side<<"; "<<"target="<<target<<"; "<<"design="<<design<<"; "<< "featuresinvicinity=";
        
    //check for multiple promoters
    if(values.size()>1){
		for(size_t j = 1; j < values.size(); ++j){
			outfile << feats.promFeatures[values[j]].genes[0]<<":"<<feats.promFeatures[values[j]].TSS ;
			if(j<(values.size()-1))
				outfile <<", ";
		}
		
		outfile <<"; ";
	}       
    else{
		outfile << "none"<<"; ";
	}
	
	outfile<< "targettss="<<feats.promFeatures[values[0]].TSS<<"; "<<"distancetotss="<<disttotss<< std::endl;
	
	int offset=0;
	if(side=="L")
		offset=1; // subtract 1 for bioio coordinates on Left side
	BaseProbe t;
	t.chr=chr;
	t.start=probestart-offset;
	t.end=probeend-offset;
	t.strand="+";
	t.feature=feats.promFeatures[values[0]].genes[0];
    
    probeList.push_back(t);
    
    return true;
}

void DesignClass::MergeAllChrOutputs(ProbeFeatureClass& Feats, PrDes::RENFileInfo& reInfo){
	
	std::string fnameAllChr = reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".AllProbes."+reInfo.REName+"."+reInfo.currTime+".gff3";    
	std::string sumFnameAllChr = reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".ProbeDesignSummary_chrAll."+reInfo.REName+"."+reInfo.currTime+".txt";
    std::string header, sumHeader;
    
    std::ofstream outfile, sumOutfile;
    outfile.open(fnameAllChr, std::fstream::out);
    sumOutfile.open(sumFnameAllChr, std::fstream::out);

	header.append("##gff-version 3.2.1");
	header.append("\n");
	header.append("##genome-build "+reInfo.genomeAssembly.substr(reInfo.genomeAssembly.find_first_of(',')+1)+" "+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')));
	header.append("\n");
    
    outfile<<header;
    
    sumHeader = std::string("Feature_Name")+"\t"+"Feature_chr"+"\t"+"Feature_coordinate"+"\t"+"Dist between Left ProbeStart and Closest RESite to fragment"+"\t"+"Dist between Right ProbeEnd and Closest RESite to fragment"+"\t"+"If Upstream(left)/Downstream(right) probes created"+"\n";
    sumOutfile<< sumHeader;

	outfile.close();
	sumOutfile.close();
	
	outfile.open(fnameAllChr, std::ios_base::binary | std::ios_base::app);
	sumOutfile.open(sumFnameAllChr, std::ios_base::binary | std::ios_base::app);
	
	for(auto &iChr : Feats.ChrNames_proms){
		
		std::string fName=reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+"_"+iChr+"."+reInfo.REName+"."+reInfo.currTime+".gff3";
		std::string sumfName=reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".ProbeDesignSummary_"+iChr+"."+reInfo.REName+"."+reInfo.currTime+".txt";
		std::ifstream readChrFiles(fName, std::ios_base::binary);
		std::ifstream readSumFiles(sumfName, std::ios_base::binary);
		readChrFiles.seekg(header.size()); //strip header
		readSumFiles.seekg(sumHeader.size()); //strip header
		outfile << readChrFiles.rdbuf();
		sumOutfile << readSumFiles.rdbuf();
		readChrFiles.close();
		readSumFiles.close();
		
		remove(fName.c_str());
		remove(sumfName.c_str());
		
	}	
}


bool DesignClass::ConstructSeq(PrDes::RENFileInfo& reInfo, bioioMod& getSeq, std::string calledChr){
	
	std::string fname = reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".ProbeSequences."+calledChr+"."+reInfo.REName+"."+reInfo.currTime+".txt";    
    std::string header;
    
	std::string probesBedFile= reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".AllProbeSequences."+reInfo.REName+"."+reInfo.currTime+".bed";
	    
    
    std::ofstream outfile, outfile2;
    outfile.open(fname, std::fstream::app);
    outfile2.open(probesBedFile, std::fstream::app);

	outfile2<<"track name=\"AllProbes_"<<reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))<<"_"<<reInfo.REName<<"_"<<reInfo.currTime<<"\""<<std::endl;
	outfile<<"TargetID"<<'\t'<<"ProbeID"<<'\t'<<"Sequence"<<'\t'<<"Replication"<<'\t'<<"Strand"<<'\t'<<"Coordinates"<<std::endl;	
	
	for(auto &probeVar : probeList){
		
		outfile2<<probeVar.chr<<'\t'<<probeVar.start <<'\t'<<probeVar.end<<'\t'<<probeVar.feature<<std::endl;
		
		std::string toFas = getSeq.GetFasta(probeVar.chr+":"+std::to_string(probeVar.start)+"-"+std::to_string(probeVar.end));
		
		if(toFas!="Error"){
			outfile<<probeVar.feature<<'\t'<< probeVar.chr<<"_"<<probeVar.start+1<<'\t'<<toFas<<'\t'<<"1"<<'\t'<<probeVar.strand<<'\t'<<probeVar.chr<<":"<<probeVar.start+1<<"-"<<probeVar.end<<std::endl;
		}
		else
			return false;	
	}
	
	return true;	
}

