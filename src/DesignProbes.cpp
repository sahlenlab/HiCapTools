#include "DesignProbes.h"
#include <fstream>
#include <cstdio>


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
        dLog<<"Error opening pipe!"<<std::endl;
        return -1;
    }
    
    while (fgets(buf, BUFSIZE, fp) != NULL) {        
        a = atof(buf);
    }
    
    if(pclose(fp))  {
        dLog<<"Command not found or exited with error status"<<std::endl;
        return -1;
    }
    return a;
    
}


bool DesignClass::overlap(RESitesClass& dpnII, Repeats repeat_trees, int& closest_re, int tss, int whichside, std::string chr, int probe_start, int probe_end, bool ifRep, bool ifMap){
    
    int *resites;
    resites = new int[2];
    bool refound = false, passed = false;
    
    int overlaprepeats;
    double mappability; 
    
    if(ifRep && !ifMap){
		overlaprepeats = repeat_trees.FindOverlaps(chr, probe_start ,probe_end);
		mappability=mapThreshold;
	}
	if(!ifRep && ifMap){
		overlaprepeats=0;
		mappability = BigWigSummary(chr, probe_start, probe_end);
	}
	if(ifRep && ifMap){
		
		overlaprepeats = repeat_trees.FindOverlaps(chr, probe_start ,probe_end);
		mappability = BigWigSummary(chr, probe_start, probe_end);
	}
    
    if(overlaprepeats > 0 || mappability < mapThreshold ){ // if overlap with repeats and low mappability, go to the next RE site
        
        refound = dpnII.GettheREPositions(chr,closest_re, resites);
        
        if (refound){
            closest_re = resites[whichside];
        }
        else
            return false;
            
        if (abs(closest_re - tss) <= MaxDistancetoTSS) {
            passed = CheckRepeatOverlaps(dpnII, chr, tss, closest_re, whichside, repeat_trees, ifRep, ifMap);
        }
        else{
            passed = false;
        }
    }
    else{
        if(abs(closest_re - tss) <= MaxDistancetoTSS)
            passed = true;
    }
    return passed;
}


bool DesignClass::CheckRepeatOverlaps(RESitesClass & dpnII, std::string chr, int tss, int& closest_re, bool rightside, Repeats repeat_trees, bool ifRep, bool ifMap){
    
    if(!ifRep && !ifMap){
		return true;
	}
	
    bool passed = false;
    
    if (rightside) {
        passed = overlap(dpnII, repeat_trees, closest_re, tss, rightside, chr, ( closest_re - ProbeLen + reRightCut -1 ), ( closest_re + reRightCut - 1 ), ifRep, ifMap); 
    }
    else{
        passed = overlap(dpnII, repeat_trees, closest_re, tss, rightside, chr, ( closest_re - reLeftCut - 1 ), (closest_re + ProbeLen - reLeftCut - 1 ), ifRep, ifMap); 
    }
    return passed;
}


int DesignClass::CheckDistanceofProbetoTSS(RESitesClass& dpnII, std::string chr, int tss, int closest_re, bool rightside){
    //If probe is less than 120 bases away from the TSS, it looks at the next RE site to design probe
    
    int *resites;
    resites = new int[2];
    bool refound = 0;

    while (abs(tss - closest_re) < ProbeLen) {
		refound = dpnII.GettheREPositions(chr,closest_re, resites);
		if (refound){
			closest_re = resites[rightside];
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
    }
    else{
        thismap[REpos].prom_indexes.push_back(promind);
    }
    return 1;
}


void DesignClass::DesignProbes(ProbeFeatureClass & Feats, RESitesClass & dpnII, Repeats repeat_trees, std::string bgwgsumbin, std::string mappabilityfilepath, std::string whichchr, int mDisttoTSS, int prlen, PrDes::RENFileInfo& reInfo, int bufSize){
    
    // There is only one design layer which contains both upstream and downstream designs. Each feature will have two probes if possible.
    ProbeLen = prlen;
    MaxDistancetoTSS=mDisttoTSS;
    bigwigsummarybinary=bgwgsumbin;
    mappabilityfile=mappabilityfilepath;
	reLeftCut = reInfo.leftOfCut;
	reRightCut = reInfo.rightOfCut;
    mapThreshold = reInfo.mappabilityThreshold;
    BUFSIZE=bufSize;
    
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
        if(pIt->second.chr == whichchr){
            it = chrIndex.find(pIt->second.chr);
            bool passed_upstream = false, passed_downstream = false;
            left_res = CheckDistanceofProbetoTSS(dpnII, pIt->second.chr, pIt->second.TSS, pIt->second.closestREsitenums[0], 0);
            right_res = CheckDistanceofProbetoTSS(dpnII, pIt->second.chr, pIt->second.TSS, pIt->second.closestREsitenums[1], 1);
            //last parameter is design direction not the promoter (RE on the left or right side of TSS)
             
            passed_upstream = CheckRepeatOverlaps(dpnII, pIt->second.chr, pIt->second.TSS, left_res, 0, repeat_trees, reInfo.ifRepeatAvail, reInfo.ifMapAvail);
            passed_downstream = CheckRepeatOverlaps(dpnII, pIt->second.chr, pIt->second.TSS, right_res, 1, repeat_trees, reInfo.ifRepeatAvail, reInfo.ifMapAvail);
            
            
            if (passed_upstream && passed_downstream) {
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, left_res,0); 
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, right_res,1); 
                summaryfile << pIt->second.genes[0] << '\t' << pIt->second.chr  << '\t' << pIt->second.TSS << '\t' << (pIt->second.closestREsitenums[0] - left_res) << '\t' << (right_res - pIt->second.closestREsitenums[1]) << '\t' << "both" << std::endl;
            }
            if((passed_upstream) && (!passed_downstream)){
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, (left_res),0); 
                summaryfile << pIt->second.genes[0] << '\t' << pIt->second.chr  << '\t' << pIt->second.TSS << '\t' << (pIt->second.closestREsitenums[0] - left_res) << '\t' << -1 << '\t' << "only_upstream" << std::endl;
            }
            if((!passed_upstream) && passed_downstream){
                createNewEntry(OneDesign.back().Layer[it->second].repmap, OneDesign.back().Layer[it->second].repmap, it->second, pIt->first, (right_res),1); 
                summaryfile << pIt->second.genes[0] << '\t' << pIt->second.chr  << '\t' << pIt->second.TSS << '\t' << -1 << '\t' << (right_res - pIt->second.closestREsitenums[1]) << '\t' << "only_downstream" << std::endl;
            }
            if (!passed_upstream && !passed_downstream) {
                summaryfile << pIt->second.genes[0]
                << '\t' << pIt->second.chr << '\t' << pIt->second.TSS << '\t' << -1 << '\t' << -1 << '\t' << "none" << std::endl;
                
            }
        }
    }
    std::unordered_map<int, PrDes::REposStruct >::iterator itr;
    for (it = chrIndex.begin(); it != chrIndex.end(); ++it) {
        for ( itr = OneDesign.back().Layer[it->second].repmap.begin(); itr != OneDesign.back().Layer[it->second].repmap.end(); ++itr)
            WritetoFile(outfile, it->first, it->second, itr->first, itr->second.prom_indexes, itr->second.whichside, reInfo.desName, Feats);
        }
    dLog << "First design layer (no filter) written to Output file " << std::endl;

    
    outfile.close();
}



int DesignClass::CheckREsiteAroundProbe(RESitesClass& dpnII, std::string chr, int probe_re, int direction){
    // It checks if for neighbouring restriction enzyme positions any probes designed. if yes, it will move it to a new design layer.
    int *resites;
    resites = new int[2];
    bool refound = 0;
    
    refound = dpnII.GettheREPositions(chr,probe_re,resites);
    
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
			target="SNP";
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
		//probestart = ( (repos+1) - reLeftCut );
		//probeend =( (repos+1) + ProbeLen - reLeftCut); //+1 to select the right fragment
		probestart = ( (repos+1) );
		probeend =( (repos+1) + ProbeLen ); //+1 to select the right fragment
		side = "L";
		disttotss = abs(probestart - feats.promFeatures[values[0]].TSS);  
	}   

    outfile  << chr  << '\t' << "." << '\t' <<"probe"<<'\t';
    
    outfile <<probestart << '\t' << probeend << '\t'<< "." << '\t'<< "." << '\t'<< "." << '\t'; // to adjust for 1-based coords
        
    outfile << "Name="<<feats.promFeatures[values[0]].genes[0]<<"; " << "side="<<side<<"; "<<"target="<<target<<"; "<<"design="<<design<<"; "<< "featuresinvicinity=";
        
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
		outfile << "none"<<"; " <<'\t';
	}
	
	outfile<< "targettss="<<feats.promFeatures[values[0]].TSS<<"; "<<"distancetotss="<<disttotss<< std::endl;
    
    return true;
}

void DesignClass::MergeAllChrOutputs(ProbeFeatureClass& Feats, PrDes::RENFileInfo& reInfo){
	
	std::string fnameAllChr = reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".AllProbes."+reInfo.REName+"."+reInfo.currTime+".gff3";    
    
    std::ofstream outfile;
    outfile.open(fnameAllChr, std::fstream::out);

    outfile<<"##gff-version 3.2.1"<<std::endl;
    outfile<<"##genome-build "<<reInfo.genomeAssembly.substr(reInfo.genomeAssembly.find_first_of(',')+1)<<" "<<reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(',')) <<std::endl;
	
	outfile.close();
	
	outfile.open(fnameAllChr, std::ios_base::binary | std::ios_base::app);
	
	for(auto &iChr : Feats.ChrNames_proms){
		
		//std::string fName=reInfo.desName+"_"+iChr+"."+reInfo.REName+"."+reInfo.currTime+".gff3";
		std::string fName=reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+"_"+iChr+"."+reInfo.REName+"."+reInfo.currTime+".gff3";
		std::ifstream readChrFiles(fName, std::ios_base::binary);
		readChrFiles.seekg(20); //strip header
		outfile << readChrFiles.rdbuf();
		readChrFiles.close();
		
		remove(fName.c_str());
		
	}
	
	
}

