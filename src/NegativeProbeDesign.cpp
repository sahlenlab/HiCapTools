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
//  NegativeProbeDesign.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include <fstream>
#include <random>
#include <algorithm>
#include <sstream>
#include "NegativeProbeDesign.h"


int NegativeProbeDesign::InitialiseDesign(ProbeFeatureClass& Features, std::string transcriptfile, std::string regRegionFile, bool ifRReg, int prlen, std::string bgwgsumbin, std::string mapfp, int maxDistToTss, PrDes::RENFileInfo& reInfo, int bufSize, std::string digestFile, int dfI, int dfReg, std::string extranscriptfile, int dfProm, int minDistToTss){
	
	//Set class variables
	MaxDistancetoTSS = maxDistToTss;
	minDistToTSS = minDistToTss;
	ProbeLen = prlen;
	bigwigsummarybinary=bgwgsumbin;
	mappabilityfile=mapfp;
	ifExistRegRegionFile=ifRReg;
	designName=reInfo.desName;
	fileName=reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".NegCtrlProbes."+reInfo.REName+"."+reInfo.currTime+".gff3";
	fasFileName = reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".NegProbeSequences."+reInfo.REName+"."+reInfo.currTime+".txt";  
	write2ProbesBedFileName =   reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".AllNegativeControlProbeSequences."+reInfo.REName+"."+reInfo.currTime+".bed";
	reLeftCut = reInfo.leftOfCut;
	reRightCut = reInfo.rightOfCut;
	mapThreshold = reInfo.mappabilityThreshold;
	repOverlapExtent = reInfo.repeatOverlapExtent;
	ifRep = reInfo.ifRepeatAvail;
	ifMap = reInfo.ifMapAvail;
	ifCountNs = reInfo.ifCountNs;
	BUFSIZE=bufSize;
	genAssem=reInfo.genomeAssembly;
	summaryFileName = reInfo.desName + "."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".NegCtrlRegions_" +reInfo.REName+"."+reInfo.currTime+".bed";
	distforbidIntergenic=dfI; 
	distforbidReg=dfReg;
	distforbidProm=dfProm;
	

	std::string temp, lineReadin;
	std::string chr, start, end, strand, intId, discardtemp;	
	std::string genename, tr_id, exonCounts, exonStarts, exonEnds;
	
	// Store gene
	struct genetemp{
		std::string chr;
		int start;
		int end;
		bool multiFlag;
	};
	
	bool flag=true;
	int istart, iend;
	
	
	std::map <std::string, genetemp> genemap;
	std::map< std::string, std::vector < Interval < std::string > > > genechrIntervals;
	std::map< std::string, std::vector < Interval < std::string > > > exonchrIntervals;
	//std::map< std::string, std::vector < Interval < std::string > > > intronchrIntervals;
	std::map< std::string, std::vector < Interval < std::string > > > regRegchrIntervals;
	std::map< std::string, std::vector < Interval < std::string > > > extraTranscriptchrIntervals;
	
	std::ifstream fileReadin;
	
	fileReadin.open(transcriptfile, std::ifstream::in);
	getline(fileReadin,temp); //discard header
		
	 while (getline(fileReadin, lineReadin)){
      		
		std::stringstream geneline ( lineReadin );
		
		//#bin	name	chrom	strand	txStart	txEnd	cdsStart	cdsEnd	exonCount	exonStarts	exonEnds	score	name2	cdsStartStat	cdsEndStat	exonFrames
		
		getline(geneline,discardtemp,'\t');
		
		getline(geneline, tr_id,'\t'); 
	
		getline(geneline, chr,'\t'); 
        	  
		getline(geneline, strand,'\t');
    
		getline(geneline, start,'\t');
		getline(geneline, end,'\t');
		
		getline(geneline,discardtemp,'\t');
		getline(geneline,discardtemp,'\t');
		
		getline(geneline, exonCounts,'\t');
		getline(geneline, exonStarts,'\t');
		getline(geneline, exonEnds,'\t');
		
		getline(geneline,discardtemp,'\t');
		
		getline(geneline,genename,'\t');
		
		getline(geneline,discardtemp,'\t');
		getline(geneline,discardtemp,'\t');
		getline(geneline,discardtemp,'\t');
		
		genetemp gt;
		gt.chr=chr;
		gt.start=std::stoi(start);
		gt.end=std::stoi(end);
		auto it =genemap.find(genename);//chr+genename
		
		if (it!=genemap.end()){		
			if(it->second.chr ==gt.chr){
				if(abs(gt.end-gt.start)>abs(it->second.end-it->second.start)){ // select longest transcript for gene
					it->second=gt;
				}
			}
			//log multiple copies of this gene were observed in different chromosomes - no probes 
		}
		else{
			genemap.emplace(genename, gt);
		}
		
		size_t eStarts=std::count(exonStarts.begin(), exonStarts.end(), ',');
		size_t eEnds=std::count(exonEnds.begin(), exonEnds.end(), ',');
		
		if(eStarts==eEnds){ // if both have same number of values
			if((exonStarts.back()!=','&&exonEnds.back()!=',')){ //add , to ends if not present
				exonStarts.push_back(',');
				exonEnds.push_back(',');
				eStarts++;
				eEnds++;
			}
			size_t pos11, pos12, pos21, pos22;
			int exS, exE;
			std::string ifProm;
						
			while(eStarts > 0 && eEnds > 0){ //while values still exist in the list
				pos11=0; pos21=0;
				pos12=exonStarts.find_first_of(','); pos22=exonEnds.find_first_of(',');    
				
				exS=std::stoi(exonStarts.substr(pos11, pos12)); //get first values
				exE=std::stoi(exonEnds.substr(pos21, pos22));
				
				//check prom tree for exons
				ifProm = FindOverlaps(Features.promIntTree, chr, exS, exE);
				
				if(ifProm=="null"){
				//exon intervals added 
					if(exonchrIntervals.find(chr)!=exonchrIntervals.end())
						exonchrIntervals[chr].push_back(Interval <std::string>(exE , exS, chr+":"+exonStarts.substr(pos11, pos12)+"-"+exonEnds.substr(pos21, pos22)));
					else{
						std::vector < Interval < std::string > > tempvector ;
						tempvector.push_back(Interval <std::string>(exE , exS, chr+":"+exonStarts.substr(pos11, pos12)+"-"+exonEnds.substr(pos21, pos22)));
						exonchrIntervals.emplace(chr, tempvector);
					}
				}
				exonStarts.erase(pos11, pos12+1); exonEnds.erase(pos21, pos22+1); //erase first values from string
				--eStarts; --eEnds;
				
				/***calculate introns
				 * 
				 * 
				if(intE<exS){
					intE=exS;
				}
				
				if(intS<intE)
				//add to intr interval intS intE
				
				intS=exE;
				intE=exE; 
				***/		
			}
			/*** Calculate introns
			 * 
			if(intE<end){
				//std::cout<<"Intron :"<<intE<< " "<<end<<std::endl;
				//Add last intron if any to interval tree 
			}
			***/
		}
	}
	
	fileReadin.close();
	
	//add gene intervals 
	for(auto it=genemap.begin(); it!=genemap.end(); ++it){
		
		if(genechrIntervals.find(it->second.chr)!=genechrIntervals.end())
			genechrIntervals[it->second.chr].push_back(Interval <std::string>(it->second.start , it->second.end, it->second.chr+":"+it->first));
		else{
			std::vector < Interval < std::string > > tempvector ;
			tempvector.push_back(Interval <std::string>(it->second.start , it->second.end, it->second.chr+":"+it->first));
			genechrIntervals.emplace(it->second.chr, tempvector);
		}	
	}
	
	//construct gene tree
	for(auto it = genechrIntervals.begin(); it != genechrIntervals.end(); ++it){
			std::vector< Interval < std::string > > temp;
            geneIntTree[it->first] = IntervalTree< std::string >(it->second);
            it->second.swap(temp);
    }
    
    //construct exon tree
    for(auto it = exonchrIntervals.begin(); it != exonchrIntervals.end(); ++it){
			std::vector< Interval < std::string > > temp;
            exonIntTree[it->first] = IntervalTree< std::string >(it->second);
            it->second.swap(temp);
    }
    
    //extra transcript checking
    
	fileReadin.open(extranscriptfile, std::ifstream::in);
	while(!fileReadin.eof()){
		getline(fileReadin, lineReadin);
		if(lineReadin.substr(0,5)!="track"){
			std::stringstream bedline(lineReadin);
			getline(bedline,chr,'\t');
			getline(bedline,start,'\t');
			istart=std::stoi(start);
			getline(bedline,end,'\t');
			iend=std::stoi(end);
			
			intId=chr+":"+start+"-"+end;
		
			if(extraTranscriptchrIntervals.find(chr)!=extraTranscriptchrIntervals.end())
				extraTranscriptchrIntervals[chr].push_back(Interval <std::string>(istart , iend, intId));
			else{
				std::vector < Interval < std::string > > tempvector ;
				tempvector.push_back(Interval <std::string>(istart ,iend, intId));
				extraTranscriptchrIntervals.emplace(chr, tempvector);
				}			
		}	
	}
		
		//construct extra transcript tree
	for(auto it = extraTranscriptchrIntervals.begin(); it != extraTranscriptchrIntervals.end(); ++it){
		std::vector< Interval < std::string > > temp;
        extraPromTree[it->first] = IntervalTree< std::string >(it->second);
        it->second.swap(temp);
	}
	
    fileReadin.close();
		
	//if user provided regulatory region file exists
    if(ifExistRegRegionFile){
		fileReadin.open(regRegionFile, std::ifstream::in);
		while(!fileReadin.eof()){
			getline(fileReadin, lineReadin);
			if(lineReadin.substr(0,5)!="track"){
				std::stringstream bedline(lineReadin);
				getline(bedline,chr,'\t');
				getline(bedline,start,'\t');
				istart=std::stoi(start);
				getline(bedline,end,'\t');
				iend=std::stoi(end);
			
				intId=chr+":"+start+"-"+end;
			
				if(regRegchrIntervals.find(chr)!=regRegchrIntervals.end())
					regRegchrIntervals[chr].push_back(Interval <std::string>(istart , iend, intId));
				else{
					std::vector < Interval < std::string > > tempvector ;
					tempvector.push_back(Interval <std::string>(istart ,iend, intId));
					regRegchrIntervals.emplace(chr, tempvector);
				}			
			}
		}
		
		//construct extra regulatory region tree
		for(auto it = regRegchrIntervals.begin(); it != regRegchrIntervals.end(); ++it){
			std::vector< Interval < std::string > > temp;
            regulRegIntTree[it->first] = IntervalTree< std::string >(it->second);
            it->second.swap(temp);
		}
	}
	ConstructPools(digestFile, Features);
	return 1;
}


void NegativeProbeDesign::ConstructPools(std::string digestFile, ProbeFeatureClass& Features){
	
	bool flag;
	int istart, iend;
	std::string temp, lineReadin, chr, start, end, ifGene, ifProm, ifRegulReg, ifExon, ifIntron, ifExtraTrans;
	std::ifstream fileReadin;
	
	fileReadin.open(digestFile, std::ifstream::in);
	
	std::getline(fileReadin,temp); //get the header row1
	std::getline(fileReadin,temp); //get the header row
	
	while(fileReadin >> chr >> start >> end >> temp >> temp >> temp >> temp){
		
		 // Read the fragments lines
		istart=std::stoi(start);
		iend=std::stoi(end);
		//Check if RE fragment overlaps with promoter regions
		ifProm = FindOverlaps(Features.promIntTree, chr, istart, iend); // Check if overlaps with forbidden region - prom tree with intervals 50 kb both up/downstream of feature 
		int transtart = istart-distforbidProm-MaxDistancetoTSS;
		if(transtart<0)
			transtart=0;
		ifExtraTrans=FindOverlaps(extraPromTree, chr, transtart, iend+distforbidProm+MaxDistancetoTSS);
		
		//Check if RE fragments overlap with other regulatory regions
		if(ifExistRegRegionFile){
			int regstart = istart-distforbidReg;
			if(regstart<0)
				regstart=0;
			ifRegulReg = FindOverlaps(regulRegIntTree, chr, regstart, iend+distforbidReg);
			if(ifRegulReg=="null"){
				flag=true;
			}
			else
				flag=false;
		}
		else{
			flag=true;
		}
		
		if(ifProm=="null" && ifExtraTrans == "null" && flag){
			
			PrDes::FeatureStruct feat;
			feat.ProbeID.push_back(chr+":"+start+"-"+end);
			feat.chr=chr;
			feat.start=istart;
			feat.end=iend;
			
			int genstart = istart-distforbidIntergenic;
			if(genstart<0)
				genstart=0;
			
			ifGene=FindOverlaps(geneIntTree, chr, genstart, iend+distforbidIntergenic); //add 50 kb around gene get value from user
			
			if(ifGene=="null"){
				if(intergenicPool.find(chr)==intergenicPool.end()){
					std::vector <PrDes::FeatureStruct> tempvector ;
					tempvector.push_back(feat);
					intergenicPool.emplace(chr, tempvector);
				}
				else
					intergenicPool[chr].push_back(feat);
			}
			else{
				ifExon=FindOverlaps(exonIntTree, chr, istart, iend);
				if(ifExon=="null"){
					ifIntron=FindOverlaps(geneIntTree, chr, istart, iend); //check only gene region					
					if(ifIntron!="null"){
							intronPool.push_back(feat);
					}
				}
				else{
					if(exonPool.find(chr)==exonPool.end()){
							std::vector <PrDes::FeatureStruct> tempvector ;
							tempvector.push_back(feat);
							exonPool.emplace(chr, tempvector);
					}
					else
						exonPool[chr].push_back(feat);
					
				}
			}
		}		
	}
}


std::string NegativeProbeDesign::FindOverlaps( std::map< std::string, IntervalTree <std::string> >& theTree ,std::string chr, unsigned long int readstart, unsigned long int readend){
    
    std::vector<Interval< std::string > > overlapResult;
    theTree[chr].findOverlapping(readstart, readend, overlapResult);
    if (overlapResult.size() > 0){ // value = probe_index
        return overlapResult[0].value;
    }
    else
        return "null";
}


int NegativeProbeDesign::ConstructNegativeControlProbes(int nCtrls,std::string nCtrlType,  Repeats& repeatTrees, PrDes::RENFileInfo& reInfo, RESitesClass& dpnII,bioioMod& getSeq){
	
	if(nCtrls<=0)
		return 1;
	
	std::string outFileName;
	std::ofstream outFile, summaryFile;
	bool flag=false;
	auto cmp = [](std::pair<std::string,float> const & a, std::pair<std::string,float> const & b) {
		return  a.second > b.second;
    };
	
	std::string indPath, chr;
	std::string temp, lineReadin, length;
	std::map<std::string, long int> totalChr;
	long int total=0;
	////////////
	if(nCtrlType!="intronic"){
		if (reInfo.fastaindexpath.empty()) { 
			indPath = reInfo.fastafilepath;
			indPath.replace(indPath.begin() + indPath.find_last_of("."), indPath.end(), ".fai");
		}
		else
			indPath=reInfo.fastaindexpath;
		
		std::ifstream indFile {indPath, std::ifstream::in};
		if (!indFile) {
			indFile.open(reInfo.fastafilepath + ".fai");
			if (!indFile) {
				dLog<<"Fasta index file not found!!!"<<std::endl;
				return 0;
			}
		}
	
		if(indFile.good()){
			while (getline(indFile, lineReadin)){  		
				std::stringstream indLine ( lineReadin );
				indLine>>chr>>length>>temp>>temp>>temp;
				if(chr.length()<6){
					totalChr.emplace(chr, std::stoi(length));
					total=total + std::stoi(length);
				}		
			}
		}
		indFile.close();
		int checkTotal=nCtrls;
		std::vector<std::pair<std::string, float>> frac;
		numProbesPerChr.clear();
		for(auto it=totalChr.begin(); it!=totalChr.end(); ++it){
			float num= (it->second/float(total))*nCtrls;
			numProbesPerChr.emplace(it->first, (int)num);
			checkTotal =checkTotal-(int)num;
			frac.push_back(std::make_pair(it->first, num-(int)num)); 
		}
		
		std::sort(frac.begin(), frac.end(), cmp);
		
		if(checkTotal>0){
			while(checkTotal>0){
				for(auto it=frac.begin(); it!=frac.end() && checkTotal>0; ++it){
					numProbesPerChr[it->first]= numProbesPerChr[it->first]+1;
					checkTotal=checkTotal-1;
				}
			}
		}
	}
		
	outFileName.append(fileName);

    std::ifstream checkFile(outFileName);
    std::ifstream checkSFile(summaryFileName);
    if(checkFile.good()&&checkSFile.good())
		flag=true;
    
    checkFile.close();
    checkSFile.close();
    
    outFile.open(outFileName, std::fstream::app); // append to existing file
	summaryFile.open(summaryFileName, std::fstream::app);
	
	if(!flag){
		outFile<<"##gff-version 3.2.1"<<std::endl;
		outFile<<"##genome-build "<<genAssem.substr(genAssem.find_first_of(',')+1)<<" "<<genAssem.substr(0, genAssem.find_first_of(',')) <<std::endl;
		summaryFile<<"track name=\""<<designName<<"\" type=bedDetail description=\"Negative control probe target fragments for "<< designName<<"\""<<std::endl;
	}
	
	size_t poolSize=0;
	
	if(nCtrlType=="exonic" && nCtrls>0){
		
		for(auto it=exonPool.begin(); it!=exonPool.end();++it){
			poolSize=poolSize+it->second.size();
		}
		if(nCtrls>poolSize){
			dLog<<"!!!Error!!! The provided number of exonic negative control probes required are greater than the number of available candidate RE fragments. No exonic negative control probes designed"<<std::endl;
			return 0;
		}
		chooseRandomProbesFromPool(nCtrls, exonPool, repeatTrees, outFile, "Exon", summaryFile, poolSize, dpnII, getSeq);
	}
	else if(nCtrlType=="intronic" && nCtrls>0){		
		poolSize=intronPool.size();
		if(nCtrls>poolSize){
			dLog<<"!!!Error!!! The provided number of intronic negative control probes required are greater than the number of available candidate RE fragments. No intronic negative control probes designed"<<std::endl;
			return 0;
		}
			chooseRandomProbesFromPool(nCtrls, intronPool, repeatTrees, outFile, "Intron", summaryFile, poolSize, dpnII, getSeq);
	}
	else if(nCtrlType=="intergenic" && nCtrls>0){
		for(auto it=intergenicPool.begin(); it!=intergenicPool.end();++it){
			poolSize=poolSize+it->second.size();
		}
		if(nCtrls>poolSize){
			dLog<<"!!!Error!!! The provided number of intergenic negative control probes required are greater than the number of available candidate RE fragments. No intergenic negative control probes designed"<<std::endl;
			return 0;
		}
		chooseRandomProbesFromPool(nCtrls, intergenicPool, repeatTrees, outFile, "Intergen", summaryFile, poolSize, dpnII, getSeq);
	}
	return 1;
}

void NegativeProbeDesign::chooseRandomProbesFromPool(int nProbesReq, std::map<std::string, std::vector<PrDes::FeatureStruct>>& whichgenPool, Repeats& repeatTrees, std::ofstream &outfile, std::string whichPool, std::ofstream &summaryfile, int poolSize, RESitesClass& dpnII, bioioMod& getSeq){
	
	bool leftEnd, rightEnd;
	int loopIndex, total=0;
	std::mt19937_64 rng;
	for(auto it=numProbesPerChr.begin(); it!=numProbesPerChr.end();++it){
		std::uniform_int_distribution<int> fragDist(0, whichgenPool[it->first].size()-1);
		std::vector<int> fragSeen;
		std::vector<int> fragSelected;
	
		rng.seed(std::random_device()());
		loopIndex=0;
		
		while(loopIndex<it->second){
		//below gives you index of negative control probe to write
			int currIndex=fragDist(rng);
			bool checkDistFlag=true;
		
			if(std::find(fragSeen.begin(), fragSeen.end(), currIndex)==fragSeen.end()){
				
				for(auto elem=fragSelected.begin(); elem < fragSelected.end(); ++elem){
					int checkDist;
					if(whichgenPool[it->first][currIndex].start<whichgenPool[it->first][*elem].start)
						checkDist=whichgenPool[it->first][*elem].start-whichgenPool[it->first][currIndex].end;
					else
						checkDist=whichgenPool[it->first][currIndex].start-whichgenPool[it->first][*elem].end;
						
					if(checkDist<5000){//hard coded value
						checkDistFlag=false;
						break;
					}
				}
		
				fragSeen.push_back(currIndex);
				if(fragSeen.size()>=whichgenPool[it->first].size()){
					dLog<<"Not enough valid candidate fragments for "<<whichPool <<" in "<<it->first<<" for provided parameters. The provided number of negative controls may not have been selected. Please try running again with/without changing parameters!!"<<std::endl;
					break;
				}
				//check repeat overlap'
				//RE fragment start
				if(checkDistFlag){
					int takeTSS=(whichgenPool[it->first][currIndex].start + whichgenPool[it->first][currIndex].end)/2;
					int left_res, right_res; 	
					bool passed_upstream = false, passed_downstream = false;
            
					left_res = CheckDistanceofProbetoTSS(dpnII, it->first, takeTSS, whichgenPool[it->first][currIndex].start, 0);
					right_res = CheckDistanceofProbetoTSS(dpnII, it->first, takeTSS, whichgenPool[it->first][currIndex].end, 1);
					
					passed_upstream = CheckRESite(dpnII, getSeq, it->first, takeTSS, left_res, 0, repeatTrees, ifRep, ifMap, ifCountNs);
					passed_downstream = CheckRESite(dpnII, getSeq, it->first, takeTSS, right_res, 1, repeatTrees, ifRep, ifMap, ifCountNs);
					
									
					if (passed_upstream && passed_downstream){
					//write to file	 - 2 probes for each fragment
						if(chrToIndex.find(whichgenPool[it->first][currIndex].chr)!=chrToIndex.end())
							chrToIndex[it->first].push_back(negProbe(left_res, right_res, "NegCtrl_"+whichPool+std::to_string(total+1)));
						else{
							std::vector<negProbe> temp;
							temp.push_back(negProbe(left_res, right_res, "NegCtrl_"+whichPool+std::to_string(total+1)));
							chrToIndex.emplace(it->first, temp);
						}
						summaryfile << it->first << '\t'  << whichgenPool[it->first][currIndex].start << '\t' << whichgenPool[it->first][currIndex].end << '\t' << "NegCtrl_"<<whichPool<<total+1 << '\t'<<"."<<'\t'<<"+"<<'\t' << "NegCtrl_"<<whichPool<<total+1 << '\t'<< "Neg_Ctrl"<< std::endl;
						++loopIndex;
						++total;
						fragSelected.push_back(currIndex);
					}				
				}
			}
		}
	}
}

void NegativeProbeDesign::chooseRandomProbesFromPool(int nProbesReq, std::vector<PrDes::FeatureStruct>& whichgenPool, Repeats& repeatTrees, std::ofstream &outfile, std::string whichPool, std::ofstream &summaryfile, int poolSize, RESitesClass& dpnII, bioioMod& getSeq){
	
	bool leftEnd, rightEnd;
	int loopIndex, total=0;
	std::mt19937_64 rng;
	
	std::uniform_int_distribution<int> fragDist(0, whichgenPool.size()-1);
	std::vector<int> fragSeen;
	std::vector<int> fragSelected;
	
	rng.seed(std::random_device()());
	loopIndex=0;
		
	while(loopIndex<nProbesReq){
	//	while(loopIndex<it->second){
		//below gives you index of negative control probe to write
		int currIndex=fragDist(rng);
		bool checkDistFlag=true;
		
		if(std::find(fragSeen.begin(), fragSeen.end(), currIndex)==fragSeen.end()){
				
			for(auto elem=fragSelected.begin(); elem < fragSelected.end(); ++elem){
				int checkDist=5000;
				if(whichgenPool[currIndex].chr==whichgenPool[*elem].chr){	
					if(whichgenPool[currIndex].start<whichgenPool[*elem].start)
						checkDist=whichgenPool[*elem].start-whichgenPool[currIndex].end;
					else
						checkDist=whichgenPool[currIndex].start-whichgenPool[*elem].end;
				}	
				if(checkDist<5000){//hard coded value
					checkDistFlag=false;
					break;
				}
			}
		
			fragSeen.push_back(currIndex);
			if(fragSeen.size()>=whichgenPool.size()){
				dLog<<"Not enough valid candidate fragments in "<<whichPool<<" Pool for provided parameters. Please try running again with/without changing parameters!!"<<std::endl;
				break;
			}
			//check repeat overlap'
			//RE fragment start
			if(checkDistFlag){
				int takeTSS=(whichgenPool[currIndex].start + whichgenPool[currIndex].end)/2;
				int left_res, right_res; 	
				bool passed_upstream = false, passed_downstream = false;
            
				left_res = CheckDistanceofProbetoTSS(dpnII, whichgenPool[currIndex].chr, takeTSS, whichgenPool[currIndex].start, 0);
				right_res = CheckDistanceofProbetoTSS(dpnII, whichgenPool[currIndex].chr, takeTSS, whichgenPool[currIndex].end, 1);
					
				passed_upstream = CheckRESite(dpnII, getSeq, whichgenPool[currIndex].chr, takeTSS, left_res, 0, repeatTrees, ifRep, ifMap, ifCountNs);
				passed_downstream = CheckRESite(dpnII, getSeq, whichgenPool[currIndex].chr, takeTSS, right_res, 1, repeatTrees, ifRep, ifMap, ifCountNs);
							
				if (passed_upstream && passed_downstream){
					//write to file	 - 2 probes for each fragment
					if(chrToIndex.find(whichgenPool[currIndex].chr)!=chrToIndex.end())
						chrToIndex[whichgenPool[currIndex].chr].push_back(negProbe(left_res, right_res, "NegCtrl_"+whichPool+std::to_string(total+1)));
					else{
						std::vector<negProbe> temp;
						temp.push_back(negProbe(left_res, right_res, "NegCtrl_"+whichPool+std::to_string(total+1)));
						chrToIndex.emplace(whichgenPool[currIndex].chr, temp);
					}
					summaryfile << whichgenPool[currIndex].chr << '\t'  << whichgenPool[currIndex].start << '\t' << whichgenPool[currIndex].end << '\t' << "NegCtrl_"<<whichPool<<total+1 << '\t'<<"."<<'\t'<<"+"<<'\t' << "NegCtrl_"<<whichPool<<total+1 << '\t'<< "Neg_Ctrl"<< std::endl;
					++loopIndex;
					++total;
					fragSelected.push_back(currIndex);
				}
			}
		}
	}
}

void NegativeProbeDesign::WritetoFile(bioioMod& getSeq, PrDes::RENFileInfo& reInfo){
	
	std::string target, side;
	int disttotss, probeStart, probeEnd;
    target ="neg_ctrl";
	
	std::ofstream outfile(fileName, std::fstream::app);
	std::ofstream outFasFile(fasFileName, std::fstream::app);
	std::ofstream write2ProbesBedFile(write2ProbesBedFileName, std::fstream::app);
	
	outFasFile<<"TargetID"<<'\t'<<"ProbeID"<<'\t'<<"Sequence"<<'\t'<<"Replication"<<'\t'<<"Strand"<<'\t'<<"Coordinates"<<std::endl;	
	write2ProbesBedFile<<"track name=\"AllNegativeControlProbes_"<<reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))<<"_"<<reInfo.REName<<"_"<<reInfo.currTime<<"\""<<std::endl;
	
	for(auto it=chrToIndex.begin(); it!=chrToIndex.end(); ++it){
		std::sort(it->second.begin(), it->second.end(), [](negProbe const &a, negProbe const &b) { return a.start < b.start; });
		for(auto ind : it->second){
			//left side
			
			probeStart = ( ind.start+1 ); /////check 1 based 0 based
			probeEnd =probeStart+ProbeLen;
			side="L";
			
			outfile  << it -> first  << '\t' << "." << '\t' <<"probe"<<'\t';
    
			outfile <<probeStart+1 << '\t' << probeEnd << '\t'<< "." << '\t'<< "." << '\t'<< "." << '\t'; // to adjust for 1-based coords
			
			outfile << "Name="<< ind.name <<"; " <<"transcriptid="<< "none"<<"; " << "side="<<side<<"; "<<"target="<<target<<"; "<<"design="<<designName<<"; "<< "featuresinvicinity=";
          
			outfile << "none"<<"; ";
	
			outfile<< "targettss="<<"none; " <<"distancetotss="<<"not applicable"<< std::endl;
			
			std::string getFas=getSeq.GetFasta(it->first+":"+std::to_string(probeStart-1)+"-"+std::to_string(probeEnd-1)); // subtract 1 for bioio coordinates on Left side
			
			outFasFile<<ind.name<<'\t'<< it->first<<"_"<<probeStart<<'\t'<<getFas<<'\t'<<"1"<<'\t'<<"+"<<'\t'<<it->first<<":"<<probeStart<<"-"<<probeEnd-1<<std::endl;
			
			write2ProbesBedFile<<it->first<<'\t'<<probeStart-1<<'\t'<<probeEnd-1<<'\t'<<ind.name<<std::endl;
			
			//right side
			probeEnd=( ind.end + reRightCut );
			probeStart=( probeEnd - ProbeLen ); 
			side="R";
			
			outfile  << it -> first  << '\t' << "." << '\t' <<"probe"<<'\t';
    
			outfile <<probeStart + 1 << '\t' << probeEnd << '\t'<< "." << '\t'<< "." << '\t'<< "." << '\t'; // to adjust for 1-based coords
        
			outfile << "Name="<< ind.name <<"; " <<"transcriptid="<< "none"<<"; " << "side="<<side<<"; "<<"target="<<target<<"; "<<"design="<<designName<<"; "<< "featuresinvicinity=";
          
			outfile << "none"<<"; " ;
	
			outfile<< "targettss="<<"none; "<<"distancetotss="<<"not applicable"<< std::endl;
			
			getFas=getSeq.GetFasta(it->first+":"+std::to_string(probeStart)+"-"+std::to_string(probeEnd));
			
			outFasFile<<ind.name<<'\t'<< it->first<<"_"<<probeStart+1<<'\t'<<getFas<<'\t'<<"1"<<'\t'<<"+"<<'\t'<<it->first<<":"<<probeStart+1<<"-"<<probeEnd<<std::endl;
			
			write2ProbesBedFile<<it->first<<'\t'<<probeStart<<'\t'<<probeEnd<<'\t'<<ind.name<<std::endl;
		}
	}
}
