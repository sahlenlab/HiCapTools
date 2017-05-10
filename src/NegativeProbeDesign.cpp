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


void NegativeProbeDesign::InitialiseDesign(ProbeFeatureClass& Features, std::string transcriptfile, std::string regRegionFile, bool ifRReg, int prlen, std::string bgwgsumbin, std::string mapfp, int minREfragLength, PrDes::RENFileInfo& reInfo, int bufSize, std::string digestFile){
	
	//Set class variables
	minREfragLen = minREfragLength;
	ProbeLen = prlen;
    bigwigsummarybinary=bgwgsumbin;
    mappabilityfile=mapfp;
	ifExistRegRegionFile=ifRReg;
	designName=reInfo.desName;
	fileName=reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".NegCtrlProbes."+reInfo.REName+"."+reInfo.currTime+".gff3";
	fasFileName = reInfo.desName+"."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".NegProbeSequences."+reInfo.REName+"."+reInfo.currTime+".txt";    
	reLeftCut = reInfo.leftOfCut;
	reRightCut = reInfo.rightOfCut;
	mapThreshold = reInfo.mappabilityThreshold;
	ifRep = reInfo.ifRepeatAvail;
	ifMap = reInfo.ifMapAvail;
	BUFSIZE=bufSize;
	genAssem=reInfo.genomeAssembly;
	summaryFileName = reInfo.desName + "."+reInfo.genomeAssembly.substr(0, reInfo.genomeAssembly.find_first_of(','))+".NegCtrlRegions_" +reInfo.REName+"."+reInfo.currTime+".bed";
	
	// Store gene
	struct genetemp{
		std::string chr;
		int start;
		int end;
		bool multiFlag;
	};
	
	bool flag=true;
	int istart, iend;
	std::string temp, lineReadin;
	std::string chr, start, end, strand, intId;	
	std::string genename, tr_id, exonCounts, exonStarts, exonEnds;
	
	std::map <std::string, genetemp> genemap;
	std::map< std::string, std::vector < Interval < std::string > > > genechrIntervals;
	std::map< std::string, std::vector < Interval < std::string > > > exonchrIntervals;
	//std::map< std::string, std::vector < Interval < std::string > > > intronchrIntervals;
	std::map< std::string, std::vector < Interval < std::string > > > regRegchrIntervals;
	
	std::ifstream fileReadin;
	
	fileReadin.open(transcriptfile, std::ifstream::in);
	getline(fileReadin,temp); //discard header
		
	 while (getline(fileReadin, lineReadin)){
      		
		std::stringstream geneline ( lineReadin );
		
		//name2	 name	chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds
		getline(geneline,genename,'\t');
		 
		getline(geneline, tr_id,'\t'); 
		
		getline(geneline, chr,'\t'); 
        	  
		getline(geneline, strand,'\t');
    
		getline(geneline, start,'\t');
		getline(geneline, end,'\t');
		
		getline(geneline, exonCounts,'\t');
		getline(geneline, exonStarts,'\t');
		getline(geneline, exonEnds,'\t');
		
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
}


void NegativeProbeDesign::ConstructPools(std::string digestFile, ProbeFeatureClass& Features){
	
	bool flag;
	int istart, iend;
	std::string temp, lineReadin, chr, start, end, ifGene, ifProm, ifRegulReg, ifExon, ifIntron;
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
		
		//Check if RE fragments overlap with other regulatory regions
		if(ifExistRegRegionFile){
			ifRegulReg = FindOverlaps(regulRegIntTree, chr, istart-distforbidReg, iend+distforbidReg);
			if(ifRegulReg=="null"){
				flag=true;
			}
			else
				flag=false;
		}
		else{
			flag=true;
		}
		
		if(ifProm=="null" && flag){
			
			PrDes::FeatureStruct feat;
			feat.ProbeID.push_back(chr+":"+start+"-"+end);
			feat.chr=chr;
			feat.start=istart;
			feat.end=iend;
			
			ifGene=FindOverlaps(geneIntTree, chr, istart-distforbidIntergenic, iend+distforbidIntergenic); //add 50 kb around gene get value from user
			
			if(ifGene=="null"){
				intergenicPool.push_back(feat);
			}
			else{
				ifExon=FindOverlaps(exonIntTree, chr, istart, iend);
				if(ifExon=="null"){
					ifIntron=FindOverlaps(geneIntTree, chr, istart, iend); //check only gene region					
					if(ifIntron!="null")						
						intronPool.push_back(feat);
				}
				else{
						exonPool.push_back(feat);
					
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


int NegativeProbeDesign::ConstructNegativeControlProbes(int nCtrls,std::string nCtrlType,  Repeats repeatTrees, int dfI, int dfReg){
	
	
	distforbidIntergenic=dfI; 
	distforbidReg=dfReg;
	
	std::string outFileName;
	std::ofstream outFile, summaryFile;
	bool flag=false;
	
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
		summaryFile<<"track name=\""<<designName<<"\" description=\"Negative control probe target fragment for "<< designName<<"\""<<std::endl;
	}

	
	if(nCtrlType=="exonic" && nCtrls>0){
		if(nCtrls>exonPool.size()){
			dLog<<"!!!Error!!! The provided number of exonic negative control probes required are greater than the number of available candidate RE fragments. No exonic negative control probes designed"<<std::endl;
			return 0;
		}
		chooseRandomProbesFromPool(nCtrls, exonPool, repeatTrees, outFile, "Exon", summaryFile);
	}
	else if(nCtrlType=="intronic" && nCtrls>0){
		if(nCtrls>intronPool.size()){
			dLog<<"!!!Error!!! The provided number of intronic negative control probes required are greater than the number of available candidate RE fragments. No intronic negative control probes designed"<<std::endl;
			return 0;
		}
		chooseRandomProbesFromPool(nCtrls, intronPool, repeatTrees, outFile, "Intron", summaryFile);
	}
	else if(nCtrlType=="intergenic" && nCtrls>0){
		if(nCtrls>intergenicPool.size()){
			dLog<<"!!!Error!!! The provided number of intergenic negative control probes required are greater than the number of available candidate RE fragments. No intergenic negative control probes designed"<<std::endl;
			return 0;
		}
		chooseRandomProbesFromPool(nCtrls, intergenicPool, repeatTrees, outFile, "Intergen", summaryFile);
	}
	return 1;
}

void NegativeProbeDesign::chooseRandomProbesFromPool(int nProbesReq, std::vector<PrDes::FeatureStruct>& whichgenPool, Repeats& repeatTrees, std::ofstream &outfile, std::string whichPool, std::ofstream &summaryfile){
	
	bool leftEnd, rightEnd;
	int loopIndex;
	std::mt19937_64 rng;
	std::uniform_int_distribution<int> fragDist(0, whichgenPool.size());
	std::vector<int> fragSeen;
	
	rng.seed(std::random_device()());
	loopIndex=0;
		
	while(loopIndex<nProbesReq){
		//below gives you index of negative control probe to write
		int currIndex=fragDist(rng);
		
		if(std::find(fragSeen.begin(), fragSeen.end(), currIndex)==fragSeen.end()){
		
			fragSeen.push_back(currIndex);
			//check repeat overlap'
			//RE fragment start
		
			if((abs(whichgenPool[currIndex].start - whichgenPool[currIndex].end)) > minREfragLen){
				
				leftEnd = CheckRepeatOverlaps(whichgenPool[currIndex].chr, whichgenPool[currIndex].start, false, repeatTrees);
		
				if(leftEnd){
					rightEnd = CheckRepeatOverlaps(whichgenPool[currIndex].chr, whichgenPool[currIndex].end, true, repeatTrees);
					if(rightEnd){
					//write to file	 - 2 probes for each fragment
						if(chrToIndex.find(whichgenPool[currIndex].chr)!=chrToIndex.end())
							chrToIndex[whichgenPool[currIndex].chr].push_back(negProbe(whichgenPool[currIndex].start, whichgenPool[currIndex].end, "NegCtrl_"+whichPool+std::to_string(loopIndex+1)));
						else{
							std::vector<negProbe> temp;
							temp.push_back(negProbe(whichgenPool[currIndex].start, whichgenPool[currIndex].end, "NegCtrl_"+whichPool+std::to_string(loopIndex+1)));
							chrToIndex.emplace(whichgenPool[currIndex].chr, temp);
						}
						summaryfile << whichgenPool[currIndex].chr << '\t'  << whichgenPool[currIndex].start << '\t' << whichgenPool[currIndex].end << '\t' << "NegCtrl_"<<whichPool<<loopIndex+1 << std::endl;
						++loopIndex;
					}				
				}	
			}
		}
	}
}


bool NegativeProbeDesign::CheckRepeatOverlaps(std::string chr, int& closest_re, bool rightside, Repeats& repeat_trees){
    
	if(!ifRep && !ifMap){
		return true;
	}
	
    int overlaprepeats;
    double mappability;  
    int probeStart, probeEnd;
    
    if (rightside) {
        
        probeStart=( closest_re - ProbeLen + reRightCut ); //add
		probeEnd=( closest_re + reRightCut );    ///////check 1 based 0 based
	
    }
    else{
		probeStart = ( closest_re - reLeftCut); /////check 1 based 0 based
		probeEnd =( closest_re + ProbeLen - reLeftCut);
	
    }
    
    if(ifRep && !ifMap){
		overlaprepeats = repeat_trees.FindOverlaps(chr, probeStart, probeEnd);
		mappability = mapThreshold;
	}
	if(!ifRep && ifMap){
		overlaprepeats = 0;
		mappability = BigWigSummary(chr, probeStart, probeEnd);
	}
	if(ifRep && ifMap){
		overlaprepeats = repeat_trees.FindOverlaps(chr, probeStart, probeEnd);
		mappability = BigWigSummary(chr, probeStart, probeEnd);
	}
	
    		
	if(overlaprepeats > 0 || mappability < mapThreshold ){
		return false;
	}
	else 
		return true;
    
    
}

void NegativeProbeDesign::WritetoFile(bioioMod& getSeq){
	
	std::string target, side;
	int disttotss, probeStart, probeEnd;
    target ="neg_ctrl";
	
	std::ofstream outfile(fileName, std::fstream::app);
	std::ofstream outFasFile(fasFileName, std::fstream::app);
	
	for(auto it=chrToIndex.begin(); it!=chrToIndex.end(); ++it){
		std::sort(it->second.begin(), it->second.end(), [](negProbe const &a, negProbe const &b) { return a.start < b.start; });
		for(auto ind : it->second){
			//left side
			
			probeStart = ( ind.start ); /////check 1 based 0 based
			probeEnd =probeStart+ProbeLen;
			side="L";
			
			outfile  << it -> first  << '\t' << "." << '\t' <<"probe"<<'\t';
    
			outfile <<probeStart << '\t' << probeEnd << '\t'<< "." << '\t'<< "." << '\t'<< "." << '\t'; // to adjust for 1-based coords
			
			outfile << "Name="<< ind.name <<"; " <<"transcriptid="<< "none"<<"; " << "side="<<side<<"; "<<"target="<<target<<"; "<<"design="<<designName<<"; "<< "featuresinvicinity=";
          
			outfile << "none"<<"; ";
	
			outfile<< "targettss="<<"none; " <<"distancetotss="<<"not applicable"<< std::endl;
			
			std::string getFas=getSeq.GetFasta(it->first+":"+std::to_string(probeStart-1)+"-"+std::to_string(probeEnd-1)); // subtract 1 for bioio coordinates on Left side
			
			outFasFile<<it->first<<'\t'<<probeStart<<'\t'<<probeEnd<<'\t'<<"+"<<'\t'<<ind.name<<'\t'<<getFas<<std::endl;
			
			//right side
			probeEnd=( ind.end + reRightCut );
			probeStart=( probeEnd - ProbeLen ); 
			side="R";
			
			outfile  << it -> first  << '\t' << "." << '\t' <<"probe"<<'\t';
    
			outfile <<probeStart << '\t' << probeEnd << '\t'<< "." << '\t'<< "." << '\t'<< "." << '\t'; // to adjust for 1-based coords
        
			outfile << "Name="<< ind.name <<"; " <<"transcriptid="<< "none"<<"; " << "side="<<side<<"; "<<"target="<<target<<"; "<<"design="<<designName<<"; "<< "featuresinvicinity=";
          
			outfile << "none"<<"; " ;
	
			outfile<< "targettss="<<"none; "<<"distancetotss="<<"not applicable"<< std::endl;
			
			getFas=getSeq.GetFasta(it->first+":"+std::to_string(probeStart)+"-"+std::to_string(probeEnd));
			
			outFasFile<<it->first<<'\t'<<probeStart<<'\t'<<probeEnd<<'\t'<<"+"<<'\t'<<ind.name<<'\t'<<getFas<<std::endl;
		}
	}
}
