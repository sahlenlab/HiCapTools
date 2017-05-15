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
//  ProbeDesignMain.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <ctime>

#include "HiCapTools.h"
#include "NegativeProbeDesign.h"
#include "PrintUsage.h"
#include "OutStream.h"
#include "bioioMod.h"

int DistanceBetweenProbes = 1000;

int HiCapTools::ProbeDesignMain(std::string whichchr, std::string extraConfig) {
	
	int ClusterPromoters  = 1200;
	int ProbeLen = 120;
	int MaxDistancetoTSS = 2500;
	int MinNegFragLen=500;
	std::string DigestedGenomeFileName;
	std::string mappabilityfile;
	std::string bigwigsummarybinary;
	std::string transcriptfile;
	std::string SNPfile;
	std::string repeatfile;
	std::string config_file_name;
    std::string line;
    std::string motif;
    std::locale l;
    int dforbidIntergen=50000;
    int dforbidProm=50000;
    int dforbidRegReg=50000;
    int exonNegCtrls, intronNegCtrls, intergenNegCtrls;
    std::string ifNeg;
    bool ifRegRegion=false;
    std::string regRegionFile;
    PrDes::RENFileInfo reFileInfo;
    bool emptyErrFlag = false;
    int featFileCount = 3; //3 - both files, 2 - transcript file only, 1 - snp file only
    int countFeatFiles = 2;
    int repeatOverlapExtent = 6;
    int BUFSIZE = 128;
    std::string fastaIndexFile;
    std::string fastaFile;
    
    
    std::time_t now_time = std::time(NULL);
    std::strftime(reFileInfo.currTime, sizeof(reFileInfo.currTime), "%H.%M.%S_%F", std::localtime(&now_time));//Date and time when run starts
	
	//Initialise
    reFileInfo.mappabilityThreshold = 0.7;
    
	//-----Output to Log and std::cout-----//
    
    std::ofstream logFile(std::string("ProbeDesignLog_")+reFileInfo.currTime+".log");
    
    OutStream log;
    log.AddStreams(&std::cout, &logFile);
    //-----------------------------------//
	
	
    
    
    log<<"READ IN INPUTS"<<std::endl;
    
    if(extraConfig.empty()){
		extraConfig="config/probeDesignConfig.txt";
	}
	
	std::ifstream config_file(extraConfig);
    
    if(config_file.good()){
		while (!config_file.eof()){
			getline(config_file, line);
			if(line.find("=")!=std::string::npos){
				if(line.substr(0, line.find('=')).find("Base File Name")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : Base File Name is empty. Will use default value: ProbesDesigned" <<std::endl;
						s = "ProbesDesigned";
						//emptyErrFlag = true;
					}		
					reFileInfo.desName=s;
				}
				if(line.substr(0, line.find('=')).find("Digested Genome File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : Digested Genome File is required" <<std::endl;
						emptyErrFlag = true;
					}
					DigestedGenomeFileName=s;
				}
				if(line.substr(0, line.find('=')).find("RE cut site motif")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : RE cut site motif is required" <<std::endl;
						emptyErrFlag = true;
					}
					motif=s;
				}
				if(line.substr(0, line.find('=')).find("Genome assembly")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : Genome assembly version is required" <<std::endl;
						emptyErrFlag = true;
					}
					reFileInfo.genomeAssembly=s;
				}
				if(line.substr(0, line.find('=')).find("Transcript List File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Warning## : Either Transcript List File or SNV List File is required" <<std::endl;
						featFileCount=featFileCount-2;
						countFeatFiles-=1;
						//emptyErrFlag = true;
					}
					transcriptfile=s;
				}
				if(line.substr(0, line.find('=')).find("SNV List File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Warning## : Either Transcript List File or SNV List File is required" <<std::endl;
						featFileCount=featFileCount-1;
						countFeatFiles-=1;
						//emptyErrFlag = true;
					}
					SNPfile=s;
				}
				if(line.substr(0, line.find('=')).find("Repeat File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Note## : Repeat File path is empty. Repeat Overlaps will not be checked" <<std::endl;
						reFileInfo.ifRepeatAvail=false;
					}
					else
						reFileInfo.ifRepeatAvail=true;
					repeatfile=s;
				}
				if(line.substr(0, line.find('=')).find("Mappability File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Note## : Mappability File path is empty. Mappability will not be checked" <<std::endl; 
						reFileInfo.ifMapAvail=false;
					}
					else
						reFileInfo.ifMapAvail=true;
					mappabilityfile=s;
				}
				if(line.substr(0, line.find('=')).find("bigWigSummary executable path")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : bigWigSummary executable path is empty and is required" <<std::endl;
						emptyErrFlag = true;
					}
					bigwigsummarybinary=s;
				}
				if(line.substr(0, line.find('=')).find("Probe Length")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						ProbeLen=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Minimum distance between Probes")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						DistanceBetweenProbes=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Maximum distance from Probe to feature start (TSS if the feature is transcript)")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						MaxDistancetoTSS=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Cluster Promoters")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						ClusterPromoters=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Extent of Repeat Overlaps")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						repeatOverlapExtent=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Mappability Threshold")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						reFileInfo.mappabilityThreshold=std::stod(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Fasta File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : Fasta File is required" <<std::endl;
						emptyErrFlag = true;
					}
					fastaFile=s;
				}
				if(line.substr(0, line.find('=')).find("Fasta Index File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Note!! : Fasta Index File will be checked in the same directory as Fasta File" <<std::endl;
					}
					fastaIndexFile=s;
				}
//--------------------------------------Negative Control Probes ---------------------------------------------------//
			
				if(line.substr(0, line.find('=')).find("Design negative probe set")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Warning## : Design negative probe set is empty and therefore set to No" <<std::endl;
						//emptyErrFlag = true;
						ifNeg="No";
					}
					else
						ifNeg=s;
				}
				if(line.substr(0, line.find('=')).find("Minimum fragment length for negative probe")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						MinNegFragLen=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Number of Intergenic Negative control fragments")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						intergenNegCtrls=std::stoi(line.substr(line.find('=')+1));
					else
						intergenNegCtrls=0;
				}
				if(line.substr(0, line.find('=')).find("Number of Exonic Negative control fragments")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						exonNegCtrls=std::stoi(line.substr(line.find('=')+1));
					else
						exonNegCtrls=0;
				}
				if(line.substr(0, line.find('=')).find("Number of Intronic Negative control fragments")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						intronNegCtrls=std::stoi(line.substr(line.find('=')+1));
					else
						intronNegCtrls=0;
				}
				if(line.substr(0, line.find('=')).find("Minimum distance to a known promoter for negative control probes")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						dforbidProm=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Minimum distance to a known gene for intergenic negative control probes")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						dforbidIntergen=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Use user provided forbidden regions?")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s=="Yes")
						ifRegRegion=true;
				}
				if(line.substr(0, line.find('=')).find("User provided forbidden regions File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(ifRegRegion && s.empty()){
						log<<"!!Error!! : Forbidden regions file is required when Use user provided regions to avoid? is set to Yes" <<std::endl;
						emptyErrFlag = true;
					}
					regRegionFile=s;
				}
				if(line.substr(0, line.find('=')).find("Minimum distance to any user provided forbidden region")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						dforbidRegReg=std::stoi(line.substr(line.find('=')+1));
				}
				
			}
		}
	}
	else{
		log<<"Provided config file not found or probeDesignConfig.txt not found in config directory."<<std::endl;
	}
	
	if(transcriptfile.empty()&&SNPfile.empty()){
		log << "Both transcript list file and SNV list file paths are not entered"<< std::endl;
		emptyErrFlag=true;
	}
	
	if(emptyErrFlag){ // a required field is empty
		return 0;
	}
	
	if(!CheckFile(DigestedGenomeFileName)){
		log<<"!!Error!! : Digested Genome File is not accessible " << std::endl;
		return 0;
	}
	if(!transcriptfile.empty() && !CheckFile(transcriptfile)){
		log<<"!!Error!! : Transcript List File is not accessible " << std::endl;
		return 0;
	}
	if(!SNPfile.empty() && !CheckFile(SNPfile)){
		log<<"!!Error!! : SNV List File is not accessible " << std::endl;
		return 0;
	}
	if(reFileInfo.ifRepeatAvail && !CheckFile(repeatfile)){
		log<<"!!Error!! : Repeats File is not accessible " << std::endl;
		return 0;
	}
	if(reFileInfo.ifMapAvail && !CheckFile(mappabilityfile)){
		log<<"!!Error!! : Mappability File is not accessible " << std::endl;
		return 0;
	}
	if(!CheckFile(fastaFile)){
		log<<"!!Error!! : Fasta File is not accessible " << std::endl;
		return 0;
	}
	if(ifRegRegion && !CheckFile(regRegionFile)){
		log<<"!!Error!! : User provided forbidden regions file is not accessible " << std::endl;
		return 0;
	}
	
            
    log << std::setw(75)<<std::left<<"Base File Name:" << reFileInfo.desName << std::endl;    
    log << std::setw(75)<<"Digested Genome File:" << DigestedGenomeFileName << std::endl;
    log << std::setw(75)<<"RE cut site motif:" << motif << std::endl;
    if(!transcriptfile.empty())
		log <<std::setw(75)<< "Transcript List File:" << transcriptfile << std::endl;
	if(!SNPfile.empty())
		log << std::setw(75)<<"SNV List File:" << SNPfile << std::endl;
	if(reFileInfo.ifRepeatAvail)	
		log << std::setw(75)<<"Repeat File:" << repeatfile << std::endl;
	if(reFileInfo.ifRepeatAvail)
		log << std::setw(75)<<"Mappability File:" << mappabilityfile << std::endl;
    log << std::setw(75)<<"bigWigSummary executable path:" << bigwigsummarybinary << std::endl;
	log << std::setw(75)<<"Probe Length:" << ProbeLen << std::endl;
	log << std::setw(75)<<"Minimum distance between Probes:" << DistanceBetweenProbes << std::endl;
	log << std::setw(75)<<"Maximum distance from Probe to TSS:"<<MaxDistancetoTSS << std::endl;
	log << std::setw(75)<<"Cluster Promoters:"<< ClusterPromoters << std::endl;
	log << std::setw(75)<<"Extent of Repeat Overlaps:"<< repeatOverlapExtent << std::endl;
	log << std::setw(75)<<"Mappability Threshold:"<< reFileInfo.mappabilityThreshold << std::endl;
	log << std::setw(75)<<"Fasta File:"<< fastaFile << std::endl;

	if(ifNeg=="Yes"){
		log << std::setw(75)<<"Design negative probe set:"<< ifNeg << std::endl;
		log << std::setw(75)<<"Minimum fragment length for negative probe:"<< MinNegFragLen << std::endl;
		log << std::setw(75)<<"Number of Intergenic Negative control probes:"<< intergenNegCtrls << std::endl;
		log << std::setw(75)<<"Number of Intronic Negative control probes:" << intronNegCtrls << std::endl;
		log << std::setw(75)<<"Number of Exonic Negative control probes:" <<exonNegCtrls << std::endl;
		log << std::setw(75)<<"Minimum distance to a known promoter for negative control probes:" << dforbidProm << std::endl;
		log << std::setw(75)<<"Minimum distance to a known gene for intergenic negative control probes:"<< dforbidIntergen << std::endl;
		if(ifRegRegion){
			log << std::setw(75)<<"Use user provided forbidden regions?:" << "Yes" << std::endl;
			log << std::setw(75)<<"User provided forbidden regions File:" <<regRegionFile << std::endl;
			log << std::setw(75)<<"Minimum distance to a user provided forbidden regions:" << dforbidRegReg << std::endl; 
			
		}
		else
			log << std::setw(75)<<"Use user provided forbidden regions?:" << "No" << std::endl;		
	}
    
    
    reFileInfo.leftOfCut = motif.substr(0, motif.find_first_of(",")).find_first_of("^");
    reFileInfo.rightOfCut = (motif.substr(0, motif.find_first_of(",")).size() -1) - reFileInfo.leftOfCut;
	reFileInfo.REName = motif.substr(motif.find_first_of(",")+1);
    
    log << "Reading Digest file and finding RE sites: Starting!" << std::endl;
    RESitesClass dpnIIsites(log);
    dpnIIsites.InitialiseVars(DigestedGenomeFileName);
    log << "Reading Digest file and finding RE sites: Done!" << std::endl;
    
    Repeats hg19repeats;
    if(reFileInfo.ifRepeatAvail){
		log << "Reading repeat file and calculating repeat overlap: Starting!" << std::endl;
		hg19repeats.ReadRepeatIntervals(repeatfile, log, repeatOverlapExtent);
		log << "Reading repeat file and calculating repeat overlap: Done!" << std::endl;
    }
//-------------------//------------------------------
    log << "Reading Feature files and annotating features: Starting!" << std::endl;
    ProbeFeatureClass Features(log);
    Features.InitialiseData(ClusterPromoters, countFeatFiles, dforbidProm); //ReadFeatureAnnotation called with 1 or 2 files
   
    switch(featFileCount){
		case 1:	
			Features.ReadFeatureAnnotation(dpnIIsites, SNPfile, "SNV");
			break;
		case 2:
			Features.ReadFeatureAnnotation(dpnIIsites, transcriptfile, "transcript");	
			break;
		case 3:
			Features.ReadFeatureAnnotation(dpnIIsites, transcriptfile, "transcript");
			Features.ReadFeatureAnnotation(dpnIIsites, SNPfile, "SNV");
			break;
		default:
			break;
			
	}
	log << "Reading Feature files and annotating features: Done!" << std::endl;
	
//-------------------//------------------------------    
    log << "Designing Probes: Starting!" << std::endl;
    DesignClass designProbes(log);
    bioioMod getSeq(log, fastaIndexFile, fastaFile);
    
    if(whichchr=="chrAll"){
		for(auto &iChr: Features.ChrNames_proms){
				designProbes.DesignProbes(Features, dpnIIsites, hg19repeats, bigwigsummarybinary, mappabilityfile, iChr, MaxDistancetoTSS, ProbeLen, reFileInfo, BUFSIZE) ;
		}
		designProbes.MergeAllChrOutputs(Features, reFileInfo);
	}
	else
		designProbes.DesignProbes(Features, dpnIIsites, hg19repeats, bigwigsummarybinary, mappabilityfile, whichchr, MaxDistancetoTSS, ProbeLen, reFileInfo, BUFSIZE) ;
	
	bool isFas = designProbes.ConstructSeq(reFileInfo, getSeq);
	
	if(!isFas){
		log << "!!Error getting sequence from fasta file!! Aborting!" << std::endl;
		return 0;
	}
	
	log << "Designing Probes: Done!" << std::endl;
    
//-------------------//------------------------------    
    
    if(ifNeg=="Yes"){
		
		log << "Designing Negative control Probes: Starting!" << std::endl;
		NegativeProbeDesign NegctrlProbes(log);
		int res;
		NegctrlProbes.InitialiseDesign(Features, transcriptfile, regRegionFile, ifRegRegion, ProbeLen, bigwigsummarybinary, mappabilityfile, MinNegFragLen, reFileInfo, BUFSIZE, DigestedGenomeFileName, dforbidIntergen, dforbidRegReg);
							
		res = NegctrlProbes.ConstructNegativeControlProbes(exonNegCtrls, "exonic",  hg19repeats);
		res = NegctrlProbes.ConstructNegativeControlProbes(intronNegCtrls,"intronic",  hg19repeats);
		res = NegctrlProbes.ConstructNegativeControlProbes(intergenNegCtrls, "intergenic", hg19repeats);
		NegctrlProbes.WritetoFile(getSeq);
		
		log << "Designing Negative control Probes: Done!" << std::endl;
	}
	
	log<<"Execution Complete...................."<<std::endl;
	
    return 1;
}

bool HiCapTools::CheckFile(std::string filename){
	std::ifstream in(filename);
	if(in.good())
		return true;
	else
		return false;
}
