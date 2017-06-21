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
//	ProxDetectMain.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//


#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cstring>
#include <fstream>
#include <iomanip>

#include "Global.h"
#include "ProcessPairs.h"
#include "Find_Interactions.h"
#include "OutStream.h"
#include "HiCapTools.h"
#include "CallHiCUP.h"


unsigned int totalNumberofPairs=0;
unsigned int NumberofPairs=0;
unsigned int NofPairs_One_on_Probe=0;
unsigned int NofPairs_Both_on_Probe=0;
unsigned int NofPairsNoAnn=0;
const int coreprom_upstream = 1000;
const int coreprom_downstream = 1000;


int HiCapTools::ProxDetectMain(std::string whichchr, std::string statsOption, std::string interactiontype, std::string extraConfig){
	
	
	struct Experiment{
		std::string filepath;
		std::string name;
		std::string designname;
	};
	int NOFEXPERIMENTS; // Number of Experiments
	int padding = 500; //For Sequence Capture Probes
	std::string DigestedGenomeFileName;
	std::string TranscriptListFileName;
	std::string SNPFile;
	std::string negCtrlRegFile;
	std::string ProbeFileName;
	std::string NegCtrlProbeFileName;
	std::string ExpFileName;
	std::string BaseFileName;
    std::vector <Experiment> Experiments; 
	Experiment Exptemp;
	std::map <std::string, std::string> probeType;
	std::locale l;
	std::ofstream statsFile;
	bool emptyErrFlag=false;
	bool generateDigest=false;
	char currTime[100];
	int featFileCount = 3; //4 - transcript and SNV files, 2 - transcript file only, 1 - snp file only
    int countFeatFiles = 2;
    int MinNumberofSupportingPairs = 2; // To be entered by the user
    int ReadLen = 80;//ReadLength
    const int ClusterPromoters  = 1200;
    int MinimumJunctionDistance = 1000; // To be entered by the user
	bool CALCULATE_P_VALUES;
	int WindowSize = 101; 
	int WindowSizeProbeProbe = 3; 
	int BinSize  = 1000; // Only To Calculate Background Interaction Frequencies
	int BinSizeProbeProbe=20000;
	
	PrDes::RENFileInfo reFileInfo;
	
	std::string line, motif, fastaFile;
   

	std::time_t now_time = std::time(NULL);
    std::strftime(reFileInfo.currTime, sizeof(currTime), "%H.%M.%S_%F", std::localtime(&now_time));//Date and time when run starts
    
    std::ofstream logFile(std::string("ProxDetLog_")+reFileInfo.currTime+".log");
    
    OutStream log;
    log.AddStreams(&std::cout, &logFile);
	
	log<<"READ IN INPUTS"<<std::endl;
	
	if(extraConfig.empty()){
		extraConfig="config/configFile.txt";
	}
	
	
	std::ifstream configFile(extraConfig);
	
	if(configFile.good()){
		while (!configFile.eof()){
			getline(configFile, line);
			if(line.find("=")!=std::string::npos){
				if(line.substr(0, line.find('=')).find("Experiment File Name Path")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"!!Error!! : Experiment File Name Path is empty. It is required!" <<std::endl;
						emptyErrFlag = true;
					}
					ExpFileName=s;
				}
				if(line.substr(0, line.find('=')).find("Minimum Number of Supporting Pairs")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						MinNumberofSupportingPairs=std::stoi(line.substr(line.find('=')+1));
					else{
						log<<"##Warning## : Minimum Number of Supporting Pairs is empty. Set to default : "<<MinNumberofSupportingPairs<<std::endl;
					}
				}
				if(line.substr(0, line.find('=')).find("Minimum Junction Distance")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						MinimumJunctionDistance=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## : Minimum Junction Distance is empty. Set to default : "<< MinimumJunctionDistance<<std::endl;
				}
				if(line.substr(0, line.find('=')).find("Padding")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						padding=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## : Padding is empty. Set to default : "<< padding<<std::endl;
				}
				if(line.substr(0, line.find('=')).find("Read Length")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						ReadLen=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## : Read Length is empty. Set to default : "<< ReadLen<<std::endl;
					
				}
				if(line.substr(0, line.find('=')).find("Base File Name")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Warning## : Base File Name is empty. Set to default: ProxDetect!" <<std::endl;
						s = "ProxDetect";
					}
					BaseFileName=s;
				}
				if(line.substr(0, line.find('=')).find("Calculate p_values")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s=="Yes" || s=="yes"|| s=="Y" || s== "y"){
						CALCULATE_P_VALUES=true;
					}
					else{
						CALCULATE_P_VALUES=false;
						if(s.empty())
							log<<"##Note## : Calculate p values is empty. Set to default: No!" <<std::endl;
					}
				}
				if(line.substr(0, line.find('=')).find("Bin Size for Probe-Distal")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						BinSize=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## :Bin Size for Probe-Distal is empty. Set to default : "<< BinSize<<std::endl;
				}
				if(line.substr(0, line.find('=')).find("Window Size for Probe-Distal")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						WindowSize=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## :Window Size for Probe-Distal is empty. Set to default : "<< WindowSize<<std::endl;
				}
				if(line.substr(0, line.find('=')).find("Bin Size for Probe-Probe")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						BinSizeProbeProbe=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## :Bin Size for Probe-Probe is empty. Set to default : "<< BinSizeProbeProbe<<std::endl;
				}
				if(line.substr(0, line.find('=')).find("Window Size for Probe-Probe")!=std::string::npos){
					if(!(line.substr(line.find('=')+1).empty()))
						WindowSizeProbeProbe=std::stoi(line.substr(line.find('=')+1));
					else
						log<<"##Warning## :Window Size for Probe-Probe is empty. Set to default : "<< WindowSizeProbeProbe<<std::endl;
				}
			}
		}
	}
	else{
		log<<"Provided config file not found or configFile.txt not found in config directory."<<std::endl;
	}
	
	//------------------------
	
	if(emptyErrFlag){
		log<< " Enter all required fields."<<std::endl;
		return 0;
	}
	
	//------------------------
	

	if(statsOption=="ComputeStatsOnly"){
		if(CALCULATE_P_VALUES){
			CALCULATE_P_VALUES=false;
			log << "Calculate p values set to No for computing Stats only"<< std::endl;
		}		
	}

	log << std::setw(35)<<std::left<<"Min Number of Supporting Pairs" << MinNumberofSupportingPairs << std::endl;
	log << std::setw(35)<<"Min Junction Distance" << MinimumJunctionDistance << std::endl;
	log << std::setw(35)<<"Padding around Probes" << padding << std::endl;
	log << std::setw(35)<<"Calculate p values" ;
	if(CALCULATE_P_VALUES)
		log<< "Yes" << std::endl;
	if(!CALCULATE_P_VALUES)
		log<< "No" << std::endl;
	
	log << std::setw(35)<<"Read Length"<<	ReadLen<<std::endl;
	log << std::setw(35)<<"Base File name"<<	BaseFileName<<std::endl;
	log << std::setw(35)<<"Bin Size for Probe-Distal"	<<BinSize<<std::endl;
	log << std::setw(35)<<"Window Size for Probe-Distal"	<<WindowSize<<std::endl;
	log << std::setw(35)<<"Bin Size for Probe-Probe"	<<BinSizeProbeProbe<<std::endl;
	log << std::setw(35)<<"Window Size for Probe-Probe"<<WindowSizeProbeProbe<<std::endl;
    
    std::ifstream ExpFile(ExpFileName.c_str());
    if(ExpFile.good()){
		while (!ExpFile.eof()){
			getline(ExpFile, line);
			if(line.find("=")!=std::string::npos){
				if(line.substr(0, line.find('=')).find("Feature Probe File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						log<<"!!Error!! : Feature Probe File Path is empty. It is required!" <<std::endl;
						emptyErrFlag=true;
					}
					ProbeFileName = s;
				}
				if(line.substr(0, line.find('=')).find("Negative control Probe File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty() && CALCULATE_P_VALUES){
						log<<"!!Error!! : Negative Control Probe File Path is empty when Calculate p_values is Yes. It is required!" <<std::endl;
						emptyErrFlag=true;
					}
					NegCtrlProbeFileName = s;
				}
				if(line.substr(0, line.find('=')).find("Digested Genome File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						log<<"##Note## : Digested Genome File is empty. The digest will be generated from the provided fasta file" <<std::endl;
						generateDigest=true;
					}
					DigestedGenomeFileName = s;
				}
				if(line.substr(0, line.find('=')).find("RE cut site motif")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Note## : RE cut site motif is required if genome restriction digest file has to be generated" <<std::endl;
						//emptyErrFlag = true;
					}
					motif=s;
				}
				if(line.substr(0, line.find('=')).find("Genome assembly")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Note## : Genome assembly version is required if genome restriction digest file has to be generated" <<std::endl;
						//emptyErrFlag = true;
					}
					reFileInfo.genomeAssembly=s;
				}
				if(line.substr(0, line.find('=')).find("Fasta File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Note## : Fasta File is required if genome restriction digest file has to be generated" <<std::endl;
						//emptyErrFlag = true;
					}
					fastaFile=s;
				}
				if(line.substr(0, line.find('=')).find("Transcript List File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						log<<"##Warning## : Either Transcript List File or SNV List File is required" <<std::endl;
						featFileCount=featFileCount-2;
						countFeatFiles-=1;
					}
					TranscriptListFileName = s;
				}
				if(line.substr(0, line.find('=')).find("SNV List File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty()){
						log<<"##Warning## : Either Transcript List File or SNV List File is required" <<std::endl;
						featFileCount=featFileCount-1;
						countFeatFiles-=1;
					}
					SNPFile=s;
				}
				if(line.substr(0, line.find('=')).find("Negative control region File")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(begin(s), end(s), [l](char ch) { return std::isspace(ch, l); }), end(s));
					if(s.empty() && CALCULATE_P_VALUES){
						log<<"!!Error!! : Negative Control region File Path is empty when Calculate p_values is Yes. It is required!" <<std::endl;
						emptyErrFlag = true;
					}
					negCtrlRegFile=s;
				}
				if(line.substr(0, line.find('=')).find("Promoters")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						s="promoter";
						log<<"##Warning## : Promoters set to default value : "<< s<<std::endl;
						
					}
					probeType["promoter"]=s;
				}
				if(line.substr(0, line.find('=')).find("SNVs")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						s="SNV";
						log<<"##Warning## : SNVs set to default value : "<< s<<std::endl;
						
					}
					probeType["SNP"]=s;
				}
				if(line.substr(0, line.find('=')).find("Negative controls")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						s="neg_ctrl";
						log<<"##Warning## : Negative Controls set to default value : "<< s<<std::endl;
						
					}
					probeType["negctrl"]=s;
				}
				if(line.substr(0, line.find('=')).find("Other")!=std::string::npos){
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					if(s.empty()){
						s="other";
						log<<"##Warning## : Other set to default value : "<< s<<std::endl;
						
					}
					probeType["other"]=s;
				}
				if(line.substr(0, line.find('=')).find("Number of Experiments")!=std::string::npos){
					NOFEXPERIMENTS=std::stoi(line.substr(line.find('=')+1));
				}
				if(line.substr(0, line.find('=')).find("Experiment BAM File Name Path")!=std::string::npos){
					
					std::string s;
					s=line.substr(line.find('=')+1);
					s.erase(std::remove_if(s.begin(), s.end(), [l](char ch) { return std::isspace(ch, l); }), s.end());
					Exptemp.filepath=s;
					
					getline(ExpFile, line);
					Exptemp.name=line.substr(line.find('=')+1);
					getline(ExpFile, line);
					Exptemp.designname=line.substr(line.find('=')+1);
					Experiments.push_back(Exptemp);
				}
			}
		}
	}
	else{
		log<<"ExperimentFile.txt not found in config directory."<<std::endl;
	}
    
    if(TranscriptListFileName.empty()&&SNPFile.empty()){
		log << "Both transcript list file and SNV list file paths are not entered"<< std::endl;
		emptyErrFlag=true;
	}
	
	if(generateDigest && motif.empty()){
		log << "Restriction enzyme Motif is required if Genome Restriction Digest file is to be generated!"<< std::endl;
		emptyErrFlag=true;
	}
	if(generateDigest && reFileInfo.genomeAssembly.empty()){
		log << "Genome Assembly information is required if Genome Restriction Digest file is to be generated!"<< std::endl;
		emptyErrFlag=true;
	}
	
	if(emptyErrFlag){ // a required field is empty
		log<< " Enter all required fields."<<std::endl;
		return 0;
	}
	
	if(!CheckFile(ExpFileName)){
		log<< "!!Error!! : Experiment File is not accessible"<<std::endl;
		return 0;
	}
	if(!CheckFile(ProbeFileName)){
		log<< "!!Error!! : Feature Probe File is not accessible"<<std::endl;
		return 0;
	}
	if(CALCULATE_P_VALUES && !CheckFile(NegCtrlProbeFileName)){
		log<< "!!Error!! : Negative control Probe File is not accessible"<<std::endl;
		return 0;
	}
	if(!DigestedGenomeFileName.empty() && !CheckFile(DigestedGenomeFileName)){
		log<< "!!Error!! : Genome Digest File is not accessible"<<std::endl;
		return 0;
	}
	if(generateDigest && !CheckFile(fastaFile)){
		log<< "!!Error!! : Fasta File is not accessible"<<std::endl;
		return 0;
	}
	if(!TranscriptListFileName.empty() && !CheckFile(TranscriptListFileName)){
		log<<"!!Error!! : Transcript List File is not accessible " << std::endl;
		return 0;
	}
	if(!SNPFile.empty() && !CheckFile(SNPFile)){
		log<<"!!Error!! : SNV List File is not accessible " << std::endl;
		return 0;
	}
	if(CALCULATE_P_VALUES && !CheckFile(negCtrlRegFile)){
		log<< "!!Error!! : Negative control Regions File is not accessible"<<std::endl;
		return 0;
	}
	
	
	
	if(generateDigest){
		log<< "Generating Restriction Digest File : Starting!"<<std::endl;
		int digestErr;
		CallHiCUP getDigest(log); 
		std::string assemblyVer =reFileInfo.genomeAssembly.substr(0, reFileInfo.genomeAssembly.find_first_of(","));
		 digestErr = getDigest.GenerateRestrictionFile("hiCUPDigester/hicup_digester",motif, assemblyVer, fastaFile, DigestedGenomeFileName);
		 if(digestErr!=1){
			 log<<"!!Digest file could not be generated. Program exiting!!"<<std::endl;
			 return 0;
		 }
		 log<< "Generating Restriction Digest File : Done!"<<std::endl;
	}
	
	
	
	log <<std::setw(35)<< "Feature Probe File"<<ProbeFileName<<std::endl;
	if(CALCULATE_P_VALUES)
		log << std::setw(35)<<"Negative Control Probe File"<<NegCtrlProbeFileName<<std::endl;
	log << std::setw(35)<<"Digested Genome File"<<DigestedGenomeFileName<<std::endl;
	if(!TranscriptListFileName.empty())
		log << std::setw(35)<<"Transcript List File"<<TranscriptListFileName<<std::endl;
	if(!SNPFile.empty())
		log << std::setw(35)<<"SNV List File"<<SNPFile<<std::endl;
	if(CALCULATE_P_VALUES)
		log << std::setw(35)<<"Negative control region File"<<negCtrlRegFile<<std::endl;
	for(auto i=probeType.begin(); i!=probeType.end();++i){
		std::string temp =  "Target-"+i->first +":";
		log << std::setw(35)<<temp<<i->second<<std::endl;
	}
	log << std::setw(35)<<"Number of Experiments"<<NOFEXPERIMENTS<<std::endl;
	for (unsigned i=0; i<Experiments.size(); ++i){ 
		log<<"EXPERIMENT "<<i+1<<std::endl;
		log << std::setw(35)<<"BAM file path"<<Experiments[i].filepath<<std::endl;
		if(!CheckFile(Experiments[i].filepath)){
			log<< "!!Error!! : BAM File "<< i +1 <<" is not accessible"<<std::endl;
			return 0;
		}
			
		log << std::setw(35)<<"Experiment Name"<<Experiments[i].name<<std::endl;
		log << std::setw(35)<<"Probe Design Name"<<Experiments[i].designname<<std::endl;
	}
    
    //------------------------------------------------------------------------------------
    RESitesClass dpnII(log);
    dpnII.InitialiseVars(DigestedGenomeFileName);
    //-------------------------------------------------------------------------------------
    
    ProcessBAM bamfile(log);
    //-------------------------------------------------------------------------------------
    DetectInteractions Interactions(log, MinNumberofSupportingPairs, CALCULATE_P_VALUES, MinimumJunctionDistance);
   //-------------------------------------------------------------------------------------
    
    std::string BAMFILENAME, ExperimentName;
    std::vector < std::string > ExperimentNames;
    int ExperimentNo = 0;
    
    FeatureClass proms(log);
   
    //-------------------------------------------------------------------------------
	log << "Reading Feature files and annotating features: Starting!" << std::endl;
    proms.InitialiseData(ClusterPromoters, countFeatFiles+CALCULATE_P_VALUES);
    
    switch(featFileCount){
		case 1:	
			proms.ReadFeatureAnnotation(dpnII, SNPFile, "SNV");
			break;
		case 2:
			proms.ReadFeatureAnnotation(dpnII, TranscriptListFileName, "transcript");	
			break;
		case 3:
			proms.ReadFeatureAnnotation(dpnII, TranscriptListFileName, "transcript");
			proms.ReadFeatureAnnotation(dpnII, SNPFile, "SNV");
			break;
		default:
			break;
			
	}
	if(CALCULATE_P_VALUES)
		proms.ReadFeatureAnnotation(dpnII, negCtrlRegFile, "neg_ctrl");
		
	log << "Reading Feature files and annotating features: Done!" << std::endl;

	//---------------------------------------------------------------------------------
	
	ProbeSet ProbeClass(log, (1+CALCULATE_P_VALUES), 0); // 1 for feature probe file, 2 if Neg ctrl probe file.
    ProbeClass.ReadProbeCoordinates(ProbeFileName, probeType, padding, false, reFileInfo);
    if(CALCULATE_P_VALUES)	
		ProbeClass.ReadProbeCoordinates(NegCtrlProbeFileName, probeType, padding, true, reFileInfo);
    
    //---------------------------------------------------------------------------------
    if(statsOption=="ComputeStatsOnly"){
		statsFile.open(BaseFileName+"."+reFileInfo.genomeAssembly.substr(0, reFileInfo.genomeAssembly.find_first_of(',')) +".computeStats."+reFileInfo.currTime+".txt");
		statsFile << "Total_Number_of_Pairs"<<'\t'<< "Total_Number_of_Pairs_on_Probes"<<'\t' << "Number_of_Pairs_Both_Reads_on_Probe" << '\t'<< "Number_of_Pairs_One_Read_on_Probe" << '\t' << "Number_of_Pairs_None_on_Probe" << '\t'<< "On_Probe_Pair_Fraction" <<std::endl;
	}
    
    //---------------------------------------------------------------------------------
    std::vector < DetermineBackgroundLevels > background;
    
    ProximityClass proximities(NOFEXPERIMENTS);
       
		
	for (unsigned i=0; i<Experiments.size(); ++i){ // Reads all the pairs in each experiment and fills the interaction maps

		Exptemp=Experiments[i];
		bamfile.Initialize(Exptemp.filepath, NOFEXPERIMENTS, padding, ReadLen);
		ExperimentNames.push_back(Exptemp.name);
		
		
		
		if(CALCULATE_P_VALUES){ //Fill NegCtrl proximities to calculate background interaction frequencies
			bamfile.ProcessSortedBamFile_NegCtrls(ProbeClass, dpnII, proximities, Exptemp.filepath, ExperimentNo, Exptemp.designname, statsOption);
		
			background.push_back(DetermineBackgroundLevels());
		
			background.back().CalculateMeanandStdRegress(Exptemp.name+".Probe_Distal", ExperimentNo, Exptemp.designname, background.back().bglevels, BinSize, "ProbeDistal", MinimumJunctionDistance, log, WindowSize);
			background.back().CalculateMeanandStdRegress(Exptemp.name+".Probe_Probe", ExperimentNo, Exptemp.designname, background.back().bglevelsProbeProbe, BinSizeProbeProbe, "ProbeProbe", MinimumJunctionDistance, log, WindowSizeProbeProbe);
			
			log << "Total_Number_of_Pairs" << '\t' << totalNumberofPairs/2 << std::endl;
			log << "Total_Number_of_Pairs on Probes" << '\t' << NumberofPairs << std::endl;
			log << "Number_of_Pairs_Both_Reads_on_Probe" << '\t' << NofPairs_Both_on_Probe << std::endl;
			log << "Number_of_Pairs_One_Read_on_Probe" << '\t' << NofPairs_One_on_Probe << std::endl;
			log << "Number_of_Pairs_None_on_Probe" << '\t' << NofPairsNoAnn << std::endl;
			log << "FractionofPairsOnProbe" << '\t' << (NumberofPairs)/double(totalNumberofPairs) << std::endl; 	    
			NumberofPairs = 0; NofPairs_Both_on_Probe = 0; NofPairs_One_on_Probe = 0; NofPairsNoAnn = 0;
		}
       	bamfile.ProcessSortedBAMFile(ProbeClass, dpnII, proximities, Exptemp.filepath, ExperimentNo, whichchr, Exptemp.designname, statsOption);
         
		log << "Total_Number_of_Pairs" << '\t' << totalNumberofPairs/2 << std::endl;
		log << "Total_Number_of_Pairs on Probes" << '\t' << NumberofPairs << std::endl;
		log << "Number_of_Pairs_Both_Reads_on_Probe" << '\t' << NofPairs_Both_on_Probe << std::endl;
		log << "Number_of_Pairs_One_Read_on_Probe" << '\t' << NofPairs_One_on_Probe << std::endl;
		log << "Number_of_Pairs_None_on_Probe" << '\t' << NofPairsNoAnn << std::endl;
		log << "On_Probe_Pair_Fraction" << '\t' << (NumberofPairs)/double(totalNumberofPairs) << std::endl; 	    
			
		if(statsOption=="ComputeStatsOnly"){
			statsFile <<  Exptemp.name << std::endl;
			statsFile << totalNumberofPairs/2 << '\t' <<NumberofPairs << '\t' << NofPairs_Both_on_Probe<< '\t' << NofPairs_One_on_Probe<< '\t' << NofPairsNoAnn << '\t' << (NumberofPairs)/double(totalNumberofPairs) << std::endl;
		}
		
        totalNumberofPairs = 0;
		NumberofPairs = 0; 
		NofPairs_Both_on_Probe = 0; 
		NofPairs_One_on_Probe = 0; 
		NofPairsNoAnn = 0;       
		++ExperimentNo;

		log << Exptemp.filepath << "     finished" << std::endl;
	}
    
	if(interactiontype=="Neg"){
			
		Interactions.CalculatePvalAndPrintInteractionsProbeDistal_NegCtrls(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSize, reFileInfo); //Print all same type of interactions
		Interactions.CalculatePvalAndPrintInteractionsProbeProbe_NegCtrls(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSizeProbeProbe, reFileInfo); //Print all same type of interactions
	}
	else if(interactiontype=="NonNeg"){
			
		Interactions.CalculatePvalAndPrintInteractionsProbeDistal(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSize, reFileInfo);//Print all same type of interactions
		Interactions.CalculatePvalAndPrintInteractionsProbeProbe(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSizeProbeProbe, reFileInfo);//Print all same type of interactions
	}
	else{	
	
		Interactions.CalculatePvalAndPrintInteractionsProbeDistal(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSize, reFileInfo);//Print all same type of interactions
		Interactions.CalculatePvalAndPrintInteractionsProbeDistal_NegCtrls(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSize, reFileInfo); //Print all same type of interactions
		Interactions.CalculatePvalAndPrintInteractionsProbeProbe(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSizeProbeProbe, reFileInfo);//Print all same type of interactions
		Interactions.CalculatePvalAndPrintInteractionsProbeProbe_NegCtrls(ProbeClass, background, BaseFileName, NOFEXPERIMENTS, ExperimentNames, whichchr, BinSizeProbeProbe, reFileInfo); //Print all same type of interactions
	}

}
