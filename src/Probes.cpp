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
//  Probes.cpp
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#include "Probes.h"
#include "Global.h"
#include <fstream>


std::map < std::string, Probe_Design > Design;
std::map < std::string, Probe_Design > Design_NegCtrl;


void ProbeSet::GetProbeFeats(std::stringstream& line, CaptureProbes& t, std::string& Name){
    
    std::string promoter, snp, negctrl, other, SOterm;
    
    std::string field, field2, temp, attributes, start, end;
    
    std::locale l;


    // gff3 field 1 chromosome
    getline(line,t.chr,'\t'); 
    
    // gff3 field 2 source - discarded
    getline(line,temp,'\t'); 
    
    // gff3 field 3 feature type - SO term
    
    getline(line,SOterm,'\t'); // SO term
    
    // gff3 field 4 start     
    getline(line,field,'\t');
    t.start = std::stoi(field);
    
    // gff3 field 5 end     
    getline(line,end,'\t');
    t.end = std::stoi(end);
    
    if(SOterm.find("probe")==std::string::npos && SOterm.find("SO:0000051")==std::string::npos){
		prLog<<"The probe with chr:"<<t.chr<<" Start:"<<t.start<< " End:"<<t.end<< " is not identified by SO term for 'probe'"<< std::endl; 
	}
     
    // gff3 field 6 score - discarded    
    getline(line,field,'\t');
    
    // gff3 field 7 strand - discarded
    getline(line,temp,'\t');
    // gff3 field 8 phase - discarded
    getline(line,temp,'\t');
    // gff3 field 9 attributes
    getline(line,attributes,'\t'); // semicolon separated tags for upstream/downstream
    
    if(attributes.back()!=';')
		attributes.push_back(';');
		
    size_t countSC=std::count(attributes.begin(), attributes.end(), ';');
    std::string parseValue;
    size_t posEnd;
    while(countSC>0){
		
		posEnd=attributes.find_first_of(';'); 
		
		parseValue=attributes.substr(0, posEnd);
		parseValue.erase(std::remove_if(parseValue.begin(), parseValue.end(), [l](char ch) { return std::isspace(ch, l); }), parseValue.end());
		
		if(parseValue.substr(0, parseValue.find('=')).find("Name")!= std::string::npos){
			Name=parseValue.substr(parseValue.find('=')+1);
		}
		
		if(parseValue.substr(0, parseValue.find('=')).find("transcriptid")!= std::string::npos){
			t.target_id=parseValue.substr(parseValue.find('=')+1); //target transcript 
		}
		
		if(parseValue.substr(0, parseValue.find('=')).find("side")!= std::string::npos){
			t.side=parseValue.substr(parseValue.find('=')+1); // Upstream (Left) or the downstream (Right) of the feature
		}
		
		if(parseValue.substr(0, parseValue.find('='))== "target"){
			std::string target = parseValue.substr(parseValue.find('=')+1); // Probe target for annotation 
			if(probeTypeMap.at("promoter").find(target)!=std::string::npos){
				t.annotated = 1;
			}
			else if(probeTypeMap.at("SNP").find(target)!=std::string::npos){
				t.annotated = 2;
			}
			else if(probeTypeMap.at("negctrl").find(target)!=std::string::npos){
				t.annotated = 3;
			}
			else if(probeTypeMap.at("other").find(target)!=std::string::npos){
				t.annotated = 4;
			}
			else{
				t.annotated = 4;
			}
		}
		
		if(parseValue.substr(0, parseValue.find('=')).find("design")!= std::string::npos){
			t.name_of_design=parseValue.substr(parseValue.find('=')+1);// Design Name 
		}
		
		attributes.erase(0, posEnd+1);
		--countSC;
	
	}
}

void ProbeSet::ProcessProbeLine(std::map< std::string, std::vector < std::string > >::iterator nameit, std::vector < int>& coords, CaptureProbes tempprobe, int& closest,  std::map<std::string, std::vector< Interval < int > > > & negctrl_probes_interval, std::map<std::string, std::vector< Interval < int > > >& probes_interval, int padding){
    
    for(int i = 0; i < nameit->second.size(); ++i){
        coords.push_back(Features[nameit->second[i]].start);
    }
    if(coords.size() == 1){
        tempprobe.feature_id = nameit->second[0];
        Features[nameit->second[0]].FeatureType = tempprobe.annotated;
    }
    else{
        closest = FindClosestFeature(tempprobe.start, coords,tempprobe.side);
        tempprobe.feature_id = nameit->second[closest];
        Features[nameit->second[closest]].FeatureType = tempprobe.annotated;
    }
    if(tempprobe.annotated == 3){
        Design_NegCtrl[tempprobe.name_of_design].Probes.push_back(tempprobe);
        //left and right
        if(tempprobe.side=="L"){
			negctrl_probes_interval[tempprobe.name_of_design].push_back(Interval<int>((tempprobe.start),(tempprobe.end +padding),(Design_NegCtrl[tempprobe.name_of_design].Probes.size()-1)));
		}
		else if(tempprobe.side=="R"){
			negctrl_probes_interval[tempprobe.name_of_design].push_back(Interval<int>((tempprobe.start -padding),(tempprobe.end),(Design_NegCtrl[tempprobe.name_of_design].Probes.size()-1)));
		}
    }
    else{
        Design[tempprobe.name_of_design].Probes.push_back(tempprobe);
        //left and right
        if(tempprobe.side=="L"){
			probes_interval[tempprobe.name_of_design].push_back(Interval<int>((tempprobe.start),(tempprobe.end + padding),(Design[tempprobe.name_of_design].Probes.size()-1)));
		}
		else if(tempprobe.side=="R"){
			probes_interval[tempprobe.name_of_design].push_back(Interval<int>((tempprobe.start -padding),(tempprobe.end),(Design[tempprobe.name_of_design].Probes.size()-1)));
		}
    }
}

void ProbeSet::ReadProbeCoordinates(std::string ProbeFileName, std::map <std::string, std::string>& probeType, int padding, bool isNeg, PrDes::RENFileInfo& reInfo){
	
	probeTypeMap=probeType;
    std::string filename,designname;
    bool found=true;
	std::locale l;

	int index = 0, i, closest;
	std::string pline, pline0, chr1, chr2, temp, Name;
    bool flag=true;
    
   
    CaptureProbes tempprobe;
    
    std::map< std::string, std::vector < std::string > >::iterator nameit; // key = name (A1BG, rs1101 etc); value = list of feature_id (e.g. chr1_199202)
    nameit = MetaFeatures.begin();
    std::vector <int> coords;

    std::map< std::string, FeatureStruct >::iterator featiter;
    
    
    filename.append(ProbeFileName);
	std::ifstream probefile(filename.c_str());
    
    //Discard headers///////
    while(flag == true){
		getline(probefile, temp);
		if(temp[0] != '#'){		// Read the first probe
			flag=false;
			pline0=temp;
		}
		else if(temp.find("genome-build ")!=std::string::npos){
			std::string d=temp.substr(temp.find("genome-build ")+13);
			reInfo.genomeAssembly=d.substr(d.find(" ")+1)+","+d.substr(0, d.find(" ")+1);
		}
	}
	
	std::stringstream probeline ( pline0 );
	GetProbeFeats(probeline, tempprobe, Name);
	chr1 = tempprobe.chr; //take the first chr outside the loop
    
    while (getline(probefile, pline)) {
        std::map< std::string, std::vector< Interval < int >  > > negctrl_probes_interval;
        std::map< std::string, std::vector< Interval < int >  > > probes_interval;
        std::map< std::string, std::vector< Interval < int >  > >::iterator it;
        
        while(index == 0){
            do{
                nameit = MetaFeatures.find(Name);
                if(nameit != MetaFeatures.end()){
                    ProcessProbeLine(nameit, coords, tempprobe, closest, negctrl_probes_interval, probes_interval, padding);
                    tempprobe.feature_id = "";
                    ++index;
                    coords.clear();
                    found = true;
                }
                else{ //Could not find the first probe in transcript file, try again
					if(pline.empty()){
						found=true;
						break;
					}
                    std::stringstream probeline ( pline );
                    GetProbeFeats(probeline, tempprobe, Name);
                    chr1 = tempprobe.chr; //take the first chr outside the loop
                    if (!getline(probefile, pline))
						break;
                    found = false;
                }
            }while(!found);
        }
        std::stringstream probeline ( pline );
		GetProbeFeats(probeline, tempprobe, Name);
		chr2 = tempprobe.chr;
        do{
            nameit = MetaFeatures.find(Name);
            if(nameit != MetaFeatures.end()){
                ProcessProbeLine(nameit, coords, tempprobe, closest, negctrl_probes_interval, probes_interval, padding);
                tempprobe.feature_id = "";
                ++index;
                coords.clear();
                found = true;
            }
            else{
                std::stringstream probeline ( pline );
                GetProbeFeats(probeline, tempprobe, Name);
                chr2 = tempprobe.chr;
                if(!getline(probefile, pline))
					break;
                found = false;
            }
        }while(!found);
        
		while(getline(probefile, pline)){
            if(chr1 == chr2){
                std::stringstream probeline ( pline );
                GetProbeFeats(probeline, tempprobe, Name);
                chr2 = tempprobe.chr;
                nameit = MetaFeatures.find(Name);
                if(nameit != MetaFeatures.end()){
                    ProcessProbeLine(nameit, coords, tempprobe, closest, negctrl_probes_interval, probes_interval, padding);
                    tempprobe.feature_id = "";
                    ++index;
                    coords.clear();
                    found = true;
                }
                else{
                    found = false;
                }
            }
            else{
                std::stringstream probeline ( pline );
                GetProbeFeats(probeline, tempprobe, Name);
                chr2 = tempprobe.chr;
                nameit = MetaFeatures.find(Name);
                if(nameit != MetaFeatures.end()){
                    ProcessProbeLine(nameit, coords, tempprobe, closest, negctrl_probes_interval, probes_interval, padding);
                    tempprobe.feature_id = "";
                    ++index;
                    coords.clear();
                    found = true;
                }
                else{
                    found = false;
                }
                break; // one chromosome finished, get out of the loop
            }
        }
        if(isNeg){
			for(auto it = negctrl_probes_interval.begin(); it != negctrl_probes_interval.end();++it){
				designname = it->first;
				AddTotheIntervalTree(probes_interval[designname], negctrl_probes_interval[designname], chr1, designname, isNeg);
			}
		}
		else{
			for(auto it = probes_interval.begin(); it != probes_interval.end();++it){
				designname = it->first;
				AddTotheIntervalTree(probes_interval[designname], negctrl_probes_interval[designname], chr1, designname, isNeg);
			}
		}
        
	
        chr1 = tempprobe.chr;
    }
	prLog << index <<  "   Probe Coordinates Read from "<< ProbeFileName<< std::endl;
    probefile.close();
    
    ++filesReadCount; 
}


int ProbeSet::FindClosestFeature(int probe_coord, std::vector<int>& isoformcoords, std::string side){
    
    int smallest_distance = 0, dist = 0, index = 0;
    if(side == "L"){ // Left probe, upstream
        smallest_distance = -10E6;
        for (int i = 0; i < isoformcoords.size(); ++i) {
            if((probe_coord - isoformcoords[i]) < 0){
                dist = (probe_coord - isoformcoords[i]);
                if (dist >= smallest_distance) {
                    smallest_distance = dist;
                    index = i;
                }
            }
        }
    }
    else{ //Right probe, downstream
        smallest_distance = 10E6;
        for (int i = 0; i < isoformcoords.size(); ++i) {
            if(probe_coord - isoformcoords[i] > 0){
                dist = abs(probe_coord - isoformcoords[i]);
                if (smallest_distance < dist) {
                    smallest_distance = dist;
                    index = i;
                }
            }
        }
    }
    return index;
}

int ProbeSet::AddTotheIntervalTree(std::vector<Interval< int > >& intervals, std::vector<Interval< int > >& intervals_negctrls, std::string chr_of_vector, std::string nameofdesign, bool isNeg){
    
    std::vector<Interval< int >  > temp, temp2;
    
    if(isNeg){
			Design_NegCtrl[nameofdesign].Probe_Tree[chr_of_vector] = IntervalTree< int >(intervals_negctrls);
			intervals_negctrls.swap(temp2); //empty the tree
	}
	else{
			Design[nameofdesign].Probe_Tree[chr_of_vector] = IntervalTree< int >(intervals);
			intervals.swap(temp); //empty the tree
	}  
    
    return 1;
}



int ProbeSet::FindOverlaps(std::string chr, unsigned long int readstart, unsigned long int readend, std::string nameofdesign){
    
    std::vector<Interval< int > > results;
    
    Design[nameofdesign].Probe_Tree[chr].findOverlapping(readstart, readend, results);

    if (results.size() > 0){ // value = probe_index
        return results[0].value;
    }
    else
        return -1;
}

int ProbeSet::FindOverlaps_NegCtrls(std::string chr, unsigned long int readstart, unsigned long int readend, std::string nameofdesign){
    
    std::vector<Interval< int > > results;
    
    Design_NegCtrl[nameofdesign].Probe_Tree[chr].findOverlapping(readstart, readend, results);
    
    if (results.size() > 0){ // value = probe_index
        return results[0].value;
    }
    else
        return -1;
}


