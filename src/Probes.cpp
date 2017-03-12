#include "Probes.h"
#include "Global.h"
#include <fstream>


std::map < std::string, Probe_Design > Design;
std::map < std::string, Probe_Design > Design_NegCtrl;


void ProbeSet::GetProbeFeats(std::stringstream& line, CaptureProbes& t, std::string& Name){
    
    std::string promoter, snp, negctrl, other, SOterm;
    
    std::string field, field2, temp, attributes;
    
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
    getline(line,field,'\t');
    t.end = std::stoi(field);
    
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
     
    size_t smpos1 = attributes.find_first_of(";");
    if(smpos1!=std::string::npos){
		std::string parsename, parseside, parsetarget, parsedesign;
		
		parsename = attributes.substr(0, smpos1);
		
		size_t smpos2 = attributes.substr(smpos1+1, std::string::npos).find_first_of(";");
		
		parseside = attributes.substr(smpos1+1, smpos2);
		
		size_t smpos3 = attributes.substr(smpos2+1, std::string::npos).find_first_of(";");
		
		parsetarget=attributes.substr(smpos2+1, smpos3);
		
		size_t smpos4 = attributes.substr(smpos3+1, std::string::npos).find_first_of(";");
		
		parsedesign=attributes.substr(smpos3+1, smpos4);
		
		if(parsename.substr(0, parsename.find('=')).find("Name")!= std::string::npos){
			parsename.erase(std::remove_if(parsename.begin(), parsename.end(), [l](char ch) { return std::isspace(ch, l); }), parsename.end());
			Name=parsename.substr(parsename.find('=')+1);
			for (auto & c: Name) c = toupper(c); // Name of Feature
		}
		
		if(parseside.substr(0, parseside.find('=')).find("side")!= std::string::npos){
			parseside.erase(std::remove_if(parseside.begin(), parseside.end(), [l](char ch) { return std::isspace(ch, l); }), parseside.end());
			t.side=parseside.substr(parseside.find('=')+1); // Upstream (Left) or the downstream (Right) of the feature
		}
		
		if(parsetarget.substr(0, parsetarget.find('=')).find("target")!= std::string::npos){
			parsetarget.erase(std::remove_if(parsetarget.begin(), parsetarget.end(), [l](char ch) { return std::isspace(ch, l); }), parsetarget.end());
			std::string target = parsetarget.substr(parsetarget.find('=')+1); // Probe target for annotation 
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
		if(parsedesign.substr(0, parsedesign.find('=')).find("design")!= -1){
			parsedesign.erase(std::remove_if(parsedesign.begin(), parsedesign.end(), [l](char ch) { return std::isspace(ch, l); }), parsedesign.end());
			t.name_of_design=parsedesign.substr(parsedesign.find('=')+1);// Design Name 
		}
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
		temp.erase(std::remove_if(temp.begin(), temp.end(), [l](char ch) { return std::isspace(ch, l); }), temp.end()); 
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
                    std::stringstream probeline ( pline );
                    GetProbeFeats(probeline, tempprobe, Name);
                    chr1 = tempprobe.chr; //take the first chr outside the loop
                    getline(probefile, pline);
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
                getline(probefile, pline);
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
        
        for(auto it = probes_interval.begin(); it != probes_interval.end();++it){
				designname = it->first;
				AddTotheIntervalTree(probes_interval[designname], negctrl_probes_interval[designname], chr1, designname, isNeg);
		}
	
        chr1 = tempprobe.chr;
    }
	prLog << index <<  "   Probe Coordinates Read from "<< ProbeFileName<< std::endl;
    probefile.close();
    
    ++filesReadCount; 
 //////////////////////////////////////////////   how to get chr?what does this do
    //if(fileReadCount==fileCount){
	//}	
/////////////////////////////////////////////////////////
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
    
    //IntervalTree< int > tree, tree_negctrls;
    std::vector<Interval< int >  > temp, temp2;
    
    //tree = IntervalTree< int >(intervals);
    //tree_negctrls = IntervalTree< int >(intervals_negctrls);
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
    
    //results.clear();
    std::vector<Interval< int > > results;
    
    Design[nameofdesign].Probe_Tree[chr].findOverlapping(readstart, readend, results);

    if (results.size() > 0){ // value = probe_index
        return results[0].value;
    }
    else
        return -1;
}

int ProbeSet::FindOverlaps_NegCtrls(std::string chr, unsigned long int readstart, unsigned long int readend, std::string nameofdesign){
    
    //results.clear();
    std::vector<Interval< int > > results;
    
    Design_NegCtrl[nameofdesign].Probe_Tree[chr].findOverlapping(readstart, readend, results);
    
    if (results.size() > 0){ // value = probe_index
        return results[0].value;
    }
    else
        return -1;
}


