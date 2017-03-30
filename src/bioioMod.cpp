/*
 
 Copyright (C) 2015 University of Oxford.
 
 Author: Daniel Cooke <dcooke@well.ox.ac.uk>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 DEALINGS IN THE SOFTWARE.  */

/*** 

	Modified and used from the bioio package
 
***/

#include "bioioMod.h"
#include <regex>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstddef>
#include <algorithm>

using GenomicRegion = std::tuple<std::string, size_t, size_t> ;

GenomicRegion bioioMod::parse_region(std::string region, const bioio::FastaIndex& index){
	
	region.erase(std::remove(region.begin(), region.end(), ','), region.end());
    
    const static std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    
    std::smatch match;
    
    if (std::regex_match(region, match, re) && match.size() == 5) {
        auto contig_name = match.str(1);
        
        if (index.count(contig_name) == 0) {
            throw std::runtime_error {"contig " + contig_name + " not found"};
        }
        
        const auto contig_size = index.at(contig_name).length;
        
        size_t begin {0}, end {0};
        
        if (match.str(2).empty()) {
            end = contig_size;
        } else {
            begin = static_cast<size_t>(std::stoull(match.str(2)));
            
            if (match.str(3).empty()) {
                end = begin + 1;
            } else if (match.str(4).empty()) {
                end = contig_size;
            } else {
                end = static_cast<size_t>(std::stoull(match.str(4)));
            }
            
            if (begin > contig_size) {
                throw std::runtime_error {"region " + region + " is larger than contig " + contig_name + ":0-" + std::to_string(contig_size)};
            }
            
            if (begin > end) {
                throw std::runtime_error {"begin position is past end position in region " + region};
            }
            
            if (end > contig_size) {
                end = contig_size;
            }
        }
        
        return GenomicRegion {std::move(contig_name), begin, end};
    }
    
    throw std::runtime_error {"could not parse region " + region};
}

std::string bioioMod::GetFasta(std::string regionToGet){
	
	  
    std::ifstream fasta {fasta_path, std::ios::binary};
    
    if (!fasta) {
        bLog << "Error: could not open fasta " << fasta_path << std::endl;
        return "Error";
    }
    
    if (index_path.empty()) { //changed
        index_path = fasta_path;
        index_path.replace(index_path.begin() + index_path.find_last_of("."), index_path.end(), ".fai");
    }
    
    std::ifstream index_file {index_path, std::ios::binary};
    
    if (!index_file) {
        index_file.open(fasta_path + ".fai");
        if (!index_file) {
            bLog << "Error: could not open index file, use samtools faidx <fasta> to make index file" << std::endl;
            return "Error";
        }
    }
    
    try {
        const auto index = bioio::read_fasta_index(index_file);
        
        const auto region  = parse_region(regionToGet, index);
        const auto& contig = std::get<0>(region);
        const auto begin   = std::get<1>(region);
        const auto length  = std::get<2>(region) - begin;
        
        const auto sequence = bioio::read_fasta_contig(fasta, index.at(contig), begin, length);
            
        return sequence;
    }
    catch (std::runtime_error& e) {
        bLog << "Error: " << e.what() << std::endl;
    }
}
