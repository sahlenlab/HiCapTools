/*  fasta.cpp
 
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
 
#ifndef HCT_INC_BIOIO_H_
#define HCT_INC_BIOIO_H_

#include <string>
#include <tuple>
#include "OutStream.h"
#include "bioio.hpp"

class bioioMod{
public :
	
	using GenomicRegion = std::tuple<std::string, size_t, size_t> ;
	
	GenomicRegion parse_region(std::string, const bioio::FastaIndex&);
	std::string GetFasta(std::string);
	bioioMod(OutStream& blog, std::string inpath, std::string faspath) : bLog (blog), index_path (inpath), fasta_path (faspath) {}
	
	
private:
	std::string index_path;
	std::string fasta_path;
	OutStream& bLog;
};

#endif // HCT_INC_BIOIO_H_
