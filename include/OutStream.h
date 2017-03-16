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
//  OutStream.h
//  HiCapTools
//
//  Created by Pelin Sahlen and Anandashankar Anil.
//

#ifndef HCT_INC_UTSTR_H_
#define HCT_INC_UTSTR_H_

#include <iostream>

typedef std::ostream& (*ostream_manipulator)(std::ostream&);

class OutStream {
    std::ostream* conOut;
    std::ostream* fileOut;

public:
    void AddStreams(std::ostream* cOut, std::ostream* fOut){
        conOut=cOut;
        fileOut=fOut;
    }

    template <typename T> OutStream &operator<<(const T &t){
		*conOut << t;
		*fileOut << t;
		
		return *this;
    }
    
    OutStream& operator<<(ostream_manipulator fp){
		*conOut << fp;
		*fileOut << fp;
		
		return *this;
   }
};

#endif //HCT_INC_UTSTR_H_
