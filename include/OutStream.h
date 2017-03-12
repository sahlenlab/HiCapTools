#ifndef PRSUITE_INC_UTSTR_H_
#define PRSUITE_INC_UTSTR_H_

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

#endif //PRSUITE_INC_UTSTR_H_
