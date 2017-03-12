#ifndef INTSUITE_INC_BIF_H_
#define INTSUITE_INC_BIF_H_
/*
Background interaction frequency is calculated 
* 
*/

#include "Probes.h"

class DetermineBackgroundLevels{
    friend class ProbeSet;
public:
	BG_signals bglevels;
	BG_signals bglevelsProbeProbe;
	
	void CalculateMeanandStdRegress(std::string, int, std::string, BG_signals&, int, std::string, int, OutStream&, int);

};

#endif  // INTSUITE_INC_BIF_H_
